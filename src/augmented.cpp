/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "tools/dag.h"
#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "augmented.h"
#include "augmented/check.h"
#include "augmented/graph.h"
#include "data.h"
#include "data/household.h"
#include "data/individual.h"
#include "data/wide.h"
#include "distribution.h"
#include "function.h"
#include "infec_profile.h"
#include "mcmc/base.h"
#include "parameter.h"
#include "registry.h"
#include "translator.h"
#include "variable.h"

namespace housepi
{
    namespace infprof
    {
        void LookUpTable::init(const Distribution *ptr, size_t n)
        {
            generic_ptr_ = ptr;
            size_ = n;
            value_.assign(n, 0);
            old_.assign(n, 0);
            materialize_();
        }
        void LookUpTable::clear() noexcept
        {
            value_.clear();
            old_.clear();
        }
        void LookUpTable::recover() noexcept
        {
            if (generic_ptr_->is_constant())
                return;
            std::swap(value_, old_);
        }
        void LookUpTable::update()
        {
            if (generic_ptr_->is_constant())
                return;
            std::swap(value_, old_);
            materialize_();
        }

        void LookUpTable::Table_::clear() noexcept
        {
            inst.clear();
            cumul.clear();
        }
        void LookUpTable::Table_::assign(size_t n, double v)
        {
            inst.assign(n, v);
            cumul.assign(n, v);
        }

        void LookUpTable::materialize_()
        {
            for (size_t i = 1; i < size_; ++i)
            {
                double x = generic_ptr_->pdf(i);
                value_.inst[i] = x;
                value_.cumul[i] = x + value_.cumul[i - 1];
            }
        }
    }

    //////////////////////////////  InfecProfile  //////////////////////////////

    InfecProfile::InfecProfile() = default;
    InfecProfile::~InfecProfile() = default;

    void InfecProfile::build(Translator &translator,
                             const str_table &input)
    {
        latent_definitions.clear();

        if (input.empty())
            return;

        // check if definition is generic or latent
        if (count_rows(input) > 1 ||
            parse::Definition(input.at("func")[0]).consumes_latent() ||
            parse::Definition(input.at("min")[0]).consumes_latent() ||
            parse::Definition(input.at("max")[0]).consumes_latent() ||
            !na::is(input.at("status")[0]))
            build_latent_(translator, input);
        else
            build_generic_(translator, input);
    }

    void InfecProfile::build_generic_(Translator &translator,
                                      const str_table &input)
    {
        mode_ = Mode::Generic;
        const bool is_discr = is_required_discrete_ || str::tobool(input.at("discr")[0]);
        generic = helper_build_(translator, input.at("func")[0], is_discr);
        lookup_.clear();
        latent_definitions.clear();

        // domain minimum
        generic->init_attr("min", translator.to_attribute(input.at("min")[0]));

        // domain maximum
        if (input.count("max"))
            generic->init_attr("max", translator.to_attribute(input.at("max")[0]));

        // discretization and domain truncation
        generic->reshape(is_discr);

        // discretized generic profile: initialize the lookup table
        if (generic->is_discretized())
        {
            mode_ = Mode::LookUp;
            size_t n = infprof::LOOKUP_SIZE;
            if (generic->max().is_callable())
            {
                // ensure enough space when 'max' is de-located
                const double loc = generic->loc().is_callable() ? generic->loc().value() : 0.;
                n = std::max(n, static_cast<size_t>(1. + generic->max().value() + loc));
            }
            lookup_.init(generic.get(), n);
        }
    }

    void InfecProfile::build_latent_(Translator &translator,
                                     const str_table &input)
    {
        // drop generic definition
        mode_ = Mode::Latent;
        generic = nullptr;
        lookup_.clear();
        latent_definitions.clear();

        // match infectivity profile and latent variable definition
        for (size_t r = 0; r < count_rows(input); ++r)
        {
            // a discrete-time model discretizes whatever profile it is handed
            bool discr = is_required_discrete_ || str::tobool(input.at("discr")[r]);

            // check profile name and discretization
            const std::string &func_name = parse::Definition(input.at("func")[r]).root.func;
            auto probe = helper_build_(translator, func_name, discr);
            if (discr)
                probe->discretize();

            latent_definitions["name"].push_back(infprof::LATENT_FUNC);
            latent_definitions["status"].push_back(input.at("status")[r]);
            latent_definitions["func"].push_back(input.at("func")[r]);
            latent_definitions["fcumul"].push_back("false");
            latent_definitions["value"].push_back("");
            latent_definitions["const"].push_back("false");
            latent_definitions["discr"].push_back(str::format(discr));
            latent_definitions["model"].push_back("");
            latent_definitions["min"].push_back(input.at("min")[r]);
            latent_definitions["max"].push_back(input.at("max")[r]);
            latent_definitions["proposal"].push_back("");
            latent_definitions["kind"].push_back("latent");
        }
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Translator &translator,
                                                              const std::string &func,
                                                              bool discr)
    {
        parse::Definition parser(func);

        std::vector<registry::Kind> kinds = {registry::Kind::InfecProfile};
        for (registry::Kind r : registry::distr_kinds(discr))
            kinds.push_back(r);

        std::unique_ptr<Distribution> output = registry::create_distribution(
            parser.root.func, ErrType::InfecProfile, kinds);

        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, translator.to_attribute(node));
            translator.mark_param_produces_transm_cache(node);
        }

        return output;
    }

    //////////////////////////////  augm::DefTable  //////////////////////////////

    namespace augm
    {
        /// Map a definition's kind label to its typed Kind.
        static Kind string_to_kind(const std::string &kstr)
        {
            if (kstr == "member")
                return augm::Kind::Member;
            if (kstr == "house")
                return augm::Kind::House;
            if (kstr == "status")
                return augm::Kind::Status;
            // if (kstr == "latent")
            return augm::Kind::Latent;
        }

        DefRow::DefRow(const str_table &input, size_t r)
            : name(input.at("name")[r]),
              kind(string_to_kind(input.at("kind")[r])),
              value(input.at("value")[r]),
              func(input.at("func")[r]),
              model(input.at("model")[r]),
              min(input.at("min")[r]),
              max(input.at("max")[r]),
              proposal(input.at("proposal")[r]),
              has_value(!na::is(input.at("value")[r])),
              has_func(!na::is(input.at("func")[r])),
              has_model(!na::is(input.at("model")[r])),
              has_proposal(!na::is(input.at("proposal")[r])),
              is_const(str::tobool(input.at("const")[r])),
              is_discr(str::tobool(input.at("discr")[r])),
              is_fcumul(str::tobool(input.at("fcumul")[r]))
        {
            // household-level and status observations carry no status partition
            const std::string &sstr = input.at("status")[r];
            if ((kind == augm::Kind::Latent || kind == augm::Kind::Member) && !na::is(sstr))
                status = static_cast<StatusLvl>(str::toi(sstr, str::format(kind) + " status"));
        }

        std::vector<const parse::Definition *> DefRow::source_cols(bool simulation) const
        {
            // any row but an override or a flat row reads every column in both runs
            if (!is_sim_override() && !is_flat())
            {
                const auto all = source_cols();
                return {all.begin(), all.end()};
            }
            // simulation reads the value, inference the model and bounds
            if (simulation)
                return {&value, &func};
            return {&model, &func, &min, &max};
        }

        void DefTable::append(const str_table &input)
        {
            for (size_t r = 0; r < count_rows(input); ++r)
                rows.emplace_back(input, r);
        }

        bool DefTable::has_row(std::string_view name, std::optional<Kind> kind) const noexcept
        {
            for (const DefRow &row : rows)
                if (row.name == name && (!kind || row.kind == *kind))
                    return true;
            return false;
        }

        bool DefTable::has_continuous_declaration(std::string_view name) const noexcept
        {
            for (const DefRow &row : rows)
                if (row.name == name && (row.has_model || row.is_sampled()) && !row.is_discr)
                    return true;
            return false;
        }

        /// Whether a parsed definition subtree references any group risk.
        static bool entry_has_group_risk(const parse::Entry &node)
        {
            if (parse::is_group_risk(node.type))
                return true;
            for (const parse::Entry &child : node.children)
                if (entry_has_group_risk(child))
                    return true;
            return false;
        }

        bool DefTable::is_pure_data(size_t k) const
        {
            for (const parse::Definition *col : rows[k].source_cols())
            {
                if (col->is_na())
                    continue;
                if (col->consumes_latent() || col->consumes_parameter() ||
                    entry_has_group_risk(col->root))
                    return false;
            }
            return true;
        }
    }

    //////////////////////////////  augm::graph  //////////////////////////////

    namespace augm::graph
    {
        dag::prod_graph producers(const DefTable &def,
                                  const std::vector<size_t> &rows,
                                  bool simulation)
        {
            dag::prod_graph pgraph;
            for (size_t k : rows)
            {
                auto &producers = pgraph[def[k].name];
                for (const parse::Definition *col : def[k].source_cols(simulation))
                    for (const std::string &astr : col->list_latent_producers())
                        producers.insert(astr);
            }
            return pgraph;
        }

        void add_data_edges(dag::prod_graph &pgraph,
                            const DefTable &def,
                            const std::vector<size_t> &rows,
                            const std::vector<std::string> &targets)
        {
            for (size_t k : rows)
            {
                auto &producers = pgraph[def[k].name];
                for (const parse::Definition *col : def[k].source_cols(true))
                    for (const std::string &dstr : col->list_data_wide_producers())
                        if (dstr != def[k].name && vect::contains(targets, dstr))
                            producers.insert(dstr);
            }
        }

        str_table model_consumers(const DefTable &def,
                                  const std::vector<size_t> &order)
        {
            const dag::prod_graph pgraph = producers(def, order, false);

            // a flat row is keyed as its own consumer: its uniform log-ratio is zero, its log-pdf constant
            std::vector<std::string> without_model;
            for (size_t k : order)
                if ((!def[k].has_model && !def[k].is_flat()) || def[k].kind == Kind::Member)
                    without_model.push_back(def[k].name);

            // walk the whole graph, so a chain through a modelless node keeps its edge
            // then keep only the models whose density a move can shift
            // an unshifted consumer cancels in the ratio, a missing one does not
            str_table output;
            for (const auto &[key, prods] : pgraph)
                if (!vect::contains(without_model, key))
                    output.emplace(
                        key,
                        vect::difference(dag::transitive_consumers(key, pgraph, ErrType::Latent),
                                         without_model));

            return output;
        }

        /// List the names a run reads.
        /// - transmission names;
        /// - observations of `order` carrying a datum;
        /// - latents a household observation reduces over.
        static std::vector<std::string> sinks(const DefTable &def,
                                              const std::vector<size_t> &order,
                                              const std::vector<std::string> &transm,
                                              const std::vector<std::string> &absent)
        {
            // always sinks: the transmission side and the profile read no datum
            std::vector<std::string> output = transm;

            // observations of the scheme, except those this member has no datum for
            for (size_t k : order)
                if (def[k].kind != Kind::Latent && !vect::contains(absent, def[k].name))
                    output.push_back(def[k].name);

            // household observations sit in no status scheme, yet the latents they reduce over are read
            for (const DefRow &row : def.rows)
                if (row.kind == Kind::House)
                    for (const parse::Definition *col : row.source_cols())
                        for (const std::string &astr : col->list_latent_producers())
                            output.push_back(astr);

            return output;
        }

        /// Whether the run reads `astr`, directly or through the definitions consuming it.
        /// A sink reaches itself: transitive_consumers returns `start` first.
        static bool reaches_sink(const std::string &astr,
                                 const dag::prod_graph &pgraph,
                                 const std::vector<std::string> &sinks)
        {
            for (const std::string &cstr : dag::transitive_consumers(astr, pgraph, ErrType::Latent))
                if (vect::contains(sinks, cstr))
                    return true;
            return false;
        }

        std::vector<std::string> dead_names(const DefTable &def,
                                            const std::vector<size_t> &order,
                                            const std::vector<std::string> &transm,
                                            const std::vector<std::string> &absent)
        {
            const dag::prod_graph pgraph = producers(def, order, false);
            const std::vector<std::string> run_sinks = sinks(def, order, transm, absent);

            std::vector<std::string> dead;
            for (size_t k : order)
                if (def[k].kind == Kind::Latent && !reaches_sink(def[k].name, pgraph, run_sinks))
                    dead.push_back(def[k].name);
            return dead;
        }

        std::vector<std::string> transm_producers(const DefTable &def,
                                                  const std::vector<size_t> &order,
                                                  const std::vector<std::string> &transm)
        {
            // whatever a node consumes, transmission consumes
            const dag::prod_graph pgraph = producers(def, order, false);
            std::vector<std::string> output;
            for (const std::string &node : transm)
                if (pgraph.count(node)) // a root needs no definition in every block
                    vect::append(output, dag::transitive_producers(node, pgraph, ErrType::Latent));
            return output;
        }
    }

    //////////////////////////////  augm::check  //////////////////////////////

    namespace augm::check
    {
        std::string undeclared_status_message(StatusLvl err_status)
        {
            if (status::is_missing(err_status))
                return "No generic latent definition was declared.\n"
                       "Provide a definition without " +
                       msg::arg("status") + " to apply it to every infected status.";

            return "No latent definition was declared for status " + std::to_string(err_status) + ".\n" +
                   "Try to associate a status-specific definition.\n" +
                   "In alternative, drop " + msg::arg("status") + " from existing definitions to make them generic.";
        }

        void producers_declared(const DefTable &def, const dag::prod_graph &pgraph, StatusLvl s)
        {
            for (const auto &node : pgraph)
                for (const std::string &astr : node.second)
                {
                    if (pgraph.count(astr))
                        continue;

                    // declared for another status: report it
                    if (def.has_row(astr))
                    {
                        if (status::is_missing(s))
                            msg::stop("No generic latent definition was found for " + msg::val(astr) + ".");
                        else
                            msg::stop("No latent definition was found for " + msg::val(astr) +
                                      " for status " + std::to_string(s) + ".");
                    }
                    msg::stop("Could not find a latent definition for " + msg::val(astr) + ".\n" +
                              "Declare it with " + msg::fn("add_latent") + ".");
                }
        }

        void status_obs_consistent(const DefTable &def, size_t row)
        {
            if (!def.is_pure_data(row))
                msg::stop(
                    "Probabilities of " + msg::fn("add_status_process") +
                    " must depend only on data or numeric constants.\n" +
                    "Inference holds the observed status fixed and does not estimate them.");
        }

        /// List the data names (source D) read by `def` that are in `names`.
        static std::vector<std::string> referred_data(const parse::Definition &def,
                                                      const std::vector<std::string> &names)
        {
            std::vector<std::string> output = {};
            if (def.is_na())
                return output;
            for (const std::string &dstr : def.list_data_wide_producers())
                if (vect::contains(names, dstr))
                    output.push_back(dstr);
            return output;
        }

        void every_member_cols(const DefTable &def,
                               const data::Wide &data_wide,
                               const std::vector<size_t> &per_infection)
        {
            for (size_t k = 0; k < def.size(); ++k)
            {
                const std::string &name = def[k].name;
                if (def[k].kind != Kind::Member ||
                    !data_wide.count(name) ||
                    !data_wide.is_required_everyone(data_wide.consult(name)))
                    continue;

                // drawn, or computed from a draw, with the infection of each member
                if (def[k].has_model || vect::contains(per_infection, k))
                    msg::stop(
                        "Data " + msg::val(name) + " needs a value for every member, escaped ones included.\n" +
                        msg::alert(str::format(def[k].kind) + " " + msg::val(name) +
                                   " has a value only for infected members.") +
                        "\nDefine it at household level with " + msg::code("level = \"house\"") +
                        ", or from data columns only.");
            }
        }

        void no_live_house_obs_cycle(const DefTable &def,
                                     const std::vector<std::string> &live_names)
        {
            if (live_names.empty())
                return;

            for (const DefRow &row : def.rows)
            {
                if (row.kind != Kind::Latent)
                    continue;

                // an override skips its model in simulation, see Binder::bind_member()
                // a min, max or model cycling through a live observation is never evaluated, nor flagged
                for (const parse::Definition *col : row.source_cols(true))
                    for (const std::string &dstr : referred_data(*col, live_names))
                        msg::stop(
                            "Found cyclic dependency in simulation.\n" +
                            msg::alert("Latent " + msg::val(row.name) +
                                       " reads generated household observation " + msg::val(dstr) + ".") +
                            "\nDefine " + msg::val(row.name) +
                            " without referencing " + msg::val(dstr) + ", or viceversa.");
            }
        }

        void no_obs_entry_collision(const DefTable &def,
                                    const std::vector<std::string> &obs_names)
        {
            if (obs_names.empty())
                return;

            // producer graph (name -> the latent names it reads), over latent + obs rows
            std::vector<size_t> rows;
            for (size_t k = 0; k < def.size(); ++k)
                if (def[k].kind == Kind::Latent || def[k].kind == Kind::Member)
                    rows.push_back(k);
            const dag::prod_graph pgraph = graph::producers(def, rows, false);

            for (size_t k : rows)
            {
                const std::string &name = def[k].name;

                // the two rules partition DefRow::source_cols() by role:
                // - model support: model, min, max;
                // - deterministic value: value, func.
                // keep them in sync with source_cols()

                // rule 1: a latent's model (distribution and support) reads no observation
                // its producers are parameters and other latents only
                // an observation's model may read another observation, drawn before it
                if (def[k].kind == Kind::Latent && def[k].has_model)
                    for (const parse::Definition *col : {&def[k].model, &def[k].min, &def[k].max})
                        for (const std::string &dstr : referred_data(*col, obs_names))
                            msg::stop(
                                "Definition collision: the model of " + msg::val(name) +
                                " reads observation " + msg::val(dstr) +
                                ".\nA latent model cannot depend on observed or generated values.");

                // rule 2: a deterministic definition may read an observation, as a downstream readout
                // unless it closes a loop: the observation is also its consumer
                for (const parse::Definition *col : {&def[k].value, &def[k].func})
                    for (const std::string &dstr : referred_data(*col, obs_names))
                    {
                        if (dstr == name)
                            continue; // observation node reading its own entry column
                        // loop iff observation `d` transitively consumes `name`
                        const auto consumers = dag::transitive_consumers(name, pgraph, ErrType::Latent);
                        if (vect::contains(consumers, dstr))
                            msg::stop(
                                "Cyclic definition: " + msg::val(name) +
                                " reads observation " + msg::val(dstr) +
                                ", but " + msg::val(dstr) + " depends on " + msg::val(name) +
                                ".\nBreak the loop by defining one of them independently of the other.");
                    }
            }
        }

        void flat_rows(const DefTable &def,
                       const Parameter &param,
                       const std::vector<std::string> &transm)
        {
            for (const DefRow &row : def.rows)
            {
                if (row.kind != Kind::Latent)
                    continue;
                const bool is_transm = vect::contains(transm, row.name);

                if (row.has_model && row.is_sampled() && is_transm)
                    msg::stop(
                        "Latent " + msg::val(row.name) +
                        " is generated by the transmission model, which already gives its probability.\n" +
                        "Drop " + msg::arg("model") + ", and keep " + msg::arg("min") + " and " +
                        msg::arg("max") + " to sample it within.");

                if (!row.is_flat())
                    continue;

                if (!is_transm)
                    msg::stop(
                        "Latent " + msg::val(row.name) +
                        " has no model, and the transmission model does not generate it.\n" +
                        "Provide a " + msg::arg("model") + ".");

                // the flat density 1 / (max - min) is left out: the bounds must not move with the draws
                for (const auto &[key, col] : {std::pair{"min", &row.min}, std::pair{"max", &row.max}})
                {
                    std::string read;
                    for (const std::string &pstr : col->list_parameter_producers())
                        if (!param.is_constant(pstr))
                            read = "parameter " + msg::val(pstr);
                    for (const std::string &astr : col->list_latent_producers())
                        read = "latent " + msg::val(astr);
                    if (read.empty() && entry_has_group_risk(col->root))
                        read = "a group risk";
                    if (!read.empty())
                        msg::stop(
                            "Bound " + msg::arg(key) + " of latent " + msg::val(row.name) + " reads " + read +
                            ".\nWithout a model, bounds can read only data and observations.");
                }
            }
        }

        void no_sim_override_obs_read(const DefTable &def,
                                      const std::vector<std::string> &obs_names)
        {
            // a simulation override (value of a sampled latent) is applied when its variable is first drawn
            // sample_obs runs after draw_sample: an observation reached from the override is still undrawn
            // an undrawn observation reads zero and leaves no trace
            // deterministic definitions resolve through: the chain is walked, not the override's own column
            // the walk stops at a sampled latent, whose read yields its own draw
            if (obs_names.empty())
                return;

            // producer graph over the column each definition resolves through in simulation
            dag::prod_graph pgraph;
            for (const DefRow &row : def.rows)
                if (row.kind == Kind::Latent)
                {
                    auto &producers = pgraph[row.name];
                    for (const std::string &astr : row.sim_source().list_latent_producers())
                        producers.insert(astr);
                }

            for (const DefRow &row : def.rows)
            {
                if (row.kind != Kind::Latent || !row.is_sim_override())
                    continue;

                for (const std::string &pstr :
                     dag::transitive_producers(row.name, pgraph, ErrType::Latent))
                    // a name spans several rows when statuses declare it apart
                    for (const DefRow &prow : def.rows)
                    {
                        if (prow.kind != Kind::Latent || prow.name != pstr)
                            continue;
                        for (const std::string &dstr : referred_data(prow.sim_source(), obs_names))
                            msg::stop(
                                "Simulation override of " + msg::val(row.name) +
                                " reads observation " + msg::val(dstr) +
                                (pstr == row.name ? "" : " through " + msg::val(pstr)) +
                                ".\nThe value of a sampled latent is applied during simulation, before " +
                                msg::val(dstr) + " is generated; depend on data, parameters, or " +
                                "latents that do not read it.");
                    }
            }
        }

        [[noreturn]] void obs_init_error(const DefTable &def,
                                         const std::vector<size_t> &rows,
                                         const Parameter &param,
                                         const std::string &context)
        {
            for (size_t k : rows)
            {
                // domain: loc + [min, max]
                std::vector<std::string> latent_producers, param_producers;
                for (const parse::Entry &child : def[k].model.root.children)
                    if (child.key == "loc")
                    {
                        child.append_latent_producers(latent_producers);
                        child.append_parameter_producers(param_producers);
                    }
                for (const parse::Definition *col : {&def[k].min, &def[k].max})
                {
                    vect::append(latent_producers, col->list_latent_producers());
                    vect::append(param_producers, col->list_parameter_producers());
                }

                if (latent_producers.size() > 0)
                    msg::stop(
                        "Exceeded attempts to initialize latent variables of " +
                        context + " due to observations.\n" +
                        "Widen latent model of " + str::format(str::unique(latent_producers), false, 3) +
                        ", revise the observation domain of " + msg::val(def[k].name) +
                        ", or check the observed individual data.");

                // non-constant parameter in domain cannot be fixed by per-individual resampling
                for (const std::string &pstr : str::unique(param_producers))
                    if (!param.is_constant(pstr))
                        msg::stop(
                            "Failed to initialize due to observations of " + context + ".\n" +
                            "Try to remove free parameter dependencies in observation domain of " + msg::val(def[k].name) + ".");
            }

            // otherwise the observed data are unreachable
            msg::stop(
                "Failed to initialize due to observations of " + context + ".\n" +
                msg::alert("The observed values lie outside of respective domain.") +
                "\nCheck the individual data or revise the domain.");
        }

        void report_missing_obs(const KeyHolder &keys,
                                const std::vector<size_t> &declared,
                                const std::vector<size_t> &missing)
        {
            for (size_t a = 0; a < missing.size(); ++a)
            {
                if (missing[a] == 0)
                    continue;

                if (missing[a] == declared[a])
                    msg::stop(
                        "Observation " + msg::val(keys.consult(a)) +
                        " is missing for every infected individual.\n"
                        "An observation model needs at least one observed value: check the data column.");

                msg::warn(
                    "Observation " + msg::val(keys.consult(a)) + " is missing for " +
                    std::to_string(missing[a]) + " of " + std::to_string(declared[a]) +
                    " infected individuals.\n" +
                    msg::alert("A missing observation does not contribute to the likelihood."));
            }
        }
    }

    //////////////////////////////  augm::Scheme  //////////////////////////////

    namespace augm
    {
        size_t Scheme::index_of(StatusLvl s) const noexcept
        {
            auto it = status_index_.find(s);
            if (it != status_index_.end())
                return it->second;
            // otherwise fall back to the generic definition, if any
            it = status_index_.find(status::MISSING);
            return it != status_index_.end() ? it->second : na::INDEX;
        }

        void SchemeBuilder::register_rows()
        {
            for (size_t k = 0; k < def_.size(); ++k)
            {
                // skip household observations: they are not per-member
                if (def_[k].kind == Kind::House)
                    continue;

                // status observation is per-member, but out of status scheme
                if (def_[k].kind == Kind::Status)
                {
                    scheme_.status_obs_row_ = k;
                    scheme_.status_obs_key_ = keys_.consult(def_[k].name);
                    continue;
                }

                // store status -> definition rows
                scheme_.status_rows_[def_[k].status].push_back(k);
            }
        }

        //////  augm::SchemeBuilder : order  //////

        std::vector<std::string> SchemeBuilder::order(const std::vector<StatusLvl> &infected,
                                                      const Parameter &param)
        {
            const auto &rows = scheme_.status_rows_;
            auto &index = scheme_.status_index_;
            Sources_ src{param};

            // skip if no member definition goes into a status scheme
            // household observations need no status scheme
            if (rows.empty())
                return {};

            // load status-specific user definitions first
            // NA, the generic status, needs no check here
            for (const auto &[s, vec] : rows)
                if (status::is_infected(s))
                    order_status_(s, s, src);

            // associate each observed status with its scheme, else with the generic one
            for (StatusLvl s : infected)
                if (!index.count(s))
                {
                    if (!index.count(status::MISSING))
                        order_status_(status::MISSING, s, src);
                    index.emplace(s, index.at(status::MISSING));
                }

            // no observed infected status or status-specific definition: load the generic definition
            if (infected.empty() ||
                (rows.size() == 1 && status::is_missing(rows.begin()->first)))
                // the generic definition serves any status left without a scheme
                // load it when declared, else the call reports it missing
                // a status drawn in simulation reports its gap then
                // see Augmented::helper_materialize_definition_order_()
                if (rows.count(status::MISSING) ||
                    !scheme_.is_any_status_obs() ||
                    index.empty())
                    order_status_(status::MISSING, status::MISSING, src);

            scheme_.n_base_ = scheme_.status_schemes_.size();
            return std::move(src.covariates);
        }

        void SchemeBuilder::order_status_(StatusLvl s,
                                          StatusLvl err_status,
                                          Sources_ &src)
        {
            auto &index = scheme_.status_index_;
            // skip a status ordered before
            if (index.count(s))
                return;
            // throw if not defined
            if (!scheme_.status_rows_.count(s))
                msg::stop(check::undeclared_status_message(err_status));

            const std::vector<size_t> rows = unsorted_rows_(s);
            dict<size_t> name_to_row; // reverse dictionary name-to-row-index
            std::vector<std::string> member_obs;
            for (size_t k : rows)
            {
                name_to_row[def_[k].name] = k;
                if (def_[k].kind == Kind::Member)
                    member_obs.push_back(def_[k].name);
            }

            // producer graph of one run, cycles reported by name
            // in simulation, a member observation is drawn before the definitions reading its data column
            // in inference, they read a fixed data cell
            const auto sorted_rows = [&](bool simulation)
            {
                dag::prod_graph pgraph = graph::producers(def_, rows, simulation);
                check::producers_declared(def_, pgraph, err_status);
                if (simulation)
                    graph::add_data_edges(pgraph, def_, rows, member_obs);
                std::vector<size_t> order;
                for (const std::string &astr : dag::order_producers_first(pgraph, ErrType::Latent))
                    order.emplace_back(name_to_row[astr]);
                return order;
            };

            // assert acyclicity and translate order from name to deftable row index
            try
            {
                // called once per status
                index.emplace(s, scheme_.status_schemes_.size());
                scheme_.status_schemes_.emplace_back();
                at_(s).order = sorted_rows(false);

                // define usage lists respecting definition order
                for (size_t k : at_(s).order)
                    list_definition_(s, k, src);

                // an override reads its value in simulation: its producers differ from inference
                try
                {
                    at_(s).sim_order = sorted_rows(true);
                }
                catch (const std::exception &e)
                {
                    msg::stop(std::string(e.what()) +
                              "\nSimulation reads the value of a sampled latent in place of its model and bounds.");
                }
                for (size_t k : at_(s).sim_order)
                {
                    const size_t a = keys_.consult(def_[k].name);
                    if (vect::contains(at_(s).sample, a))
                        at_(s).sim_sample.push_back(a);
                }
            }
            catch (const std::exception &e)
            {
                std::string sms = e.what();
                if (status::is_missing(err_status))
                    sms += "\nCheck your generic definitions.";
                else
                    sms += "\nCheck status-specific definitions for status " + std::to_string(err_status) + ".";
                msg::stop(sms);
            }
        }

        std::vector<size_t> SchemeBuilder::unsorted_rows_(StatusLvl s) const
        {
            const auto &rows = scheme_.status_rows_;

            // load the list of status-specific definition rows
            std::vector<size_t> allrows = rows.at(s);

            // extract names
            std::vector<std::string> name_to_status = {};
            for (size_t k : allrows)
                name_to_status.push_back(def_[k].name);

            // add generic definitions to the list, skipping names a status-specific one declares
            if (status::is_infected(s) && rows.count(status::MISSING))
                for (size_t k : rows.at(status::MISSING))
                    if (!vect::contains(name_to_status, def_[k].name))
                        allrows.push_back(k);

            return allrows;
        }

        void SchemeBuilder::list_definition_(StatusLvl s,
                                             size_t k,
                                             Sources_ &src)
        {
            const std::string &astr = def_[k].name;
            size_t a = keys_.consult(astr);

            bool is_to_export = true;
            bool is_to_draw = true;

            if (def_[k].is_const)
            {
                is_to_draw = false;
            }
            // deterministic definition: a function, or a value neither sampled nor constant
            // a value paired with a model or a proposal is a simulation-only override, see Binder::bind_member()
            else if (!def_[k].is_sampled())
            {
                if (!def_[k].has_func && !def_[k].has_value)
                    msg::stop(
                        "Impossible to initialize latent " + msg::val(astr) + " without a value or model");
                std::tie(is_to_export, is_to_draw) = list_function_(s, k, src);
            }
            else if (!def_[k].has_proposal)
            {
                msg::stop(
                    "No proposal was provided for non-constant definition of latent " +
                    msg::val(astr));
            }
            else
            {
                scheme_.is_any_sample_ = true;
                at_(s).sample.emplace_back(a);

                if (parse::is_adaptive_proposal(def_[k].proposal))
                {
                    scheme_.is_any_adaptive_ = true;
                    at_(s).adaptive.emplace_back(a);
                }
            }

            // observation exports are listed separately
            // the profile latent is an engine artefact, reported through the profile itself
            if (def_[k].kind == Kind::Member || astr == infprof::LATENT_FUNC)
                is_to_export = false;

            // nothing to do if not to export
            if (!is_to_export)
                return;

            if (!vect::contains(scheme_.list_latent_export, a))
                scheme_.list_latent_export.emplace_back(a);

            if (is_to_draw)
                at_(s).draws.emplace_back(a);
        }

        std::pair<bool, bool> SchemeBuilder::list_function_(StatusLvl s,
                                                            size_t k,
                                                            Sources_ &src)
        {
            bool is_to_export = true;
            bool is_to_draw = true;

            // a value definition has priority over a function
            const parse::Definition &parser = def_[k].has_value ? def_[k].value : def_[k].func;

            // raw data are equivalence-class covariates, unless produced
            for (const std::string &dstr : parser.list_data_wide_producers())
                if (!def_.has_row(dstr, Kind::Member) && !def_.has_row(dstr, Kind::House))
                    src.covariates.push_back(dstr);

            // do not export if not callable: it would crash
            if (!parser.is_callable_implicitly("x"))
                is_to_export = false;

            bool all_consts = true; // start assuming all constants
            for (const std::string &pstr : parser.list_parameter_producers())
                all_consts &= src.param.is_constant(pstr);

            // inspect status-specific latent definitions: draw unless constant
            // relies on the definition order
            for (const std::string &bstr : parser.list_latent_producers())
                if (vect::contains(at_(s).draws, keys_.consult(bstr)))
                    all_consts = false;

            // to draw if not all constant dependencies
            is_to_draw = !all_consts;

            return {is_to_export, is_to_draw};
        }

        //////  augm::SchemeBuilder : surrogates and nodes  //////

        void SchemeBuilder::init_surrogates()
        {
            // look for surrogate column definitions: deterministic definitions standing in for a data column
            dict<size_t> name_to_row;
            for (size_t k = 0; k < def_.size(); ++k)
                // household reduction or member-level definition carrying no model
                if (def_[k].kind == Kind::House ||
                    (def_[k].kind == Kind::Member && !def_[k].has_model))
                    name_to_row[def_[k].name] = k;

            if (name_to_row.empty())
                return;

            // order producer-before-consumer surrogate columns
            std::vector<size_t> rows;
            std::vector<std::string> names;
            for (const auto &[nm, k] : name_to_row)
                rows.push_back(k), names.push_back(nm);
            dag::prod_graph pgraph;
            graph::add_data_edges(pgraph, def_, rows, names);

            // extract member-level live observations (with a model)
            std::unordered_set<std::string> live_names;
            for (const DefRow &row : def_.rows)
                if (row.kind == Kind::Member && row.has_model)
                    live_names.insert(row.name);

            // build list of surrogate data columns in order
            for (const std::string &nm : dag::order_producers_first(pgraph, ErrType::Observation))
            {
                Surrogate surr;
                surr.row = name_to_row.at(nm);
                surr.is_house = def_[surr.row].kind == Kind::House;
                surr.is_static = def_.is_pure_data(surr.row);

                // update household observations
                if (surr.is_house)
                {
                    surr.slot = scheme_.list_obs_house.size();
                    scheme_.list_obs_house.push_back(surr.row);
                }

                // live if it depends on a live column
                if (surr.is_static)
                    for (const parse::Definition *col : def_[surr.row].source_cols())
                        for (const std::string &dstr : col->list_data_wide_producers())
                            if (live_names.count(dstr))
                                surr.is_static = false;
                if (!surr.is_static)
                    live_names.insert(nm);

                scheme_.list_surrogate.push_back(std::move(surr));
            }
        }

        void SchemeBuilder::init_observations()
        {
            for (StatusScheme &ss : scheme_.status_schemes_)
            {
                ss.obs_consumers.resize(keys_.size());

                // an observation reaches its sampled producers through deterministic nodes
                // the consumer map follows the transitive closure, as the latent-model map does
                // otherwise moving such a producer would skip the observation ratio
                dag::prod_graph pgraph = graph::producers(def_, ss.order, false);

                // simulation order: an observation is drawn after the observations it reads
                for (size_t k : ss.sim_order)
                {
                    if (def_[k].kind != Kind::Member)
                        continue;

                    // a surrogate column is data, not drawn
                    bool skip = false;
                    for (const Surrogate &surr : scheme_.list_surrogate)
                        if (surr.row == k)
                            skip = surr.is_static;
                    if (skip)
                        continue;

                    const size_t obs_key = keys_.consult(def_[k].name);

                    // unique observation keys to return as simulated columns
                    if (!vect::contains(scheme_.list_obs_member, obs_key))
                        scheme_.list_obs_member.push_back(obs_key);

                    // a deterministic observation node carries no distribution
                    if (!def_[k].has_model)
                        continue;

                    scheme_.is_any_obs_model_ = true;
                    ss.obs_models.push_back(obs_key);

                    // list the observation under each of its producers, so moving any of them recomputes it
                    for (const std::string &astr : dag::transitive_producers(def_[k].name, pgraph, ErrType::Latent))
                        if (astr != def_[k].name) // transitive_producers returns 'start' last
                            ss.obs_consumers[keys_.consult(astr)].push_back(obs_key);
                }
            }
        }

        void SchemeBuilder::init_transm(const std::vector<std::string> &transm)
        {
            for (StatusScheme &ss : scheme_.status_schemes_)
            {
                ss.produces_transm.assign(keys_.size(), false);
                // declaring no transmission key skips the computation
                if (transm.empty())
                    continue;

                for (const std::string &pstr : graph::transm_producers(def_, ss.order, transm))
                    ss.produces_transm[keys_.consult(pstr)] = true;
            }
        }

        //////  augm::SchemeBuilder : evaluation  //////

        void SchemeBuilder::init_models()
        {
            for (StatusScheme &ss : scheme_.status_schemes_)
            {
                ss.consumers.resize(keys_.size());

                // index the consumer lists by their producer key
                const str_table model_consumers = graph::model_consumers(def_, ss.order);
                for (const auto &[astr, consumers] : model_consumers)
                    ss.consumers[keys_.consult(astr)] = keys_.consult(consumers);

                // a flat row adds no density to the latent log-likelihood
                std::vector<std::string> with_model = extract_keys(model_consumers), flat;
                for (size_t k : ss.order)
                    if (def_[k].is_flat())
                        flat.push_back(def_[k].name);
                vect::remove(with_model, flat);
                ss.latent_models = keys_.consult(with_model);
                scheme_.is_any_latent_model_ |= !ss.latent_models.empty();
            }
        }

        size_t SchemeBuilder::make_variant(size_t base,
                                           const std::vector<size_t> &absent,
                                           const std::vector<std::string> &transm)
        {
            StatusScheme ss = scheme_.status_schemes_[base]; // copy: the store reallocates below
            const std::vector<size_t> dead = keys_.consult(
                graph::dead_names(def_, ss.order, transm, keys_.consult(absent)));

            vect::remove(ss.obs_models, absent);
            vect::remove(ss.draws, dead);
            vect::remove(ss.sample, dead);
            vect::remove(ss.adaptive, dead);
            vect::remove(ss.latent_models, dead);

            std::vector<size_t> order;
            for (size_t k : ss.order)
                if (!vect::contains(dead, keys_.consult(def_[k].name)))
                    order.push_back(k);
            ss.order = std::move(order);

            // a dead key shifts no density and an absent observation reads none
            for (std::vector<size_t> &vec : ss.consumers)
                vect::remove(vec, dead);
            for (std::vector<size_t> &vec : ss.obs_consumers)
                vect::remove(vec, absent);
            for (size_t a : dead)
                ss.consumers[a].clear(), ss.obs_consumers[a].clear();

            scheme_.status_schemes_.push_back(std::move(ss));
            return scheme_.status_schemes_.size() - 1;
        }

        void SchemeBuilder::init_param_producers(const KeyHolder &param)
        {
            // mark parameters producing latent or observation models
            scheme_.param_produces_latent_ll_.assign(param.size(), false);
            scheme_.param_produces_obs_ll_.assign(param.size(), false);

            std::vector<bool> latent_produces_obs(keys_.size(), false);
            for (const StatusScheme &ss : scheme_.status_schemes_)
                for (size_t a = 0; a < ss.obs_consumers.size(); ++a)
                    if (!ss.obs_consumers[a].empty())
                        latent_produces_obs[a] = true;

            for (const DefRow &row : def_.rows)
            {
                if (row.name == infprof::LATENT_FUNC ||
                    row.kind == Kind::House ||
                    row.kind == Kind::Status)
                    continue;

                const bool is_obs = row.kind == Kind::Member; // else Kind::Latent
                // a deterministic latent forwards its parameters to the observation it produces
                const bool to_obs = is_obs ||
                                    (!row.is_sampled() &&
                                     latent_produces_obs[keys_.consult(row.name)]);

                for (const parse::Definition *col : row.source_cols(false))
                    for (size_t p : param.consult(col->list_parameter_producers()))
                    {
                        if (!is_obs)
                            scheme_.param_produces_latent_ll_[p] = true;
                        if (to_obs)
                            scheme_.param_produces_obs_ll_[p] = true;
                    }
            }
        }
    }

    //////////////////////////////  augm::State  //////////////////////////////

    namespace augm
    {
        void State::allocate(size_t n_indiv,
                             const KeyHolder &keys,
                             const DefTable &def)
        {
            blank.assign(keys.size(), Variable());
            for (size_t a = 0; a < blank.size(); ++a)
                if (def.has_row(keys.consult(a), Kind::Member))
                    blank[a].name(keys.consult(a), ErrType::Observation);
                else
                    blank[a].name(keys.consult(a), ErrType::Latent);

            varia.resize(n_indiv);
            for (IndivState &st : varia)
                st.vec = blank;
        }

        double commit_latent_log_lik(const StatusScheme &ss,
                                     IndivState &st)
        {
            if (ss.latent_models.empty())
                return 0.;
            double new_ll = 0.0;
            // guarded: a moved parameter feeding loc/min/max may push a value out of its shifted support
            // that value reads -Inf, not a stale finite density
            for (size_t a : ss.latent_models)
                new_ll += st.vec[a].guarded_log_pdf();
            st.latent_log_lik.set(new_ll);
            return new_ll;
        }

        double commit_latent_log_lik_ratio(const StatusScheme &ss,
                                           IndivState &st,
                                           size_t a)
        {
            double delta_ll = 0;
            // consumers[a] includes `a`: the model of `a` is counted
            for (size_t b : ss.consumers[a])
                delta_ll += st.vec[b].log_ratio();
            st.latent_log_lik.increment(delta_ll);
            return delta_ll;
        }

        double sum_guarded_latent_log_lik(const StatusScheme &ss,
                                          const IndivState &st,
                                          size_t a)
        {
            double ll = 0.;
            // consumers[a] includes `a`: the model of `a` is counted
            for (size_t b : ss.consumers[a])
                ll += st.vec[b].guarded_log_pdf();
            return ll;
        }

        double commit_latent_log_lik_ratio_with_consumers(const StatusScheme &ss,
                                                          IndivState &st,
                                                          size_t a,
                                                          double L_old)
        {
            // with `a` proposed, its model consumers read the moved value live
            // compare their guarded model sum against the pre-move snapshot `L_old`
            const double delta_ll = sum_guarded_latent_log_lik(ss, st, a) - L_old;
            st.latent_log_lik.increment(delta_ll);
            return delta_ll;
        }

        double eval_log_obs(const StatusScheme &ss,
                            const IndivState &st)
        {
            double ll = 0.0;
            for (size_t a : ss.obs_models)
                ll += st.vec[a].guarded_log_pdf();
            return ll;
        }

        double commit_obs_log_lik(const StatusScheme &ss,
                                  IndivState &st)
        {
            if (ss.obs_models.empty())
                return 0.;
            const double new_ll = eval_log_obs(ss, st);
            st.obs_log_lik.set(new_ll);
            return new_ll;
        }

        double commit_obs_log_lik_ratio(const StatusScheme &ss,
                                        IndivState &st,
                                        size_t a)
        {
            // skip unless the moved latent `a` parameterizes an observation node
            if (ss.obs_consumers[a].empty())
                return 0.;

            // recompute the log-proba: an observation value is a datum
            const double new_ll = eval_log_obs(ss, st);
            const double delta_ll = new_ll - st.obs_log_lik.value();
            st.obs_log_lik.set(new_ll);
            return delta_ll;
        }

        void restart(const StatusScheme &ss, IndivState &st, rnd_gen &prng)
        {
            for (size_t a : ss.sample)
                st.vec[a].restart(prng);
        }
        void draw_sample(const std::vector<size_t> &keys, IndivState &st, rnd_gen &prng)
        {
            for (size_t a : keys)
                st.vec[a].draw_value(prng);
        }
    }

    //////////////////////////////  augm::Binder  //////////////////////////////

    namespace augm
    {
        void Binder::bind_variable(Variable &var,
                                   size_t k,
                                   const Individual &anchor,
                                   const Household *house) const
        {
            const DefRow &row = def_[k];
            var.name(row.name, kind_to_err_type(row.kind));

            // a deterministic definition is driven by its value, anything else by its function
            // the value of a sampled row is a simulation-only override, applied before this
            const parse::Definition &parser = row.is_sampled() || !row.has_value ? row.func : row.value;

            // a household definition binds no member: its per-member cells fold over them
            const Individual *indiv_ptr = row.kind == Kind::House ? nullptr : &anchor;

            // pointer mode: a bare reference points at what it names
            if (parser.is_terminal())
            {
                var.transfer(translator_.to_attribute(parser.root, indiv_ptr, house));
                return;
            }

            // function mode: compose it once
            std::unique_ptr<Function> func = registry::create_function(
                parser.root.func,
                ErrType::Function,
                registry::func_distr_kinds(row.is_discr),
                row.is_fcumul);
            if (translator_.is_house_fold(parser.root, house))
                func->set_eval_as_fold();
            var.init_function(std::move(func));

            // a member surrogate materializes a data column for every member
            // where the cells it reads are empty, it has nothing to stand for
            const bool is_surrogate = row.kind == Kind::Member && !row.has_model;

            bool is_defined = true;
            if (row.kind == Kind::House)
                is_defined = translator_.bind_arguments(var, parser.root, nullptr, house);
            else
                for (const parse::Entry &child : parser.root.children)
                {
                    const Attribute &attr = translator_.to_attribute(child, &anchor);
                    if (is_surrogate && attr.is_constant() && attr.is_na())
                    {
                        is_defined = false;
                        break;
                    }
                    var.init_func_attr(child.key, attr);
                    translator_.mark_param_produces_transm_cache(child);
                }

            // no operand carries a value: the definition yields NA here
            if (!is_defined)
            {
                var.init_numeric(na::DOUBLE);
                return;
            }

            // domain and reshaping
            var.init_func_attr("min", translator_.to_attribute(row.min.root, &anchor));
            var.init_func_attr("max", translator_.to_attribute(row.max.root, &anchor));
            var.reshape_function(row.is_discr);
        }

        void Binder::bind_member(Variable &var, size_t k,
                                 const Individual &indiv,
                                 bool simulation) const
        {
            const DefRow &row = def_[k];

            try
            {
                if (row.is_const)
                {
                    // placeholder for observation node if generated by simulation
                    if (simulation && (row.kind == Kind::Member || row.kind == Kind::Status))
                        var.init_numeric(0.);
                    else
                        var.init_numeric(translator_.to_numeric(row.value.root, &indiv).value());
                    var.enforce_user_const(); // enforce user constant choice
                }
                // simulation-only override: a sampled latent carrying a value is pinned to it in simulation
                // the inference model is a nuisance law, unfit to generate
                // to_attribute binds the pin lazily: data, parameters and earlier latents resolve at read time
                // a constant collapses back to a numeric
                // no domain is attached, and the early return skips the model: the pin survives every draw
                else if (simulation && row.is_sim_override())
                {
                    var.transfer(translator_.to_attribute(row.value.root, &indiv));
                    return;
                }
                else
                    bind_variable(var, k, indiv, nullptr);

                // initialize model attributes
                bind_distr_(var, k, indiv);

                // assign a proposal if the variable is sampled: its arguments fold into the proposal call
                // simulation draws from the model and never proposes: data-bound args need not resolve
                // per-member context lets proposal args fold data cells
                if (!simulation && !var.is_constant())
                    var.init_proposal(row.proposal.root.func, translator_.to_numeric_args(row.proposal.root, &indiv));

                // check discreteness and positivity
                keys_.assert_requests(keys_.consult(row.name), var);
            }
            catch (const std::exception &e)
            {
                msg::stop(std::string(e.what()) +
                          "\nCheck " + var.str() + " for " + data_.wide().str(indiv) +
                          " (status " + std::to_string(indiv.status()) + ").");
            }
        }

        void Binder::bind_distr_(Variable &var,
                                 size_t k,
                                 const Individual &indiv) const
        {
            const DefRow &row = def_[k];
            // a flat row is uniform within its bounds, its density left out of the latent log-likelihood
            var.init_distr(row.is_flat() ? std::string_view("uniform") : std::string_view(row.model.root.func),
                           row.is_discr);

            if (var.has_distr())
                translator_.bind_distr(var, row.model.root, row.min.root, row.max.root, row.is_discr, &indiv, false);
        }

        void Binder::bind_surrogate(State &state,
                                    const Surrogate &surr,
                                    const Household &house,
                                    bool for_evaluation) const
        {
            const std::string &name = def_[surr.row].name;
            const size_t d = data_.wide().consult(name);
            const auto build = [&](Variable &var,
                                   const Individual &anchor,
                                   const Household *hptr)
            {
                try
                {
                    bind_variable(var, surr.row, anchor, hptr);
                }
                catch (const std::exception &e)
                {
                    msg::stop(std::string(e.what()) + "\nCheck " + var.str() + " for " +
                              (hptr ? data_.wide().str(*hptr) : data_.wide().str(anchor)) + ".");
                }
            };

            // household level
            if (surr.is_house)
            {
                Variable &hvar = state.house_varia[Index::of(house)][surr.slot];
                build(hvar, house.at(0), &house);
                if (for_evaluation)
                {
                    assert_datum_requests_(d, hvar, data_.wide().str(house));
                    hvar.finalize(); // enforce inference-finalization
                }
                for (const Individual &indiv : house)
                    data_.state.wire_readcell(d, indiv, &hvar);
                return;
            }

            // member level
            const size_t a = keys_.consult(name);
            for (const Individual &indiv : house)
            {
                Variable &avar = state.of(indiv).vec[a];
                build(avar, indiv, nullptr);
                if (for_evaluation)
                {
                    assert_datum_requests_(d, avar, data_.wide().str(indiv));
                    avar.finalize(); // enforce inference-finalization
                }
                data_.state.wire_readcell(d, indiv, &avar);
            }
        }

        void Binder::assert_datum_requests_(size_t d,
                                            const Variable &var,
                                            const std::string &context) const
        {
            // a column with nothing to derive from carries no verdict
            // its absence is reported at the run wall, not as a violated requirement
            if (na::is(var.value()))
                return;
            try
            {
                data_.wide().assert_requests(d, var);
            }
            catch (const std::exception &e)
            {
                msg::stop(std::string(e.what()) + "\nCheck " + context + ".");
            }
        }
    }

    //////////////////////////////  Augmented  //////////////////////////////

    namespace
    {
        /// Infected status of a member whose status observation draws category `k`.
        StatusLvl infected_level(double k) noexcept
        {
            return static_cast<StatusLvl>(status::INFECTED + static_cast<StatusLvl>(k));
        }
    }

    Augmented::Augmented(MCMCBase *mcmc, std::vector<std::string> keys)
        : KeyHolderWithRequests(std::move(keys)), mcmc_(mcmc) {}
    Augmented::~Augmented() = default;

    void Augmented::request_transmission(std::vector<size_t> vec)
    {
        vect::append_unique(user_transm_vec_, vec);
    }

    void Augmented::build(InfecProfile &infec_profile,
                          const str_table &input)
    {
        build_restore_internal_(); // reset internal state
        if (!input.empty())        // DefTable::append indexes its columns
            def_.append(input);
        build_append_profile_(infec_profile); // add infectivity profile

        // skip if no latent definitions
        if (def_.empty())
            return;

        // register member definition names, household observations excluded
        for (size_t k = 0; k < def_.size(); ++k)
            if (def_[k].kind != augm::Kind::House)
                add(def_[k].name);

        augm::SchemeBuilder builder(scheme_, def_, *this);
        // status routing
        builder.register_rows();
        // per-member storage
        state.allocate(mcmc_->data.wide().report_indiv.total, *this, def_);

        // infected statuses the data carries
        std::unordered_set<StatusLvl> infected;
        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house)
                if (indiv.is_infected_status())
                    infected.emplace(indiv.status());

        // definition order per status, and equivalence-class covariates
        mcmc_->data.wide().set_eq_class_data_names(
            builder.order({infected.begin(), infected.end()}, mcmc_->param));

        // surrogate columns
        builder.init_surrogates();
        // observation nodes
        builder.init_observations();
        // keys reaching the transmission side
        builder.init_transm(user_transm_vec_.empty() ? std::vector<std::string>{}
                                                     : transm_names_());

        // assert status observation
        if (scheme_.is_any_status_obs())
            augm::check::status_obs_consistent(def_, scheme_.status_obs_row());
        // no definition reads an obs entry column, no simulation override reads an observation
        augm::check::no_obs_entry_collision(def_, consult(scheme_.list_obs_member));
        augm::check::no_sim_override_obs_read(def_, consult(scheme_.list_obs_member));
        // a flat latent is generated by transmission, a generated latent has no model
        augm::check::flat_rows(def_, mcmc_->param, consult(user_transm_vec_));
    }

    void Augmented::build_restore_internal_()
    {
        def_.clear();
        scheme_ = augm::Scheme();
        clear_optional();
        state.clear();
    }

    void Augmented::build_append_profile_(InfecProfile &infec_profile)
    {
        if (!infec_profile.is_latent())
            return;

        def_.append(infec_profile.latent_definitions);
        add(infprof::LATENT_FUNC);
        infec_profile.latent_key = consult(infprof::LATENT_FUNC);
        infec_profile.bind_augm(this); // for runtime per-member profile resolution
    }

    //////  Augmented : export columns  //////

    bool Augmented::is_allocated() const noexcept
    {
        return state.varia.size() == mcmc_->data.wide().report_indiv.total;
    }

    size_t Augmented::n_export_cols() const noexcept
    {
        size_t n = scheme_.list_latent_export.size() + scheme_.list_obs_member.size();
        // add number of surrogate columns to be exported
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            n += surr.is_own_column();
        return n;
    }

    bool Augmented::is_col_discrete_(size_t i, bool is_house) const
    {
        // a household column carries no status partition, hence no declaration to check
        if (!is_house && def_.has_continuous_declaration(consult(i)))
            return false;

        bool seen = false;
        const size_t n = is_house ? state.house_varia.size() : state.varia.size();
        for (size_t j = 0; j < n; ++j)
        {
            const Variable &var = is_house ? state.house_varia[j][i] : state.varia[j].vec[i];
            if (var.is_na())
                continue;
            if (!var.is_discrete())
                return false;
            seen = true;
        }
        return seen;
    }

    std::vector<bool> Augmented::discrete_export_cols() const
    {
        std::vector<bool> out;
        out.reserve(n_export_cols());

        for (size_t a : vect::concatenate(scheme_.list_latent_export, scheme_.list_obs_member))
            out.push_back(is_col_discrete_(a, false));

        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            if (surr.is_own_column())
                out.push_back(surr.is_house ? is_col_discrete_(surr.slot, true)
                                            : is_col_discrete_(consult(def_[surr.row].name), false));

        return out;
    }

    //////  Augmented : initialize latent  //////

    void Augmented::init_for_evaluation()
    {
        if (empty())
            return;

        assert_every_member_cols_();
        materialize_surrogate_columns(false);
        augm::SchemeBuilder builder(scheme_, def_, *this);
        builder.init_models();
        const augm::Binder binder = binder_();

        // define latent variables for all latent individuals in order
        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
                for (size_t k : helper_materialize_definition_order_(indiv, false))
                    binder.bind_member(state.of(indiv).vec[consult(def_[k].name)], k, indiv, false);

        // observations the data provide a value for, and what still reaches a sink
        init_live_schemes_(builder);

        // parameters that may move latent/observation models
        builder.init_param_producers(mcmc_->param);
    }

    void Augmented::init_live_schemes_(augm::SchemeBuilder &builder)
    {
        // a previous run left its variants behind
        // the definition-order pass points every individual back at the base scheme of its status
        builder.clear_variants();
        const std::vector<std::string> transm = transm_names_();

        // per observation key: members declaring it, and those without a datum
        std::vector<size_t> declared(size(), 0), missing(size(), 0);
        // base scheme and observations left out -> the variant standing for them
        std::vector<std::pair<std::vector<size_t>, size_t>> variants;

        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
            {
                std::vector<size_t> absent;
                for (size_t a : scheme_of_(indiv).obs_models)
                {
                    ++declared[a];
                    if (na::is(state.of(indiv).vec[a].value()))
                        ++missing[a], absent.push_back(a);
                }
                const std::vector<size_t> key =
                    vect::concatenate(std::vector<size_t>({state.of(indiv).idx}), absent);
                const auto it = std::find_if(variants.begin(), variants.end(),
                                             [&key](const auto &v)
                                             { return v.first == key; });
                state.of(indiv).idx =
                    it != variants.end()
                        ? it->second
                        : variants.emplace_back(key, builder.make_variant(key.front(), absent, transm)).second;
            }

        augm::check::report_missing_obs(*this, declared, missing);
    }

    std::vector<std::string> Augmented::transm_names_() const
    {
        // transmission reads model-declared keys, and the latent infectivity profile if any
        std::vector<std::string> transm = consult(user_transm_vec_);
        transm.push_back(infprof::LATENT_FUNC);
        return transm;
    }

    void Augmented::assert_every_member_cols_() const
    {
        // a member surrogate drawn with the members it reads has a value per infection
        std::vector<size_t> per_infection;
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            if (!surr.is_own_column())
                per_infection.push_back(surr.row);
        augm::check::every_member_cols(def_, mcmc_->data.wide(), per_infection);
    }

    //////  Augmented : surrogate columns  //////

    std::vector<std::string> Augmented::get_produced_cols()
    {
        std::vector<std::string> produced_cols = consult(scheme_.list_obs_member);
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            produced_cols.push_back(def_[surr.row].name);
        return produced_cols;
    }

    std::vector<std::string> Augmented::get_live_cols() const
    {
        std::vector<std::string> live_cols = consult(scheme_.list_obs_member);
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            if (!surr.is_static)
                live_cols.push_back(def_[surr.row].name);
        return live_cols;
    }

    void Augmented::init_house_for_simulation()
    {
        data::Wide &data_wide = mcmc_->data.wide();

        // generated household observations
        std::vector<std::string> live_house;
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            if (surr.is_house && !surr.is_static)
                live_house.push_back(def_[surr.row].name);
        augm::check::no_live_house_obs_cycle(def_, live_house);
        assert_every_member_cols_();
        mcmc_->data.reset_readcells(get_produced_cols());

        // reserve household-level observations
        state.house_varia.clear();
        if (!scheme_.list_obs_house.empty())
            state.house_varia.assign(data_wide.report_house.total,
                                     std::vector<Variable>(scheme_.list_obs_house.size()));

        const augm::Binder binder = binder_();
        for (const Household &house : data_wide)
        {
            // wire stochastic member observations to their data read-cells
            for (const Individual &indiv : house)
                for (size_t a : scheme_.list_obs_member)
                    mcmc_->data.state.wire_readcell(data_wide.consult(consult(a)), indiv,
                                                    &state.of(indiv).vec[a]);

            // build and wire surrogate columns in dependency order
            for (const augm::Surrogate &surr : scheme_.list_surrogate)
                // defer live member observation (drawn with the member)
                if (surr.is_house || surr.is_static)
                    binder.bind_surrogate(state, surr, house, false);
        }
    }

    void Augmented::init_members_for_simulation()
    {
        if (empty() || !is_allocated())
            return;

        // member variables bound for the whole run: status observation, static surrogates
        sim_fixed_.assign(size(), false);
        if (scheme_.is_any_status_obs())
            sim_fixed_[scheme_.status_obs_key()] = true;
        for (const augm::Surrogate &surr : scheme_.list_surrogate)
            if (surr.is_static && !surr.is_house)
                sim_fixed_[consult(def_[surr.row].name)] = true;

        // bind one block per reachable status scheme, swapped in by the household loop
        state.sim.assign(state.varia.size(), {});
        const augm::Binder binder = binder_();
        msg::WarningBuffer warnings;
        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house)
                bind_simulation_blocks_(binder, indiv);
        warnings.flush();
    }

    //////  Augmented : data-bound surrogate columns  //////

    void Augmented::materialize_surrogate_columns(bool for_export)
    {
        if (scheme_.list_surrogate.empty())
            return;

        data::Wide &data_wide = mcmc_->data.wide();
        const augm::Binder binder = binder_();

        // build the surrogate Variables, shared with the simulation init
        // point each standing-in data column's read-surface at its Variable
        // data readers then resolve the surrogate value without a data copy
        state.house_varia.assign(data_wide.report_house.total,
                                 std::vector<Variable>(scheme_.list_obs_house.size()));

        for (const augm::Surrogate &surr : scheme_.list_surrogate)
        {
            const std::string &astr = def_[surr.row].name;
            // a surrogate column must resolve from data
            // a latent, parameter or group-risk source has no fixed value in inference
            // it would impose an indicator constraint
            if (!def_.is_pure_data(surr.row))
            {
                // exporting the data view is not inferring: the column keeps its datum, having no derived value
                // the run wall still refuses it
                if (for_export)
                    continue;
                msg::stop(
                    "Cannot infer from " + str::format(def_[surr.row].kind) + " " +
                    msg::val(astr) + ": it must resolve from data alone.\n" +
                    "A definition reading a latent variable or a parameter is only meaningful in simulation.");
            }

            for (Household &house : data_wide)
                binder.bind_surrogate(state, surr, house, true);
        }
    }

    //////  Augmented : init missing attribute  //////

    void Augmented::clear_all()
    {
        if (empty() || !is_allocated())
            return;

        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
            {
                augm::IndivState &st = state.of(indiv);
                augm::SimMember &sm = state.sim[Index::of(indiv)];
                for (size_t a : sm.swapped)
                    st.vec[a].swap(sm.blocks[st.idx][a]);
                sm.swapped.clear();

                for (size_t a = 0; a < st.vec.size(); ++a)
                    if (!sim_fixed_[a])
                        st.vec[a].clear();
            }
    }

    void Augmented::release_simulation()
    {
        for (augm::IndivState &st : state.varia)
            for (Variable &avar : st.vec)
                avar.clear();
        state.sim.clear();
    }

    void Augmented::init_definitions_for_simulation(const Individual &indiv)
    {
        const std::vector<size_t> &order = helper_materialize_definition_order_(indiv, true);
        augm::IndivState &st = state.of(indiv);
        augm::SimMember &sm = state.sim[Index::of(indiv)];
        if (sm.errors[st.idx])
            std::rethrow_exception(sm.errors[st.idx]);
        std::vector<Variable> &block = sm.blocks[st.idx];

        for (size_t k : order)
        {
            const size_t a = consult(def_[k].name);
            Variable &avar = st.vec[a];
            // set by transmission, or bound for the whole run
            if (avar.has_definition())
                continue;

            // only transmission sets a flat row, which has no block entry
            if (!block[a].has_definition())
                msg::stop(
                    "No " + msg::arg("value") + " to simulate latent " + msg::val(def_[k].name) +
                    ": it has no model, and transmission did not set it.\nProvide " +
                    msg::arg("value") + " to set it in simulation." +
                    "\nCheck " + avar.str() + " for " + mcmc_->data.wide().str(indiv) +
                    " (status " + std::to_string(indiv.status()) + ").");

            avar.swap(block[a]);
            sm.swapped.push_back(a);
        }
    }

    void Augmented::bind_simulation_blocks_(const augm::Binder &binder, const Individual &indiv)
    {
        augm::IndivState &st = state.of(indiv);
        augm::SimMember &sm = state.sim[Index::of(indiv)];

        std::vector<size_t> reached;
        const auto reach = [&](StatusLvl s)
        {
            const size_t i = scheme_.index_of(s);
            // an unmapped status reports its gap when drawn
            if (!na::is(i) && !vect::contains(reached, i))
                reached.push_back(i);
        };

        // reachable statuses:
        // - the observed one, kept by a model that does not redraw it;
        // - the infected level, shifted by the status observation draw.
        if (indiv.is_infected_status())
            reach(indiv.status());
        double k_min = 0., k_max = 0.;
        if (scheme_.is_any_status_obs())
        {
            Variable &svar = st.vec[scheme_.status_obs_key()];
            try
            {
                binder.bind_member(svar, scheme_.status_obs_row(), indiv, true);
                svar.release_user_const(); // readers follow the draw, not the placeholder
                const Codomain c = svar.codomain();
                k_min = c.inf;
                k_max = c.sup;
            }
            // raised if the member is ever drawn infected
            catch (const std::exception &)
            {
                sm.status_error = std::current_exception();
                svar.clear();
                k_max = -1.;
            }
        }
        for (double k = k_min; k <= k_max; ++k)
            reach(infected_level(k));

        for (size_t i : reached)
        {
            if (sm.blocks.size() <= i)
            {
                sm.blocks.resize(i + 1);
                sm.errors.resize(i + 1);
            }
            st.idx = i;
            std::vector<Variable> &block = sm.blocks[i];
            block = state.blank;

            // bind into the member's variables, the context the definitions read
            // a latent transmission generates binds into the block
            // its readers bind to a placeholder in the member's variable
            // they read what transmission or the block sets there
            std::vector<size_t> bound, generated;
            try
            {
                for (size_t k : scheme_.at(i).sim_order)
                {
                    const augm::DefRow &row = def_[k];
                    const size_t a = consult(row.name);
                    if (st.vec[a].has_definition()) // bound for the whole run
                        continue;
                    if (vect::contains(user_transm_vec_, a))
                    {
                        generated.push_back(a);
                        st.vec[a].set(0.);
                        // flat: only transmission sets it
                        if (!row.is_flat() || row.is_sim_override())
                            binder.bind_member(block[a], k, indiv, true);
                        continue;
                    }
                    bound.push_back(a);
                    binder.bind_member(st.vec[a], k, indiv, true);
                }
            }
            // raised if the member reaches this status scheme
            catch (const std::exception &)
            {
                sm.errors[i] = std::current_exception();
            }

            // move bound variables to the block, leaving the member's variables empty
            for (size_t a : generated)
                st.vec[a].clear();
            for (size_t a : bound)
                block[a].swap(st.vec[a]);
        }
    }

    const std::vector<size_t> &Augmented::helper_materialize_definition_order_(const Individual &indiv,
                                                                               bool simulation)
    {
        StatusLvl s = indiv.status();

        // its own status scheme, otherwise the generic one if it exists
        const size_t i = scheme_.index_of(s);
        if (!na::is(i))
        {
            state.of(indiv).idx = i;
            return simulation ? scheme_.at(i).sim_order : scheme_.at(i).order;
        }

        if (indiv.is_infected_status())
            msg::stop(
                "No latent definitions were found for " +
                mcmc_->data.wide().str(indiv) + " with status " + std::to_string(s) + ".");

        msg::stop(
            "Generic latent definition for non-infectious " +
            mcmc_->data.wide().str(indiv) + " was missing.");
    }

    //////  Augmented : user reset  //////

    void Augmented::simulate_infected_status(Household &house,
                                             Individual &indiv,
                                             rnd_gen &prng)
    {
        double k = 0.;
        if (scheme_.is_any_status_obs())
        {
            if (const std::exception_ptr &err = state.sim[Index::of(indiv)].status_error)
                std::rethrow_exception(err);
            Variable &svar = state.of(indiv).vec[scheme_.status_obs_key()];
            svar.draw_value(prng);
            k = svar.value();
        }

        house.set_new_status_(indiv, infected_level(k));
    }

    void Augmented::restart(rnd_gen &prng)
    {
        if (!scheme_.is_any_sample())
            return;

        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
                restart(indiv, prng);
    }

    void Augmented::draw_sample(rnd_gen &prng)
    {
        if (!scheme_.is_any_sample())
            return;

        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
                augm::draw_sample(scheme_of_(indiv).sample, state.of(indiv), prng);
    }

    void Augmented::finalize()
    {
        if (empty())
            return;

        for (const Household &house : mcmc_->data.wide())
            for (const Individual &indiv : house.infected_indivs())
                for (Variable &avar : state.of(indiv).vec)
                    avar.finalize();
    }

    [[noreturn]] void Augmented::throw_obs_init_error(const Individual &indiv) const
    {
        // observations of the member whose datum lies outside the domain
        std::vector<size_t> rows;
        for (size_t k : scheme_of_(indiv).order)
        {
            if (def_[k].kind != augm::Kind::Member)
                continue;
            const size_t a = consult(def_[k].name);
            if (vect::contains(list_obs_models(indiv), a) &&
                !std::isfinite(state.of(indiv).vec[a].guarded_log_pdf()))
                rows.push_back(k);
        }
        augm::check::obs_init_error(def_, rows, mcmc_->param, mcmc_->data.wide().str(indiv));
    }

    augm::Binder Augmented::binder_()
    {
        return augm::Binder(def_, *this, mcmc_->translator, mcmc_->data);
    }
}
