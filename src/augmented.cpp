/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <string>
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
#include "data/dataset.h"
#include "data/household.h"
#include "data/individual.h"
#include "distribution.h"
#include "function.h"
#include "infec_profile.h"
#include "parameter.h"
#include "registry/functions.h"
#include "rules.h"
#include "status_process.h"
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

    InfecProfile::InfecProfile(const Augmented &augm) : augm_(augm) {}
    InfecProfile::~InfecProfile() = default;

    void InfecProfile::build(Translator &translator,
                             const str_table &input)
    {
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
        const bool is_discr = is_discr_(input, 0);
        generic = helper_build_(translator, input.at("func")[0], is_discr);
        lookup_.clear();

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

        // check profile name and discretization
        for (size_t r = 0; r < count_rows(input); ++r)
        {
            const bool discr = is_discr_(input, r);
            auto probe = helper_build_(translator, parse::Definition(input.at("func")[r]).root.func, discr);
            if (discr)
                probe->discretize();
        }
    }

    str_table InfecProfile::list_latent_definitions(const str_table &input) const
    {
        str_table output;
        if (!is_latent())
            return output;

        for (size_t r = 0; r < count_rows(input); ++r)
        {
            output["name"].push_back(infprof::LATENT_FUNC);
            output["status"].push_back(input.at("status")[r]);
            output["func"].push_back(input.at("func")[r]);
            output["fcumul"].push_back("false");
            output["const"].push_back("");
            output["discr"].push_back(str::format(is_discr_(input, r)));
            output["model"].push_back("");
            output["min"].push_back(input.at("min")[r]);
            output["max"].push_back(input.at("max")[r]);
            output["proposal"].push_back("");
            output["kind"].push_back("latent");
        }
        return output;
    }

    bool InfecProfile::is_discr_(const str_table &input, size_t r) const
    {
        return is_required_discrete_ || str::tobool(input.at("discr")[r]);
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Translator &translator,
                                                              const std::string &func,
                                                              bool discr)
    {
        parse::Definition parser(func);

        std::unique_ptr<Distribution> output = registry::create_distribution(
            parser.root.func, ErrType::InfecProfile, registry::infec_profile_search_order(discr));

        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, translator.to_attribute(node));
            translator.mark_param_produces_transm_cache(node);
        }

        return output;
    }

    void InfecProfile::assert_zero_at_origin(const Function &f) const
    {
        const bool is_discr = is_required_discrete_;
        // skip for hand-written infectivity profiles
        if (!is_discr && registry::instance(registry::Kind::InfecProfile).contains(f.name()))
            return;

        const Domain &domain = f.domain();
        // relocated domain, not support(): the dense rule of support() would relax a floor at zero
        const Codomain reach = domain.relocated();
        // the bounds keep every attainable draw past the origin
        if (reach.is_positive())
            return;

        std::string sms = "Infectivity profile " + msg::val(f.name());

        if (na::is(domain.delocate(0.)))
        {
            // the origin reads an undrawn value the bounds do not keep away from it
            sms += " is not kept away from generation " + std::string(is_discr ? "day" : "time") + " 0.\n" +
                   msg::alert("Its support is " + reach.str() + ", so it can put infectivity on the " +
                              (is_discr ? "day" : "instant") + " of infection.");
        }
        else
        {
            const double y = f.pdf(0.);
            if (std::abs(y) < 1.e-10 || na::is(y))
                return;

            sms += " is non-zero at generation ";
            if (is_discr)
                sms += "day 0 (F[0] = " + str::format(y) + ").\n" +
                       "A discrete-time profile admits no transmission on the day of infection.";
            else
                sms += "time 0 (f(0) = " + str::format(y) + ").\n" +
                       "A continuous-time profile admits no transmission at the instant of infection.";
        }

        msg::stop(
            sms +
            "\nConsider truncating the distribution with " + msg::arg("min") +
            ", or shifting with " + msg::arg("loc") + ".");
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
            // if (kstr == "latent")
            return augm::Kind::Latent;
        }

        DefRow::DefRow(const str_table &input, size_t r)
            : name(input.at("name")[r]),
              kind(string_to_kind(input.at("kind")[r])),
              value(input.at("const")[r]),
              func(input.at("func")[r]),
              model(input.at("model")[r]),
              min(input.at("min")[r]),
              max(input.at("max")[r]),
              proposal(input.at("proposal")[r]),
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
            if (!is_bare())
            {
                const auto all = source_cols();
                return {all.begin(), all.end()};
            }
            if (simulation)
                return {};
            return {&min, &max};
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
                if (row.name == name && (row.has_model() || row.is_sampled()) && !row.is_discr)
                    return true;
            return false;
        }

        /// Whether a parsed definition subtree references any group rate.
        static bool entry_has_group_rate(const parse::Entry &node)
        {
            if (parse::is_group_rate(node.type))
                return true;
            for (const parse::Entry &child : node.children)
                if (entry_has_group_rate(child))
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
                    entry_has_group_rate(col->root))
                    return false;
            }
            return true;
        }

        std::vector<std::string> DefTable::list_eq_class_keys() const
        {
            std::vector<std::string> output;
            for (const DefRow &row : rows)
            {
                // a modelled observation reads its datum in inference
                if (row.kind == Kind::Member && row.has_model())
                    output.push_back(row.name);
                // a produced column is told apart by the sources of its own row
                for (const parse::Definition *col : row.source_cols())
                    if (!col->is_na())
                        for (const std::string &dstr : col->list_dataset_producers())
                            if (!has_row(dstr, Kind::Member) && !has_row(dstr, Kind::House))
                                output.push_back(dstr);
            }
            vect::sort_unique(output);
            return output;
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
                    for (const std::string &dstr : col->list_dataset_producers())
                        if (dstr != def[k].name && vect::contains(targets, dstr))
                            producers.insert(dstr);
            }
        }

        str_table model_consumers(const DefTable &def,
                                  const std::vector<size_t> &order)
        {
            const dag::prod_graph pgraph = producers(def, order, false);

            // a bare row is keyed as its own consumer: its uniform log-ratio is zero, its log-pdf constant
            std::vector<std::string> without_model;
            for (size_t k : order)
                if ((!def[k].has_model() && !def[k].is_bare()) || def[k].kind == Kind::Member)
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
        /// * transmission names;
        /// * observations of `order` carrying a datum;
        /// * latents a household observation reduces over.
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
                // a root needs no definition in every block
                if (pgraph.count(node))
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

        /// List the data names (source D) read by `def` that are in `names`.
        static std::vector<std::string> referred_data(const parse::Definition &def,
                                                      const std::vector<std::string> &names)
        {
            std::vector<std::string> output = {};
            if (def.is_na())
                return output;
            for (const std::string &dstr : def.list_dataset_producers())
                if (vect::contains(names, dstr))
                    output.push_back(dstr);
            return output;
        }

        void every_member_cols(const DefTable &def,
                               const data::Dataset &dataset,
                               const std::vector<size_t> &per_infection)
        {
            for (size_t k = 0; k < def.size(); ++k)
            {
                const std::string &name = def[k].name;
                if (def[k].kind != Kind::Member ||
                    !dataset.count(name) ||
                    !dataset.is_required_everyone(dataset.consult(name)))
                    continue;

                // drawn, or computed from a draw, with the infection of each member
                if (def[k].has_model() || vect::contains(per_infection, k))
                    msg::stop(
                        "Data " + msg::val(name) + " needs a value for every member, escaped ones included.\n" +
                        msg::alert(str::format(def[k].kind) + " " + msg::val(name) +
                                   " has a value only for infected members.") +
                        "\nDefine it at household level with " + msg::code("level = \"house\"") +
                        ", or from data columns only.");
            }
        }

        void no_drawn_house_obs_cycle(const DefTable &def,
                                      const std::vector<std::string> &drawn_names)
        {
            if (drawn_names.empty())
                return;

            for (const DefRow &row : def.rows)
            {
                if (row.kind != Kind::Latent)
                    continue;

                // the bounds of a bare row are never evaluated in simulation, nor flagged
                for (const parse::Definition *col : row.source_cols(true))
                    for (const std::string &dstr : referred_data(*col, drawn_names))
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
                // * model support: model, min, max;
                // * deterministic value: value, func.
                // keep them in sync with source_cols()

                // rule 1: a latent's model (distribution and support) reads no observation
                // its producers are parameters and other latents only
                // an observation's model may read another observation, drawn before it
                if (def[k].kind == Kind::Latent && def[k].has_model())
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
                        // observation node reading its own entry column
                        if (dstr == name)
                            continue;

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

        void bare_rows(const DefTable &def,
                       const Parameter &param,
                       const std::vector<std::string> &transm)
        {
            for (const DefRow &row : def.rows)
            {
                if (row.kind != Kind::Latent)
                    continue;
                const bool is_transm = vect::contains(transm, row.name);

                if (row.has_model() && row.is_sampled() && is_transm)
                    msg::stop(
                        "Latent " + msg::val(row.name) +
                        " is generated by the transmission model, which already gives its probability.\n" +
                        "Drop " + msg::arg("model") + ", and keep " + msg::arg("min") + " and " +
                        msg::arg("max") + " to sample it within.");

                if (!row.is_bare())
                    continue;

                if (!is_transm)
                    msg::stop(
                        "Latent " + msg::val(row.name) +
                        " has no model, and the transmission model does not generate it.\n" +
                        "Provide a " + msg::arg("model") + ".");

                // the uniform density 1 / (max - min) is left out: the bounds must not move with the draws
                for (const auto &[key, col] : {std::pair{"min", &row.min}, std::pair{"max", &row.max}})
                {
                    std::string read;
                    for (const std::string &pstr : col->list_parameter_producers())
                        if (!param.is_constant(pstr))
                            read = "parameter " + msg::val(pstr);
                    for (const std::string &astr : col->list_latent_producers())
                        read = "latent " + msg::val(astr);
                    if (read.empty() && entry_has_group_rate(col->root))
                        read = "a group rate";
                    if (!read.empty())
                        msg::stop(
                            "Bound " + msg::arg(key) + " of latent " + msg::val(row.name) + " reads " + read +
                            ".\nWithout a model, bounds can read only data and observations.");
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
                        "Match the latent model of " + str::format(str::unique(latent_producers), false, 3) +
                        " to the values observation " + msg::val(def[k].name) +
                        " accepts, revise its domain, or check the observed individual data.");

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

                // store status -> definition rows
                scheme_.status_rows_[def_[k].status].push_back(k);
            }
        }

        //////  augm::SchemeBuilder : order  //////

        void SchemeBuilder::order(const std::vector<StatusLvl> &infected,
                                  const Parameter &param,
                                  bool draws_status)
        {
            const auto &rows = scheme_.status_rows_;
            auto &index = scheme_.status_index_;
            Sources_ src{param};

            // skip if no member definition goes into a status scheme
            // household observations need no status scheme
            if (rows.empty())
                return;

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
                    !draws_status ||
                    index.empty())
                    order_status_(status::MISSING, status::MISSING, src);

            scheme_.n_base_ = scheme_.status_schemes_.size();
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
            dict<size_t> name_to_row;
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

                // simulation draws a member observation before its readers: its producers differ from inference
                at_(s).sim_order = sorted_rows(true);
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

            if (def_[k].is_const())
            {
                is_to_draw = false;
            }
            // deterministic definition: a function
            else if (def_[k].has_func())
            {
                // do not export if not callable: it would crash
                is_to_export = def_[k].func.is_callable_implicitly("x");
                is_to_draw = list_function_(s, k, src);
            }
            else
            {
                scheme_.is_any_sample_ = true;
                at_(s).sample.emplace_back(a);

                if (is_declared_adaptive(def_[k].proposal))
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

        bool SchemeBuilder::list_function_(StatusLvl s, size_t k, Sources_ &src)
        {
            const parse::Definition &parser = def_[k].func;

            bool all_consts = true;
            for (const std::string &pstr : parser.list_parameter_producers())
                all_consts &= src.param.is_constant(pstr);

            // inspect status-specific latent definitions: draw unless constant
            // relies on the definition order
            for (const std::string &bstr : parser.list_latent_producers())
                if (vect::contains(at_(s).draws, keys_.consult(bstr)))
                    all_consts = false;

            // to draw if not all constant dependencies
            return !all_consts;
        }

        //////  augm::SchemeBuilder : deterministic observations and nodes  //////

        void SchemeBuilder::init_deterministic_obs()
        {
            // look for deterministic observations, each producing a data column
            dict<size_t> name_to_row;
            for (size_t k = 0; k < def_.size(); ++k)
                // household reduction or member-level definition carrying no model
                if (def_[k].kind == Kind::House ||
                    (def_[k].kind == Kind::Member && !def_[k].has_model()))
                    name_to_row[def_[k].name] = k;

            if (name_to_row.empty())
                return;

            // order deterministic observations, producer before consumer
            std::vector<size_t> rows;
            std::vector<std::string> names;
            for (const auto &[nm, k] : name_to_row)
                rows.push_back(k), names.push_back(nm);
            dag::prod_graph pgraph;
            graph::add_data_edges(pgraph, def_, rows, names);

            // member observations carrying a model are drawn
            std::unordered_set<std::string> drawn_names;
            for (const DefRow &row : def_.rows)
                if (row.kind == Kind::Member && row.has_model())
                    drawn_names.insert(row.name);

            // list deterministic observations in order
            for (const std::string &nm : dag::order_producers_first(pgraph, ErrType::Observation))
            {
                DeterministicObs dobs;
                dobs.row = name_to_row.at(nm);
                dobs.is_house = def_[dobs.row].kind == Kind::House;
                dobs.is_static = def_.is_pure_data(dobs.row);

                // update household observations
                if (dobs.is_house)
                {
                    dobs.slot = scheme_.list_obs_house.size();
                    scheme_.list_obs_house.push_back(dobs.row);
                }

                // drawn if it reads a drawn column
                if (dobs.is_static)
                    for (const parse::Definition *col : def_[dobs.row].source_cols())
                        for (const std::string &dstr : col->list_dataset_producers())
                            if (drawn_names.count(dstr))
                                dobs.is_static = false;
                if (!dobs.is_static)
                    drawn_names.insert(nm);

                scheme_.list_deterministic_obs.push_back(std::move(dobs));
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

                    // a static deterministic observation is data, not drawn
                    bool skip = false;
                    for (const DeterministicObs &dobs : scheme_.list_deterministic_obs)
                        if (dobs.row == k)
                            skip = dobs.is_static;
                    if (skip)
                        continue;

                    const size_t obs_key = keys_.consult(def_[k].name);

                    // unique observation keys to return as simulated columns
                    if (!vect::contains(scheme_.list_obs_member, obs_key))
                        scheme_.list_obs_member.push_back(obs_key);

                    // a deterministic observation node carries no distribution
                    if (!def_[k].has_model())
                        continue;

                    scheme_.is_any_obs_model_ = true;
                    ss.obs_models.push_back(obs_key);

                    // list the observation under each of its producers, so moving any of them recomputes it
                    for (const std::string &astr : dag::transitive_producers(def_[k].name, pgraph, ErrType::Latent))
                        // transitive_producers returns 'start' last
                        if (astr != def_[k].name)
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

                // a bare row adds no density to the latent log-likelihood
                std::vector<std::string> with_model = extract_keys(model_consumers), bare;
                for (size_t k : ss.order)
                    if (def_[k].is_bare())
                        bare.push_back(def_[k].name);
                vect::remove(with_model, bare);
                ss.latent_models = keys_.consult(with_model);
                scheme_.is_any_latent_model_ |= !ss.latent_models.empty();
            }
        }

        size_t SchemeBuilder::make_variant(size_t base,
                                           const std::vector<size_t> &absent,
                                           const std::vector<std::string> &transm)
        {
            // copy: the store reallocates below
            StatusScheme ss = scheme_.status_schemes_[base];
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
                    row.kind == Kind::House)
                    continue;

                // Kind::Member or Kind::Latent
                const bool is_obs = row.kind == Kind::Member;
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

        void SchemeBuilder::init_latent_block(size_t i, std::vector<size_t> keys)
        {
            StatusScheme &ss = scheme_.status_schemes_[i];
            if (keys.size() < 2)
                keys.clear();
            ss.sample_alone = vect::difference(ss.sample, keys);

            // member densities recomputed whole, transmission if a key produces it
            ss.block_terms = {!ss.obs_models.empty(), false, !ss.latent_models.empty()};
            for (size_t a : keys)
                ss.block_terms.transm = ss.block_terms.transm || ss.produces_transm[a];
            ss.block = std::move(keys);
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
    }

    //////////////////////////////  StatusProcess  //////////////////////////////

    void StatusProcess::build(const data::Dataset &dataset,
                              const str_table &input,
                              size_t n_levels)
    {
        // every column but weight and rule group is a covariate
        std::vector<std::string> covariates = extract_keys(input);
        vect::remove(covariates, {"weight", "rule_group"});
        const rules::Member member(input, dataset, covariates);
        std::vector<double> row_weight(member.rules.n_rows);
        for (size_t r = 0; r < member.rules.n_rows; ++r)
            row_weight[r] = str::tod(input.at("weight")[r]);

        // one block per covariate class
        const rules::CovariateClasses classes(dataset, member.keys);
        block_of_ = classes.of_indiv;
        blocks_.assign(classes.count, {});
        totals_.assign(classes.count, 0.);

        // weight per infected level: product of the rule groups, on the first member of each class
        for (const Household &house : dataset)
            for (const Individual &indiv : house)
            {
                const size_t b = block_of_[Index::of(indiv)];
                std::vector<double> &weights = blocks_[b];
                if (!weights.empty())
                    continue;
                weights.assign(n_levels, 1.);
                for (size_t c = 0; c < n_levels; ++c)
                    for (const auto &[gid, r] : member.rules.resolve<size_t>(
                             member.probe(dataset, indiv, static_cast<double>(status::INFECTED + c)),
                             [&](size_t row)
                             { return member.rules.row_gid[row]; }))
                        weights[c] *= row_weight[r];
                totals_[b] = std::accumulate(weights.begin(), weights.end(), 0.);
            }
    }

    StatusLvl StatusProcess::draw(const Individual &indiv, RndGen &prng) const
    {
        if (!is_defined())
            return status::INFECTED;

        // inverse cumulative weight of one uniform draw
        const size_t b = block_of_[Index::of(indiv)];
        const std::vector<double> &weights = blocks_[b];
        const double target = standard_uniform(prng) * totals_[b];
        size_t c = 0;
        double acc = weights[0];
        while (c + 1 < weights.size() && acc < target)
            acc += weights[++c];
        return static_cast<StatusLvl>(status::INFECTED + c);
    }

    std::vector<StatusLvl> StatusProcess::drawable(const Individual &indiv) const
    {
        if (!is_defined())
            return {status::INFECTED};
        std::vector<StatusLvl> out;
        const std::vector<double> &weights = blocks_[block_of_[Index::of(indiv)]];
        for (size_t c = 0; c < weights.size(); ++c)
            if (weights[c] > 0.)
                out.push_back(static_cast<StatusLvl>(status::INFECTED + c));
        return out;
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

            // a constant is driven by its value, anything else by its function
            const parse::Definition &parser = row.is_const() ? row.value : row.func;

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
                registry::func_distr_search_order(row.is_discr),
                row.is_fcumul);
            if (translator_.is_house_fold(parser.root, house))
                func->set_eval_as_fold();
            var.init_function(std::move(func));

            // a member deterministic observation materializes a data column for every member
            // where the cells it reads are empty, it has nothing to stand for
            const bool is_deterministic_obs = row.kind == Kind::Member && !row.has_model();

            bool is_defined = true;
            if (row.kind == Kind::House)
                is_defined = translator_.bind_arguments(var, parser.root, nullptr, house);
            else
                for (const parse::Entry &child : parser.root.children)
                {
                    const Attribute &attr = translator_.to_attribute(child, &anchor);
                    if (is_deterministic_obs && attr.is_constant() && attr.is_na())
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
                if (row.is_const())
                {
                    // drawn observation: non-constant placeholder, so finalize keeps its readers bound
                    if (simulation && row.kind == Kind::Member)
                    {
                        var.init_numeric(0.);
                        var.release_user_const();
                    }
                    else
                    {
                        var.init_numeric(translator_.to_numeric(row.value.root, &indiv).value());
                        var.enforce_user_const();
                    }
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
                          "\nCheck " + var.str() + " for " + data_.dataset().str(indiv) +
                          " (status " + std::to_string(indiv.status()) + ").");
            }
        }

        void Binder::bind_distr_(Variable &var,
                                 size_t k,
                                 const Individual &indiv) const
        {
            const DefRow &row = def_[k];
            // a bare row is uniform within its bounds, its density left out of the latent log-likelihood
            var.init_distr(row.is_bare() ? std::string_view("uniform") : std::string_view(row.model.root.func),
                           row.is_discr);

            if (var.has_distr())
                translator_.bind_distr(var, row.model.root, row.min.root, row.max.root, row.is_discr, &indiv, false);
        }

        void Binder::bind_deterministic_obs(State &state,
                                            const DeterministicObs &dobs,
                                            const Household &house,
                                            bool for_evaluation) const
        {
            const std::string &name = def_[dobs.row].name;
            const size_t d = data_.dataset().consult(name);
            const auto build = [&](Variable &var,
                                   const Individual &anchor,
                                   const Household *hptr)
            {
                try
                {
                    bind_variable(var, dobs.row, anchor, hptr);
                }
                catch (const std::exception &e)
                {
                    msg::stop(std::string(e.what()) + "\nCheck " + var.str() + " for " +
                              (hptr ? data_.dataset().str(*hptr) : data_.dataset().str(anchor)) + ".");
                }
            };

            // household level
            if (dobs.is_house)
            {
                Variable &hvar = state.house_varia[Index::of(house)][dobs.slot];
                build(hvar, house.at(0), &house);
                if (for_evaluation)
                {
                    assert_datum_requests_(d, hvar, data_.dataset().str(house));
                    hvar.finalize();
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
                    assert_datum_requests_(d, avar, data_.dataset().str(indiv));
                    avar.finalize();
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
                data_.dataset().assert_requests(d, var);
            }
            catch (const std::exception &e)
            {
                msg::stop(std::string(e.what()) + "\nCheck " + context + ".");
            }
        }
    }

    //////////////////////////////  Augmented  //////////////////////////////

    Augmented::Augmented(std::vector<std::string> keys,
                         const augm::DefTable &def,
                         Data &data,
                         const Parameter &param,
                         InfecProfile &infec_profile,
                         const StatusProcess &status_process,
                         const Translator &translator)
        : KeyHolderWithRequests(std::move(keys)),
          data_(data),
          param_(param),
          infec_profile_(infec_profile),
          status_process_(status_process),
          translator_(translator),
          def_(def) {}
    Augmented::~Augmented() = default;

    void Augmented::request_transmission(std::vector<size_t> vec)
    {
        vect::append_unique(user_transm_vec_, vec);
    }

    void Augmented::build()
    {
        build_restore_internal_();
        build_register_profile_();

        // data columns the definitions read tell equivalence classes apart
        data_.set_eq_class_keys(Data::Source::Latent, def_.list_eq_class_keys());
        // skip if no latent definitions
        if (def_.empty())
            return;

        // register member definition names, household observations excluded
        for (size_t k = 0; k < def_.size(); ++k)
            if (def_[k].kind != augm::Kind::House)
                add_(def_[k].name);

        augm::SchemeBuilder builder(scheme_, def_, *this);
        // status routing
        builder.register_rows();
        // per-member storage
        state.allocate(data_.dataset().report_indiv.total, *this, def_);

        // infected statuses the data carries
        std::unordered_set<StatusLvl> infected;
        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house)
                if (indiv.is_infected_status())
                    infected.emplace(indiv.status());

        // definition order per status
        builder.order({infected.begin(), infected.end()}, param_, status_process_.is_defined());

        // deterministic observations
        builder.init_deterministic_obs();
        // observation nodes
        builder.init_observations();
        // keys reaching the transmission side
        builder.init_transm(user_transm_vec_.empty() ? std::vector<std::string>{}
                                                     : transm_names_());

        // no definition reads an obs entry column
        augm::check::no_obs_entry_collision(def_, consult(scheme_.list_obs_member));
        // a bare latent is generated by transmission, a generated latent has no model
        augm::check::bare_rows(def_, param_, consult(user_transm_vec_));
    }

    void Augmented::build_restore_internal_()
    {
        scheme_ = augm::Scheme();
        clear_optional_();
        state.clear();
        pins_.clear();
    }

    void Augmented::build_register_profile_()
    {
        if (!infec_profile_.is_latent())
            return;

        add_(infprof::LATENT_FUNC);
        infec_profile_.latent_key = consult(infprof::LATENT_FUNC);
    }

    //////  Augmented : export columns  //////

    bool Augmented::is_allocated() const noexcept
    {
        return state.varia.size() == data_.dataset().report_indiv.total;
    }

    size_t Augmented::n_export_cols() const noexcept
    {
        size_t n = scheme_.list_latent_export.size() + scheme_.list_obs_member.size();
        // add number of deterministic observations to be exported
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            n += dobs.is_own_column();
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

        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            if (dobs.is_own_column())
                out.push_back(dobs.is_house ? is_col_discrete_(dobs.slot, true)
                                            : is_col_discrete_(consult(def_[dobs.row].name), false));

        return out;
    }

    //////  Augmented : initialize latent  //////

    void Augmented::init_for_evaluation()
    {
        // no definition -> no latent block, the lookup sized for every member
        if (empty())
            return state.block_factor_of.assign(data_.dataset().report_indiv.total, na::INDEX);

        assert_every_member_cols_();
        materialize_deterministic_obs(false);
        augm::SchemeBuilder builder(scheme_, def_, *this);
        builder.init_models();
        const augm::Binder binder = binder_();

        // bind latents of every infected member in order, each distinct warning raised once
        msg::WarningBuffer warnings;
        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house.infected_indivs())
                for (size_t k : helper_materialize_definition_order_(indiv, false))
                    binder.bind_member(state.of(indiv).vec[consult(def_[k].name)], k, indiv, false);
        warnings.flush();

        // observations the data provide a value for, and what still reaches a sink
        init_live_schemes_(builder);

        // parameters that may move latent/observation models
        builder.init_param_producers(param_);
        init_latent_blocks_(builder);
    }

    void Augmented::init_latent_blocks_(augm::SchemeBuilder &builder)
    {
        state.block_factors.clear();
        state.block_factor_of.assign(data_.dataset().report_indiv.total, na::INDEX);
        std::vector<bool> is_built = {};
        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house.infected_indivs())
            {
                // members of one status scheme share its latent block
                const size_t i = state.of(indiv).idx;
                if (i >= is_built.size())
                    is_built.resize(i + 1, false);
                if (!is_built[i])
                {
                    builder.init_latent_block(i, list_block_keys_(indiv));
                    is_built[i] = true;
                }
                if (list_block(indiv).empty())
                    continue;
                state.block_factor_of[Index::of(indiv)] = state.block_factors.size();
                state.block_factors.emplace_back();
            }
    }

    std::vector<size_t> Augmented::list_block_keys_(const Individual &indiv) const
    {
        const augm::StatusScheme &ss = scheme_of_(indiv);
        std::vector<size_t> candidates = {};
        for (size_t a : ss.sample)
            if (at(a, indiv).can_join_block())
                candidates.push_back(a);

        // candidates each value moves with: itself, through functions, none for any other row
        dict<std::vector<size_t>> moved_by;
        const auto list_moved_by = [&](const std::vector<std::string> &names)
        {
            std::vector<size_t> out;
            for (const std::string &astr : names)
                if (const auto it = moved_by.find(astr); it != moved_by.end())
                    vect::append_unique(out, it->second);
            return out;
        };

        // a latent whose bounds or location shift move with another candidate moves alone
        std::vector<size_t> keys = {};
        for (size_t k : ss.order)
        {
            const augm::DefRow &row = def_[k];
            if (row.has_func())
                moved_by[row.name] = list_moved_by(row.func.list_latent_producers());
            const size_t a = consult(row.name);
            if (!vect::contains(candidates, a))
                continue;
            std::vector<std::string> support = vect::concatenate(row.min.list_latent_producers(),
                                                                 row.max.list_latent_producers());
            for (const parse::Entry &child : row.model.root.children)
                if (child.key == "loc")
                    child.append_latent_producers(support);
            if (list_moved_by(support).empty())
                keys.push_back(a);
            moved_by[row.name] = {a};
        }
        return keys;
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

        for (const Household &house : data_.dataset())
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
        // a drawn member deterministic observation has a value per infection
        std::vector<size_t> per_infection;
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            if (!dobs.is_own_column())
                per_infection.push_back(dobs.row);
        augm::check::every_member_cols(def_, data_.dataset(), per_infection);
    }

    //////  Augmented : deterministic observations  //////

    std::vector<std::string> Augmented::get_produced_cols()
    {
        std::vector<std::string> produced_cols = consult(scheme_.list_obs_member);
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            produced_cols.push_back(def_[dobs.row].name);
        return produced_cols;
    }

    std::vector<std::string> Augmented::get_drawn_cols() const
    {
        std::vector<std::string> drawn_cols = consult(scheme_.list_obs_member);
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            if (!dobs.is_static)
                drawn_cols.push_back(def_[dobs.row].name);
        return drawn_cols;
    }

    void Augmented::init_house_for_simulation()
    {
        data::Dataset &dataset = data_.dataset();

        // generated household observations
        std::vector<std::string> drawn_house;
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            if (dobs.is_house && !dobs.is_static)
                drawn_house.push_back(def_[dobs.row].name);
        augm::check::no_drawn_house_obs_cycle(def_, drawn_house);
        assert_every_member_cols_();
        data_.reset_readcells(get_produced_cols());

        // reserve household-level observations
        state.house_varia.clear();
        if (!scheme_.list_obs_house.empty())
            state.house_varia.assign(dataset.report_house.total,
                                     std::vector<Variable>(scheme_.list_obs_house.size()));

        const augm::Binder binder = binder_();
        for (const Household &house : dataset)
        {
            // wire stochastic member observations to their data read-cells
            for (const Individual &indiv : house)
                for (size_t a : scheme_.list_obs_member)
                    data_.state.wire_readcell(dataset.consult(consult(a)), indiv,
                                              &state.of(indiv).vec[a]);

            // build and wire deterministic observations in dependency order
            for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
                // defer drawn member observation, built with the member
                if (dobs.is_house || dobs.is_static)
                    binder.bind_deterministic_obs(state, dobs, house, false);
        }
    }

    void Augmented::init_members_for_simulation()
    {
        if (empty() || !is_allocated())
            return;

        // member variables bound for the whole run: static deterministic observations
        sim_fixed_.assign(size(), false);
        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
            if (dobs.is_static && !dobs.is_house)
                sim_fixed_[consult(def_[dobs.row].name)] = true;

        // bind one block per reachable status scheme, swapped in by the household loop
        state.sim.assign(state.varia.size(), {});
        const augm::Binder binder = binder_();
        msg::WarningBuffer warnings;
        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house)
                bind_simulation_blocks_(binder, indiv);
        warnings.flush();
    }

    //////  Augmented : data-bound deterministic observations  //////

    void Augmented::materialize_deterministic_obs(bool for_export)
    {
        if (scheme_.list_deterministic_obs.empty())
            return;

        data::Dataset &dataset = data_.dataset();
        const augm::Binder binder = binder_();

        // build Variables of deterministic observations, shared with the simulation init
        // point each standing-in data column's read-surface at its Variable
        // data readers then resolve the computed value without a data copy
        state.house_varia.assign(dataset.report_house.total,
                                 std::vector<Variable>(scheme_.list_obs_house.size()));

        for (const augm::DeterministicObs &dobs : scheme_.list_deterministic_obs)
        {
            const std::string &astr = def_[dobs.row].name;
            // a deterministic observation must resolve from data
            // a latent, parameter or group-rate source has no fixed value in inference
            // it would impose an indicator constraint
            if (!def_.is_pure_data(dobs.row))
            {
                // exporting the data view is not inferring: the column keeps its datum, having no computed value
                // the run wall still refuses it
                if (for_export)
                    continue;
                msg::stop(
                    "Cannot infer from " + str::format(def_[dobs.row].kind) + " " +
                    msg::val(astr) + ": it must resolve from data alone.\n" +
                    "A definition reading a latent variable or a parameter is only meaningful in simulation.");
            }

            for (Household &house : dataset)
                binder.bind_deterministic_obs(state, dobs, house, true);
        }
    }

    //////  Augmented : init missing attribute  //////

    void Augmented::clear_all()
    {
        if (empty() || !is_allocated())
            return;

        for (const Household &house : data_.dataset())
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
        // * the observed one, kept by a model that does not redraw it;
        // * every status the status process can draw.
        if (indiv.is_infected_status())
            reach(indiv.status());
        for (StatusLvl s : status_process_.drawable(indiv))
            reach(s);

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
                    // bound for the whole run
                    if (st.vec[a].has_definition())
                        continue;
                    if (vect::contains(user_transm_vec_, a))
                    {
                        generated.push_back(a);
                        st.vec[a].set(0.);
                        // bare: only transmission sets it
                        if (!row.is_bare())
                            binder.bind_member(block[a], k, indiv, true);
                        continue;
                    }
                    bound.push_back(a);
                    binder.bind_member(st.vec[a], k, indiv, true);
                }

                // the profile reads undrawn values: judged on the intervals they reach
                if (infec_profile_.is_latent())
                    if (const Function *f = st.vec[infec_profile_.latent_key].function())
                        infec_profile_.assert_zero_at_origin(*f);

                // finalize while the member's variables hold this scheme's definitions
                for (size_t a : bound)
                    st.vec[a].finalize();
                for (size_t a : generated)
                    block[a].finalize();
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
                data_.dataset().str(indiv) + " with status " + std::to_string(s) + ".");

        msg::stop(
            "Generic latent definition for non-infectious " +
            data_.dataset().str(indiv) + " was missing.");
    }

    //////  Augmented : user reset  //////

    void Augmented::restart(RndGen &prng)
    {
        if (!scheme_.is_any_sample())
            return;

        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house.infected_indivs())
                restart(indiv, prng);
    }

    void Augmented::draw_sample(RndGen &prng)
    {
        if (!scheme_.is_any_sample())
            return;

        for (const Household &house : data_.dataset())
            for (const Individual &indiv : house.infected_indivs())
            {
                augm::IndivState &st = state.of(indiv);
                // definition order: a latent reading a pinned one follows it
                for (size_t a : scheme_of_(indiv).sample)
                {
                    const double x = pins_.empty() || !pins_[a] ? na::DOUBLE : (*pins_[a])[Index::of(indiv)];
                    if (na::is(x))
                        st.vec[a].draw_value(prng);
                    else
                        st.vec[a].set(x);
                }
            }
    }

    void Augmented::pin_values(const dbl_table &table)
    {
        pins_.assign(size(), nullptr);
        for (const auto &[name, values] : table)
            if (const size_t a = find(name); !na::is(a))
                pins_[a] = &values;
    }

    void Augmented::finalize()
    {
        if (empty())
            return;

        for (const Household &house : data_.dataset())
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
        augm::check::obs_init_error(def_, rows, param_, data_.dataset().str(indiv));
    }

    [[noreturn]] void Augmented::throw_latent_init_error(const Individual &indiv) const
    {
        // latents whose log-likelihood is not finite at their value
        std::vector<std::string> names;
        for (size_t a : list_latent_models(indiv))
            if (!std::isfinite(state.of(indiv).vec[a].guarded_log_pdf()))
                names.push_back(consult(a));
        msg::stop(
            "The latent model of " + str::format(names, false, 3) + " is not finite at the initial values of " +
            data_.dataset().str(indiv) + ".\n" +
            "Check the parameters the latent model reads and the values it accepts.");
    }

    augm::Binder Augmented::binder_()
    {
        return augm::Binder(def_, *this, translator_, data_);
    }
}
