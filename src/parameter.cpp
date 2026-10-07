/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tools/dag.h"
#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "data/household.h"
#include "data/individual.h"
#include "data/wide.h"
#include "distribution.h"
#include "function.h"
#include "parameter.h"
#include "proposal.h"
#include "registry.h"
#include "risks.h"
#include "risks/definition.h"
#include "translator.h"
#include "variable.h"

namespace housepi
{
    //////////////////////////////  Variable  //////////////////////////////

    Variable::Variable() = default;
    Variable::~Variable() = default;

    Variable::Variable(const Variable &other)
        : Attribute(other),
          valid_mode_(other.valid_mode_),
          ratio_mode_(other.ratio_mode_),
          old_(other.old_),
          pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr) {}
    Variable::Variable(Variable &&) = default;

    Variable &Variable::operator=(Variable other)
    {
        swap(other);
        return *this;
    }

    void Variable::swap(Variable &other) noexcept
    {
        Attribute::swap(other);
        std::swap(valid_mode_, other.valid_mode_);
        std::swap(ratio_mode_, other.ratio_mode_);
        std::swap(old_, other.old_);
        std::swap(pDistr_, other.pDistr_);
        std::swap(pProp_, other.pProp_);
    }

    void Variable::init_distr(std::string_view opt, bool discr)
    {
        pDistr_ = registry::create_distribution(
            opt, ErrType::Distribution, registry::distr_kinds(discr));
    }

    void Variable::reshape_distribution(bool is_discr)
    {
        pDistr_->reshape(is_discr);

        const Domain &domain = pDistr_->domain();
        domain.assert_consistent(str());

        // nothing to do on a non-degenerate domain
        if (!domain.is_degenerate())
            return;

        // a live bound or location has no point to pin yet
        if (!domain.min.is_constant() || !domain.max.is_constant() || is_user_const_ ||
            (!domain.loc.is_empty() && !domain.loc.is_constant()))
            return;

        const std::string sms = "Definitions of " + msg::arg("min") + " as " + domain.min.str() +
                                " and " + msg::arg("max") + " as " + domain.max.str() +
                                " produced a degenerate domain";
        // parameter
        if (err_type_ == ErrType::Parameter)
            msg::stop(sms + " for " + str() + ".");

        // latent: pinned to the point, its degenerate model replaced by a null one
        msg::warn("Collapsed (at least) one " + str() + " to constant.\n" +
                  msg::alert(sms + "."));
        const double x0 = domain.relocate(domain.min.value());
        init_numeric(x0);
        pDistr_ = std::make_unique<NullDistribution>(x0);
        helper_init_null_proposal_();
    }

    void Variable::helper_init_null_proposal_()
    {
        // proposal that never moves
        pProp_ = std::make_unique<NullProposal>();

        // no ratio, no validity check
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    //////  Variable : proposal  //////

    void Variable::init_proposal(const std::string &proposal,
                                 const dict<Attribute> &args)
    {
        if (mode_ != Mode::Empty && mode_ != Mode::Numeric)
            return;

        if (mode_ == Mode::Empty)
            init_numeric(na::DOUBLE);
        is_user_const_ = false;

        try
        {
            // "forward" or a distribution name: independence sampler
            // "rw", "rw_log", "rw_unif": random walk tuned by a "scale" argument
            if (proposal == "forward")
                helper_init_distr_proposal_();
            else if (proposal == "rw" || proposal == "rw_log" || proposal == "rw_unif")
            {
                const auto it = args.find("scale");
                helper_init_random_walk_proposal_(proposal, it != args.end() ? it->second.value() : na::DOUBLE);
            }
            else
                helper_init_independence_proposal_(proposal, args);
        }
        catch (const std::exception &e)
        {
            std::string sm = e.what();
            if (!sm.empty() && sm.back() == '.')
                sm.pop_back();
            msg::stop(sm + " for " + str() + ".");
        }
    }

    void Variable::helper_init_distr_proposal_()
    {
        pProp_ = std::make_unique<ForwardProposal>(pDistr_.get());
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    void Variable::helper_init_random_walk_proposal_(const std::string &proposal,
                                                     double scale)
    {
        const Distribution *prior = pDistr_.get();
        if (proposal == "rw_unif")
            pProp_ = std::make_unique<RndWalkUnifProposal>(scale, prior);
        else if (proposal == "rw")
            pProp_ = std::make_unique<RndWalkProposal>(scale, prior);
        else // "rw_log"
            pProp_ = std::make_unique<RndWalkLogProposal>(scale, prior);

        ratio_mode_ = Ratio::Compute;
        valid_mode_ = Validity::Check;
    }

    void Variable::helper_init_independence_proposal_(const std::string &proposal,
                                                      const dict<Attribute> &args)
    {
        // valid by construction: no per-step validity check
        pProp_ = std::make_unique<IndependenceProposal>(proposal, args, pDistr_.get());
        ratio_mode_ = Ratio::Compute;
        valid_mode_ = Validity::Always;
    }

    void Variable::restart(rnd_gen &prng)
    {
        // @todo check whether this guard is reachable
        if (is_constant() || !has_distr())
            return;

        // independence proposal: absolute, data-informed start, e.g. near the observed datum
        // otherwise: a prior draw
        if (pProp_->is_independent())
            set(pProp_->propose(prng, value()));
        else
            set(pDistr_->rand(prng));
    }

    void Variable::finalize()
    {
        Attribute::finalize();
        if (pDistr_)
            pDistr_->finalize();
        if (pProp_)
            pProp_->finalize(str());
    }

    void Variable::clear() noexcept
    {
        mode_ = Mode::Empty;
        value_ = na::DOUBLE;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        pDistr_ = nullptr;
        pProp_ = nullptr;
    }

    //////////////////////////////  param::DefTable  //////////////////////////////

    namespace param
    {
        DefRow::DefRow(const str_table &input, size_t r)
            : name(input.at("name")[r]),
              value(input.at("value")[r]),
              prior(input.at("prior")[r]),
              func(input.at("func")[r]),
              min(input.at("min")[r]),
              max(input.at("max")[r]),
              proposal(input.at("proposal")[r]),
              has_func(!na::is(input.at("func")[r])),
              is_const(str::tobool(input.at("const")[r])),
              is_discr(str::tobool(input.at("discr")[r]))
        {
            // a function parameter may only read other parameters
            if (!has_func)
                return;
            if (func.consumes_latent())
                msg::stop(
                    "Parameter function " + msg::val(name) +
                    " cannot depend on latent variables.");
            if (func.consumes_data_wide())
                msg::stop(
                    "Parameter function " + msg::val(name) +
                    " cannot depend on data variables.");
        }

        void DefTable::build(const str_table &input)
        {
            rows.clear();
            for (size_t r = 0; r < count_rows(input); ++r)
                rows.emplace_back(input, r);
        }

        std::vector<size_t> DefTable::order() const
        {
            // producers first when function parameters exist, row order otherwise
            std::vector<size_t> out(rows.size());
            if (std::none_of(rows.begin(), rows.end(), [](const DefRow &row)
                             { return row.has_func; }))
            {
                std::iota(out.begin(), out.end(), size_t{0});
                return out;
            }

            dict<size_t> name_to_row;
            for (size_t k = 0; k < rows.size(); ++k)
                name_to_row[rows[k].name] = k;

            out.clear();
            dag::prod_graph producer_graph;
            for (const DefRow &row : rows)
            {
                auto &producers = producer_graph[row.name];
                for (const parse::Definition *col : {&row.func, &row.min, &row.max})
                    for (const std::string &pstr : col->list_parameter_producers())
                        producers.insert(pstr);
            }
            for (const std::string &pstr : dag::order_producers_first(producer_graph, ErrType::Parameter))
                out.push_back(name_to_row.at(pstr));
            return out;
        }
    }

    //////////////////////////////  param::Scheme  //////////////////////////////

    namespace param
    {
        void Scheme::build(const DefTable &def,
                           const KeyHolder &khold,
                           std::vector<size_t> rows)
        {
            order = std::move(rows);
            constant.clear();
            func.clear();
            sample.clear();
            adaptive.clear();
            sample_cached.clear();

            for (size_t k : order)
            {
                const DefRow &row = def[k];
                size_t p = khold.from_string(row.name);
                if (row.has_func)
                {
                    constant.push_back(p);
                    func.push_back(p);
                }
                else if (row.is_const)
                {
                    constant.push_back(p);
                }
                else
                {
                    sample.push_back(p);
                    if (parse::is_adaptive_proposal(row.proposal))
                        adaptive.push_back(p);
                }
            }

            // no parameter feeds the transmission cache until marked
            sample_standard = sample;
            // one slot per parameter
            simulation_values.assign(khold.size(), {});
        }

        void Scheme::pin(size_t p, std::vector<double> values)
        {
            vect::remove(sample, {p});
            vect::remove(sample_standard, {p});
            vect::remove(sample_cached, {p});
            vect::remove(adaptive, {p});
            if (!vect::contains(constant, p))
                constant.push_back(p);
            simulation_values[p] = std::move(values);
        }

        void Scheme::mark_produces_transm_cache(size_t p)
        {
            // register as cached, once
            if (!vect::contains(sample_cached, p))
                sample_cached.push_back(p);
            vect::remove(sample_standard, {p});
        }
    }

    //////////////////////////////  param::State  //////////////////////////////

    namespace param
    {
        void State::allocate(const KeyHolder &khold)
        {
            vec.clear();
            vec.resize(khold.size());
            for (size_t p = 0; p < vec.size(); ++p)
                vec[p].name(khold.as_string(p), ErrType::Parameter);
        }

        void restart(const Scheme &scheme, State &state, rnd_gen &prng)
        {
            for (size_t p : scheme.sample)
                state.vec[p].restart(prng);
        }

        void load_sample_values(const Scheme &scheme, State &state, size_t r)
        {
            for (size_t p = 0; p < scheme.simulation_values.size(); ++p)
                if (const std::vector<double> &vec = scheme.simulation_values[p]; !vec.empty())
                    state.vec[p].set(vec[r]);
        }

        void finalize(const Scheme &scheme, State &state)
        {
            for (size_t p : scheme.constant)
                state.vec[p].finalize();
            for (size_t p : scheme.sample)
                state.vec[p].finalize();
        }
    }

    //////////////////////////////  Parameter  //////////////////////////////

    void Parameter::build(const str_table &input,
                          const Translator &translator)
    {
        khold_.clear_optional();

        if (input.empty())
            return;

        def_.build(input);
        std::vector<size_t> order = def_.order();

        // user parameters take keys after model-declared ones, producer first
        for (size_t k : order)
            if (!khold_.count(def_[k].name))
                khold_.add(def_[k].name);

        scheme_.build(def_, khold_, std::move(order));
        state.allocate(khold_);
        bind_(translator);
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        // classification from build() is kept: a parameter with no table column is drawn
        scheme_.simulation_values.assign(khold_.size(), {});
        for (const std::string &pstr : khold_.all())
            if (auto found = input.find(pstr); found != input.end())
            {
                // a value per repetition: readers follow it
                at(khold(pstr)).release_user_const();
                scheme_.pin(khold(pstr), found->second);
            }
    }

    void Parameter::bind_(const Translator &translator)
    {
        // producers bound before consumers
        for (size_t k : scheme_.order)
        {
            const param::DefRow &row = def_[k];
            size_t p = khold(row.name);
            Variable &var = at(p);

            if (row.has_func)
                bind_function_(var, row, translator);
            else if (row.is_const)
                bind_constant_(var, row, translator);
            else
                bind_to_sample_(var, row, translator);

            // check discreteness and positivity
            assert_requests(p, var);
        }
    }

    void Parameter::bind_constant_(Variable &var,
                                   const param::DefRow &row,
                                   const Translator &translator)
    {
        var.transfer(translator.to_numeric(row.value.root));
    }

    void Parameter::bind_to_sample_(Variable &var,
                                    const param::DefRow &row,
                                    const Translator &translator)
    {
        // initial value
        var.transfer(translator.to_numeric(row.value.root));

        // prior distribution within the parameter domain, its arguments constant
        var.init_distr(row.prior.root.func, row.is_discr);
        translator.bind_distr(var, row.prior.root, row.min.root, row.max.root, row.is_discr, nullptr, true);

        // proposal, its arguments folded into the call
        var.init_proposal(row.proposal.root.func, translator.to_numeric_args(row.proposal.root));
    }

    void Parameter::bind_function_(Variable &var,
                                   const param::DefRow &row,
                                   const Translator &translator)
    {
        // pointer mode: a terminal is an alias for another parameter
        if (row.func.is_terminal())
        {
            var.transfer(translator.to_attribute(row.func.root));
            return;
        }

        // function mode
        var.init_function(registry::create_function(
            row.func.root.func, ErrType::Function, registry::func_kinds()));

        // bind function arguments
        for (const auto &child : row.func.root.children)
            var.init_func_attr(child.key, translator.to_attribute(child));

        // domain
        var.init_func_attr("min", translator.to_attribute(row.min.root));
        var.init_func_attr("max", translator.to_attribute(row.max.root));

        // discretization
        var.reshape_function(row.is_discr);
    }

    //////////////////////////////  risks::Anchor  //////////////////////////////

    namespace risks
    {
        bool ScaleAnchor::is_free(const Parameter &param) const
        {
            if (anchored)
                return false;
            for (size_t p : scaling_params)
                if (param.is_pinned(p))
                    return false;
            return true;
        }

        size_t count_scale_free(const std::vector<ScaleAnchor> &anchors,
                                const Parameter &param)
        {
            size_t n = 0;
            for (const ScaleAnchor &anchor : anchors)
                n += anchor.is_free(param);
            return n;
        }

        //////////////////////////////  risks::Store  //////////////////////////////

        void Store::finalize()
        {
            for (Attribute &risk : flat_)
                risk.finalize();
        }

        void Store::reset_store_(size_t n_owners)
        {
            flat_.clear();
            base_.assign(n_owners, 0);
        }

        //////////////////////////////  risks::Group  //////////////////////////////

        void Group::build(data::Wide &data_wide,
                          const Parameter &param,
                          const str_table &input)
        {
            // drop the previous build
            khold_.clear_optional();
            reset_store_(data_wide.report_indiv.total);
            scale_anchors_.clear();
            slot_base_.fill(0);
            n_slots_ = 1;

            // covariate headers, or skip when nothing is declared nor defined
            std::vector<std::string> spec_headers = {};
            if (!input.empty())
                spec_headers = helper_process_definitions_(data_wide, input);
            else if (khold_.empty())
                return;
            scale_anchors_.assign(khold_.size(), {});

            // group-risk name column, if any
            std::vector<std::string> g_input_vec = {};
            if (input.count("group_risk"))
                g_input_vec = input.at("group_risk");
            assert_status_stratified_(input, g_input_vec);

            // a storage slot of its own for every status level a rule keys on
            StatusTable<StatusSlot> slot_of = {};
            if (input.count(STATUS_KEY))
                slot_of = helper_init_status_slots_(input.at(STATUS_KEY));

            const GroupDefinitions defs(input, data_wide, param, spec_headers, slot_of);

            // scale anchoring per risk name and rule group
            for (const auto &[key, anchor] : read_scale_anchors<GroupKey>(
                     defs.rules, [&](size_t r)
                     { return GroupKey{g_input_vec[r], defs.rules.row_gid[r]}; }))
                scale_anchors_[khold(key.first)].push_back(anchor);

            // unique definitions: parameter keys per status slot and risk, the reference left out
            const size_t n_keys = khold_.size();
            std::map<std::vector<std::vector<size_t>>, RiskOffset> unique_defin_map;
            for (Household &house : data_wide)
                for (Individual &indiv : house)
                {
                    std::vector<std::vector<size_t>> cell_params(n_slots_ * n_keys);
                    for (size_t s = 0; s < n_slots_; ++s)
                        for (auto &[key, r] : defs.resolve(data_wide, indiv, static_cast<StatusSlot>(s)))
                            if (const size_t p = defs.rules.row_param[r]; !na::is(p))
                                cell_params[s * n_keys + khold(key.first)].push_back(p);

                    auto [it, is_new] = unique_defin_map.try_emplace(
                        std::move(cell_params), static_cast<RiskOffset>(flat_.size()));
                    if (is_new)
                        helper_init_unique_risk_(param, it->first);
                    base_[Index::of(indiv)] = it->second;
                }
        }

        std::vector<std::string> Group::helper_process_definitions_(data::Wide &data_wide,
                                                                    const str_table &input)
        {
            const std::vector<std::string> file_headers = {"group_risk", "param_name", "rule_group"};

            // every column but the template headers is a covariate
            std::vector<std::string> data_headers = extract_keys(input);
            vect::remove(data_headers, file_headers);

            // register the group-risk names
            for (const std::string &gstr : str::unique(input.at("group_risk")))
                khold_.add(gstr);

            // the status is carried by the member: no data key, the equivalence-class hash mixes it in
            data_wide.set_eq_class_data_names(
                vect::difference(data_headers, {STATUS_KEY}));

            return data_headers;
        }

        void Group::assert_status_stratified_(const str_table &input,
                                              const std::vector<std::string> &g_input_vec) const
        {
            if (!input.count(STATUS_KEY))
                return;
            const std::vector<std::string> &status_col = input.at(STATUS_KEY);
            for (size_t r = 0; r < status_col.size(); ++r)
            {
                if (na::is(status_col[r]))
                    continue;
                const size_t g = khold(g_input_vec[r]);
                if (vect::contains(status_stratified_, g))
                    continue;
                msg::stop(
                    "Group risk " + msg::val(g_input_vec[r]) +
                    " cannot be stratified by the infected status.\n"
                    "It is read on members that carry no infected status yet.");
            }
        }

        StatusTable<StatusSlot> Group::helper_init_status_slots_(const std::vector<std::string> &status_col)
        {
            // levels are narrow: the table is the set of keyed levels
            StatusTable<bool> is_keyed = {};
            for (const std::string &sval : status_col)
                if (!na::is(sval))
                    is_keyed[static_cast<StatusLvl>(str::toi(sval))] = true;

            // one slot per keyed level
            // other levels share slot 0 with escaped and missing ones, resolving catch-all rules alone
            StatusTable<StatusSlot> slot_of = {};
            for (size_t level = 0; level < slot_of.size(); ++level)
            {
                if (is_keyed[level])
                    slot_of[level] = static_cast<StatusSlot>(n_slots_++);
                slot_base_[level] =
                    static_cast<SlotOffset>(slot_of[level] * khold_.size());
            }
            return slot_of;
        }

        void Group::helper_init_unique_risk_(const Parameter &param,
                                             const std::vector<std::vector<size_t>> &cell_params)
        {
            // every cell starts at the reference: status rules never cover escaped members
            const std::vector<std::string> names = khold_.all();
            for (size_t c = 0; c < cell_params.size(); ++c)
            {
                flat_.emplace_back(VALUE_REFERENCE);
                flat_.back().name(names[c % names.size()], ErrType::GroupRisk);

                // product of the parameters the cell resolves to
                std::vector<const Attribute *> param_ptrs = {};
                for (size_t p : cell_params[c])
                    param_ptrs.push_back(param.at(p).pointer());
                if (!param_ptrs.empty())
                    flat_.back().init_product(param_ptrs);
            }
        }

        //////////////////////////////  risks::Contact  //////////////////////////////

        void Contact::build(data::Wide &data_wide,
                            const Parameter &param,
                            const str_table &input)
        {
            reset_store_(data_wide.report_house.total);
            scale_anchors_.clear();
            if (input.empty())
                return;

            ContactDefinitions defs(input, data_wide, param);

            data_wide.set_eq_class_data_names(defs.data);

            // scale anchoring per rule group
            for (const auto &[gid, anchor] : read_scale_anchors<size_t>(
                     defs.rules, [&](size_t r)
                     { return defs.rules.row_gid[r]; }))
                scale_anchors_.push_back(anchor);

            // pairs of covariate classes resolve on first encounter only
            const CovariateClasses classes(data_wide, defs.keys);
            std::vector<PairRisk> pair_risks(classes.count * classes.count);

            // every household's matrix, row-major, into the one flat store
            size_t n_cells = 0;
            for (size_t H = 0; H < data_wide.report_house.total; ++H)
            {
                base_[H] = static_cast<RiskOffset>(n_cells);
                const size_t house_size = data_wide.at(H).size();
                n_cells += house_size * house_size;
            }
            flat_.assign(n_cells, VALUE_REFERENCE);

            for (size_t H = 0; H < data_wide.report_house.total; ++H)
            {
                Household &house = data_wide.at(H);
                size_t house_size = house.size();

                for (size_t from = 0; from < house_size; ++from)
                    for (size_t to = 0; to < house_size; ++to)
                    {
                        PairRisk &pair =
                            pair_risks[classes.of_indiv[Index::of(house.at(from))] * classes.count +
                                       classes.of_indiv[Index::of(house.at(to))]];
                        if (!pair.resolved)
                        {
                            pair.param_ptrs =
                                defs.associate_param_ptrs(data_wide, param, house.at(from), house.at(to));
                            pair.resolved = true;
                        }

                        if (!pair.param_ptrs.empty())
                            flat_[base_[H] + from * house_size + to].init_product(pair.param_ptrs);
                    }
            }
        }
    }

    //////////////////////////////  Risks  //////////////////////////////

    bool Risks::same_partition(const Household &house,
                               MemberIdx i,
                               MemberIdx j) const
    {
        return group_.same_group(house.at(i), house.at(j)) &&
               contact_.same_contacts(house, i, j);
    }

    void Risks::assert_identifiability(const Parameter &param) const
    {
        for (const risks::Anchor &anchor : anchors_)
        {
            size_t n_scale_free = 0;
            std::vector<std::string> free_names = {};

            const size_t n_contact = anchor.contacts ? contact_.n_scale_free(param) : 0;
            if (n_contact > 0)
            {
                n_scale_free += n_contact;
                free_names.emplace_back("contact");
            }

            for (size_t g : anchor.risks)
            {
                const size_t n_group = group_.n_scale_free(g, param);
                if (n_group == 0)
                    continue;

                n_scale_free += n_group;
                free_names.emplace_back(group_.khold(g));
            }

            // the R front end holds baselines at the reference for absolute risks
            if (n_scale_free > 1)
                msg::stop(
                    "The product of " + std::to_string(n_scale_free) + " absolute rule groups introduces non-identifiability: " +
                    str::format(free_names) +
                    "\nRisks are absolute unless a " + msg::code("reference") + " is given: make all but one relative," +
                    " naming their reference group or asking for " + msg::val(".auto") + ".");
        }
    }
}
