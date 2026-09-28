/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <map>
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
          pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
          pProp_(other.pProp_ ? other.pProp_->clone() : nullptr)
    {
        if (pProp_)
            pProp_->rebind_prior(pDistr_.get());
    }
    Variable::Variable(Variable &&) = default;

    Variable &Variable::operator=(Variable other)
    {
        swap(other);
        return *this;
    }

    void Variable::swap(Variable &other) noexcept
    {
        Attribute::swap(other);
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

        if (!domain.min.is_constant() || !domain.max.is_constant() || is_user_const_)
            return;

        const std::string sms = "Definitions of " + msg::arg("min") + " as " + domain.min.str() +
                                " and " + msg::arg("max") + " as " + domain.max.str() +
                                " produced a degenerate domain";
        // parameter
        if (err_type_ == ErrType::Parameter)
            msg::stop(sms + " for " + str() + ".");

        // latent
        msg::warn("Collapsed (at least) one " + str() + " to constant.\n" +
                  msg::alert(sms + "."));
        init_numeric(domain.min.value());
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
                double scale = na::DOUBLE;
                for (const auto &[key, attr] : args)
                    if (key == "scale")
                        scale = attr.value();
                helper_init_random_walk_proposal_(proposal, scale);
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

        const std::vector<std::string> options = {"rw_unif", "rw", "rw_log"};
        if (proposal == options[0])
            pProp_ = std::make_unique<RndWalkUnifProposal>(scale, prior);
        else if (proposal == options[1])
            pProp_ = std::make_unique<RndWalkProposal>(scale, prior);
        else if (proposal == options[2])
            pProp_ = std::make_unique<RndWalkLogProposal>(scale, prior);
        else
            msg::stop(
                "Unknown proposal " + msg::val(proposal) + "." + str::find_closest(proposal, options));

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

    double Variable::operator()(double x) const
    {
        switch (mode_)
        {
        case Mode::Function:
            return pFunc_->evaluate(x);
        default: // cases Mode::Empty, Mode::Numeric, Mode::Pointer
            msg::stop(
                "Called undefined evaluate operator for non-function " + str() + ".");
        }
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
            dict<size_t> name_to_row;
            for (size_t k = 0; k < rows.size(); ++k)
                name_to_row[rows[k].name] = k;

            // producers first when function parameters exist, row order otherwise
            std::vector<size_t> out;
            out.reserve(rows.size());
            if (std::none_of(rows.begin(), rows.end(), [](const DefRow &row)
                             { return row.has_func; }))
            {
                for (const DefRow &row : rows)
                    out.push_back(name_to_row.at(row.name));
                return out;
            }

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
                scheme_.pin(khold(pstr), found->second);
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

        // prior distribution and its arguments
        var.init_distr(row.prior.root.func, row.is_discr);
        for (const auto &child : row.prior.root.children)
            var.init_distr_attr(child.key, translator.to_numeric(child));
        // parameter domain
        var.init_distr_attr("min", translator.to_numeric(row.min.root));
        var.init_distr_attr("max", translator.to_numeric(row.max.root));
        var.reshape_distribution(row.is_discr);

        // proposal, its arguments folded into the call
        dict<Attribute> prop_args;
        for (const auto &child : row.proposal.root.children)
            prop_args.emplace(child.key, translator.to_numeric(child));
        var.init_proposal(row.proposal.root.func, prop_args);
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
            row.func.root.func,
            ErrType::Function,
            {registry::Kind::Custom,
             registry::Kind::Transcendental,
             registry::Kind::Arithmetic}));

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
        void ScaleAnchor::update(std::string_view pstr,
                                 const Parameter &param)
        {
            if (na::is(pstr) || pstr == REFERENCE)
                anchored = true;
            else if (size_t p = param.khold(pstr); !vect::contains(scaling_params, p))
                scaling_params.push_back(p);
        }

        bool ScaleAnchor::is_free(const Parameter &param) const
        {
            if (anchored || scaling_params.empty())
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

            // group-risk and parameter name columns, if any
            std::vector<std::string> g_input_vec = {};
            std::vector<std::string> p_input_vec = {};
            if (input.count("group_risk"))
            {
                g_input_vec = input.at("group_risk");
                p_input_vec = input.at("param_name");
            }
            assert_status_stratified_(input, g_input_vec);

            // risks no rule defines read the reference
            str_table undefined_risks = {};
            for (const std::string &gstr : khold_.all())
                if (!vect::contains(g_input_vec, gstr))
                    undefined_risks[gstr].push_back(REFERENCE);

            // a storage slot of its own for every status level a rule keys on
            StatusTable<StatusSlot> slot_of = {};
            if (input.count(STATUS_KEY))
                slot_of = helper_init_status_slots_(input.at(STATUS_KEY));

            const GroupDefinitions defs(input, data_wide, spec_headers, slot_of);

            // scale anchoring per risk name and rule group, resolved member by member
            std::map<GroupKey, ScaleAnchor> anchors;
            for (size_t r = 0; r < defs.rules.n_rows; ++r)
                anchors.emplace(GroupKey{g_input_vec[r], defs.rules.row_gid[r]}, ScaleAnchor{});

            // unique definitions, one resolution per status slot
            std::map<std::vector<str_table>, RiskOffset> unique_defin_map;
            for (Household &house : data_wide)
                for (Individual &indiv : house)
                {
                    // status still to be drawn: anchors on every slot; settled: on its own
                    StatusSlot own_slot = 0;
                    const bool any_slot = indiv.is_missing_status();
                    if (!any_slot)
                        own_slot = slot_of[indiv.status()];

                    std::vector<str_table> slot_risks(n_slots_, undefined_risks);
                    for (size_t s = 0; s < n_slots_; ++s)
                    {
                        const auto resolved =
                            defs.resolve(data_wide, indiv, static_cast<StatusSlot>(s));

                        for (auto &[key, r] : resolved)
                            slot_risks[s][key.first].push_back(p_input_vec[r]);

                        if (any_slot || s == own_slot)
                            update_anchors(anchors, resolved, p_input_vec, param);
                    }

                    auto it = unique_defin_map.find(slot_risks);
                    if (it == unique_defin_map.end())
                    {
                        it = unique_defin_map.emplace(
                                                 slot_risks, static_cast<RiskOffset>(flat_.size()))
                                 .first;
                        helper_init_unique_risk_(param, slot_risks);
                    }
                    base_[Index::of(indiv)] = it->second;
                }

            for (const auto &[key, anchor] : anchors)
                scale_anchors_[khold(key.first)].push_back(anchor);
        }

        std::vector<std::string> Group::helper_process_definitions_(data::Wide &data_wide,
                                                                    const str_table &input)
        {
            const std::vector<std::string> file_headers = {"group_risk", "param_name", "rule_group"};

            // every column but the template headers is a covariate
            std::vector<std::string> data_headers = extract_keys(input);
            vect::remove(data_headers, file_headers);

            // register the group-risk names
            dict<size_t> definition_counts = {};
            for (const std::string &gstr : str::unique(input.at("group_risk")))
            {
                khold_.add(gstr);
                definition_counts.emplace(gstr, 0);
            }

            // reserve the cartesian product of distinct definitions per group risk
            for (const std::string &gstr : input.at("group_risk"))
                ++definition_counts[gstr];
            size_t max_unique = 1;
            for (auto &[gstr, count] : definition_counts)
                max_unique *= count;
            flat_.reserve(max_unique * khold_.size());

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
            // other levels share slot 0 with escaped and missing ones, resolving wildcard rules alone
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
                                             const std::vector<str_table> &slot_risks)
        {
            // every cell starts at the reference: status rules never cover escaped members
            const size_t base = flat_.size();
            for (size_t s = 0; s < n_slots_; ++s)
                for (const std::string &gstr : khold_.all())
                {
                    flat_.emplace_back(VALUE_REFERENCE);
                    flat_.back().name(gstr, ErrType::GroupRisk);
                }

            // product of the non-reference parameters each cell resolves to
            for (size_t s = 0; s < n_slots_; ++s)
                for (const auto &[gstr, pstr_vec] : slot_risks[s])
                {
                    std::vector<const Attribute *> param_pointers = {};
                    for (const std::string &pstr : pstr_vec)
                        if (pstr != REFERENCE)
                            param_pointers.emplace_back(
                                param.at(param.khold(pstr)).pointer());

                    if (!param_pointers.empty())
                        flat_[base + s * khold_.size() + khold(gstr)].init_product(param_pointers);
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

            ContactDefinitions defs(input, data_wide);

            data_wide.set_eq_class_data_names(defs.data);

            // scale anchoring per rule group, resolved pair by pair
            std::map<size_t, ScaleAnchor> anchors;
            for (size_t r = 0; r < defs.rules.n_rows; ++r)
                anchors.emplace(defs.rules.row_gid[r], ScaleAnchor{});

            // pairs of covariate classes resolve on first encounter only:
            // ScaleAnchor::update is idempotent, anchors end as one pass over every pair leaves them
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
                                defs.associate_param_ptrs(data_wide, param, house.at(from),
                                                          house.at(to), anchors);
                            pair.resolved = true;
                        }

                        if (!pair.param_ptrs.empty())
                            flat_[base_[H] + from * house_size + to].init_product(pair.param_ptrs);
                    }
            }

            for (const auto &[gid, anchor] : anchors)
                scale_anchors_.push_back(anchor);
        }
    }

    //////////////////////////////  Risks  //////////////////////////////

    bool Risks::same_partition(const Household &house,
                               MemberIdx i,
                               MemberIdx j) const
    {
        if (i == j)
            return true;
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
                    "The product of " + msg::qty(free_names) + " risk{?s} introduces non-identifiability: " + str::format(free_names) +
                    "\nRisks are absolute unless a " + msg::code("reference") + " is given: make all but one relative," +
                    " naming their reference group or asking for " + msg::val(".auto") + ".");
        }
    }
}
