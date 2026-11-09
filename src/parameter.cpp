/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <map>
#include <numeric>
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
#include "data/dataset.h"
#include "data/household.h"
#include "data/individual.h"
#include "distribution.h"
#include "function.h"
#include "parameter.h"
#include "proposal.h"
#include "rates.h"
#include "rates/definition.h"
#include "registry.h"
#include "translator.h"
#include "variable.h"

namespace housepi
{
    //////////////////////////////  Variable  //////////////////////////////

    Variable::Variable() = default;
    Variable::~Variable() = default;

    Variable::Variable(const Variable &other)
        : Attribute(other),
          old_(other.old_),
          distr_(other.distr_ ? other.distr_->clone_distrib() : nullptr) {}
    Variable::Variable(Variable &&) noexcept = default;

    Variable &Variable::operator=(Variable other)
    {
        swap(other);
        return *this;
    }

    void Variable::swap(Variable &other) noexcept
    {
        Attribute::swap(other);
        std::swap(move_mode_, other.move_mode_);
        std::swap(old_, other.old_);
        std::swap(distr_, other.distr_);
        std::swap(prop_, other.prop_);
    }

    void Variable::init_distr(std::string_view opt, bool discr)
    {
        distr_ = registry::create_distribution(
            opt, ErrType::Distribution, registry::distr_search_order(discr));
    }

    void Variable::reshape_distribution(bool is_discr)
    {
        distr_->reshape(is_discr);

        const Domain &domain = distr_->domain();
        domain.assert_consistent(str());

        // nothing to pin: a non-degenerate domain, or a live bound or location
        if (!domain.is_degenerate() || !domain.min.is_constant() || !domain.max.is_constant() ||
            (!domain.loc.is_empty() && !domain.loc.is_constant()))
            return;

        // bounds named by source, identical across members so the warning is raised once
        const std::string detail = msg::alert(msg::arg("min") + " (" + domain.min.str_source() + ") equals " +
                                              msg::arg("max") + " (" + domain.max.str_source() + ").");
        const double x0 = domain.relocate(domain.min.value());
        const Codomain point{x0, x0, Interval::Closed, Interval::Closed, distr_->support().lattice};

        // observation or user constant: its model replaced by a point mass, its value kept
        if (err_type_ == ErrType::Observation || is_user_const_)
        {
            msg::warn(str() + " has a single possible value for (at least) one member.\n" + detail);
            return pin_to(point);
        }

        // parameter
        if (err_type_ == ErrType::Parameter)
            msg::stop("Found a degenerate domain for " + str() + ".\n" + detail);

        // latent: pinned to the point, its degenerate model replaced by a pinned one
        msg::warn("Collapsed (at least) one " + str() + " to constant.\n" + detail);
        init_numeric(x0);
        pin_to(point);
    }

    void Variable::pin_to(Codomain range)
    {
        // pinned distribution reporting the range
        distr_ = std::make_unique<PinnedDistribution>(range);
        // proposal that never moves
        prop_ = std::make_unique<NullProposal>();
        // no ratio, no validity check
        move_mode_ = Move::Forward;
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
            // none declared: the default random walk
            // "forward", or distribution name: independence sampler
            // "rw", "rw_log", "rw_unif": random walk tuned by a "scale" argument
            if (proposal.empty())
                helper_init_default_proposal_();
            else if (proposal == "forward")
                helper_init_forward_proposal_();
            else if (is_random_walk(proposal))
            {
                const auto it = args.find("scale");
                helper_init_random_walk_proposal_(proposal, it != args.end() ? it->second.value() : na::DOUBLE);
            }
            else
                helper_init_independence_proposal_(proposal, args);
        }
        catch (const std::exception &e)
        {
            msg::stop(msg::join(e.what(), "for " + str()));
        }
    }

    void Variable::helper_init_forward_proposal_()
    {
        prop_ = std::make_unique<ForwardProposal>(distr_.get());
        move_mode_ = Move::Forward;
    }

    void Variable::helper_init_random_walk_proposal_(const std::string &proposal,
                                                     double scale)
    {
        const Distribution *prior = distr_.get();
        if (proposal == RndWalkUnifProposal::NAME)
            prop_ = std::make_unique<RndWalkUnifProposal>(scale, prior);
        else if (proposal == RndWalkProposal::NAME)
            prop_ = std::make_unique<RndWalkProposal>(scale, prior, RndWalkProposal::Coord::Identity);
        // RndWalkProposal::LOG_NAME
        else
            prop_ = std::make_unique<RndWalkProposal>(scale, prior, RndWalkProposal::Coord::Log);

        move_mode_ = Move::RandomWalk;
    }

    void Variable::helper_init_default_proposal_()
    {
        // the proposal refuses a value rounded onto a bound, or a lattice value outside the domain
        prop_ = std::make_unique<RndWalkProposal>(distr_.get());
        move_mode_ = Move::ProposalChecked;
    }

    void Variable::helper_init_independence_proposal_(const std::string &proposal,
                                                      const dict<Attribute> &args)
    {
        // unknown name -> suggest the random walks and `forward` alongside the distributions
        const registry::SearchOrder order = registry::distr_search_order(distr_->is_discretized());
        if (!order.contains(proposal))
            msg::stop(
                "Unknown proposal " + msg::val(proposal) + "." +
                str::find_closest(proposal, vect::concatenate(kernel_names(), order.names()),
                                  "See the " + msg::arg("proposal") + " argument in " + msg::help("add_parameter") + "."));

        // valid by construction: no per-move validity check
        prop_ = std::make_unique<IndependenceProposal>(proposal, args, distr_.get());
        move_mode_ = Move::Independence;
    }

    void Variable::restart(RndGen &prng)
    {
        // @todo check whether this guard is reachable
        if (is_constant() || !has_distr())
            return;

        // independence proposal: absolute, data-informed start, e.g. near the observed datum
        // otherwise: a prior draw
        if (prop_->is_independent())
            set(prop_->propose(prng, value()));
        else
            set(distr_->rand(prng));
    }

    void Variable::finalize()
    {
        Attribute::finalize();
        if (distr_)
            distr_->finalize();
        if (prop_)
            prop_->finalize(str());
    }

    void Variable::clear() noexcept
    {
        mode_ = Mode::Empty;
        value_ = na::DOUBLE;
        attr_ = nullptr;
        func_ = nullptr;
        distr_ = nullptr;
        prop_ = nullptr;
        move_mode_ = Move::Forward;
    }

    //////////////////////////////  param::DefTable  //////////////////////////////

    namespace param
    {
        DefRow::DefRow(const str_table &input, size_t r)
            : name(input.at("name")[r]),
              value(input.at("const")[r]),
              prior(input.at("prior")[r]),
              func(input.at("func")[r]),
              min(input.at("min")[r]),
              max(input.at("max")[r]),
              proposal(input.at("proposal")[r]),
              is_discr(str::tobool(input.at("discr")[r]))
        {
            // a function parameter may only read other parameters
            if (!has_func())
                return;
            if (func.consumes_latent())
                msg::stop(
                    "Parameter function " + msg::val(name) +
                    " cannot depend on latent variables.");
            if (func.consumes_dataset())
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
                             { return row.has_func(); }))
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
                           const KeyHolder &keys,
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
                size_t p = keys.consult(row.name);
                if (row.has_func())
                {
                    constant.push_back(p);
                    func.push_back(p);
                }
                else if (row.is_const())
                {
                    constant.push_back(p);
                }
                else
                {
                    sample.push_back(p);
                    if (is_declared_adaptive(row.proposal))
                        adaptive.push_back(p);
                }
            }

            // no parameter feeds the transmission cache until marked
            sample_standard = sample;
            // one slot per parameter
            simulation_values.assign(keys.size(), nullptr);
        }

        void Scheme::pin(size_t p, const std::vector<double> &values)
        {
            vect::remove(sample, {p});
            vect::remove(sample_standard, {p});
            vect::remove(sample_cached, {p});
            vect::remove(adaptive, {p});
            if (!vect::contains(constant, p))
                constant.push_back(p);
            simulation_values[p] = &values;
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
        void State::allocate(const KeyHolder &keys)
        {
            vec.clear();
            vec.resize(keys.size());
            for (size_t p = 0; p < vec.size(); ++p)
                vec[p].name(keys.consult(p), ErrType::Parameter);
        }

        void restart(const Scheme &scheme, State &state, RndGen &prng)
        {
            for (size_t p : scheme.sample)
                state.vec[p].restart(prng);
        }

        void load_sample_values(const Scheme &scheme, State &state, size_t r)
        {
            for (size_t p = 0; p < scheme.simulation_values.size(); ++p)
                if (const std::vector<double> *vec = scheme.simulation_values[p])
                    state.vec[p].set((*vec)[r]);
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

    void Parameter::build(const Translator &translator)
    {
        clear_optional_();
        std::vector<size_t> order = def_.order();

        // user parameters take keys after model-declared ones, producer first
        for (size_t k : order)
            if (!count(def_[k].name))
                add_(def_[k].name);

        scheme_.build(def_, *this, std::move(order));
        state.allocate(*this);
        bind_(translator);
    }

    namespace
    {
        /// Closed range of the non-NA `values`, on the integers when every one is whole.
        /// @note No non-NA value -> empty range, its infimum above its supremum.
        Codomain closed_range_(const std::vector<double> &values)
        {
            Codomain c{inf::DOUBLE, -inf::DOUBLE, Interval::Closed, Interval::Closed, Lattice::Integers};
            for (double x : values)
            {
                if (na::is(x))
                    continue;
                c.inf = std::min(c.inf, x);
                c.sup = std::max(c.sup, x);
                if (x != std::floor(x))
                    c.lattice = Lattice::Doubles;
            }
            return c;
        }
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        // classification from build() is kept: a parameter with no table column is drawn
        scheme_.simulation_values.assign(size(), nullptr);
        for (const std::string &pstr : all())
            if (auto found = input.find(pstr); found != input.end())
            {
                const size_t p = consult(pstr);
                Variable &var = at(p);
                // a value per repetition: readers follow it, judged on the range of the column
                var.release_user_const();
                if (const Codomain range = closed_range_(found->second); range.inf <= range.sup)
                {
                    var.pin_to(range);
                    assert_requests(p, var);
                }
                scheme_.pin(p, found->second);
            }
    }

    bool Parameter::is_rate_only(size_t p) const
    {
        if (is_declared(p) || vect::contains(scheme_.sample_cached, p))
            return false;
        const std::string &pstr = consult(p);
        for (const param::DefRow &row : def_.rows)
            if (vect::contains(row.func.list_parameter_producers(), pstr))
                return false;
        return true;
    }

    void Parameter::bind_(const Translator &translator)
    {
        // producers bound before consumers
        for (size_t k : scheme_.order)
        {
            const param::DefRow &row = def_[k];
            size_t p = consult(row.name);
            Variable &var = at(p);

            if (row.has_func())
                bind_function_(var, row, translator);
            else if (row.is_const())
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
            row.func.root.func, ErrType::Function, registry::func_search_order()));

        // bind function arguments
        for (const auto &child : row.func.root.children)
            var.init_func_attr(child.key, translator.to_attribute(child));

        // domain
        var.init_func_attr("min", translator.to_attribute(row.min.root));
        var.init_func_attr("max", translator.to_attribute(row.max.root));

        // discretization
        var.reshape_function(row.is_discr);
    }

    //////////////////////////////  rates::Anchor  //////////////////////////////

    namespace rates
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

        //////////////////////////////  rates::Store  //////////////////////////////

        void Store::finalize()
        {
            for (Attribute &rate : flat_)
                rate.finalize();
        }

        void Store::reset_store_(size_t n_owners)
        {
            flat_.clear();
            base_.assign(n_owners, 0);
            houses_of_param_.clear();
        }

        //////////////////////////////  rates::Group  //////////////////////////////

        std::vector<std::string> Group::build(data::Dataset &dataset,
                                              const Parameter &param,
                                              const str_table &input)
        {
            // drop the previous build
            clear_optional_();
            reset_store_(dataset.report_indiv.total);
            scale_anchors_.clear();
            slot_base_.fill(0);
            n_slots_ = 1;

            // covariate headers, or skip when nothing is declared nor defined
            std::vector<std::string> spec_headers = {};
            if (!input.empty())
                spec_headers = helper_process_definitions_(input);
            else if (empty())
                return {};
            scale_anchors_.assign(size(), {});

            // group-rate name column, if any
            std::vector<std::string> g_input_vec = {};
            if (input.count("group_rate"))
                g_input_vec = input.at("group_rate");
            assert_status_stratified_(input, g_input_vec);

            // a storage slot of its own for every status level a rule keys on
            StatusTable<StatusSlot> slot_of = {};
            if (input.count(rules::STATUS_KEY))
                slot_of = helper_init_status_slots_(input.at(rules::STATUS_KEY));

            const GroupDefinitions defs(input, dataset, param, spec_headers, slot_of);

            // scale anchoring per rate name and rule group
            for (const auto &[key, anchor] : read_scale_anchors<GroupKey>(
                     defs.row_param, read_row_units(input), [&](size_t r)
                     { return GroupKey{g_input_vec[r], defs.rules.row_gid[r]}; }))
                scale_anchors_[consult(key.first)].push_back(anchor);

            // unique definitions: parameter keys per status slot and rate, the reference left out
            const size_t n_keys = size();
            std::map<std::vector<std::vector<size_t>>, RateOffset> unique_defin_map;
            for (Household &house : dataset)
                for (Individual &indiv : house)
                {
                    std::vector<std::vector<size_t>> cell_params(n_slots_ * n_keys);
                    for (size_t s = 0; s < n_slots_; ++s)
                        for (auto &[key, r] : defs.resolve(dataset, indiv, static_cast<StatusSlot>(s)))
                            if (const size_t p = defs.row_param[r]; !na::is(p))
                            {
                                cell_params[s * n_keys + consult(key.first)].push_back(p);
                                mark_reads_(p, Index::of(house));
                            }

                    auto [it, is_new] = unique_defin_map.try_emplace(
                        std::move(cell_params), static_cast<RateOffset>(flat_.size()));
                    if (is_new)
                        helper_init_unique_rate_(param, it->first);
                    base_[Index::of(indiv)] = it->second;
                }

            // the status is carried by the member:
            // no data key, the equivalence-class key holds it
            return vect::difference(spec_headers, {rules::STATUS_KEY});
        }

        std::vector<std::string> Group::helper_process_definitions_(const str_table &input)
        {
            const std::vector<std::string> file_headers = {"group_rate", "param_name", "n", "rule_group"};

            // every column but the template headers is a covariate
            std::vector<std::string> data_headers = extract_keys(input);
            vect::remove(data_headers, file_headers);

            // register the group-rate names
            for (const std::string &gstr : str::unique(input.at("group_rate")))
                add_(gstr);

            return data_headers;
        }

        void Group::assert_status_stratified_(const str_table &input,
                                              const std::vector<std::string> &g_input_vec) const
        {
            if (!input.count(rules::STATUS_KEY))
                return;
            const std::vector<std::string> &status_col = input.at(rules::STATUS_KEY);
            for (size_t r = 0; r < status_col.size(); ++r)
            {
                if (na::is(status_col[r]))
                    continue;
                const size_t g = consult(g_input_vec[r]);
                if (vect::contains(status_stratified_, g))
                    continue;
                msg::stop(
                    "Group rate " + msg::val(g_input_vec[r]) +
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
                    static_cast<SlotOffset>(slot_of[level] * size());
            }
            return slot_of;
        }

        void Group::helper_init_unique_rate_(const Parameter &param,
                                             const std::vector<std::vector<size_t>> &cell_params)
        {
            // every cell starts at the reference: status rules never cover escaped members
            const std::vector<std::string> names = all();
            for (size_t c = 0; c < cell_params.size(); ++c)
            {
                flat_.emplace_back(VALUE_REFERENCE);
                flat_.back().name(names[c % names.size()], ErrType::GroupRate);

                // product of the parameters the cell resolves to
                std::vector<const Attribute *> param_ptrs = {};
                for (size_t p : cell_params[c])
                    param_ptrs.push_back(param.at(p).pointer());
                if (!param_ptrs.empty())
                    flat_.back().init_product(param_ptrs);
            }
        }

        //////////////////////////////  rates::Contact  //////////////////////////////

        std::vector<std::string> Contact::build(data::Dataset &dataset,
                                                const Parameter &param,
                                                const str_table &input)
        {
            reset_store_(dataset.report_house.total);
            scale_anchors_.clear();
            if (input.empty())
                return {};

            ContactDefinitions defs(input, dataset, param);

            // scale anchoring per rule group
            for (const auto &[gid, anchor] : read_scale_anchors<size_t>(
                     defs.row_param, read_row_units(input), [&](size_t r)
                     { return defs.rules.row_gid[r]; }))
                scale_anchors_.push_back(anchor);

            // pairs of covariate classes resolve on first encounter only
            const rules::CovariateClasses classes(dataset, defs.keys);
            std::vector<PairRate> pair_rates(classes.count * classes.count);

            // every household's matrix, row-major, into the one flat store
            size_t n_cells = 0;
            for (size_t H = 0; H < dataset.report_house.total; ++H)
            {
                base_[H] = static_cast<RateOffset>(n_cells);
                const size_t house_size = dataset.at(H).size();
                n_cells += house_size * house_size;
            }
            flat_.assign(n_cells, VALUE_REFERENCE);

            for (size_t H = 0; H < dataset.report_house.total; ++H)
            {
                Household &house = dataset.at(H);
                size_t house_size = house.size();

                for (size_t from = 0; from < house_size; ++from)
                    for (size_t to = 0; to < house_size; ++to)
                    {
                        PairRate &pair =
                            pair_rates[classes.of_indiv[Index::of(house.at(from))] * classes.count +
                                       classes.of_indiv[Index::of(house.at(to))]];
                        if (!pair.resolved)
                        {
                            pair.params = defs.associate_params(dataset, house.at(from), house.at(to));
                            for (size_t p : pair.params)
                                pair.param_ptrs.push_back(param.at(p).pointer());
                            pair.resolved = true;
                        }

                        if (!pair.param_ptrs.empty())
                            flat_[base_[H] + from * house_size + to].init_product(pair.param_ptrs);
                        // transmission reads a contact from an infected member only
                        if (house.at(from).is_infected_status())
                            for (size_t p : pair.params)
                                mark_reads_(p, H);
                    }
            }
            return defs.data;
        }
    }

    //////////////////////////////  Rates  //////////////////////////////

    bool Rates::same_partition(const Household &house,
                               MemberIdx i,
                               MemberIdx j) const
    {
        return group_.same_group(house.at(i), house.at(j)) &&
               contact_.same_contacts(house, i, j);
    }

    std::vector<std::vector<HouseIdx>> Rates::list_houses_of_param(size_t n_params) const
    {
        std::vector<std::vector<HouseIdx>> out(n_params);
        const rates::Store *group = &group_, *contact = &contact_;
        for (const rates::Store *store : {group, contact})
            for (size_t p = 0; p < store->houses_of_param().size(); ++p)
                vect::append(out[p], store->houses_of_param()[p]);
        for (std::vector<HouseIdx> &houses : out)
            vect::sort_unique(houses);
        return out;
    }

    std::vector<std::vector<size_t>> Rates::list_param_rw_blocks(const Parameter &param) const
    {
        const auto is_movable = [&](size_t p)
        { return vect::contains(param.scheme().sample, p) && param.at(p).can_join_block(); };

        // candidate members: each sampled baseline with its heavy rates, then each requested set
        std::vector<std::vector<size_t>> sets = {};
        for (const rates::Anchor &anchor : anchors_)
        {
            if (!is_movable(anchor.baseline))
                continue;
            std::vector<size_t> set = {anchor.baseline};
            for (size_t g : anchor.rates)
                for (const rates::ScaleAnchor &scaling : group_.scale_anchors(g))
                    vect::append(set, scaling.list_heavy());
            if (anchor.contacts)
                for (const rates::ScaleAnchor &scaling : contact_.scale_anchors())
                    vect::append(set, scaling.list_heavy());
            sets.push_back(std::move(set));
        }
        vect::append(sets, param.rw_blocks());

        std::vector<std::vector<size_t>> blocks = {};
        for (const std::vector<size_t> &set : sets)
        {
            std::vector<size_t> block = {};
            for (size_t p : set)
                if (is_movable(p) && !vect::contains(block, p))
                    block.push_back(p);

            // blocks sharing a member merge
            for (size_t b = 0; b < blocks.size();)
            {
                const bool is_shared = std::any_of(blocks[b].begin(), blocks[b].end(),
                                                   [&](size_t p)
                                                   { return vect::contains(block, p); });
                if (!is_shared)
                {
                    ++b;
                    continue;
                }
                vect::append_unique(block, blocks[b]);
                blocks.erase(blocks.begin() + static_cast<std::ptrdiff_t>(b));
            }
            blocks.push_back(std::move(block));
        }

        // a single member moves alone
        blocks.erase(std::remove_if(blocks.begin(), blocks.end(),
                                    [](const std::vector<size_t> &block)
                                    { return block.size() < 2; }),
                     blocks.end());
        return blocks;
    }

    void Rates::assert_identifiability(const Parameter &param) const
    {
        for (const rates::Anchor &anchor : anchors_)
        {
            size_t n_scale_free = 0;
            std::vector<std::string> free_names = {};

            const size_t n_contact = anchor.contacts ? contact_.n_scale_free(param) : 0;
            if (n_contact > 0)
            {
                n_scale_free += n_contact;
                free_names.emplace_back("contact");
            }

            for (size_t g : anchor.rates)
            {
                const size_t n_group = group_.n_scale_free(g, param);
                if (n_group == 0)
                    continue;

                n_scale_free += n_group;
                free_names.emplace_back(group_.consult(g));
            }

            // the R front end holds baselines at the reference for absolute rates
            if (n_scale_free > 1)
                msg::stop(
                    "The product of " + std::to_string(n_scale_free) + " absolute rule groups introduces non-identifiability: " +
                    str::format(free_names) +
                    "\nRates are absolute unless a " + msg::code("reference") + " is given: make all but one relative," +
                    " naming their reference group or asking for " + msg::val(".auto") + ".");
        }
    }
}
