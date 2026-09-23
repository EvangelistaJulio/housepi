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
#include <vector>

#include "tools/dependency.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/rule_table.h"
#include "tools/string.h"
#include "tools/types.h"
#include "tools/vector.h"

#include "contact_risk.h"
#include "data_wide.h"
#include "distribution.h"
#include "function.h"
#include "group_risk.h"
#include "household.h"
#include "individual.h"
#include "parameter.h"
#include "proposal.h"
#include "registry.h"
#include "variable.h"
#include "translator.h"

namespace housepi
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** Variable : class ***  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

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

        // skip if the domain is not degenerate
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
        // assign the pointer to the proposal
        pProp_ = std::make_unique<NullProposal>();

        // setup sample proposal functional methods
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Variable : Setter (Proposal)  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

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
            // "forward" and any distribution name are independence samplers;
            // "rw"/"rw_log"/"rw_unif" are random-walk kernels tuned by a "scale" argument.
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

    void Variable::helper_init_independence_proposal_(
        const std::string &proposal, const dict<Attribute> &args)
    {
        // construction guarantees move never needs a per-step validity check
        pProp_ = std::make_unique<IndependenceProposal>(proposal, args, pDistr_.get());
        ratio_mode_ = Ratio::Compute;
        valid_mode_ = Validity::Always;
    }

    void Variable::restart(rnd_gen &prng)
    {
        // UPDATE: this should be useless
        if (is_constant() || !has_distr())
            return;

        // an independence proposal draws an absolute, data-informed value (e.g.
        // centered near the observed datum), seeding an observation-feasible start;
        // otherwise fall back to the prior
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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Parameter::build(const str_table &input, const Translator &translator)
    {
        khold_.clear_optional();

        if (input.empty())
            return;

        // check entries
        const std::vector<std::string> file_headers = {"name", "value", "const", "discr",
                                                       "func", "prior", "min", "max",
                                                       "proposal"};
        const std::vector<std::string> headers = read::extract_keys(input);

        // build the producer graph over parameters and validate function parameters
        bool has_func_params = false;
        dict<size_t> name_to_row;
        for (size_t k = 0; k < input.at("name").size(); ++k)
        {
            const std::string &pstr = input.at("name")[k];
            name_to_row[pstr] = k;

            // validate: a function parameter may only read other parameters
            if (na::is(input.at("func")[k]))
                continue;

            has_func_params = true;
            parse::Definition parser(input.at("func")[k]);
            if (parser.consumes_latent())
                msg::stop(
                    "Parameter function " + msg::val(pstr) +
                    " cannot depend on latent variables.");
            if (parser.consumes_data_wide())
                msg::stop(
                    "Parameter function " + msg::val(pstr) +
                    " cannot depend on data variables.");
        }

        // topological sort only when parameter-functions exist
        // preserves input order otherwise
        std::vector<std::string> ordered;
        if (has_func_params)
        {
            dpdcy::prod_graph producer_graph;
            for (size_t k = 0; k < input.at("name").size(); ++k)
            {
                auto &producers = producer_graph[input.at("name")[k]];
                for (const auto col : {"func", "min", "max"})
                    for (const std::string &pstr : parse::Definition(input.at(col)[k]).list_parameter_producers())
                        producers.insert(pstr);
            }
            ordered = dpdcy::order_producers_first(producer_graph, ErrType::Parameter);
        }
        else
        {
            ordered = input.at("name");
        }

        init_param_vec_();
        clear_lists_();
        // load parameters producer before consumer
        for (const std::string &pstr : ordered)
        {
            size_t k = name_to_row.at(pstr);

            if (!khold_.count(pstr))
            {
                khold_.add(pstr);
                param_vec_.emplace_back();
                param_vec_.back().name(pstr, ErrType::Parameter);
            }

            if (!na::is(input.at("func")[k]))
                helper_build_function_(input, k, translator);
            else if (str::tobool(input.at("const")[k]))
                helper_build_constant_(input, k, translator);
            else
                helper_build_to_sample_(input, k, translator);

            // check discreteness and positivity
            size_t p = khold(pstr);
            assert_requests(p, at(p));
        }

        // initialize all parameters to standard, i.e. none feeds transmission cache
        list_sample_standard = list_sample;
        for (size_t p : list_sample)
            if (at(p).is_adaptive())
                list_adaptive.push_back(p);

        // reserve one slot per parameter
        simulation_values.assign(khold_.size(), {});
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        const std::vector<std::string> headers = read::extract_keys(input);

        // keep the list_const/list_func/list_sample classification from build()
        // a parameter with no table column is drawn instead
        simulation_values.assign(khold_.size(), {});
        for (const std::string &pstr : khold_.all())
        {
            if (!vect::contains(headers, pstr))
                continue;

            size_t p = khold(pstr);
            switch_to_const_(p);
            simulation_values[p] = input.at(pstr);
        }
    }

    void Parameter::init_param_vec_()
    {
        size_t n = khold_.size();
        param_vec_.clear();
        param_vec_.resize(n);
        for (size_t p = 0; p < n; ++p)
            at(p).name(khold(p), ErrType::Parameter);
    }

    void Parameter::clear_lists_()
    {
        list_const.clear();
        list_func.clear();
        list_sample.clear();
        list_sample_standard.clear();
        list_adaptive.clear();
        list_sample_cached.clear();

        simulation_values.clear();
    }

    void Parameter::helper_build_constant_(const str_table &input,
                                           size_t k,
                                           const Translator &translator)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

        // assign constant parameter to the set of variables
        size_t p = khold(pstr);
        at(p).transfer(translator.to_numeric(input.at("value")[k]));

        // update the lists of constant parameters
        list_const.push_back(p);
    }

    void Parameter::helper_build_to_sample_(const str_table &input,
                                            size_t k,
                                            const Translator &translator)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

        // create a new parameter Variable
        size_t p = khold(pstr);
        Variable &xparam = at(p);
        xparam.transfer(translator.to_numeric(input.at("value")[k]));

        // assign a prior distribution to the parameter
        const bool is_discr = str::tobool(input.at("discr")[k]);
        parse::Definition parser(input.at("prior")[k]);
        xparam.init_distr(parser.root.func, is_discr);

        // load prior arguments
        for (const auto &child : parser.root.children)
            xparam.init_distr_attr(child.key, translator.to_numeric(child));
        // load parameter domain
        for (const auto key : {"min", "max"})
            xparam.init_distr_attr(key, translator.to_numeric(input.at(key)[k]));
        xparam.reshape_distribution(is_discr);

        // assign proposal (arguments are folded into the proposal call, if any)
        parse::Definition prop(input.at("proposal")[k]);
        dict<Attribute> prop_args;
        for (const auto &child : prop.root.children)
            prop_args.emplace(child.key, translator.to_numeric(child));
        xparam.init_proposal(prop.root.func, prop_args);

        // update the lists of parameters to sample
        list_sample.push_back(p);
    }

    void Parameter::helper_build_function_(const str_table &input,
                                           size_t k,
                                           const Translator &translator)
    {
        const std::string &pstr = input.at("name")[k];
        size_t p = khold(pstr);
        Variable &xparam = at(p);

        parse::Definition parser(input.at("func")[k]);

        // pointer mode (terminal = alias for another parameter)
        if (parser.is_terminal())
        {
            xparam.transfer(translator.to_attribute(parser.root));
            list_const.push_back(p);
            list_func.push_back(p);
            return;
        }

        // function mode
        xparam.init_function(registry::create_function(
            parser.root.func,
            ErrType::Function,
            {registry::Kind::Custom,
             registry::Kind::Transcendental,
             registry::Kind::Arithmetic}));

        // bind function arguments
        for (const auto &child : parser.root.children)
            xparam.init_func_attr(child.key, translator.to_attribute(child));

        // domain
        for (const auto key : {"min", "max"})
            xparam.init_func_attr(key, translator.to_attribute(input.at(key)[k]));

        // discretization
        xparam.reshape_function(str::tobool(input.at("discr")[k]));

        // classify: deterministic, not sampled
        list_const.push_back(p);
        list_func.push_back(p);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Parameter : risk identifiability  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Parameter::request_risk_identifiability(size_t p,
                                                 std::vector<size_t> risks,
                                                 bool contacts)
    {
        risk_anchors_.push_back({p, std::move(risks), contacts});
    }

    void Parameter::switch_to_const_(size_t p)
    {
        vect::remove(list_sample, {p});
        vect::remove(list_sample_standard, {p});
        vect::remove(list_sample_cached, {p});
        vect::remove(list_adaptive, {p});
        if (!vect::contains(list_const, p))
            list_const.push_back(p);
    }

    void Parameter::assert_risk_identifiability(const GroupRisk &group_risk,
                                                const ContactRisk &contact_risk)
    {
        for (const RiskAnchor &anchor : risk_anchors_)
        {
            size_t n_scale_free = 0;
            std::vector<std::string> free_names = {};

            const size_t n_contact = anchor.contacts ? contact_risk.n_scale_free(*this) : 0;
            if (n_contact > 0)
            {
                n_scale_free += n_contact;
                free_names.emplace_back("contact");
            }

            for (size_t g : anchor.risks)
            {
                const size_t n_group = group_risk.n_scale_free(g, *this);
                if (n_group == 0)
                    continue;

                n_scale_free += n_group;
                free_names.emplace_back(group_risk.khold(g));
            }

            // R front end already holds baselines at reference for absolute risks
            if (n_scale_free > 1)
                msg::stop(
                    "The product of " + msg::qty(free_names) + " risk{?s} introduces non-identifiability: " + str::format(free_names) +
                    "\nRisks are absolute unless a " + msg::code("reference") + " is given: make all but one relative," +
                    " naming their reference group or asking for " + msg::val(".auto") + ".");
        }
    }

    void Parameter::mark_produces_transm_cache(size_t p)
    {
        // add the parameter as cache if not already present
        if (!vect::contains(list_sample_cached, p))
            list_sample_cached.push_back(p);

        // if a parameter was originally defined as standard, remove it
        if (!vect::contains(list_sample_standard, p))
            return;
        vect::remove(list_sample_standard, {p});
    }

    void Parameter::restart(rnd_gen &prng)
    {
        for (size_t p : list_sample)
            at(p).restart(prng);
    }

    void Parameter::finalize()
    {
        for (size_t p : list_const)
            at(p).finalize();
        for (size_t p : list_sample)
            at(p).finalize();
    }

    void Parameter::load_sample_values(size_t r)
    {
        for (size_t p = 0; p < khold_.size(); ++p)
        {
            const std::vector<double> &vec = simulation_values[p];
            if (!vec.empty())
                at(p).set(vec[r]);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Risk : shared mechanisms ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void ScaleAnchor::update(std::string_view pstr, const Parameter &param)
    {
        if (na::is(pstr) || pstr == risks::REFERENCE)
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

    void RiskStorage::finalize()
    {
        for (Attribute &risk : flat_)
            risk.finalize();
    }

    void RiskStorage::reset_store_(size_t n_owners)
    {
        flat_.clear();
        base_.assign(n_owners, 0);
    }

    size_t RiskStorage::count_scale_free_(const std::vector<ScaleAnchor> &anchors,
                                          const Parameter &param)
    {
        size_t n = 0;
        for (const ScaleAnchor &anchor : anchors)
            n += anchor.is_free(param);
        return n;
    }

    namespace
    {
        /// Register on every rule group the risk one resolution read,
        /// 'resolved' holding the rule row each group matched, if any.
        template <class Key>
        void update_anchors_(std::map<Key, ScaleAnchor> &anchors,
                             const std::map<Key, size_t> &resolved,
                             const std::vector<std::string> &param_name,
                             const Parameter &param)
        {
            for (auto &[key, anchor] : anchors)
            {
                // a group no rule of which matched reads the reference value
                std::string_view pstr = {};
                if (auto found = resolved.find(key); found != resolved.end())
                    pstr = param_name[found->second];
                anchor.update(pstr, param);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** GroupRisk : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        using idx_vec = std::vector<size_t>;

        /// A group risk rule is resolved per risk name and rule group.
        using GroupKey_ = std::pair<std::string, size_t>;

        /// Auxiliary helper to extract group risk definitions.
        struct GroupDefinitions_
        {
        public:
            std::vector<std::string> data; /// covariates held by the data
            idx_vec keys;                  /// data-wide keys
            RuleTable rules;               /// rule_group / specificity

            /// Extract the definitions from input headers,
            /// 'spec_headers' listing every covariate a rule may key on.
            GroupDefinitions_(const str_table &input,
                              const DataWide &data_wide,
                              const std::vector<std::string> &spec_headers,
                              const StatusTable<StatusSlot> &slot_of)
                : data(vect::difference(spec_headers, {risks::STATUS_KEY})),
                  keys(data_wide.khold(data)),
                  rules(input, spec_headers),
                  slot_of_(slot_of)
            {
                // the status is carried by the member, every other covariate by the data
                if (input.count(risks::STATUS_KEY))
                    status_col_ = &input.at(risks::STATUS_KEY);
            }

            /// Associate the most-specific matching rule to each (risk name, rule group),
            /// reading the member covariates at status slot 'slot'.
            std::map<GroupKey_, size_t> resolve(const DataWide &data_wide,
                                                const Individual &indiv,
                                                StatusSlot slot) const
            {
                return rules.resolve<GroupKey_>(
                    [&](size_t r)
                    { return matches_(data_wide, r, indiv, slot); },
                    [&](size_t r)
                    { return GroupKey_{rules.input.at("group_risk")[r], rules.row_gid[r]}; },
                    [&](const GroupKey_ &key)
                    {
                        return "group risk " + msg::val(key.first) + " for " +
                               data_wide.str(indiv) + " (group " +
                               std::to_string(key.second) + ")";
                    });
            }

        private:
            const std::vector<std::string> *status_col_ = nullptr; /// status the rules key on, none if absent
            const StatusTable<StatusSlot> &slot_of_;               /// status level to the slot it reads

            /// Whether the member matches the rule at row 'r', read at status slot 'slot'.
            bool matches_(const DataWide &data_wide,
                          size_t r,
                          const Individual &indiv,
                          StatusSlot slot) const
            {
                if (status_col_)
                {
                    const std::string &sval = (*status_col_)[r];
                    if (!na::is(sval) && slot != slot_of_[str::toi(sval)])
                        return false;
                }
                for (size_t i = 0; i < data.size(); ++i)
                {
                    const std::string &dval = rules.input.at(data[i])[r];
                    if (!na::is(dval) && data_wide.data(keys[i], indiv) != str::tod(dval))
                        return false;
                }
                return true;
            }
        };
    }

    void GroupRisk::build(DataWide &data_wide,
                          const Parameter &param,
                          const str_table &input)
    {
        // reset group risk storage and unbind individuals
        khold_.clear_optional();
        reset_store_(data_wide.report_indiv.total);
        scale_anchors_.clear();
        slot_base_.fill(0);
        n_slots_ = 1;

        // load required covariate headers or skip if no definitions provided
        std::vector<std::string> spec_headers = {};
        if (!input.empty())
            spec_headers = helper_process_definitions_(data_wide, input);
        else if (khold_.empty())
            return;
        scale_anchors_.assign(khold_.size(), {});

        // load columns with group risk and parameter names, if any
        std::vector<std::string> g_input_vec = {};
        std::vector<std::string> p_input_vec = {};
        if (input.count("group_risk"))
        {
            g_input_vec = input.at("group_risk");
            p_input_vec = input.at("param_name");
        }
        assert_status_stratified_(input, g_input_vec);

        // check for undefined group risks, assumed to be equal to reference
        str_table undefined_risks = {};
        for (const std::string &gstr : khold_.all())
            if (!vect::contains(g_input_vec, gstr))
                undefined_risks[gstr].push_back(risks::REFERENCE);

        // give every status level a rule keys on a storage slot of its own
        StatusTable<StatusSlot> slot_of = {};
        if (input.count(risks::STATUS_KEY))
            slot_of = helper_init_status_slots_(input.at(risks::STATUS_KEY));

        // precompute rule_group / specificity and resolve column indices once
        const GroupDefinitions_ defs(input, data_wide, spec_headers, slot_of);

        // per-(risk_name, rule_group) scale anchoring, resolved member by member below
        std::map<GroupKey_, ScaleAnchor> anchors;
        for (size_t r = 0; r < defs.rules.n_rows; ++r)
            anchors.emplace(GroupKey_{g_input_vec[r], defs.rules.row_gid[r]}, ScaleAnchor{});

        // unique definition map, one resolution per status slot
        std::map<std::vector<str_table>, RiskOffset> unique_defin_map;
        for (Household &house : data_wide)
            for (Individual &indiv : house)
            {
                // a member whose status is still to be drawn may take any slot,
                // one whose status is settled anchors on the slot it reads
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
                        update_anchors_(anchors, resolved, p_input_vec, param);
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

        // store the scale anchors of every rule group
        for (const auto &[key, anchor] : anchors)
            scale_anchors_[khold(key.first)].push_back(anchor);
    }

    std::vector<std::string> GroupRisk::helper_process_definitions_(DataWide &data_wide,
                                                                    const str_table &input)
    {
        const std::vector<std::string> file_headers = {"group_risk", "param_name", "rule_group"};

        // check that tempalte headers are present
        std::vector<std::string> data_headers = read::extract_keys(input);
        // remove template headers to keep only data
        vect::remove(data_headers, file_headers);

        // assert the group risk names are allowed
        dict<size_t> definition_counts = {};
        for (const std::string &gstr : str::unique(input.at("group_risk")))
        {
            khold_.add(gstr);
            definition_counts.emplace(gstr, 0);
        }

        // allocate unique group risk definitions
        // NOTE: reserve cartesian product of distinct definitions per group risk
        for (const std::string &gstr : input.at("group_risk"))
            ++definition_counts[gstr];
        size_t max_unique = 1;
        for (auto &[gstr, count] : definition_counts)
            max_unique *= count;
        flat_.reserve(max_unique * khold_.size());

        // the status is carried by the member: it is no data key, and the
        // equivalence-class hash already mixes it in
        data_wide.set_eq_class_data_names(
            vect::difference(data_headers, {risks::STATUS_KEY}));

        return data_headers;
    }

    void GroupRisk::assert_status_stratified_(const str_table &input,
                                              const std::vector<std::string> &g_input_vec) const
    {
        if (!input.count(risks::STATUS_KEY))
            return;
        const std::vector<std::string> &status_col = input.at(risks::STATUS_KEY);
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

    StatusTable<StatusSlot> GroupRisk::helper_init_status_slots_(const std::vector<std::string> &status_col)
    {
        // levels are narrow enough for the table to be the set of the keyed ones
        StatusTable<bool> is_keyed = {};
        for (const std::string &sval : status_col)
            if (!na::is(sval))
                is_keyed[static_cast<StatusLvl>(str::toi(sval))] = true;

        // one slot per status a rule keys on. Every level no rule distinguishes
        // shares slot 0 with the escaped and the missing ones: they all resolve
        // the wildcard rules alone, hence read the same risks.
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

    void GroupRisk::helper_init_unique_risk_(const Parameter &param, const std::vector<str_table> &slot_risks)
    {
        const size_t base = flat_.size();
        // a risk the member resolved no rule for stays at the reference:
        // status rules, in particular, never cover the escaped members
        for (size_t s = 0; s < n_slots_; ++s)
            for (const std::string &gstr : khold_.all())
            {
                flat_.emplace_back(risks::VALUE_REFERENCE);
                flat_.back().name(gstr, ErrType::GroupRisk);
            }

        // assign custom definitions
        for (size_t s = 0; s < n_slots_; ++s)
            for (const auto &[gstr, pstr_vec] : slot_risks[s])
            {
                // extract non-reference parameter pointers
                std::vector<const Attribute *> param_pointers = {};
                for (const std::string &pstr : pstr_vec)
                    if (pstr != risks::REFERENCE)
                        param_pointers.emplace_back(
                            param.at(param.khold(pstr)).pointer());

                if (!param_pointers.empty())
                    flat_[base + s * khold_.size() + khold(gstr)].init_product(param_pointers);
            }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** ContactRisk : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        /// Auxiliary helper to extract contact risk definition.
        struct ContactDefinitions_
        {
        public:
            std::vector<std::string> from; /// 'from'-prefixed covariates
            std::vector<std::string> to;   /// 'to'-prefixed covariates
            std::vector<std::string> data; /// covariates
            idx_vec keys;                  /// data-wide keys
            RuleTable rules;               /// rule_group / specificity

            /// Extract from/to definitions from input headers.
            // NOTE: 'from' and 'to' expected to have the same size
            ContactDefinitions_(const str_table &input,
                                const DataWide &data_wide)
                : from(extract_from_(input)),
                  to(extract_to_(input)),
                  data(strip_from_prefix_(from)),
                  keys(data_wide.khold(data)),
                  rules(input, vect::concatenate(from, to)) {}

            /// Associate parameter pointers to pair. Accumulate the scale anchoring of each group.
            std::vector<const Attribute *> associate_param_ptrs(
                const DataWide &data_wide,
                const Parameter &param,
                const Individual &indiv_from,
                const Individual &indiv_to,
                std::map<size_t, ScaleAnchor> &anchors) const
            {
                // per-group: keep the most-specific matching rule
                const auto resolved = rules.resolve<size_t>(
                    [&](size_t r)
                    { return matches_(data_wide, r, indiv_from, indiv_to); },
                    [&](size_t r)
                    { return rules.row_gid[r]; },
                    [&](size_t)
                    {
                        return "contact risk from " + data_wide.str(indiv_from) +
                               " to " + data_wide.str(indiv_to);
                    });

                // across groups: multiply the resolved risks (skip reference)
                std::vector<const Attribute *> param_ptrs;
                param_ptrs.reserve(resolved.size());
                for (auto &[gid, r] : resolved)
                {
                    const auto &pstr = rules.input.at("param_name")[r];
                    if (!na::is(pstr) && pstr != risks::REFERENCE)
                        param_ptrs.push_back(param.at(param.khold(pstr)).pointer());
                }

                update_anchors_(anchors, resolved, rules.input.at("param_name"), param);
                return param_ptrs;
            }

        private:
            static std::vector<std::string> extract_from_(const str_table &input)
            {
                std::vector<std::string> out;
                for (const auto &key : read::extract_keys(input))
                    if (str::starts_with(key, "from_"))
                        out.push_back(key);
                return out;
            }
            static std::vector<std::string> extract_to_(const str_table &input)
            {
                std::vector<std::string> out;
                for (const auto &key : read::extract_keys(input))
                    if (str::starts_with(key, "to_"))
                        out.push_back(key);
                return out;
            }
            static std::vector<std::string> strip_from_prefix_(
                const std::vector<std::string> &from)
            {
                std::vector<std::string> out;
                out.reserve(from.size());
                for (const auto &k : from)
                    out.push_back(k.substr(5));
                return out;
            }

            /// Whether the pair of members match the contact rule at row 'r'.
            bool matches_(const DataWide &data_wide,
                          size_t r,
                          const Individual &indiv_from,
                          const Individual &indiv_to) const
            {
                for (size_t i = 0; i < from.size(); ++i)
                {
                    const std::string &fval = rules.input.at(from[i])[r];
                    if (!na::is(fval) && data_wide.data(keys[i], indiv_from) != str::tod(fval))
                        return false;

                    const std::string &tval = rules.input.at(to[i])[r];
                    if (!na::is(tval) && data_wide.data(keys[i], indiv_to) != str::tod(tval))
                        return false;
                }
                return true;
            };
        };

        /// Order covariate vectors of equal length.
        struct CovariateLess_
        {
            bool operator()(const std::vector<double> &a,
                            const std::vector<double> &b) const
            {
                for (size_t i = 0; i < a.size(); ++i)
                {
                    // a missing cell matches no rule value: it is a class of its own
                    const bool na_a = na::is(a[i]), na_b = na::is(b[i]);
                    if (na_a != na_b)
                        return na_b;
                    if (!na_a && a[i] != b[i])
                        return a[i] < b[i];
                }
                return false;
            }
        };

        /// A rule reads the covariates of the two members and nothing else (see
        /// ContactDefinitions_::matches_), so every member carrying the same
        /// covariates resolves the same way: class them, and a pair of classes
        /// resolves once for all the households repeating it.
        struct CovariateClasses_
        {
        public:
            idx_vec of_indiv; /// class every member belongs to
            size_t count = 0; /// distinct classes

            CovariateClasses_(const DataWide &data_wide, const idx_vec &keys)
                : of_indiv(data_wide.report_indiv.total, 0)
            {
                std::map<std::vector<double>, size_t, CovariateLess_> class_of;
                std::vector<double> covariates(keys.size());
                for (const Household &house : data_wide)
                    for (const Individual &indiv : house)
                    {
                        for (size_t i = 0; i < keys.size(); ++i)
                            covariates[i] = data_wide.data(keys[i], indiv);
                        of_indiv[Index::of(indiv)] =
                            class_of.emplace(covariates, class_of.size()).first->second;
                    }
                count = class_of.size();
            }
        };

        /// Risks an ordered pair of covariate classes multiplies, resolved on first encounter.
        struct PairRisk_
        {
            std::vector<const Attribute *> param_ptrs = {}; /// risks the pair applies
            bool resolved = false;                          /// whether the pair was resolved
        };
    }

    void ContactRisk::build(DataWide &data_wide,
                            const Parameter &param,
                            const str_table &input)
    {
        // reset storage
        reset_store_(data_wide.report_house.total);
        scale_anchors_.clear();
        if (input.empty())
            return;

        // load contact risk definitions
        ContactDefinitions_ defs(input, data_wide);

        data_wide.set_eq_class_data_names(defs.data);

        // per-rule_group scale anchoring, resolved pair by pair below
        std::map<size_t, ScaleAnchor> anchors;
        for (size_t r = 0; r < defs.rules.n_rows; ++r)
            anchors.emplace(defs.rules.row_gid[r], ScaleAnchor{});

        // resolved pair by pair, on first encounter only: ScaleAnchor::update is
        // idempotent, so the anchors accumulate exactly as one pass would leave them
        const CovariateClasses_ classes(data_wide, defs.keys);
        std::vector<PairRisk_> pair_risks(classes.count * classes.count);

        // lay every household's matrix row-major into the one flat store
        size_t n_cells = 0;
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
        {
            base_[H] = static_cast<RiskOffset>(n_cells);
            const size_t house_size = data_wide.at(H).size();
            n_cells += house_size * house_size;
        }
        flat_.assign(n_cells, risks::VALUE_REFERENCE);

        // loop over all households
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
        {
            Household &house = data_wide.at(H);
            size_t house_size = house.size();

            // NOTE: a priori contact risks are non-symmetric
            for (size_t from = 0; from < house_size; ++from)
                for (size_t to = 0; to < house_size; ++to)
                {
                    PairRisk_ &pair =
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

        // store the scale anchors of every rule group
        for (const auto &[gid, anchor] : anchors)
            scale_anchors_.push_back(anchor);
    }
}
