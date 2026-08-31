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
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/rule_table.h"
#include "tools/string.h"
#include "tools/types.h"

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
        helper_assert_non_degenerate_range_();
    }

    void Variable::helper_assert_non_degenerate_range_()
    {
        const Attribute &lower = pDistr_->lower();
        const Attribute &upper = pDistr_->upper();

        // skip if one of the range extremes is not provided
        if (!lower.is_callable() || !upper.is_callable())
            return;

        // skip if range is not degenerate
        const double sx = lower.value();
        const double dx = upper.value();
        if (sx < dx)
            return;

        // case empty interval, lower and upper are "open"
        if (sx > dx || !pDistr_->is_valid(sx) || !pDistr_->is_valid(dx))
            msg::stop(
                "Definition produced an inconsistent range " + pDistr_->str_range() +
                " for " + err_name() + ".");

        if (!lower.is_constant() || !upper.is_constant())
            return;

        if (is_user_const_)
            return;

        const std::string sms = "Definitions of " + msg::arg("min") + " as " + lower.err_name() +
                                " and " + msg::arg("max") + " as " + upper.err_name() +
                                " produced a degenerate range";
        // parameter
        if (err_type_ == ErrType::Parameter)
            msg::stop(sms + " for " + err_name() + ".");

        // augmented
        msg::warn("Collapsed (at least) one " + err_name() + " to constant.\n" + sms + ".");
        init_numeric(sx);
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

    void Variable::init_proposal(const std::string &proposal, const dict<Attribute> &args)
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
            msg::stop(sm + " for " + err_name() + ".");
        }
    }

    void Variable::helper_init_distr_proposal_()
    {
        pProp_ = std::make_unique<ForwardProposal>(pDistr_.get());
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    void Variable::helper_init_random_walk_proposal_(
        const std::string &proposal, double scale)
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

    void Variable::draw_value(rnd_gen &prng)
    {
        if (has_distr())
            set(pDistr_->rand(prng));
    }

    void Variable::finalize()
    {
        Attribute::finalize();
        if (pDistr_)
            pDistr_->finalize();
        if (pProp_)
            pProp_->finalize();
    }

    double Variable::operator()(double x) const
    {
        switch (mode_)
        {
        case Mode::Function:
            return pFunc_->evaluate(x);
        default: // cases Mode::Empty, Mode::Numeric, Mode::Pointer
            msg::stop(
                "Called undefined evaluate operator for non-function " + err_name() + ".");
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
            if (parser.consumes_augmented())
                msg::stop(
                    "Parameter function " + msg::val(pstr) +
                    " cannot depend on augmented variables.");
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

            // check if discrete as requested
            size_t p = khold(pstr);
            if (is_required_discrete(p) && !at(p).is_discrete())
                msg::stop(
                    "Definition of " + at(p).err_name() +
                    " did not return a discrete object as requested by model.");
        }

        // initialize all parameters to standard, i.e. none feeds transmission cache
        list_sample_standard = list_sample;
        for (size_t p : list_sample)
            if (at(p).is_adaptive())
                list_adaptive.push_back(p);
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        const std::vector<std::string> headers = read::extract_keys(input);

        clear_lists_();
        simulation_values.clear();
        // update list of constant parameters to all parameters
        list_const = khold(khold_.all());

        // reset definition dictionary to all constants NAs.
        for (const std::string &pstr : khold_.all())
            if (!contains(headers, pstr))
                simulation_values.emplace_back();
            else
                simulation_values.emplace_back(input.at(pstr));
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

    void Parameter::helper_build_constant_(const str_table &input, size_t k, const Translator &translator)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

        // assign constant parameter to the set of variables
        size_t p = khold(pstr);
        at(p).transfer(translator.to_numeric(input.at("value")[k]));

        // update the lists of constant parameters
        list_const.push_back(p);
    }

    void Parameter::helper_build_to_sample_(const str_table &input, size_t k, const Translator &translator)
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
        // load parameter range
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

    void Parameter::helper_build_function_(const str_table &input, size_t k, const Translator &translator)
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

        // range
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

    void Parameter::request_risk_identifiability(size_t p, std::vector<size_t> risks, bool contacts)
    {
        risk_anchors_.push_back({p, std::move(risks), contacts});
    }

    void Parameter::helper_pin_to_reference_(size_t p)
    {
        Variable &xparam = at(p);
        xparam.clear();
        xparam.init_numeric(risks::VALUE_REFERENCE);

        // update parameter sampling lists
        const std::vector<size_t> vec({p});
        list_sample = list_sample - vec;
        list_sample_standard = list_sample_standard - vec;
        list_sample_cached = list_sample_cached - vec;
        list_adaptive = list_adaptive - vec;
        list_const.push_back(p);
    }

    void Parameter::assert_risk_identifiability(const GroupRisk &group_risk, const ContactRisk &contact_risk)
    {
        for (const RiskAnchor &anchor : risk_anchors_)
        {
            size_t n_scale_free = 0;
            std::vector<std::string> free_names = {};

            if (anchor.contacts && contact_risk.n_scale_free() > 0)
            {
                n_scale_free += contact_risk.n_scale_free();
                free_names.emplace_back("contact");
            }

            for (size_t g : anchor.risks)
            {
                if (group_risk.n_scale_free(g) == 0)
                    continue;

                n_scale_free += group_risk.n_scale_free(g);
                free_names.emplace_back(group_risk.khold(g));
            }

            Variable &xparam = at(anchor.baseline);
            if (n_scale_free == 1 && !at(anchor.baseline).is_constant())
            {
                // UPDATE: when baseline is a function, it's a user choice
                // NOTE: at the moment is not even possible, but it might be in the future
                if (xparam.is_function())
                    continue;

                // if both baseline and group risk are free, pins the baseline
                msg::warn(
                    "No reference provided for " + msg::qty(free_names) + " risk{?s:} " + str::format(free_names) +
                    "\nBaseline " + xparam.err_name() + " was turned into constant reference.");
                helper_pin_to_reference_(anchor.baseline);
            }
            else if (n_scale_free > 1)
                msg::stop(
                    "The product of " + msg::qty(free_names) + " risk{?s} introduces non-identifiability: " + str::format(free_names) +
                    "\nRisks are absolute unless a " + msg::code("reference") + " is given: make all but one relative," +
                    " naming their reference group or asking for " + msg::val(".auto") + ".");
        }
    }

    void Parameter::mark_produces_transm_cache(size_t p)
    {
        // add the parameter as cache if not already present
        if (!contains(list_sample_cached, p))
            list_sample_cached.push_back(p);

        // if a parameter was originally defined as standard, remove it
        if (!contains(list_sample_standard, p))
            return;
        list_sample_standard = list_sample_standard - std::vector<size_t>({p});
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
    ///////////////////////////////////  *** GroupRisk : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        /// Scale anchoring of a single group (joint/independent sub-groups).
        struct ScaleAnchor_
        {
            bool anchored = false;   /// some sub-group has reference
            bool free_param = false; /// some winning parameter can move

            /// Register the parameter winning one cell. An empty name means no rule matched.
            void update(std::string_view pstr, const Parameter &param)
            {
                if (na::is(pstr) || pstr == risks::REFERENCE)
                    anchored = true;
                else if (!param.at(param.khold(pstr)).is_constant())
                    free_param = true;
            }

            /// Whether the overall scale of the rule group is left free.
            bool is_free() const noexcept { return !anchored && free_param; }
        };
    }

    void GroupRisk::build(DataWide &data_wide, const Parameter &param, const str_table &input)
    {
        // reset group risk storage and unbind individuals
        khold_.clear_optional();
        unique_groups_.clear();
        group_index_.assign(data_wide.report_indiv.total, 0);
        n_scale_free_.clear();

        // load required data headers or skip if no definitions provided
        std::vector<std::string> data_headers = {};
        if (!input.empty())
            data_headers = helper_process_definitions_(data_wide, input);
        else if (khold_.empty())
            return;
        n_scale_free_.assign(khold_.size(), 0);

        // load column with group risk names, if any
        std::vector<std::string> g_input_vec = {};
        if (input.count("group_risk"))
            g_input_vec = input.at("group_risk");

        // check for undefined group risks, assumed to be equal to reference
        str_table undefined_risks = {};
        for (const std::string &gstr : khold_.all())
            if (!contains(g_input_vec, gstr))
                undefined_risks[gstr].push_back(risks::REFERENCE);

        // precompute rule_group / specificity and resolve column indices once
        const RuleTable rules(input, data_headers);
        const auto data_ix = data_wide.khold(data_headers);

        // predicate: does rule row r match this member's feature vector?
        auto matches_with = [&](const Individual &mb)
        {
            return [&, mb_ptr = &mb](size_t r)
            {
                for (size_t d = 0; d < data_headers.size(); ++d)
                {
                    const auto &fv = input.at(data_headers[d])[r];
                    if (!na::is(fv) && data_wide.data(data_ix[d], *mb_ptr) != str::tod(fv))
                        return false;
                }
                return true;
            };
        };

        // per-(risk_name, rule_group) scale anchoring, resolved member by member below
        using Key = std::pair<std::string, size_t>;
        std::map<Key, ScaleAnchor_> anchors;
        for (size_t r = 0; r < rules.n_rows; ++r)
            anchors.emplace(Key{g_input_vec[r], rules.row_gid[r]}, ScaleAnchor_{});

        // unique definition map
        std::map<str_table, size_t> unique_defin_map;
        for (Household &house : data_wide)
            for (Individual &indiv : house)
            {
                str_table these_risks = undefined_risks;

                // per-(risk_name, rule_group): keep the most-specific matching rule
                const auto resolved = rules.resolve<Key>(
                    matches_with(indiv),
                    [&](size_t r)
                    { return Key{g_input_vec[r], rules.row_gid[r]}; },
                    [&](const Key &key)
                    {
                        return "group risk " + msg::val(key.first) + " for " +
                               data_wide.err_name(indiv) + " (group " +
                               std::to_string(key.second) + ")";
                    });

                for (auto &[key, r] : resolved)
                    these_risks[key.first].push_back(input.at("param_name")[r]);

                // a member matching no rule of the group reads the reference value
                for (auto &[key, anchor] : anchors)
                {
                    auto found = resolved.find(key);
                    anchor.update(found == resolved.end()
                                      ? std::string_view{}
                                      : std::string_view(input.at("param_name")[found->second]),
                                  param);
                }

                auto it = unique_defin_map.find(these_risks);
                if (it == unique_defin_map.end())
                {
                    it = unique_defin_map.emplace(these_risks, unique_groups_.size()).first;
                    helper_init_unique_risk_(data_wide, param, indiv, these_risks);
                }
                group_index_[Index::of(indiv)] = it->second;
            }

        // finalize the number of free scale parameters per group risk
        for (const auto &[key, anchor] : anchors)
            if (anchor.is_free())
                ++n_scale_free_[khold(key.first)];
    }

    void GroupRisk::finalize()
    {
        for (IndivGroup_ &group : unique_groups_)
            for (Attribute &risk : group.vec)
                risk.finalize();
    }

    std::vector<std::string> GroupRisk::helper_process_definitions_(DataWide &data_wide, const str_table &input)
    {
        const std::vector<std::string> file_headers = {"group_risk", "param_name", "rule_group"};

        // check that tempalte headers are present
        std::vector<std::string> data_headers = read::extract_keys(input);
        // remove template headers to keep only data
        data_headers = data_headers - file_headers;

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
        unique_groups_.reserve(max_unique);

        data_wide.set_eq_class_data_names(data_headers);

        return data_headers;
    }

    void GroupRisk::helper_init_unique_risk_(const DataWide &data_wide, const Parameter &param, const Individual &indiv, const str_table &these_risks)
    {
        unique_groups_.emplace_back();
        IndivGroup_ &igr = unique_groups_.back();
        igr.vec.reserve(khold_.size());
        for (const std::string &gstr : khold_.all())
        {
            if (!these_risks.count(gstr))
                msg::stop(
                    "A group risk definition " + msg::val(gstr) +
                    " was missing for " + data_wide.err_name(indiv) + ".");

            igr.vec.emplace_back(risks::VALUE_REFERENCE);
            igr.vec.back().name(gstr, ErrType::GroupRisk);
        }

        // assign custom definitions
        for (const auto &[gstr, pstr_vec] : these_risks)
        {
            // extract non-reference parameter pointers
            std::vector<const Attribute *> param_pointers = {};
            for (const std::string &pstr : pstr_vec)
                if (pstr != risks::REFERENCE)
                    param_pointers.emplace_back(
                        param.at(param.khold(pstr)).pointer());

            if (!param_pointers.empty())
                igr.vec[khold(gstr)].init_product(param_pointers);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** ContactRisk : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        using idx_vec = std::vector<size_t>;

        /// Concatenate 'from_X' and 'to_X' headers into the spec-column list.
        std::vector<std::string> concat_(
            const std::vector<std::string> &a,
            const std::vector<std::string> &b)
        {
            std::vector<std::string> out;
            out.reserve(a.size() + b.size());
            out.insert(out.end(), a.begin(), a.end());
            out.insert(out.end(), b.begin(), b.end());
            return out;
        }

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
            ContactDefinitions_(
                const str_table &input,
                const DataWide &data_wide)
                : from(extract_from_(input)),
                  to(extract_to_(input)),
                  data(strip_from_prefix_(from)),
                  keys(data_wide.khold(data)),
                  rules(input, concat_(from, to)) {}

            /// Associate parameter pointers to pair. Accumulate the scale anchoring of each group.
            std::vector<const Attribute *> associate_param_ptrs(
                const DataWide &data_wide,
                const Parameter &param,
                const Individual &indiv_from,
                const Individual &indiv_to,
                std::map<size_t, ScaleAnchor_> &anchors) const
            {
                // per-group: keep the most-specific matching rule
                const auto resolved = rules.resolve<size_t>(
                    [&](size_t r)
                    { return matches_(data_wide, r, indiv_from, indiv_to); },
                    [&](size_t r)
                    { return rules.row_gid[r]; },
                    [&](size_t)
                    {
                        return "contact risk from " + data_wide.err_name(indiv_from) +
                               " to " + data_wide.err_name(indiv_to);
                    });

                // across groups: multiply the winning params (skip reference)
                std::vector<const Attribute *> param_ptrs;
                param_ptrs.reserve(resolved.size());
                for (auto &[gid, r] : resolved)
                {
                    const auto &pstr = rules.input.at("param_name")[r];
                    if (!na::is(pstr) && pstr != risks::REFERENCE)
                        param_ptrs.push_back(param.at(param.khold(pstr)).pointer());
                }

                // a pair matching no rule of the group reads the reference value
                for (auto &[gid, anchor] : anchors)
                {
                    auto found = resolved.find(gid);
                    anchor.update(found == resolved.end()
                                      ? std::string_view{}
                                      : std::string_view(rules.input.at("param_name")[found->second]),
                                  param);
                }
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
            bool matches_(
                const DataWide &data_wide,
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
    }

    void ContactRisk::build(
        DataWide &data_wide,
        const Parameter &param,
        const str_table &input)
    {
        // reset storage
        house_contact_risks_.clear();
        n_scale_free_ = 0;
        if (input.empty())
            return;

        // load contact risk definitions
        ContactDefinitions_ defs(input, data_wide);

        data_wide.set_eq_class_data_names(defs.data);

        // per-rule_group scale anchoring, resolved pair by pair below
        std::map<size_t, ScaleAnchor_> anchors;
        for (size_t r = 0; r < defs.rules.n_rows; ++r)
            anchors.emplace(defs.rules.row_gid[r], ScaleAnchor_{});

        // loop over all households
        // UPDATE: this should be just over eq_class_head households.
        house_contact_risks_.resize(data_wide.report_house.total);
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
        {
            Household &house = data_wide.at(H);
            size_t house_size = house.size();

            house_contact_risks_[H].mtx.assign(
                house_size, std::vector<Attribute>(house_size, risks::VALUE_REFERENCE));

            // NOTE: a priori contact risks are non-symmetric
            for (size_t from = 0; from < house_size; ++from)
                for (size_t to = 0; to < house_size; ++to)
                {
                    std::vector<const Attribute *> param_ptrs =
                        defs.associate_param_ptrs(data_wide, param, house.at(from),
                                                  house.at(to), anchors);

                    if (!param_ptrs.empty())
                        house_contact_risks_[H].mtx[from][to].init_product(param_ptrs);
                }
        }

        // finalize the number of free scale parameters per contact risk
        for (const auto &[gid, anchor] : anchors)
            if (anchor.is_free())
                ++n_scale_free_;
    }

    void ContactRisk::finalize()
    {
        for (IndivContact_ &contact : house_contact_risks_)
            for (std::vector<Attribute> &row : contact.mtx)
                for (Attribute &risk : row)
                    risk.finalize();
    }
}
