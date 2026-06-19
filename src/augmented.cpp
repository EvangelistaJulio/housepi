/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dependency.h"
#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "augmented.h"
#include "mcmc_base.h"
#include "data_wide.h"
#include "distribution.h"
#include "function.h"
#include "household.h"
#include "individual.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"
#include "translator.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** InfecProfile : struct ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // UPDATE: this could be just forward declared maybe
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

    InfecProfile::InfecProfile() = default;
    InfecProfile::~InfecProfile() = default;

    void InfecProfile::init_template(Translator &translator, str_table input)
    {
        // if not previous build, proceed with template
        // if no template, assume no need for an infectivity profile
        if (generic || !augm_definitions.empty() || input.empty())
            return;

        try
        {
#ifndef NDEBUG
            if (!input.count("func"))
                throw std::invalid_argument(
                    "A function name was missing for infectivity profile template.");
#endif

            size_t n_rows = read::csv_n_rows(input);
            const auto file_headers = {"func", "status", "discr", "max"};

            for (const auto &key : file_headers)
                if (!input.count(key))
                    input.emplace(key, std::vector<std::string>(n_rows, ""));

            // discretization
            input["discr"] = std::vector<std::string>(n_rows, "false");
            if (is_discrete_requested_)
                input["discr"] = std::vector<std::string>(n_rows, "true");

            else if (input.count("discr"))
                for (size_t r = 0; r < n_rows; ++r)
                    if (!na::is(input.at("discr")[r]))
                        input["discr"][r] = input.at("discr")[r];

            build(translator, input); // run build as if external input

            // forbid augmented infectivity profile for template
            if (is_augmented())
                throw std::invalid_argument(
                    "Impossible to define infectivity profile template with augmented dependencies.");

            // UPDATE: maybe also not depending on parameters?
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While building infectivity profile template:\n" + std::string(e.what()));
        }
    }

    void InfecProfile::init_template(Translator &translator, dict<std::string> input)
    {
        str_table new_input = {};
        for (const auto &[key, val] : input)
            new_input[key] = {val};
        init_template(translator, std::move(new_input));
    }

    void InfecProfile::build(Translator &translator, const str_table &input)
    {
        augm_definitions.clear();

        if (input.empty())
            return;

#ifndef NDEBUG
        if (!input.count("func"))
            throw std::invalid_argument(
                "A function name choice was missing in infectivity profile definition.");
        for (const std::string &in_name : input.at("func"))
            if (parse::Definition(in_name).is_na())
                throw std::invalid_argument(
                    "Provided NA as a function name in infectivity profile definition.");
#endif

        // check if definition is generic or augmented
        if (read::csv_n_rows(input) > 1 ||
            parse::Definition(input.at("func")[0]).is_augmented_depen() ||
            parse::Definition(input.at("max")[0]).is_augmented_depen() ||
            !na::is(input.at("status")[0]))
            build_augmented_(translator, input);
        else
            build_generic_(translator, input);
    }

    void InfecProfile::build_generic_(Translator &translator, const str_table &input)
    {
        mode_ = Mode::Generic;
        generic = helper_build_(translator,
                                {{"func", {input.at("func")[0]}},
                                 {"discr", {input.count("discr") ? input.at("discr")[0] : "false"}}});
        lookup_.clear();
        augm_definitions.clear();

        // range minimum
        std::string min_expr = generic->is_discrete() ? "1|N" : "0|N"; // default
        if (input.count("min"))
            min_expr = input.at("min")[0];
        generic->init_attr("min", translator.to_attribute(min_expr));

        // range maximum
        if (input.count("max"))
            generic->init_attr("max", translator.to_attribute(input.at("max")[0]));

        // range truncation
        generic->init_truncation();

        // if discrete and generic, initializes lookup table
        if (generic->is_discrete())
        {
            mode_ = Mode::LookUp;
            size_t n = infprof::LOOKUP_SIZE;
            if (!generic->upper().is_empty())
                n = 1 + generic->upper().value();
            lookup_.init(generic.get(), n);
        }
    }

    void InfecProfile::build_augmented_(Translator &translator, const str_table &input)
    {
        // drop generic definition
        mode_ = Mode::Augmented;
        generic = nullptr;
        lookup_.clear();
        augm_definitions.clear();

        // use helper to check profile name and discretization
        for (size_t r = 0; r < read::csv_n_rows(input); ++r)
            helper_build_(translator, {{"func", {parse::Definition(input.at("func")[r]).root.func}},
                                       {"discr", {input.at("discr")[r]}}});

        // match infectivity profile and augmented variable definition
        for (size_t r = 0; r < read::csv_n_rows(input); ++r)
        {
            bool is_discrete = str::tobool(input.at("discr")[r]);

            augm_definitions["name"].push_back(infprof::AUGM_FUNC);
            augm_definitions["status"].push_back(input.at("status")[r]);
            augm_definitions["func"].push_back(input.at("func")[r]);
            augm_definitions["fcumul"].push_back("false");
            augm_definitions["value"].push_back("");
            augm_definitions["const"].push_back("false");
            augm_definitions["discr"].push_back(str::format(is_discrete));
            augm_definitions["prior"].push_back("");
            augm_definitions["min"].push_back(is_discrete ? "1|N" : "0|N");
            augm_definitions["max"].push_back(input.at("max")[r]);
            augm_definitions["proposal"].push_back("");
            augm_definitions["proposal_scale"].push_back("");
            augm_definitions["kind"].push_back("augm");
        }
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Translator &translator, const str_table &input)
    {
        parse::Definition parser(parse::Definition(input.at("func")[0]));
        auto output = compose_distr_registries(parser.root.func,
                                               "infectivity profile",
                                               {"infecprofile", "continuos", "discrete"});
        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, translator.to_attribute(node));
            translator.make_cache_dependent(node);
        }

        if (is_discrete_requested_ ||
            (input.count("discr") && str::tobool(input.at("discr")[0])))
            output->init_discretization();
        return output;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Augmented : class ***  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Augmented::Augmented(MCMCBase *mcmc)
        : mcmc_(mcmc) {}
    Augmented::~Augmented() = default;

    void Augmented::build(InfecProfile &infec_profile, const str_table &input)
    {
        // drop previous definitions
        is_any_sample_ = false;
        is_any_prior_ = false;
        is_any_obs_ = false;
        prior_param_depen_.clear();
        obs_param_depen_.clear();

        khold_.clear_optional();
        map_status_to_def_row_.clear();
        map_status_to_block_.clear();
        status_blocks_.clear();
        list_export.clear();
        list_obs_export.clear();
        list_house_export.clear();
        varia.clear();
        house_varia.clear();

        // skip if no augmented definitions
        if (input.empty() && !infec_profile.is_augmented())
            return;

        helper_build_assert_augm_input_(input);            // check entries
        helper_build_append_infec_profile_(infec_profile); // add infectivity profile
        helper_build_definition_by_status_();              // process definition by status
        helper_build_set_allocation_order_();              // allocation and order
        helper_init_observation_list_();                   // observation nodes
        helper_init_house_obs_list_();                         // household-level observations
    }

    /// Check augmented names and assert no reserved names are used.
    void Augmented::helper_build_assert_augm_input_(const str_table &input)
    {
        deftable_.clear();

#ifndef NDEBUG
        const std::vector<std::string> file_headers = {"name", "status",
                                                       "func", "fcumul",
                                                       "value", "const", "discr",
                                                       "prior", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = read::extract_keys(input);
#endif

        for (const std::string &astr : input.at("name"))
        {
#ifndef NDEBUG
            if (na::is(astr))
                throw std::invalid_argument(
                    "Found an augmented variable with undefined name.");
#endif
            if (astr == infprof::AUGM_FUNC)
                throw std::invalid_argument(
                    "Tried to use the reserved augmented name '" +
                    astr + "'.\nPlease, choose a different one.");
        }

        deftable_.append(input);
    }

    void Augmented::set_prior_param_depen_()
    {
        // initialize list of parameter dependent augmented/observation
        prior_param_depen_.assign(mcmc_->param.khold().size(), false);
        obs_param_depen_.assign(mcmc_->param.khold().size(), false);

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_name[k] == infprof::AUGM_FUNC ||
                deftable_.c_kind[k] == "house-obs")
                continue;

            std::vector<bool> &vec =
                deftable_.c_kind[k] == "obs" ? obs_param_depen_ : prior_param_depen_;
            for (const std::vector<std::string> *col : deftable_.source_cols())
                for (const std::string &pstr : parse::Definition((*col)[k]).list_parameter_depen())
                    vec[mcmc_->param.khold(pstr)] = true;
        }
    }

    void Augmented::DefTable_::clear()
    {
        for (auto col : {&c_name, &c_status, &c_func, &c_fcumul,
                         &c_value, &c_const, &c_discr, &c_prior,
                         &c_min, &c_max, &c_proposal, &c_scale, &c_kind})
            col->clear();
    }

    void Augmented::DefTable_::append(const str_table &input)
    {
        c_name += input.at("name");
        c_status += input.at("status");
        c_func += input.at("func");
        c_fcumul += input.at("fcumul");
        c_value += input.at("value");
        c_const += input.at("const");
        c_discr += input.at("discr");
        c_prior += input.at("prior");
        c_min += input.at("min");
        c_max += input.at("max");
        c_proposal += input.at("proposal");
        c_scale += input.at("proposal_scale");
        c_kind += input.at("kind");
    }

    void Augmented::helper_build_append_infec_profile_(InfecProfile &infec_profile)
    {
        if (!infec_profile.is_augmented())
            return;

        deftable_.append(infec_profile.augm_definitions);
        khold_.add(infprof::AUGM_FUNC);
        infec_profile.augm_key = khold(infprof::AUGM_FUNC);
        infec_profile.bind_augm(this); // for runtime per-member profile resolution

        // append discrete infectivity profile definitions if any
        if (infec_profile.is_discrete_requested())
            user_discrete_set_.emplace(infec_profile.augm_key);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  Augmented : assert univocal definitions  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_definition_by_status_()
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            // household-level observations are not per-member augmented variables:
            // keep them out of the key-holder and the status machinery (see
            // helper_init_house_obs_list_ / init_house_for_simulation).
            if (deftable_.c_kind[k] == "house-obs")
                continue;

            // fill registry with optional augmented definitions
            const std::string astr = deftable_.c_name[k];
            khold_.add(astr);

            // check associated status
            int status = status::MISSING;
            const std::string &s = deftable_.c_status[k];
            if (!na::is(s))
            {
                status = str::toi(s, "augmented status");
#ifndef NDEBUG
                if (!status::is_infected(status))
                    throw std::invalid_argument(
                        "Impossible to associate augmented '" + astr +
                        "' for non infectious status " + s + ".");
#endif
            }
            map_status_to_def_row_[status].push_back(k); // store status to definition row index
        }

        // assert that name and status are univocal
        for (const auto &[status, idxs] : map_status_to_def_row_)
        {
            std::vector<std::string> name_to_status = {};
            for (size_t k : idxs)
                name_to_status.push_back(deftable_.c_name[k]);

            const std::vector<std::string> degenerate = name_to_status - str::unique(name_to_status);
#ifndef NDEBUG
            if (!degenerate.empty())
                throw std::invalid_argument(
                    "Incompatible definitions found for augmented " + str::format(degenerate));
#endif
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : build definition order  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_set_allocation_order_()
    {
        std::unordered_set<int> all_infected_statuses = {};

        // allocate storage
        varia.resize(mcmc_->data_wide().report_indiv.total);

        size_t n = khold_.size();
        // UPDATE: loop on households it's just used for statuses: remove
        for (Household &house : mcmc_->data_wide())
            for (Individual &member : house)
            {
                IndivBlock &iaugm = varia[Index::of(member)];
                iaugm.vec.resize(n);
                for (size_t a = 0; a < n; ++a)
                    iaugm.vec[a].name(khold(a), "augmented");

                if (member.is_infected_status())
                    all_infected_statuses.emplace(member.status());
            }

        // first, load status specific user definitions
        // NOTE: NA (i. e. generic) does not necessarily need checking
        for (auto &[status, vec] : map_status_to_def_row_)
            if (status::is_infected(status))
                set_definition_dependency_order_for_(status, status);

        // associate observed status to corresponding status block
        // if not defined specifically, then associate to generic
        for (int status : all_infected_statuses)
            if (!map_status_to_block_.count(status))
            {
                if (!map_status_to_block_.count(status::MISSING))
                    set_definition_dependency_order_for_(status::MISSING, status);
                map_status_to_block_.emplace(status, map_status_to_block_.at(status::MISSING));
            }

        // if no observed infected status or status-specific definitions
        // then load generic definitition
        if (all_infected_statuses.empty() ||
            (map_status_to_def_row_.size() == 1 &&
             status::is_missing(map_status_to_def_row_.begin()->first)))
            set_definition_dependency_order_for_(status::MISSING, status::MISSING);
    }

    void Augmented::set_definition_dependency_order_for_(int status, int err_status)
    {
        // skip if already passed
        if (map_status_to_block_.count(status))
            return;
        // throw if not defined
        if (!map_status_to_def_row_.count(status))
            throw std::invalid_argument(
                "No definitions were compatible with augmented status " +
                std::to_string(status) + ".");

        str_table tree;           // populate augmented dependency tree
        dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
        for (size_t k : get_unsorted_definition_rows_(status))
        {
            std::vector<std::string> xdepen = {};
            for (const std::vector<std::string> *col : deftable_.source_cols())
                xdepen += parse::Definition((*col)[k]).list_augmented_depen();
            tree[deftable_.c_name[k]] = str::unique(xdepen);
            name_to_row[deftable_.c_name[k]] = k;
        }

        // check dependecies and translate order from name to deftable row index
        try
        {
            // NOTE: this method is just called once per status
            map_status_to_block_.emplace(status, status_blocks_.size());
            status_blocks_.emplace_back(status);

            // Ignore spurios private definition if no public definition depends on it
            // NOTE: but ensure private infectivity is stored.
            for (const std::string &astr : dependency::order_dependencies(tree))
                if (astr == infprof::AUGM_FUNC ||
                    !str::is_private(astr) ||
                    !str::drop_private(dependency::order_conditioned_on(astr, tree)).empty())
                    status_blocks_.back().order.emplace_back(name_to_row[astr]);

            // define usage lists respecting definition order
            set_sample_list_(status);
        }
        catch (const std::exception &e)
        {
            std::string sms = e.what();
            if (status::is_missing(err_status))
                sms += "\nCheck your generic definitions.";
            else
                sms += " for status " + std::to_string(err_status) + ".";
            throw std::invalid_argument(sms);
        }
    }

    std::vector<size_t> Augmented::get_unsorted_definition_rows_(int status) const
    {
        // load the list of status-specific definition rows
        std::vector<size_t> allrows = map_status_to_def_row_.at(status);

        // extract names
        std::vector<std::string> name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.c_name[k]);

        // add generic definitions to current list (if not generic already)
        // skip generic if a status-specific already exist
        if (status::is_infected(status) &&
            map_status_to_def_row_.count(status::MISSING))
            for (size_t k : map_status_to_def_row_.at(status::MISSING))
                if (!contains(name_to_status, deftable_.c_name[k]))
                    allrows.push_back(k);

        return allrows;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Augmented : sample lists  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::set_sample_list_(int status)
    {
        StatusBlock_ &block = status_blocks_[map_status_to_block_.at(status)];
        for (size_t k : block.order)
            set_sample_list_for_single_defininition_(status, k);

        // NOTE: by construction of draws list, I am forced to drop private at the end
        list_export = khold(str::drop_private(khold(list_export)));
        block.draws = khold(str::drop_private(khold(block.draws)));
    }

    void Augmented::set_sample_list_for_single_defininition_(int status, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        bool is_to_export = true;
        bool is_to_draw = true;

        if (str::tobool(deftable_.c_const[k]))
        {
            is_to_draw = false;
        }
        else if (!(na::is(deftable_.c_func[k]) && na::is(deftable_.c_value[k])))
        {
            auto [is_to_export_, is_to_draw_] = set_sample_list_for_function_(status, k);
            is_to_export = is_to_export_;
            is_to_draw = is_to_draw_;
        }
        else if (na::is(deftable_.c_value[k]) && na::is(deftable_.c_prior[k]))
        {
            throw std::invalid_argument(
                "Impossible to initialize augmented '" + astr + "'  without a value or prior");
        }
        else if (na::is(deftable_.c_proposal[k]))
        {
            throw std::runtime_error(
                "No proposal was provided for non-constant definition of augmented '" +
                astr + "'");
        }
        else
        {
            is_any_sample_ = true;
            status_blocks_[map_status_to_block_.at(status)].sample.emplace_back(a);

            // adaptive if proposal is parametric and scale is NA
            // NOTE: user defined proposal_scale freezes the scale
            if (deftable_.c_proposal[k] != "sample" &&
                na::is(deftable_.c_scale[k]))
            {
                is_any_adaptive_ = true;
                status_blocks_[map_status_to_block_.at(status)].adaptive.emplace_back(a);
            }
        }

        // NOTE: observation exports are listed separately
        if (deftable_.c_kind[k] == "obs")
            is_to_export = false;

        // nothing to do if not to export
        if (!is_to_export)
            return;

        if (!contains(list_export, a))
            list_export.emplace_back(a);

        if (is_to_draw)
            status_blocks_[map_status_to_block_.at(status)].draws.emplace_back(a);
    }

    std::pair<bool, bool> Augmented::set_sample_list_for_function_(int status, size_t k)
    {
        bool is_to_export = true;
        bool is_to_draw = true;

        // NOTE: value definition has priority over function
        const std::string &str_val = deftable_.c_value[k];
        const parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

        mcmc_->data_wide().set_eq_class_data_names(parser.list_data_wide_depen());

        // NOTE: do not export if not callable (otherwise crashes)
        if (!parser.is_callable_implicitly("x"))
            is_to_export = false;

        bool all_consts = true; // start assuming all constants
        for (const std::string &pstr : parser.list_parameter_depen())
        {
            size_t p = mcmc_->param.khold().from_string(pstr);
            all_consts &= mcmc_->param.at(p).is_constant();
        }

        // inspect status specific augmented definitions and check if to draw (i. e. not const)
        // NOTE: taking advatange of the order!
        for (const std::string &bstr : parser.list_augmented_depen())
            if (contains(status_blocks_[map_status_to_block_.at(status)].draws, khold(bstr)))
                all_consts = false;

        // to draw if not all constant dependencies
        is_to_draw = !all_consts;

        return {is_to_export, is_to_draw};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : initialize augmented  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_for_evaluation()
    {
        if (khold_.empty())
            return;

        // raise if deterministic observation node for evaluation
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == "obs" && na::is(deftable_.c_prior[k]))
                throw std::invalid_argument(
                    "Deterministic observation node '" + deftable_.c_name[k] +
                    "' cannot be used for evaluation, only simulation.");

        helper_init_prior_list_();

        // define augmented variables for all augmented individuals in ordefr
        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                for (size_t k : helper_materialize_definition_order_(member))
                    helper_init_attributes_(member, k);
            }

        // parameters that may move the augmented log-prior (params now registered)
        set_prior_param_depen_();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init augmented priors  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_prior_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.depen.resize(khold_.size());

            // associate blocks for ordered prior lists
            const str_table &ordered_priors = get_ordered_priors_(block);
            if (!ordered_priors.empty())
            {
                is_any_prior_ = true;
                block.priors = khold(read::extract_keys(ordered_priors));

                for (const auto &[astr, others] : ordered_priors)
                    block.depen[khold(astr)] = khold(others);
            }
        }
    }

    str_table Augmented::get_ordered_priors_(const StatusBlock_ &block) const
    {
        str_table tree;                         // populate augmented dependency tree
        str_table ordered_priors;               // store hierarchical order of priors
        std::vector<std::string> without_prior; // list of augmented without priors
        for (size_t k : block.order)
        {
            const std::string &aumgstr = deftable_.c_name[k];
            std::vector<std::string> xdepen = {};
            for (const std::vector<std::string> *col : deftable_.source_cols())
                xdepen += parse::Definition((*col)[k]).list_augmented_depen();

            // NOTE: observations hold distributions in c_prior but are NOT priors
            if (na::is(deftable_.c_prior[k]) || deftable_.c_kind[k] == "obs")
                without_prior.push_back(aumgstr);
            else
                tree[aumgstr] = str::unique(xdepen);
        }

        // prune dependency tree from augmented without prior
        for (auto &[key, vec] : tree)
            vec = vec - without_prior;

        // assign order of dependents for augmented prior computation
        for (const auto &[key, vec] : tree)
            ordered_priors.emplace(key, dependency::order_conditioned_on(key, tree));

        return ordered_priors;
    }

    void Augmented::helper_init_observation_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.obs_depen.resize(khold_.size());

            for (size_t k : block.order)
            {
                if (deftable_.c_kind[k] != "obs")
                    continue;

                const size_t obs_key = khold(deftable_.c_name[k]);

                // unique observation keys to return as simulated columns
                if (!contains(list_obs_export, obs_key))
                    list_obs_export.push_back(obs_key);

                // a deterministic observation node carries no distribution
                if (na::is(deftable_.c_prior[k]))
                    continue;

                is_any_obs_ = true;
                block.obs.push_back(obs_key);

                // store observations that depend on augmented
                std::vector<std::string> depen = {};
                for (const std::vector<std::string> *col : deftable_.source_cols())
                    depen += parse::Definition((*col)[k]).list_augmented_depen();
                for (const std::string &lstr : str::unique(depen))
                    block.obs_depen[khold(lstr)].push_back(obs_key);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : household observations  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_house_obs_list_()
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == "house-obs")
                list_house_export.push_back(k);
    }

    void Augmented::init_house_for_simulation()
    {
        house_varia.clear();
        if (list_house_export.empty())
            return;

        const size_t n_houses = mcmc_->data_wide().report_house.total;
        house_varia.resize(n_houses);
        for (size_t H = 0; H < n_houses; ++H)
        {
            Household &house = mcmc_->data_wide().at(H);
            house_varia[H].resize(list_house_export.size());
            for (size_t b = 0; b < list_house_export.size(); ++b)
                helper_init_house_attributes_(house, b);
        }
    }

    void Augmented::helper_init_house_attributes_(Household &house, size_t j)
    {
        const size_t k = list_house_export[j];
        const std::string &astr = deftable_.c_name[k];
        Variable &hvar = house_varia[Index::of(house)][j];
        hvar.name(astr, "household observation");

        try
        {
            const std::string &str_val = deftable_.c_value[k];
            parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

            // bare reference (no reduction): bind to the first member, which for a
            // singleton household is the household value itself.
            if (parser.is_terminal())
            {
                hvar.transfer(mcmc_->translator.to_attribute(parser.root, &house.member(0)));
                return;
            }

            // reduction: compose the (variadic) function once
            hvar.init_function(compose_func_registries(
                parser.root.func, "function",
                {"custom", "transcendental", "arithmetic",
                 "continuos", "discrete", "infecprofile"},
                str::tobool(deftable_.c_fcumul[k])));

            // bind operands: a data (@D) source is spread one operand per member,
            // skipping NA values (na.rm at build time, since NA constant args are
            // rejected by the function); any other source is bound once at household
            // level. The first operand seeds the reduction ('x'), the rest are
            // appended arguments.
            size_t op = 0;
            const auto append_operand = [&](const Attribute &attr)
            {
                const std::string key = op == 0 ? "x" : "x" + std::to_string(op);
                ++op;
                hvar.init_func_attr(key, attr);
            };
            for (const auto &child : parser.root.children)
                if (parse::is_data_wide(child.type))
                    for (const Individual &member : house)
                    {
                        const Attribute &attr = mcmc_->translator.to_attribute(child, &member);
                        if (!na::is(attr.value()))
                            append_operand(attr);
                    }
                else
                    append_operand(mcmc_->translator.to_attribute(child, &house.member(0)));

            // every operand was NA: the reduction is empty, yield NA (na.rm convention)
            if (op == 0)
            {
                hvar.init_numeric(na::DOUBLE);
                return;
            }

            // range and reshaping
            hvar.init_func_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &house.member(0)));
            hvar.init_func_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &house.member(0)));
            hvar.init_func_reshaping(str::tobool(deftable_.c_discr[k]));
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While initializing household observation '" + astr + "' for " +
                house.err_name() + ":\n" + std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init missing attribute  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::clear_all()
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
                clear_definitions(house.member(i));
    }

    void Augmented::clear_definitions(Individual &member)
    {
        for (Variable &avar : varia[Index::of(member)].vec)
            avar.clear();
    }

    void Augmented::init_definitions_for_simulation(Individual &member)
    {
        for (size_t k : helper_materialize_definition_order_(member))
            if (varia[Index::of(member)].vec[khold(deftable_.c_name[k])].is_empty())
                helper_init_attributes_(member, k, true);
    }

    std::vector<size_t> Augmented::helper_materialize_definition_order_(Individual &member)
    {
        int status = member.status();

        // if known status, then load its definition
        const auto it = map_status_to_block_.find(status);
        if (it != map_status_to_block_.end())
        {
            const size_t i = it->second;
            varia[Index::of(member)].idx = i;
            return status_blocks_[i].order;
        }

        // otherwise if it exists fallback to generic definition
        const auto generic_it = map_status_to_block_.find(status::MISSING);
        if (generic_it != map_status_to_block_.end())
        {
            const size_t i = generic_it->second;
            map_status_to_block_.emplace(status, i);
            varia[Index::of(member)].idx = i;
            return status_blocks_[i].order;
        }

        if (member.is_infected_status())
            throw std::invalid_argument(
                "No augmented definitions were found for " +
                member.err_name() + " with status " + std::to_string(status) + ".");

        throw std::invalid_argument(
            "Generic augmented definition for non-infectious " +
            member.err_name() + " was missing.");
    }

    void Augmented::helper_init_attributes_(Individual &member, size_t k, bool simulation)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        try
        {
            Variable &avar = varia[Index::of(member)].vec[a];

            if (str::tobool(deftable_.c_const[k]))
            {
                // placeholder for observation node if generated by simulation
                if (simulation && deftable_.c_kind[k] == "obs")
                    avar.init_numeric(0.);
                else
                    avar.init_numeric(
                        mcmc_->translator.to_numeric(deftable_.c_value[k], &member).value());
                avar.enforce_user_const(); // enforce user constant choice
            }
            else
                helper_init_func_attributes_(member, k);

            // initialize prior
            helper_init_prior_attributes_(member, k);

            // assign proposal if variable to be sampled
            if (!avar.is_constant())
                avar.init_proposal(deftable_.c_proposal[k], str::tod(deftable_.c_scale[k]));

            // check if satisfies discrete requirement, if any
            if (user_discrete_set_.count(a) && !avar.is_discrete())
                throw std::invalid_argument(
                    "Definition of " + avar.err_name() +
                    " did not return a discrete object as requested by prior.");
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While initializing augmented '" + astr + "' for " +
                member.err_name() + " (status " + std::to_string(member.status()) + "):\n" +
                std::string(e.what()));
        }
    }

    void Augmented::helper_init_prior_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        Variable &avar = varia[Index::of(member)].vec[khold(astr)];
        parse::Definition parser(deftable_.c_prior[k]);

        avar.init_prior(parser.root.func);

        if (!avar.has_prior())
            return;
        for (const auto &child : parser.root.children)
            avar.init_prior_attr(child.key, mcmc_->translator.to_attribute(child, &member));
        avar.init_prior_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &member));
        avar.init_prior_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &member));
        avar.init_prior_reshaping(str::tobool(deftable_.c_discr[k]));
    }

    void Augmented::helper_init_func_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        Variable &avar = varia[Index::of(member)].vec[khold(astr)];
        std::string str_val = deftable_.c_value[k];
        parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

        // pointer mode
        if (parser.is_terminal())
        {
            avar.transfer(mcmc_->translator.to_attribute(parser.root, &member));
            return;
        }

        // function mode
        bool is_primitive = str::tobool(deftable_.c_fcumul[k]);
        avar.init_function(compose_func_registries(
            parser.root.func, "function",
            {"custom", "transcendental", "arithmetic",
             "continuos", "discrete", "infecprofile"},
            is_primitive));

        // attribute
        for (const auto &child : parser.root.children)
        {
            avar.init_func_attr(child.key, mcmc_->translator.to_attribute(child, &member));
            mcmc_->translator.make_cache_dependent(child);
        }

        // range
        avar.init_func_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &member));
        avar.init_func_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &member));

        // discretization and truncation
        avar.init_func_reshaping(str::tobool(deftable_.c_discr[k]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::restart(rnd_gen &prng, Individual &member)
    {
        for (size_t a : list_sample(member))
            varia[Index::of(member)].vec[a].restart(prng);
    }

    void Augmented::sample_obs(rnd_gen &prng, Individual &member)
    {
        for (size_t a : list_obs(member))
            varia[Index::of(member)].vec[a].draw_value(prng);
    }

    void Augmented::restart(rnd_gen &prng)
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
                restart(prng, house.member(i));
    }

    void Augmented::optimize()
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                for (Variable &avar : varia[Index::of(member)].vec)
                    avar.optimize();
            }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    const std::vector<size_t> &Augmented::list_sample(const Individual &member) const
    {
        return status_blocks_[varia[Index::of(member)].idx].sample;
    }
    const std::vector<size_t> &Augmented::list_adaptive(const Individual &member) const
    {
        return status_blocks_[varia[Index::of(member)].idx].adaptive;
    }
    const std::vector<size_t> &Augmented::list_draws(const Individual &member) const
    {
        return status_blocks_[varia[Index::of(member)].idx].draws;
    }
    const std::vector<size_t> &Augmented::list_priors(const Individual &member) const
    {
        return status_blocks_[varia[Index::of(member)].idx].priors;
    }
    const std::vector<size_t> &Augmented::list_obs(const Individual &member) const
    {
        return status_blocks_[varia[Index::of(member)].idx].obs;
    }

    std::string format_exception_for(Individual &member, const std::exception &e)
    {
        std::string sms = e.what();
        if (!sms.empty() && sms.back() == '.')
            sms.pop_back();
        return sms + " for " + member.err_name() + ".";
    }

    double Augmented::commit_augm_log_prior(Individual &member)
    {
        try
        {
            const std::vector<size_t> &priors = list_priors(member);

            if (priors.empty())
                return 0.;
            double new_ll = 0.0;
            for (size_t a : priors)
                new_ll += varia[Index::of(member)].vec[a].get_log_prior();
            varia[Index::of(member)].log_prior.set(new_ll);
            return new_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    double Augmented::commit_augm_log_prior_ratio(Individual &member, size_t a)
    {
        try
        {
            const StatusBlock_ &block = status_blocks_[varia[Index::of(member)].idx];
            if (block.depen[a].empty())
                return 0.;
            double delta_ll = 0;
            // NOTE: vector of dependencies 'depen' includes also 'a'
            for (size_t b : block.depen[a])
                delta_ll += varia[Index::of(member)].vec[b].get_prior_log_ratio();
            varia[Index::of(member)].log_prior.increment(delta_ll);
            return delta_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    double Augmented::commit_log_lik_obs(Individual &member)
    {
        try
        {
            const std::vector<size_t> &obs = list_obs(member);

            if (obs.empty())
                return 0.;
            double new_ll = 0.0;
            for (size_t a : obs)
                new_ll += varia[Index::of(member)].vec[a].get_log_prior();
            varia[Index::of(member)].log_obs.set(new_ll);
            return new_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    double Augmented::commit_log_lik_obs_ratio(Individual &member, size_t a)
    {
        try
        {
            const StatusBlock_ &block = status_blocks_[varia[Index::of(member)].idx];
            // skip unless the moved augmented 'a' parameterizes some observation node
            if (block.obs_depen[a].empty())
                return 0.;

            // log-proba must be recomputed (observation value is just a datum)
            double new_ll = 0.0;
            for (size_t b : block.obs)
                new_ll += varia[Index::of(member)].vec[b].get_log_prior();
            const double delta_ll = new_ll - varia[Index::of(member)].log_obs.value();
            varia[Index::of(member)].log_obs.set(new_ll);
            return delta_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Augmented : sample draws  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<const Individual *> Augmented::init_individual_sample(
        const std::vector<std::string> &id_draws) const
    {
        std::vector<const Individual *> output;

        if (id_draws.empty())
            return output;

        if (!is_any_sample())
        {
            msg::warning(
                "Ignored id draws. No augmented variables to sample were defined.");
            return output;
        }

        output.reserve(id_draws.size());
        std::vector<std::string> existing_ids = {};
        for (const Household &house : mcmc_->data_wide())
            for (const Individual &member : house)
            {
                if (!contains(id_draws, member.id()))
                    continue;

                // unique IDs are granted by class DataWide load
                existing_ids.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    output.emplace_back(&member);
                else
                    msg::warning(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }
        output.shrink_to_fit();

        return output;
    }
}
