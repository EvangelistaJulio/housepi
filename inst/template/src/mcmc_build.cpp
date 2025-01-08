/*
 * mcmc_build.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// declaration for type in data file
const std::string TYPE_PARAM = "P";   // Parameter type declaration
const std::string TYPE_DATA = "D";    // Data type declaration
const std::string TYPE_SCALAR = "S";  // Scalar type declaration
const std::string TYPE_COMMAND = "C"; // Command type declaration

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//===== CONSTRUCTOR =====//
Monte_Carlo_Build::Monte_Carlo_Build() : _uniform_distr_(0., 1.)
{
#ifdef _OPENMP
    _max_n_available_threads = omp_get_max_threads();
#endif
}
Monte_Carlo_Build::~Monte_Carlo_Build() = default;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Auxiliary function to call percentuage. */
inline const std::string _aux_percent_(int num, int den)
{
    return "(" + tools::percentage(static_cast<double>(num) / den, 1) + ")";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Monte_Carlo_Build : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double Monte_Carlo_Build::type_to_double(const std::string &value,
                                         const std::string &type,
                                         const Household &house, indiv i)
{
    double output = UNDEFINED_DOUBLE;

    try
    {
        if (!(tools::is_na(value) && tools::is_na(type)))
        {
            if (!tools::is_na(value) && tools::is_na(type))
                output = tools::stod(value, "argument without type");

            else if (tools::is_na(value) && !tools::is_na(type))
                tools::warning(
                    "Ignoring value type without a paired value.");

            else if (type == TYPE_SCALAR)
                output = tools::stod(value, "argument with scalar type");

            else if (type == TYPE_DATA)
                output = house.data(i, value);

            else if (type == TYPE_PARAM)
                // UPDATE: here it won't be straightforward to drop maps
                if (param_map.at(value).is_constant())
                    output = get_param(value);
                else
                    throw std::invalid_argument(
                        "Parameter to be inferred '" + value + "' assigned as constant argument.");                    

            else if (type == TYPE_COMMAND)
                if (value == "default")
                    output = UNDEFINED_DOUBLE;
                else
                    throw std::invalid_argument(
                        "Unknown command type '" + value + "'.");
            else
                throw std::invalid_argument(
                    "Unknown variable type specifier '" + type + "'.");
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void Monte_Carlo_Build::assign_augm_type_args(Household &house, indiv i, augmID a,
                                              const std::string &tkn_key,
                                              const std::string &tkn_val,
                                              const std::string &tkn_typ)
{

    // parse the string with arguments
    str_vec keys = tools::parse_list_args(tkn_key);
    str_vec vals = tools::parse_list_args(tkn_val);
    str_vec typs = tools::parse_list_args(tkn_typ);

    // check all's good
    if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
        throw std::invalid_argument(
            "Mismatching definitions in arguments.");

    // loop over all arguments
    for (size_t x = 0; x < keys.size(); ++x)
    {
        std::string key = keys[x];
        std::string value = vals[x];
        std::string type = typs[x];

        if (!(tools::is_na(value) && tools::is_na(type)))
        {
            if (!tools::is_na(value) && tools::is_na(type))
                throw std::invalid_argument(
                    "Impossible to interpret default value of  '" + value + "'. Please provide a type.");

            else if (tools::is_na(value) && !tools::is_na(type))
                tools::warning(
                    "Found a type without a paired value. Falling back to default.");

            else if (type == TYPE_DATA)
                house.assign_augm_model_scalar_arg(i, a, key, house.data(i, value));

            else if (type == TYPE_PARAM) // UPDATE: we an improve for constant parameters
                house.assign_augm_model_pointer_arg(i, a, key, get_param_pointer(value));

            else if (type == TYPE_SCALAR)
                house.assign_augm_model_scalar_arg(i, a, key, tools::stod(value, "argument with scalar type"));

            else
                throw std::invalid_argument(
                    "Unsopported type '" + type + "' found while loading '" + value + "' for pointer value.");
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Emplace a new household. */
void Monte_Carlo_Build::assign_house(const std::string &ID)
{
    house_vector.emplace_back(ID, n_households);
    ++n_households;
}

/* Use the loaded data to fill-in the respective report. */
void Monte_Carlo_Build::update_report_data()
{
    _report_data_ = {{"n.cases", 0},
                     {"n.miss.stat", 0},
                     {"no.HH.case", 0},
                     {"lst.1.case", 0},
                     {"lst.2.case", 0},
                     {"max.n.HH.cases", 0},
                     {"n.HH.with.max.cases", 0},
                     {"min.HH.size", 1000},
                     {"n.HH.with.min.size", 0},
                     {"max.HH.size", 0},
                     {"n.HH.with.max.size", 0}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (const Household &house : house_vector)
    {
        ushort house_size = house.size();
        ushort n_infec = house.n_infectious();

        _report_data_["n.cases"] += n_infec;
        _report_data_["n.miss.stat"] += house.n_missing_status();

        _report_data_["no.HH.case"] += (n_infec == 0);
        _report_data_["lst.1.case"] += (n_infec > 0);
        _report_data_["lst.2.case"] += (n_infec > 1);

        // max number of house cases
        if (_report_data_["max.n.HH.cases"] < n_infec)
        {
            _report_data_["n.HH.with.max.cases"] = 1;
            _report_data_["max.n.HH.cases"] = n_infec;
        }
        else if (_report_data_["max.n.HH.cases"] == n_infec)
        {
            _report_data_["n.HH.with.max.cases"] += 1;
        }

        // minimum house number of members
        if (_report_data_["min.HH.size"] > house_size)
        {
            _report_data_["n.HH.with.min.size"] = 1;
            _report_data_["min.HH.size"] = house_size;
        }
        else if (_report_data_["min.HH.size"] == house_size)
        {
            _report_data_["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (_report_data_["max.HH.size"] < house_size)
        {
            _report_data_["n.HH.with.max.size"] = 1;
            _report_data_["max.HH.size"] = house_size;
        }
        else if (_report_data_["max.HH.size"] == house_size)
        {
            _report_data_["n.HH.with.max.size"] += 1;
        }
    }
}

/* Return a value from report data. */
size_t Monte_Carlo_Build::report_data(const std::string &key) const
{
    size_t output = 0;
    if (_report_data_.find(key) != _report_data_.end())
        output = _report_data_.at(key);

    else
        throw std::invalid_argument(
            "Unknown summary entry '" + key + "'.");

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Monte_Carlo_Build : build augmented  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void Monte_Carlo_Build::read_augm(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_table(pathtofile, _augm_deftable_, augm_file_headers);
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While reading augmented definitions: " + std::string(e.what()));
    }

    // use loaded definitions to build
    build_augm();
    // store the file name in the report
    _report_pathtofile_["augm"] = pathtofile;
}

/* Assign augmented variables to all individuals. */
void Monte_Carlo_Build::build_augm()
{
    try
    {
        if (!_augm_deftable_.empty())
        {
            // n. of definitions
            size_t table_length = _augm_deftable_.at("name").size();

            // auxiliary variable to load definitions : a -> status lists
            dict_str_vec name_status_check = {};

            //===== Check Augmented Variables lists of usage =====//
            for (size_t k = 0; k < table_length; ++k)
            {
                std::string a = _augm_deftable_.at("name")[k];
                std::string status = _augm_deftable_.at("status")[k];
                bool is_init_rand = false;

                // check that augmID and status are correct
                if (tools::is_na(a))
                    throw std::invalid_argument(
                        "Augmented variable name can not be undefined.");

                else if (name_status_check.find(a) == name_status_check.end())
                    // new augmented variable status
                    name_status_check.emplace(a, str_vec{status});

                else
                {
                    // check no multiple augmID with same status
                    const str_vec &mylist = name_status_check.at(a);
                    if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
                        name_status_check.at(a).push_back(status);

                    else
                        throw std::invalid_argument(
                            "Incompatible definitions for augmented variable '" + a + "'.");
                }

                // option for random initialization
                if (!tools::is_na(_augm_deftable_.at("init_rand")[k]))
                    is_init_rand = tools::stobool(_augm_deftable_.at("init_rand")[k]);

                if ((tools::is_na(_augm_deftable_.at("model")[k])) && (is_init_rand))
                    throw std::invalid_argument(
                        "Impossible to initialize '" + a + "' without a model.");
            }

            //===== Drop previous definitions =====//
            augm_list_full.clear();
            augm_list_infer.clear();
            augm_list_model.clear();
            augm_status_map_infer.clear();
            augm_status_map_model.clear();
            is_augm_infer = false;
            is_augm_model = false;

            //===== Check Univocal Definitions and Allocate New Lists =====//
            for (const auto &pair : name_status_check)
            {
                const std::string &a = pair.first;
                const str_vec &mylist = pair.second;

                for (const std::string &status : mylist)
                {
                    if (tools::is_na(status))
                    {
                        if (mylist.size() > 1)
                            throw std::invalid_argument(
                                "Conflicting status-specific and general definitions for augmented '" + a + "'.");
                    }
                    else
                    {
                        int s = tools::stoi(status, "augmented status");
                        if (s < INFECTIOUS)
                        {
                            throw std::invalid_argument(
                                "Impossible to associate augmented '" + a + "' for non infectious status.");
                        }
                        else
                        {
                            augm_status_map_infer.emplace(s, std::unordered_set<std::string>{});
                            augm_status_map_model.emplace(s, std::unordered_set<std::string>{});
                        }
                    }

                    augm_list_full.emplace(a);
                }
            }

            //===== Case of Only Generic Definitions =====//
            bool is_only_generic = (augm_list_full.size() > 0) && (augm_status_map_infer.empty());
            if (is_only_generic)
                _func_augm_model_by_definition_ = &Monte_Carlo_Build::_augm_model_generic_;

            //===== Assign new definitions to lists =====//
            // note : I can't I collapse this with the previous because of generic status
            for (size_t k = 0; k < table_length; ++k)
            {
                std::string a = _augm_deftable_.at("name")[k];
                std::string status = _augm_deftable_.at("status")[k];

                //----- Augmented Variables to be Inferred -----//
                if (!tools::is_na(_augm_deftable_.at("rule")[k]))
                {
                    is_augm_infer = true;

                    if (tools::is_na(status)) // generic definition
                        if (is_only_generic)
                            augm_list_infer.emplace(a);

                        else
                            for (const auto &pair : augm_status_map_infer)
                                augm_status_map_infer.at(pair.first).emplace(a);

                    else // status-specific definition
                        augm_status_map_infer.at(tools::stoi(status)).emplace(a);
                }

                //----- Augmented Variables with a Model -----//
                if (!tools::is_na(_augm_deftable_.at("model")[k]))
                {
                    is_augm_model = true;

                    if (tools::is_na(status)) // generic definition
                        if (is_only_generic)
                            augm_list_model.emplace(a);

                        else
                            for (const auto &pair : augm_status_map_model)
                                augm_status_map_model.at(pair.first).emplace(a);

                    else // status-specific definition
                        augm_status_map_model.at(tools::stoi(status)).emplace(a);
                }
            }

            //===== Allocate memory for augmented variables to each individual =====//
            for (Household &house : house_vector)
                for (indiv i = 0; i < house.size(); ++i)
                    clear_all_augm(house, i);

            //===== Functional method for augmented model contribution to loglikelihood ====//
            if (is_augm_model)
                _func_augm_model_loglike_ = &Monte_Carlo_Build::_augm_model_loglike_;
            else
                _func_augm_model_loglike_ = &Monte_Carlo_Build::_no_augm_model_loglike_;
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While building augmented variable definitions: " + std::string(e.what()));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Monte_Carlo_Build : assign augmented variables  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void Monte_Carlo_Build::assign_all_augm(Household &house, indiv i)
{
    try
    {
        for (size_t k = 0; k < _augm_deftable_.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = _augm_deftable_.at("status")[k];

            // assign augmented variable of the corresponding status or undefined status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                //----- Index -----//
                // |====> UPDATE: to be implemented <===|
                augmID a = _augm_deftable_.at("name")[k];

                //----- Name -----//
                const std::string &name = _augm_deftable_.at("name")[k];

                //----- is Discrete -----//
                bool is_discrete = tools::stobool(_augm_deftable_.at("discrete")[k]);

                //----- Value -----//
                double value = type_to_double(_augm_deftable_.at("value")[k],
                                              _augm_deftable_.at("value_type")[k],
                                              house, i);

                //----- Model -----//
                const std::string &model_name = _augm_deftable_.at("model")[k];

                //----- Range -----//
                // minimum value
                double lower = type_to_double(_augm_deftable_.at("min")[k],
                                              _augm_deftable_.at("min_type")[k],
                                              house, i);

                // maximum value
                double upper = type_to_double(_augm_deftable_.at("max")[k],
                                              _augm_deftable_.at("max_type")[k],
                                              house, i);

                //----- Proposer Rule -----//
                const std::string &rule_name = _augm_deftable_.at("rule")[k];
                double rule_scale = type_to_double(_augm_deftable_.at("rule_scale")[k],
                                                   _augm_deftable_.at("rule_scale_type")[k],
                                                   house, i);
                double rule_from = type_to_double(_augm_deftable_.at("rule_from")[k],
                                                  _augm_deftable_.at("rule_from_type")[k],
                                                  house, i);

                //----- Assign Augmented Variable -----//
                house.assign_augm(i, a, name, value,
                                  is_discrete, model_name,
                                  rule_name, rule_scale, rule_from);

                assign_augm_type_args(house, i, a,
                                      _augm_deftable_.at("arg_name")[k],
                                      _augm_deftable_.at("arg_value")[k],
                                      _augm_deftable_.at("arg_type")[k]);

                house.assign_augm_model_range(i, a, lower, upper);

                // random initialization at assignment
                if (tools::stobool(_augm_deftable_.at("init_rand")[k]))
                    house.augm_model_rand(i, name, get_prng());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While initializing augmented variables: " + std::string(e.what()));
    }
}

/* Clear all augmented definition of individual. */
void Monte_Carlo_Build::clear_all_augm(Household &house, indiv i)
{
    for (const augmID &a : augm_list_full)
        house.clear_augm(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void Monte_Carlo_Build::init_diagnostics(bool is_power_post)
{
    ulong sample_size = setup_c_sample_size * setup_n_chains;

    dgn_count = 0;
    dgn_glob_loglike = dbl_vec(sample_size);

    if (!is_power_post)
    {
        //----- Initialize diagnostic vectors -----//
        // FLAG: delicate point
        for (const parID &p : param_to_be_inferred)
            dgn_param_value[p] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());

        // diagnostics over households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].init_diagnostics(sample_size);
    }
}

/* Update Bayesian diagnostics */
void Monte_Carlo_Build::update_diagnostics(double glob_loglike, bool is_power_post)
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count] = glob_loglike;

    if (!is_power_post)
    {
        // diagnostics for all parameters
        for (const parID &p : param_to_be_inferred)
            dgn_param_value.at(p)[dgn_count] = get_param(p);

        // diagnostics for all households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].update_diagnostics();
    }

    ++dgn_count;
}

/* Disable all assert methods for user definitions. */
void Monte_Carlo_Build::disable_assert()
{
    // disable parameter assertion for all parameters
    _func_get_param_ = &Monte_Carlo_Build::_get_param_default_;
    _func_get_param_pointer_ = &Monte_Carlo_Build::_get_param_pointer_default_;

    // disable augmented variable assertion for all individuals
    for (Household &house : house_vector)
        house.disable_assert();
}

/* Compute the loglikelihood term associated to augmented variables. */
double Monte_Carlo_Build::compute_augm_model_loglike(const Household &house) const
{
    return (this->*_func_augm_model_loglike_)(house);
}

/* No model contributions. */
double Monte_Carlo_Build::_no_augm_model_loglike_(const Household &house) const
{
    UNUSED(house);
    return 0.0;
}

/* House loglikelihood including augmented modl contribution. */
double Monte_Carlo_Build::_augm_model_loglike_(const Household &house) const
{
    double hh_loglike = 0.;
    // Contribution of model associated to augmented variables
    for (const indiv &i : house.infectious_set())
        (this->*_func_augm_model_by_definition_)(house, i, hh_loglike);

    return hh_loglike;
}

void Monte_Carlo_Build::_augm_model_by_status_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, a);
}

void Monte_Carlo_Build::_augm_model_generic_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_list_model)
        delta += house.augm_model_logpdf(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void Monte_Carlo_Build::show_summary() const
{
    std::ostringstream oss;

    oss << "n. of individuals: " << n_individuals
        << "\n    * n. of cases: " << _report_data_.at("n.cases")
        << " " << _aux_percent_(_report_data_.at("n.cases"), n_individuals);

    if (_report_data_.at("n.miss.stat") > 0)
        oss << "\n    * n. of missing status : " << _report_data_.at("n.miss.stat")
            << " " << _aux_percent_(_report_data_.at("n.miss.stat"), n_individuals);

    oss << "\nn. of households: " << n_households;

    if (_report_data_.at("no.HH.case") > 0)
        oss << "\n    * with no cases: " << _report_data_.at("no.HH.case")
            << " " << _aux_percent_(_report_data_.at("no.HH.case"), n_households);

    oss << "\n    * with at least 1 case: " << _report_data_.at("lst.1.case")
        << " " << _aux_percent_(_report_data_.at("lst.1.case"), n_households);

    if (_report_data_.at("lst.1.case") > 0)
        oss << "\n    * with at least 2 cases: " << _report_data_.at("lst.2.case")
            << " " << _aux_percent_(_report_data_.at("lst.2.case"), n_households);

    oss << "\nmax n. house cases: " << _report_data_.at("max.n.HH.cases")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.cases"), n_households)
        << "\nmin n. house members: " << _report_data_.at("min.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.min.size"), n_households)
        << "\nmax n. house members: " << _report_data_.at("max.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.size"), n_households);

    tools::print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : show full build  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Monte_Carlo_Build::show_full_build() const
{
    //----- SETUP -----//
    std::string add_file_info = UNDEFINED_STRING;
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("setup")))
            add_file_info = "# import: " + _report_pathtofile_.at("setup");

    tools::print_decorator("Setup", add_file_info);
    show_setup();

    //----- PARAMETERS -----//
    add_file_info = UNDEFINED_STRING;
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("params")))
            add_file_info = "# import: " + _report_pathtofile_.at("params");

    tools::print_decorator("Parameters", add_file_info);
    show_all_params();

    //----- DATABASE -----//
    add_file_info = UNDEFINED_STRING;
    if (_report_pathtofile_.find("data") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("data")))
            add_file_info = "# import: " + _report_pathtofile_.at("data");

    tools::print_decorator("Database", add_file_info);
    show_summary();

    //----- AUGMENTED -----//
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("augm")))
            add_file_info = "# import: " + _report_pathtofile_.at("augm");

    tools::print_decorator("Augmented", add_file_info);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : print augm variables  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void Monte_Carlo_Build::print_augm_variables(const std::string &pathtofile) const
{
    if (!tools::is_na(pathtofile))
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            //===== Print column headers =====//
            // default headers
            augm_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

            // data headers
            for (const dataID &d : data_headers)
                augm_ofstream << C_DELIM << d;

            // augmented variables headers
            for (const augmID &a : augm_list_full)
                augm_ofstream << C_DELIM << a;

            augm_ofstream << std::endl;

            //===== Print data and augmented variables =====//
            for (const Household &house : house_vector)
            {
                for (indiv i = 0; i < house.size(); ++i)
                {
                    augm_ofstream << house.ID(i)
                                  << C_DELIM << house.house_ID()
                                  << C_DELIM << house.status(i);
                    for (const dataID &d : data_headers)
                        augm_ofstream << C_DELIM << house.data(i, d);

                    for (const augmID &a : augm_list_full)
                        if (house.status(i) == INFECTIOUS)
                            augm_ofstream << C_DELIM << house.augm(i, a);
                        else
                            augm_ofstream << C_DELIM << std::numeric_limits<double>::quiet_NaN();

                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            throw std::runtime_error(
                "Unable to open file: " + pathtofile);
        }
    }
}
