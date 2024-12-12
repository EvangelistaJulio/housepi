/*
 * build.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "build.h"

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

/* Auxiliary function to call percentuage. */
inline const std::string _aux_percent_(int num, int den)
{
    return "(" + tools::percentage(static_cast<double>(num) / den, 1) + ")";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Monte_Carlo_Build : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double Monte_Carlo_Build::type_to_double(const std::string value, const std::string type, const Household &house, indiv i)
{
    double output = INF_DOUBLE;

    try
    {
        if (tools::is_na(value) && tools::is_na(type))
            ;

        else if (!tools::is_na(value) && tools::is_na(type))
            output = tools::stod(value, "argument without type");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Ignoring value type without a paired value.");

        else if (type == TYPE_SCALAR)
            output = tools::stod(value, "argument with scalar type");

        else if (type == TYPE_DATA)
            output = house.data(i, value);

        else if (type == TYPE_PARAM)
            if (get_infer(value))
                tools::error(
                    "Non constant parameter '" + value + "' assigned as immutable argument.");
            else
                output = get_param(value);

        else if (type == TYPE_COMMAND)
            if (value == "default")
                output = INF_DOUBLE;
            else
                throw std::invalid_argument(
                    "Unknown command type '" + value + "'.");
        else
            throw std::invalid_argument(
                "Unknown variable type specifier '" + type + "'.");
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a type-std::string to a double-pointer */
const double *Monte_Carlo_Build::type_to_doubleptr(const std::string value, const std::string type, const Household &house, indiv i)
{
    const double *output = nullptr;

    try
    {
        if (tools::is_na(value) && tools::is_na(type))
        {
            output = nullptr; // default value
        }
        else if (!tools::is_na(value) && tools::is_na(type))
        {
            throw std::invalid_argument(
                "Impossible to interpret default value of  '" + value + "'. Please provide a type.");
        }
        else if (tools::is_na(value) && !tools::is_na(type))
        {
            tools::warning("Found a type without a paired value. Falling back to default.");

            output = nullptr; // default value
        }
        else
        {
            if (type == TYPE_DATA)
                output = house.data_pointer(i, value);

            else if (type == TYPE_PARAM)
                output = get_param_pointer(value);

            else if (type == TYPE_SCALAR)
                output = assign_scalar_pointer(tools::stod(value, "argument with scalar type"));

            else
                throw std::invalid_argument(
                    "Unsopported type '" + type + "' found while loading '" + value + "' for pointer value.");
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
dictargs Monte_Carlo_Build::typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i)
{
    std::vector<std::string> keys = tools::parse_list_args(tkn_key);
    std::vector<std::string> vals = tools::parse_list_args(tkn_val);
    std::vector<std::string> typs = tools::parse_list_args(tkn_typ);

    dictargs output = {};

    if ((keys.size() == vals.size()) && (typs.size() == vals.size()))
        for (size_t x = 0; x < keys.size(); ++x)
            output.emplace(keys[x], type_to_double(vals[x], typs[x], house, i));
    else
        throw std::invalid_argument(
            "Mismatching definitions in arguments.");

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double-pointer */
dictargs_ptr Monte_Carlo_Build::typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i)
{
    std::vector<std::string> keys = tools::parse_list_args(tkn_key);
    std::vector<std::string> vals = tools::parse_list_args(tkn_val);
    std::vector<std::string> typs = tools::parse_list_args(tkn_typ);

    dictargs_ptr output = {};

    if ((keys.size() == vals.size()) && (typs.size() == vals.size()))
        for (size_t x = 0; x < keys.size(); ++x)
            output.emplace(keys[x], type_to_doubleptr(vals[x], typs[x], house, i));
    else
        throw std::invalid_argument(
            "Mismatching definitions in arguments.");

    return output;
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
        tools::error("Unknown summary entry '" + key + "'.");

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
    const std::vector<std::string> augm_file_headers = {
        "name", "status", "value", "valuetype", "discrete",
        "model", "args", "argsval", "argstype", "rule",
        "rule_scale", "rule_scaletype", "rule_from", "rule_fromtype",
        "min", "mintype", "max", "maxtype", "init_rand"};

    tools::read_csv_table(pathtofile, _augm_deftable_, augm_file_headers);
    build_augm();
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

            // auxiliary variable to load definitions : augmID -> status lists
            std::map<std::string, std::vector<std::string>> name_status_check = {};

            //===== Check Augmented Variables lists of usage =====//
            for (size_t a = 0; a < table_length; ++a)
            {
                std::string augmID = _augm_deftable_.at("name")[a];
                std::string status = _augm_deftable_.at("status")[a];
                bool is_init_rand = false;

                // check that augmID and status are correct
                if (tools::is_na(augmID))
                    throw std::invalid_argument(
                        "Augmented variable name can not be undefined.");

                else if (name_status_check.find(augmID) == name_status_check.end())
                    // new augmented variable status
                    name_status_check.emplace(augmID, std::vector<std::string>{status});

                else
                {
                    // check no multiple augmID with same status
                    const std::vector<std::string> &mylist = name_status_check.at(augmID);
                    if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
                        name_status_check.at(augmID).push_back(status);

                    else
                        throw std::invalid_argument(
                            "Incompatible definitions for augmented variable '" + augmID + "'.");
                }

                // option for random initialization
                if (!tools::is_na(_augm_deftable_.at("init_rand")[a]))
                    is_init_rand = tools::stobool(_augm_deftable_.at("init_rand")[a]);

                if ((tools::is_na(_augm_deftable_.at("model")[a])) && (is_init_rand))
                    throw std::invalid_argument(
                        "Impossible to initialize '" + augmID + "' without a model.");
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
                const std::string &augmID = pair.first;
                const std::vector<std::string> &mylist = pair.second;

                for (const std::string &status : mylist)
                {
                    if (tools::is_na(status))
                    {
                        if (mylist.size() > 1)
                            throw std::invalid_argument(
                                "Conflicting status-specific and general definitions for augmented '" + augmID + "'.");
                    }
                    else
                    {
                        int s = tools::stoi(status, "augmented status");
                        if (s < INFECTIOUS)
                        {
                            throw std::invalid_argument(
                                "Impossible to associate augmented '" + augmID + "' for non infectious status.");
                        }
                        else
                        {
                            augm_status_map_infer.emplace(s, std::unordered_set<std::string>{});
                            augm_status_map_model.emplace(s, std::unordered_set<std::string>{});
                        }
                    }

                    augm_list_full.emplace(augmID);
                }
            }

            //===== Case of Only Generic Definitions =====//
            bool is_only_generic = (augm_list_full.size() > 0) && (augm_status_map_infer.empty());
            if (is_only_generic)
                _func_augm_model_by_definition_ = &Monte_Carlo_Build::_augm_model_generic_;

            //===== Assign new definitions to lists =====//
            // note : I can't I collapse this with the previous because of generic status
            for (size_t a = 0; a < table_length; ++a)
            {
                std::string augmID = _augm_deftable_.at("name")[a];
                std::string status = _augm_deftable_.at("status")[a];

                //----- Augmented Variables to be Inferred -----//
                if (!tools::is_na(_augm_deftable_.at("rule")[a]))
                {
                    is_augm_infer = true;

                    if (tools::is_na(status)) // generic definition
                        if (is_only_generic)
                            augm_list_infer.emplace(augmID);

                        else
                            for (const auto &pair : augm_status_map_infer)
                                augm_status_map_infer.at(pair.first).emplace(augmID);

                    else // status-specific definition
                        augm_status_map_infer.at(tools::stoi(status)).emplace(augmID);
                }

                //----- Augmented Variables with a Model -----//
                if (!tools::is_na(_augm_deftable_.at("model")[a]))
                {
                    is_augm_model = true;

                    if (tools::is_na(status)) // generic definition
                        if (is_only_generic)
                            augm_list_model.emplace(augmID);

                        else
                            for (const auto &pair : augm_status_map_model)
                                augm_status_map_model.at(pair.first).emplace(augmID);

                    else // status-specific definition
                        augm_status_map_model.at(tools::stoi(status)).emplace(augmID);
                }
            }

            //===== Assign augmented variables to each individual =====//
            for (Household &house : house_vector)
            {
                for (indiv i = 0; i < house.size(); ++i)
                {
                    // allocate memory for all augmented for all individuals
                    clear_all_augm(house, i);
                    assign_all_augm(house, i);
                }
            }

            //===== Functional method for augmented model contribution to loglikelihood ====//
            if (is_augm_model)
                _func_augm_model_loglike_ = &Monte_Carlo_Build::_augm_model_loglike_;
            else
                _func_augm_model_loglike_ = &Monte_Carlo_Build::_no_augm_model_loglike_;
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
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
        for (size_t a = 0; a < _augm_deftable_.at("name").size(); ++a)
        {
            // assign augmented variable of the corresponding status
            std::string augm_status = _augm_deftable_.at("status")[a];

            // assign augmented variable of the corresponding status or undefined status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                //----- Name -----//
                const std::string &augmID = _augm_deftable_.at("name")[a];

                //----- is Discrete -----//
                bool is_discr = tools::stobool(_augm_deftable_.at("discrete")[a]);

                //----- Value -----//
                double value = type_to_double(_augm_deftable_.at("value")[a],
                                              _augm_deftable_.at("valuetype")[a],
                                              house, i);

                //----- Model -----//
                const std::string &model_name = _augm_deftable_.at("model")[a];
                dictargs_ptr model_args = typedict_to_sdargs_ptr(_augm_deftable_.at("args")[a],
                                                                 _augm_deftable_.at("argsval")[a],
                                                                 _augm_deftable_.at("argstype")[a],
                                                                 house, i);

                //----- Range -----//
                dictargs range = {};
                // minimum value
                double aux_min = type_to_double(_augm_deftable_.at("min")[a],
                                                _augm_deftable_.at("mintype")[a],
                                                house, i);
                if (aux_min != INF_DOUBLE)
                    range.emplace("min", aux_min);

                // maximum value
                double aux_max = type_to_double(_augm_deftable_.at("max")[a],
                                                _augm_deftable_.at("maxtype")[a],
                                                house, i);
                if (aux_max != INF_DOUBLE)
                    range.emplace("max", aux_max);

                //----- Proposer Rule -----//
                const std::string &rule_name = _augm_deftable_.at("rule")[a];
                double rule_scale = type_to_double(_augm_deftable_.at("rule_scale")[a],
                                                   _augm_deftable_.at("rule_scaletype")[a],
                                                   house, i);
                double rule_from = type_to_double(_augm_deftable_.at("rule_from")[a],
                                                  _augm_deftable_.at("rule_fromtype")[a],
                                                  house, i);

                //----- Assign Augmented Variable -----//
                house.assign_augm(i, augmID, value, is_discr,
                                  model_name, model_args, range,
                                  rule_name, rule_scale, rule_from);

                // random initialization at assignment
                if (tools::stobool(_augm_deftable_.at("init_rand")[a]))
                    house.augm_model_rand(i, augmID, get_prng());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

/* Clear all augmented definition of individual. */
void Monte_Carlo_Build::clear_all_augm(Household &house, indiv i)
{
    for (const std::string &augmID : augm_list_full)
        house.clear_augm(i, augmID);
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
    dgn_glob_loglike = std::vector<double>(sample_size);

    if (!is_power_post)
    {
        //----- Initialize diagnostic vectors -----//
        for (const std::string &parID : param_to_be_inferred)
            param_map.at(parID).init_diagnostics(sample_size);

        // diagnostics over households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].init_diagnostics(sample_size);
    }
}

/* Update Bayesian diagnostics */
void Monte_Carlo_Build::update_diagnostics(double glob_loglike, bool is_power_post)
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count++] = glob_loglike;

    if (!is_power_post)
    {
        // diagnostics for all parameters
        for (const std::string &parID : param_to_be_inferred)
            param_map.at(parID).update_diagnostics();

        // diagnostics for all households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].update_diagnostics();
    }
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
    for (const std::string &augmID : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, augmID);
}

void Monte_Carlo_Build::_augm_model_generic_(const Household &house, indiv i, double &delta) const
{
    for (const std::string &augmID : augm_list_model)
        delta += house.augm_model_logpdf(i, augmID);
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
    std::string add_file_info = UNDEFINED;
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("setup")))
            add_file_info = "# import: " + _report_pathtofile_.at("setup");

    tools::print_decorator("Setup", add_file_info);
    show_setup();

    //----- PARAMETERS -----//
    add_file_info = UNDEFINED;
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("params")))
            add_file_info = "# import: " + _report_pathtofile_.at("params");

    tools::print_decorator("Parameters", add_file_info);
    show_all_params();

    //----- DATABASE -----//
    add_file_info = UNDEFINED;
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
            for (const std::string &dataID : data_headers)
                augm_ofstream << C_DELIM << dataID;

            // augmented variables headers
            for (const std::string &augmID : augm_list_full)
                augm_ofstream << C_DELIM << augmID;

            augm_ofstream << std::endl;

            //===== Print data and augmented variables =====//
            for (const Household &house : house_vector)
            {
                for (indiv i = 0; i < house.size(); ++i)
                {
                    augm_ofstream << house.ID(i)
                                  << C_DELIM << house.house_ID()
                                  << C_DELIM << house.status(i);
                    for (const std::string &dataID : data_headers)
                        augm_ofstream << C_DELIM << house.data(i, dataID);

                    for (const std::string &augmID : augm_list_full)
                        if (house.status(i) == INFECTIOUS)
                            augm_ofstream << C_DELIM << house.augm(i, augmID);
                        else
                            augm_ofstream << C_DELIM << std::numeric_limits<double>::quiet_NaN();

                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            tools::error("Unable to open file: " + pathtofile);
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame then run build. */
void Monte_Carlo_Build::import_augm(const Rcpp::DataFrame &augm)
{
    tools::rdataframe_to_map(augm, _augm_deftable_);
    build_augm();
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
        _report_pathtofile_.erase("augm");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Monte_Carlo_Build : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame Monte_Carlo_Build::export_data() const
{
    using namespace Rcpp;

    CharacterVector df_index(n_individuals);
    CharacterVector df_house_index(n_individuals);
    IntegerVector df_status(n_individuals);
    std::vector<NumericVector> df_dbl;
    for (size_t s = 0; s < data_headers.size() + augm_list_full.size(); ++s)
        df_dbl.emplace_back(n_individuals);

    size_t r = 0, c = 0;
    for (const Household &house : house_vector)
    {
        for (indiv i = 0; i < house.size(); ++i)
        {
            df_index[r] = house.ID(i);
            df_house_index[r] = house.house_ID();
            df_status[r] = house.status(i);

            c = 0;
            for (const std::string &dataID : data_headers)
                df_dbl[c++][r] = house.data(i, dataID);

            for (const std::string &augmID : augm_list_full)
                df_dbl[c++][r] = house.augm(i, augmID);

            ++r;
        }
    }

    DataFrame output;
    output = DataFrame::create(_["id"] = df_index,
                               _["house_id"] = df_house_index,
                               _["status"] = df_status);
    c = 0;
    for (const std::string &dataID : data_headers)
        output.push_back(df_dbl[c++], dataID);

    for (const std::string &augmID : augm_list_full)
        output.push_back(df_dbl[c++], augmID);

    return output;
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix Monte_Carlo_Build::export_loglike_matrix() const
{
    using namespace Rcpp;

    // FIXME : this should be prevented if is_power_post

    // Create a NumericMatrix with n_house columns
    ulong sample_size = setup_n_chains * setup_c_sample_size;
    NumericMatrix output(sample_size, n_households);

    // Assign each NumericVector to a column of the matrix
    for (size_t H = 0; H < n_households; ++H)
    {
        NumericVector column = wrap(house_vector[H].get_dgn_loglike());
        output(_, H) = column;
    }

    return output;
}
#endif
