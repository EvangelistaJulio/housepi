/*
 * build.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "build.h"
#include "core.h"

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

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Monte_Carlo_Build : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Monte_Carlo_Build::Monte_Carlo_Build() : pCore(std::make_unique<Monte_Carlo_Core>()) {}
Monte_Carlo_Build::~Monte_Carlo_Build() = default;

#ifdef USE_RCPP
//===== Rcpp IMPORTS =====//
void Monte_Carlo_Build::import_setup(Rcpp::List x)
{
    pCore->import_setup(x);
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
    {
        _report_pathtofile_.erase("setup");
    }
}
void Monte_Carlo_Build::import_params(Rcpp::DataFrame x)
{
    pCore->import_params(x);
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
    {
        _report_pathtofile_.erase("params");
    }
}

//===== Rcpp EXPORTS =====//
Rcpp::List Monte_Carlo_Build::export_setup() const { return pCore->export_setup(); }
Rcpp::DataFrame Monte_Carlo_Build::export_params() const { return pCore->export_params(); }
Rcpp::DataFrame Monte_Carlo_Build::export_traces() const { return pCore->export_traces(); }
#else
//===== Cpp Show/Print =====//
void Monte_Carlo_Build::show_results() const { pCore->show_results(); }
void Monte_Carlo_Build::print_infer_state(const std::string &pathtofile) const { pCore->print_infer_state(pathtofile); }
#endif

void Monte_Carlo_Build::read_setup(const std::string &pathtofile)
{
    pCore->read_setup(pathtofile);
    _report_pathtofile_["setup"] = pathtofile;
}

/* Read parameter definitions and run assignment. */
void Monte_Carlo_Build::read_params(const std::string &pathtofile)
{
    pCore->read_params(pathtofile);
    _report_pathtofile_["params"] = pathtofile;
}

bool Monte_Carlo_Build::get_rnd_init_params() const { return pCore->setup_rnd_init_params; }
ulong Monte_Carlo_Build::get_num_chains() const { return pCore->setup_num_chains; }
ulong Monte_Carlo_Build::get_num_iter() const { return pCore->setup_num_iter; }
ulong Monte_Carlo_Build::get_num_adapt() const { return pCore->setup_num_adapt; }
ulong Monte_Carlo_Build::get_thinning() const { return pCore->setup_thinning; }
ulong Monte_Carlo_Build::get_burn_in() const { return pCore->setup_burn_in; }
double Monte_Carlo_Build::get_inv_T() const { return pCore->setup_inv_T; }
unsigned Monte_Carlo_Build::get_num_threads() const { return pCore->setup_num_threads; }
ulong Monte_Carlo_Build::get_c_num_sweeps() const { return pCore->setup_c_num_sweeps; }

void Monte_Carlo_Build::set_rnd_init_params(bool x) { pCore->setup_rnd_init_params = x; }
void Monte_Carlo_Build::set_inv_T(double x) { pCore->setup_inv_T = x; }

sys_rnd_gen &Monte_Carlo_Build::get_prng() { return pCore->get_prng(); }
double Monte_Carlo_Build::uniform() { return pCore->uniform(); }

void Monte_Carlo_Build::init_rand_params() { pCore->init_rand_params(); }
void Monte_Carlo_Build::assert_param(const std::string &parID) const { pCore->assert_param(parID); }
bool Monte_Carlo_Build::get_infer(const std::string &parID) { return pCore->get_infer(parID); }
double Monte_Carlo_Build::get_param(const std::string &parID) { return pCore->get_param(parID); }
const double *Monte_Carlo_Build::get_param_pointer(const std::string &parID) const { return pCore->get_param_pointer(parID); }
void Monte_Carlo_Build::dgn_glob_ll_stats(double &avg, double &var) const { pCore->dgn_glob_ll_stats(avg, var); }

std::vector<std::string> Monte_Carlo_Build::get_param_to_be_inferred() const { return pCore->param_to_be_inferred; }
bool Monte_Carlo_Build::propose_valid_param_move(const std::string &parID) { return pCore->propose_valid_param_move(parID); }
void Monte_Carlo_Build::restore_param_value(const std::string &parID) { pCore->restore_param_value(parID); }
void Monte_Carlo_Build::update_param_accept_rate(const std::string &parID, bool is_accepted) { pCore->update_param_accept_rate(parID, is_accepted); }
double Monte_Carlo_Build::param_move_logratio(const std::string &parID) const { return pCore->param_move_logratio(parID); }
void Monte_Carlo_Build::adapt_param_rule(const std::string &parID, double ratio) { pCore->adapt_param_rule(parID, ratio); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Monte_Carlo_Build : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double Monte_Carlo_Build::type_to_double(const std::string val, const std::string typ, const Individual &this_indiv)
{
    double output = INF_DOUBLE;

    try
    {
        if (tools::is_na(val) && tools::is_na(typ))
        {
            output = INF_DOUBLE; // default value
        }
        else if (!tools::is_na(val) && tools::is_na(typ))
        {
            output = tools::stod(val, "argument without type");
        }
        else if (tools::is_na(val) && !tools::is_na(typ))
        {
            tools::warning("Ignoring value type without a paired value.");

            output = INF_DOUBLE; // FIXME default value
        }
        else
        {
            if (typ == TYPE_SCALAR)
            {
                output = tools::stod(val, "argument with scalar type");
            }
            else if (typ == TYPE_DATA)
            {
                output = this_indiv.get_data(val);
            }
            else if (typ == TYPE_PARAM)
            {
                output = get_param(val);
                if (get_infer(val))
                {
                    std::ostringstream oss;
                    // FIXME : this warning is not clear as exposed here
                    oss << "Non constant parameter '" << val << "' assigned as fixed argument."
                        << " Changes of the parameter will not affect the argument.";
                    tools::warning(oss.str());
                }
            }
            else if (typ == TYPE_COMMAND)
            {
                if (val == "default")
                {
                    output = INF_DOUBLE;
                }
                else
                {
                    throw std::invalid_argument("Unknown command type '" + val + "'.");
                }
            }
            else
            {
                throw std::invalid_argument("Unknown variable type specifier '" + typ + "'.");
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a type-std::string to a double-pointer */
const double *Monte_Carlo_Build::type_to_doubleptr(const std::string val, const std::string typ, const Individual &this_indiv)
{
    const double *output = nullptr;

    try
    {
        if (tools::is_na(val) && tools::is_na(typ))
        {
            output = nullptr; // default value
        }
        else if (!tools::is_na(val) && tools::is_na(typ))
        {
            throw std::invalid_argument("Impossible to interpret default value of  '" + val + "'. Please provide a type.");
        }
        else if (tools::is_na(val) && !tools::is_na(typ))
        {
            tools::warning("Found a type without a paired value. Falling back to default.");

            output = nullptr; // default value
        }
        else
        {
            if (typ == TYPE_DATA)
            {
                output = this_indiv.get_data_pointer(val);
            }
            else if (typ == TYPE_PARAM)
            {
                output = get_param_pointer(val);
            }
            else if (typ == TYPE_SCALAR)
            {
                double new_val = tools::stod(val, "argument with scalar type");
                if (_scalar_type_params_.find(new_val) == _scalar_type_params_.end())
                {
                    Parameter new_param(UNDEFINED, new_val, UNDEFINED, {}, UNDEFINED, {}, {});
                    _scalar_type_params_.emplace(new_val, new_param);
                }
                /* Add a new scalar for augmented typedict FIXME : temporary. */
                output = _scalar_type_params_.at(new_val).get_pointer();
            }
            else
            {
                throw std::invalid_argument("Unsopported type '" + typ + "' found loading '" + val + "' for pointer value.");
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
dictargs Monte_Carlo_Build::typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv)
{
    std::vector<std::string> keys = tools::parse_list_args(tkn_key);
    std::vector<std::string> vals = tools::parse_list_args(tkn_val);
    std::vector<std::string> typs = tools::parse_list_args(tkn_typ);

    dictargs output = {};

    if (!((keys.size() == vals.size()) && (typs.size() == vals.size())))
    {
        // inform augmented variable about the error
        output.emplace(IO_MISMATCH, INF_DOUBLE);
    }
    else
    {
        for (size_t x = 0; x < keys.size(); ++x)
        {
            output.emplace(keys[x], type_to_double(vals[x], typs[x], this_indiv));
        }
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double-pointer */
dictargs_ptr Monte_Carlo_Build::typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv)
{
    std::vector<std::string> keys = tools::parse_list_args(tkn_key);
    std::vector<std::string> vals = tools::parse_list_args(tkn_val);
    std::vector<std::string> typs = tools::parse_list_args(tkn_typ);

    dictargs_ptr output = {};

    if (!((keys.size() == vals.size()) && (typs.size() == vals.size())))
    {
        // inform augmented variable about the error
        output.emplace(IO_MISMATCH, nullptr);
    }
    else
    {
        for (size_t x = 0; x < keys.size(); ++x)
        {
            output.emplace(keys[x], type_to_doubleptr(vals[x], typs[x], this_indiv));
        }
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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
    for (const Household &this_house : house)
    {
        int this_house_size = this_house.size();
        int num_infec = this_house.num_infectious;

        _report_data_["n.cases"] += num_infec;
        _report_data_["n.miss.stat"] += this_house.num_missing_status;

        _report_data_["no.HH.case"] += (num_infec == 0);
        _report_data_["lst.1.case"] += (num_infec > 0);
        _report_data_["lst.2.case"] += (num_infec > 1);

        // max number of house cases
        if (_report_data_["max.n.HH.cases"] < num_infec)
        {
            _report_data_["n.HH.with.max.cases"] = 1;
            _report_data_["max.n.HH.cases"] = num_infec;
        }
        else if (_report_data_["max.n.HH.cases"] == num_infec)
        {
            _report_data_["n.HH.with.max.cases"] += 1;
        }

        // minimum house number of members
        if (_report_data_["min.HH.size"] > this_house_size)
        {
            _report_data_["n.HH.with.min.size"] = 1;
            _report_data_["min.HH.size"] = this_house_size;
        }
        else if (_report_data_["min.HH.size"] == this_house_size)
        {
            _report_data_["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (_report_data_["max.HH.size"] < this_house_size)
        {
            _report_data_["n.HH.with.max.size"] = 1;
            _report_data_["max.HH.size"] = this_house_size;
        }
        else if (_report_data_["max.HH.size"] == this_house_size)
        {
            _report_data_["n.HH.with.max.size"] += 1;
        }
    }
}

/* Return a value from report data. */
int Monte_Carlo_Build::report_data(const std::string &key) const
{
    int output = 0;
    if (_report_data_.find(key) != _report_data_.end())
    {
        output = _report_data_.at(key);
    }
    else
    {
        std::ostringstream oss;
        oss << "Unknown summary entry `" << key << "`.";
        tools::error(oss.str());
    }
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
        "name", "status", "value", "valuetype",
        "model", "args", "argsval", "argstype", "rule",
        "rule_args", "rule_argsval", "rule_argstype", "rule_from", "rule_fromtype",
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
            size_t table_length = _augm_deftable_.at("name").size();

            //===== Drop previous definitions =====//
            std::map<std::string, std::vector<std::string>> augm_name_status_check = {};

            augm_full_set.clear();
            augm_to_be_inferred.clear();
            augm_with_model_set.clear();

            //===== Augmented variables lists of usage =====//
            for (size_t a = 0; a < table_length; ++a)
            {
                std::string augmID = _augm_deftable_.at("name")[a];
                std::string status = _augm_deftable_.at("status")[a];
                bool is_init_rand = false;

                if (tools::is_na(augmID))
                {
                    std::ostringstream oss;
                    oss << "Augmented variable name can not be undefined.";
                    throw std::invalid_argument(oss.str());
                }
                else if (augm_full_set.find(augmID) == augm_full_set.end())
                {
                    augm_full_set.emplace(augmID);
                    augm_name_status_check[augmID].push_back(status);
                }
                else
                {
                    const std::vector<std::string> &mylist = augm_name_status_check.at(augmID);
                    if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
                    {
                        augm_name_status_check.at(augmID).push_back(status);
                    }
                    else
                    {
                        std::ostringstream oss;
                        oss << "Multiple definitions for augmented variable '" << augmID << "'."
                            << " Check your augmented variable file.";
                        throw std::invalid_argument(oss.str());
                    }
                }

                // option for random initialization
                if (!tools::is_na(_augm_deftable_.at("init_rand")[a]))
                {
                    is_init_rand = tools::stobool(_augm_deftable_.at("init_rand")[a]);
                }

                // list of augmented variables to be inferred
                if (!tools::is_na(_augm_deftable_.at("rule")[a]))
                {
                    augm_to_be_inferred.emplace(augmID);
                }

                // list of augmented variables with a model
                // note: no contribution of constant variables to Metropolis
                if (tools::is_na(_augm_deftable_.at("model")[a]))
                {
                    if (is_init_rand)
                    {
                        throw std::invalid_argument("Impossible to initialize '" + augmID + "' without a model.");
                    }
                }
                else
                {
                    augm_with_model_set.emplace(augmID);
                }
            }

            //===== Check Univocal Definitions =====//
            for (const auto &pair : augm_name_status_check)
            {
                const std::vector<std::string> &mylist = pair.second;

                for (const std::string &z : mylist)
                {
                    if (tools::is_na(z))
                    {
                        if (mylist.size() > 1)
                        {
                            std::ostringstream oss;
                            oss << "Forbidden default and status definitions for augmented variable '"
                                << pair.first << "'.";
                            throw std::invalid_argument(oss.str());
                        }
                    }
                    else if (tools::stoi(z, "augmented status") < INFECTIOUS)
                    {
                        std::ostringstream oss;
                        oss << "Impossible to associate augmented '" << pair.first << "' for non infectious status.";
                        throw std::invalid_argument(oss.str());
                    }
                }
            }

            //===== Set internal state variable =====//
            is_move_augm = !augm_to_be_inferred.empty();

            //===== Assign augmented variables to each individual =====//
            for (Household &this_house : house)
            {
                for (Individual &this_indiv : this_house.indiv)
                {
                    assign_augm(this_indiv);
                }
            }
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
void Monte_Carlo_Build::assign_augm(Individual &this_indiv)
{
    try
    {
        for (size_t a = 0; a < _augm_deftable_.at("name").size(); ++a)
        {
            // assign augmented variable of the corresponding status
            std::string augm_status = _augm_deftable_.at("status")[a];

            if ((tools::is_na(augm_status)) || (this_indiv.get_status() == tools::stoi(augm_status, "augmented status")))
            {
                std::string augmID = _augm_deftable_.at("name")[a];

                //----- load range -----//
                dictargs range = {};
                // min value
                double aux_min = type_to_double(_augm_deftable_.at("min")[a],
                                                _augm_deftable_.at("mintype")[a], this_indiv);
                if (aux_min != INF_DOUBLE)
                {
                    range.emplace("min", aux_min);
                }

                // max value
                double aux_max = type_to_double(_augm_deftable_.at("max")[a],
                                                _augm_deftable_.at("maxtype")[a], this_indiv);
                if (aux_max != INF_DOUBLE)
                {
                    range.emplace("max", aux_max);
                }

                // assign augmented variable
                Augmented new_augm(augmID,
                                   type_to_double(_augm_deftable_.at("value")[a],
                                                  _augm_deftable_.at("valuetype")[a], this_indiv),
                                   _augm_deftable_.at("model")[a],
                                   typedict_to_sdargs_ptr(_augm_deftable_.at("args")[a],
                                                          _augm_deftable_.at("argsval")[a],
                                                          _augm_deftable_.at("argstype")[a], this_indiv),
                                   _augm_deftable_.at("rule")[a],
                                   typedict_to_sdargs(_augm_deftable_.at("rule_args")[a],
                                                      _augm_deftable_.at("rule_argsval")[a],
                                                      _augm_deftable_.at("rule_argstype")[a], this_indiv),
                                   range,
                                   type_to_double(_augm_deftable_.at("rule_from")[a],
                                                  _augm_deftable_.at("rule_fromtype")[a], this_indiv));
                this_indiv.assign_augm(new_augm);

                // random initialization at assignment
                if (tools::stobool(_augm_deftable_.at("init_rand")[a]))
                {
                    this_indiv.augm_model_rand(augmID, get_prng());
                }
            }
            else if (this_indiv.get_status() < INFECTIOUS)
            {
                // placeholders for non infectious
                clear_augm(this_indiv);
            }
        }

        // assert that individual has all definitions
        for (const std::string &augmID : augm_full_set)
        {
            this_indiv.assert_augm(augmID);
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

/* Remove augmented definition of individual. */
void Monte_Carlo_Build::clear_augm(Individual &this_indiv)
{
    for (const std::string &augmID : augm_full_set)
    {
        Augmented new_augm(augmID);
        this_indiv.assign_augm(new_augm);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void Monte_Carlo_Build::show_summary() const
{
    using tools::percent;
    //===== PRINT =====//

    std::ostringstream oss;

    oss << "n. of individuals: " << num_indiv
        << "\n    * n. of cases: " << _report_data_.at("n.cases")
        << " " << percent(_report_data_.at("n.cases"), num_indiv);
    if (_report_data_.at("n.miss.stat") > 0)
    {
        oss << "\n    * n. of missing status : " << _report_data_.at("n.miss.stat")
            << " " << percent(_report_data_.at("n.miss.stat"), num_indiv);
    }

    oss << "\nn. of households: " << num_house;
    if (_report_data_.at("no.HH.case") > 0)
    {
        oss << "\n    * with no cases: " << _report_data_.at("no.HH.case")
            << " " << percent(_report_data_.at("no.HH.case"), num_house);
    }

    oss << "\n    * with at least 1 case: " << _report_data_.at("lst.1.case")
        << " " << percent(_report_data_.at("lst.1.case"), num_house);
    if (_report_data_.at("lst.1.case") > 0)
    {
        oss << "\n    * with at least 2 cases: " << _report_data_.at("lst.2.case")
            << " " << percent(_report_data_.at("lst.2.case"), num_house);
    }

    oss << "\nmax n. house cases: " << _report_data_.at("max.n.HH.cases")
        << " " << percent(_report_data_.at("n.HH.with.max.cases"), num_house)
        << "\nmin n. house members: " << _report_data_.at("min.HH.size")
        << " " << percent(_report_data_.at("n.HH.with.min.size"), num_house)
        << "\nmax n. house members: " << _report_data_.at("max.HH.size")
        << " " << percent(_report_data_.at("n.HH.with.max.size"), num_house);

    tools::print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Monte_Carlo_Build : show group classes  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show model group classes. */
void Monte_Carlo_Build::show_group_classes() const
{
    std::ostringstream oss;

    // format group classes
    for (const auto &entry : _report_group_classes_)
    {
        if ((entry.first != CONTACT_RATE) && (entry.first != EMPTY_RATE))
        {
            oss << entry.first << ":";
            for (const auto &x : entry.second)
            {
                oss << "\n    * " << x.first << ": " << x.second
                    << " " << tools::percent(x.second, num_indiv);
            }
            oss << "\n";
        }
    }

    // format contact rates
    if (_report_group_classes_.at(CONTACT_RATE).find(EMPTY_RATE) != _report_group_classes_.at(CONTACT_RATE).end())
    {
        oss << CONTACT_RATE << ":";
        for (const auto &x : _report_group_classes_.at(CONTACT_RATE))
        {
            oss << "\n    * " << x.first << ": " << x.second
                << " " << tools::percent(x.second, _num_contacts_);
        }
        oss << "\n";
    }

    // print only if there's something to print
    if (oss.str().size() > 0)
    {
        tools::print_stdout(oss.str());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Set accept rates to zero. */
void Monte_Carlo_Build::init_move_counts()
{
    // ... for parameters
    pCore->init_param_move_counts();

    // ... for augmented variables
    for (Household &this_house : house)
    {
        for (const std::string &augmID : augm_to_be_inferred)
        {
            this_house.reset_augm_accept_rate(augmID);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void Monte_Carlo_Build::init_diagnostics(bool is_power_post)
{
    ulong sample_size = pCore->setup_c_sample_size * pCore->setup_num_chains;

    pCore->dgn_count = 0;
    pCore->dgn_glob_loglike = std::vector<double>(sample_size);

    if (!is_power_post)
    {
        // diagnostics for all parameters
        pCore->init_param_diagnostics(sample_size);

        // diagnostics over households
        for (size_t H = 0; H < num_house; ++H)
        {
            house[H].init_diagnostics(sample_size);
        }
    }
}

/* Update Bayesian diagnostics */
void Monte_Carlo_Build::update_diagnostics(double glob_loglike, bool is_power_post)
{
    // diagnostics for global loglike
    pCore->dgn_glob_loglike[pCore->dgn_count++] = glob_loglike;

    if (!is_power_post)
    {
        // diagnostics for all parameters
        pCore->update_param_diagnostics();

        // diagnostics over households
        for (size_t H = 0; H < num_house; ++H)
        {
            house[H].update_diagnostics();
        }
    }
}

/* Disable all assert methods for user definitions. */
void Monte_Carlo_Build::disable_assert()
{
    pCore->disable_assert_param();

    // disable augmented variable assertion for all individuals
    for (Household &this_house : house)
    {
        for (Individual &this_indiv : this_house.indiv)
        {
            this_indiv.disable_assert();
        }
    }
}

/* House loglikelihood including augmented modl contribution. */
double Monte_Carlo_Build::compute_augm_loglike(const Household &this_house) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const std::string &augmID : augm_with_model_set)
    {
        for (int i : this_house.infectious_set)
        {
            hh_loglike += this_house.indiv[i].augm_model_logpdf(augmID);
        }
    }

    return hh_loglike;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef USE_RCPP

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame then run build. */
void Monte_Carlo_Build::import_augm(Rcpp::DataFrame augm)
{
    tools::rdataframe_to_map(augm, _augm_deftable_);
    build_augm();
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
    {
        _report_pathtofile_.erase("augm");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Monte_Carlo_Build : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame Monte_Carlo_Build::export_data() const
{
    using namespace Rcpp;

    // pre-allocate memories
    // initialize the dataframe as a list of named std::vectors
    CharacterVector df_index(num_indiv);
    CharacterVector df_house_index(num_indiv);
    IntegerVector df_status(num_indiv);
    std::map<std::string, NumericVector> df_dbl;

    for (const std::string &dataID : data_headers)
    {
        df_dbl.insert({dataID, NumericVector(num_indiv)});
    }
    for (const std::string &augmID : augm_full_set)
    {
        df_dbl.insert({augmID, NumericVector(num_indiv)});
    }

    size_t c = 0;
    for (const Household &this_house : house)
    {
        for (const Individual &this_indiv : this_house.indiv)
        {
            df_index[c] = this_indiv.get_ID();
            df_house_index[c] = this_house.get_ID();
            df_status[c] = this_indiv.get_status();

            for (const std::string &dataID : data_headers)
            {
                df_dbl.at(dataID)[c] = this_indiv.get_data(dataID);
            }
            for (const std::string &augmID : augm_full_set)
            {
                df_dbl.at(augmID)[c] = this_indiv.get_augm(augmID);
            }

            ++c;
        }
    }

    DataFrame output;
    output = DataFrame::create(_["id"] = df_index,
                               _["house_id"] = df_house_index,
                               _["status"] = df_status);
    for (const std::string &dataID : data_headers)
    {
        output.push_back(df_dbl.at(dataID), dataID);
    }
    for (const std::string &augmID : augm_full_set)
    {
        output.push_back(df_dbl.at(augmID), augmID);
    }

    return output;
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix Monte_Carlo_Build::export_loglike_matrix() const
{
    using namespace Rcpp;

    // FIXME : this should be prevented if is_power_post

    // Create a NumericMatrix with num_house columns
    ulong sample_size = pCore->setup_num_chains * pCore->setup_c_sample_size;
    NumericMatrix output(sample_size, num_house);

    // Assign each NumericVector to a column of the matrix
    for (size_t H = 0; H < num_house; ++H)
    {
        NumericVector column = wrap(house[H].get_dgn_loglike());
        output(_, H) = column;
    }

    return output;
}

#else

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : show full build  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Monte_Carlo_Build::show_full_build() const
{
    //----- SETUP -----//
    std::string add_file_info = UNDEFINED;
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
    {
        if (!tools::is_na(_report_pathtofile_.at("setup")))
        {
            add_file_info = "# import: " + _report_pathtofile_.at("setup");
        }
    }
    tools::print_decorator("Setup", add_file_info);
    pCore->show_setup();

    //----- PARAMETERS -----//
    add_file_info = UNDEFINED;
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
    {
        if (!tools::is_na(_report_pathtofile_.at("params")))
        {
            add_file_info = "# import: " + _report_pathtofile_.at("params");
        }
    }
    tools::print_decorator("Parameters", add_file_info);
    pCore->show_all_params();

    //----- DATABASE -----//
    add_file_info = UNDEFINED;
    if (_report_pathtofile_.find("data") != _report_pathtofile_.end())
    {
        if (!tools::is_na(_report_pathtofile_.at("data")))
        {
            add_file_info = "# import: " + _report_pathtofile_.at("data");
        }
    }
    tools::print_decorator("Database", add_file_info);
    show_summary();
    show_group_classes();

    //----- AUGMENTED -----//
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
    {
        if (!tools::is_na(_report_pathtofile_.at("augm")))
        {
            add_file_info = "# import: " + _report_pathtofile_.at("augm");
        }
    }
    tools::print_decorator("Augmented", add_file_info);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Monte_Carlo_Build : print augm variables  ///////////////////////////////////
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
            {
                augm_ofstream << C_DELIM << dataID;
            }

            // augmented variables headers
            for (const std::string &augmID : augm_full_set)
            {
                augm_ofstream << C_DELIM << augmID;
            }

            augm_ofstream << std::endl;

            //===== Print data and augmented variables =====//
            for (const Household &this_house : house)
            {
                for (const Individual &this_indiv : this_house.indiv)
                {
                    augm_ofstream << this_indiv.get_ID()
                                  << C_DELIM << this_house.get_ID()
                                  << C_DELIM << this_indiv.get_status();
                    for (const std::string &dataID : data_headers)
                    {
                        augm_ofstream << C_DELIM << this_indiv.get_data(dataID);
                    }
                    for (const std::string &augmID : augm_full_set)
                    {
                        if (this_indiv.get_status() == INFECTIOUS)
                        {
                            augm_ofstream << C_DELIM << this_indiv.get_augm(augmID);
                        }
                        else
                        {
                            augm_ofstream << C_DELIM << std::numeric_limits<double>::quiet_NaN();
                        }
                    }
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

#endif