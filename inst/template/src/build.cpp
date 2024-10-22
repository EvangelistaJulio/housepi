/*
 * build.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

using basics::is_na;
using basics::raise_warning;

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Param_Infer_Fixed_Move_Warning = false;
bool Non_Specified_Type_Warning = false;
bool Non_Specified_Value_Warning = false;
bool Fix_Init_Mult_Chains_Warning = false;
bool Unknown_Setup_Entry_Warning = false;
bool Unknown_Household_Data_Structure_Warning = false;
bool Unknown_Individual_Data_Warning = false;

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

void MCMC::show_full_build() const
{
    //----- SETUP -----//
    std::string add_file_info = UNDEFINED;
    if (report_pathtofile.find("setup") != report_pathtofile.end())
    {
        if (!is_na(report_pathtofile.at("setup")))
        {
            add_file_info = "# import: " + report_pathtofile.at("setup");
        }
    }
    basics::print_decorator("Setup", add_file_info);
    show_setup();

    //----- PARAMETERS -----//
    add_file_info = UNDEFINED;
    if (report_pathtofile.find("params") != report_pathtofile.end())
    {
        if (!is_na(report_pathtofile.at("params")))
        {
            add_file_info = "# import: " + report_pathtofile.at("params");
        }
    }
    basics::print_decorator("Parameters", add_file_info);
    show_params();

    //----- DATABASE -----//
    add_file_info = UNDEFINED;
    if (report_pathtofile.find("data") != report_pathtofile.end())
    {
        if (!is_na(report_pathtofile.at("data")))
        {
            add_file_info = "# import: " + report_pathtofile.at("data");
        }
    }
    basics::print_decorator("Database", add_file_info);
    show_summary();
    show_group_classes();

    //----- AUGMENTED -----//
    if (report_pathtofile.find("augm") != report_pathtofile.end())
    {
        if (!is_na(report_pathtofile.at("augm")))
        {
            add_file_info = "# import: " + report_pathtofile.at("augm");
        }
    }
    basics::print_decorator("Augmented", add_file_info);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Setup file and produce a definition list. */
void MCMC::read_setup(const std::string &pathtofile)
{
    read_csv_list(pathtofile, _setup_deflist_);
    build_setup();
    report_pathtofile["setup"] = pathtofile;
}

#ifdef USE_RCPP
void MCMC::import_setup(Rcpp::List setup)
{
    rlist_to_map(setup, _setup_deflist_);
    build_setup();
    if (report_pathtofile.find("setup") == report_pathtofile.end())
    {
        report_pathtofile.erase("setup");
    }
}
#endif

/* Assign setup. */
void MCMC::build_setup()
{
    try
    {
        // default number of threads
        unsigned num_threads = 1;

        if (!_setup_deflist_.empty())
        {
            //===== Check that all provided setups are correct =====//
            std::vector<std::string> file_entries = {"num_threads", "seed", "num_chains",
                                                     "num_iter", "num_adapt", "thinning",
                                                     "burn_in", "rnd_init_params",
                                                     "num_rungs", "inv_temper"};
            std::vector<std::string> wrong_entries = {};
            for (const auto &pair : _setup_deflist_)
            {
                if (find(file_entries.begin(), file_entries.end(), pair.first) == file_entries.end())
                {
                    wrong_entries.push_back(pair.first);
                }
            }
            if (wrong_entries.size() > 0)
            {
                std::ostringstream oss;
                oss << "Unknown setup option(s): '" << wrong_entries[0] << "'";
                for (size_t k = 1; k < wrong_entries.size(); ++k)
                {
                    oss << ", '" << wrong_entries[k] << "'";
                }
                oss << ". Available setup options are: '" << file_entries[0] << "'";
                for (size_t k = 1; k < file_entries.size(); ++k)
                {
                    oss << ", '" << file_entries[k] << "'";
                }
                raise_warning(oss.str() + ".", Unknown_Setup_Entry_Warning);
            }

            //===== ASSIGN SETUP CONFIGURATIONS =====//

            //===== NUMBER OF THREADS =====//
            if (_setup_deflist_.find("num_threads") != _setup_deflist_.end())
            {
                num_threads = basics::stoi(_setup_deflist_.at("num_threads"), "number of threads");
            }

            //===== PRNG SEED =====//
            if ((_setup_deflist_.find("seed") != _setup_deflist_.end()) && (!is_na(_setup_deflist_.at("seed"))))
            {
                setup_seed = basics::stoull(_setup_deflist_.at("seed"), "random seed");
            }

            //===== NUMBER OF CHAINS =====//
            if (_setup_deflist_.find("num_chains") != _setup_deflist_.end())
            {
                setup_num_chains = basics::stoul(_setup_deflist_.at("num_chains"), "number of chains");
                if (setup_num_iter < 1)
                {
                    throw std::out_of_range("Number of chains must be greater than 0.");
                }
            }

            //===== NUMBER OF ITERATIONS =====//
            if (_setup_deflist_.find("num_iter") != _setup_deflist_.end())
            {
                setup_num_iter = basics::stoul(_setup_deflist_.at("num_iter"), "number of iterations");
                if (setup_num_iter < 2)
                {
                    throw std::out_of_range("Number of iterations must be greater than 1.");
                }
            }

            //===== NUMBER OF (FINITE) ADAPTIVE ITERATIONS =====//
            if (_setup_deflist_.find("num_adapt") != _setup_deflist_.end())
            {
                setup_num_adapt = basics::stoul(_setup_deflist_.at("num_adapt"), "number of adaptive steps");
            }

            //===== THINNING =====//
            if (_setup_deflist_.find("thinning") != _setup_deflist_.end())
            {
                setup_thinning = basics::stoul(_setup_deflist_.at("thinning"), "thinning");
                if (setup_thinning >= setup_num_iter)
                {
                    throw std::out_of_range("Thinning must be smaller than number of iterations. ");
                }
                else if (setup_thinning == 0)
                {
                    throw std::out_of_range("Thinning must be greater than 0.");
                }
            }

            //===== BURN-IN =====//
            if (_setup_deflist_.find("burn_in") != _setup_deflist_.end())
            {
                setup_burn_in = basics::stoul(_setup_deflist_.at("burn_in"), "burn-in");
                if (setup_burn_in >= setup_num_iter)
                {
                    throw std::out_of_range("Burn-in must be smaller than number of iterations.");
                }
            }

            //===== RANDOM INITIALIZATION OF PARAMETERS =====//
            if (_setup_deflist_.find("rnd_init_params") != _setup_deflist_.end())
            {
                setup_rnd_init_params = stobool(_setup_deflist_.at("rnd_init_params"));
                if ((setup_num_chains > 1) && (setup_rnd_init_params == false))
                {
                    raise_warning("Unable to fix initial parameters when running multiple chains.",
                                  Fix_Init_Mult_Chains_Warning);
                    setup_rnd_init_params = true;
                }
            }

            //===== INVERSE TEMPERATURE =====//
            if (_setup_deflist_.find("inv_temper") != _setup_deflist_.end())
            {
                double inv_T = basics::stod(_setup_deflist_.at("inv_temper"), "inverse temperature");
                if ((inv_T < 0.) || (inv_T > 1.))
                {
                    std::ostringstream oss;
                    oss << "Inverse temperature " << inv_T << " out of range [0, 1].";
                    throw std::out_of_range(oss.str());
                }
                else
                {
                    set_inver_temper(inv_T);
                }
            }
        }

        //===== Define Number of Threads (Monte_Carlo_Core) =====//
        set_num_threads(num_threads);

        //===== Pseudo-Random Number Generator Seeding =====//
        set_prng_seed(setup_seed);

        //===== Number of MCMC sweeps per single chain =====//
        // note : Euclidean divisions!
        _thin_sweeps_ = setup_num_adapt / setup_thinning;
        _thin_sweeps_ += setup_burn_in / setup_thinning;
        _thin_sweeps_ += (setup_num_iter - setup_burn_in) / setup_thinning;
        _thin_sweeps_ *= setup_num_chains;
    }
    catch (const std::exception &e)
    {
        basics::raise_error(e.what());
    }
}

/* Set inverse temperature and its relative functional methods. */
void MCMC::set_inver_temper(double inv_T)
{
    if ((inv_T < 0.) || (inv_T > 1.))
    {
        std::ostringstream oss;
        oss << "Inverse temperature " << inv_T << " out of range [0, 1].";
        throw std::range_error(oss.str());
    }
    else if (inv_T == 1.0)
    {
        setup_inv_T = 1.0;
        _param_metropolis_method_ = &MCMC::_param_metropolis_1_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_1_;
    }
    else if (inv_T == 0.0)
    {
        setup_inv_T = 0.0;
        _param_metropolis_method_ = &MCMC::_param_metropolis_0_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_0_;
    }
    else
    {
        setup_inv_T = inv_T;
        _param_metropolis_method_ = &MCMC::_param_metropolis_T_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_T_;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : build parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read parameter definitions and run assignment. */
void MCMC::read_params(const std::string &pathtofile)
{
    read_csv_table(pathtofile, _param_deftable_, param_file_headers);
    build_params();
    report_pathtofile["params"] = pathtofile;
}

#ifdef USE_RCPP
void MCMC::import_params(Rcpp::DataFrame params)
{
    rdataframe_to_map(params, _param_deftable_);
    build_params();
    if (report_pathtofile.find("params") == report_pathtofile.end())
    {
        report_pathtofile.erase("params");
    }
}
#endif

/* Assign parameter definitions to dictionary. */
void MCMC::build_params()
{
    try
    {
        if (_param_deftable_.empty())
        {
            // error if nothing was previously defined
            if (param_to_be_static.empty() && param_to_be_inferred.empty())
            {
                throw std::invalid_argument("Undefined parameter definitions.");
            }
        }
        else
        {
            //===== Read Parameter Definitions =====//
            size_t table_length = _param_deftable_.at("infer").size();

            //===== Empty current list of parameters =====//
            drop_all_params();

            //===== Load parameters from the row with options =====//
            for (size_t p = 0; p < table_length; ++p)
            {
                std::string parID = _param_deftable_.at("name")[p];
                double value = INF_DOUBLE;

                // add entry to shared parameter dictionary
                if (stobool(_param_deftable_.at("infer")[p]))
                {
                    // load value
                    if (is_na(_param_deftable_.at("value")[p]))
                    {
                        if (!setup_rnd_init_params)
                        {
                            std::ostringstream oss;
                            oss << "Undefined value without random initialization for parameter to infer '" << parID << "'.";
                            throw std::invalid_argument(oss.str());
                        }
                    }
                    else
                    {
                        value = basics::stod(_param_deftable_.at("value")[p], "parameter");
                    }

                    // load range
                    dictargs range = {};
                    if (!is_na(_param_deftable_.at("min")[p]))
                    {
                        range.emplace("min", basics::stod(_param_deftable_.at("min")[p], "parameter minimum"));
                    }
                    if (!is_na(_param_deftable_.at("max")[p]))
                    {
                        range.emplace("max", basics::stod(_param_deftable_.at("max")[p], "parameter maximum"));
                    }

                    // assign parameter to the set of variables
                    Parameter new_param(parID, value,
                                        _param_deftable_.at("prior")[p],
                                        merge_list_args(_param_deftable_.at("args")[p],
                                                        _param_deftable_.at("argsval")[p]),
                                        _param_deftable_.at("rule")[p],
                                        merge_list_args(_param_deftable_.at("rule_args")[p],
                                                        _param_deftable_.at("rule_argsval")[p]),
                                        range);
                    assign_param(parID, new_param);
                }
                else
                {
                    // load value
                    if (is_na(_param_deftable_.at("value")[p]))
                    {
                        throw std::invalid_argument("Undefined value for immutable parameter '" + parID + "'.");
                    }
                    else
                    {
                        value = basics::stod(_param_deftable_.at("value")[p], "parameter");
                    }

                    // assign immutable parameter to the set of variables
                    Parameter new_param(parID, value);
                    assign_param(parID, new_param);
                }
            }
        }
    }
    catch (const std::exception &e)
    {
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : build data  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read data file and then run assignment. */
void MCMC::read_data(const std::string &pathtofile)
{
    std::map<std::string, std::vector<std::string>> datatable = {};
    read_csv_table(pathtofile, datatable);
    build_data(datatable);

    report_pathtofile["data"] = pathtofile;
}

#ifdef USE_RCPP
/* Import data from R DataFrame and then run assignment. */
void MCMC::import_data(Rcpp::DataFrame data)
{
    std::map<std::string, std::vector<std::string>> datatable = {};
    rdataframe_to_map(data, datatable);
    build_data(datatable);
    if (report_pathtofile.find("data") == report_pathtofile.end())
    {
        report_pathtofile.erase("data");
    }
}
#endif

/* Assign data after pre-processing. */
void MCMC::build_data(const std::map<std::string, std::vector<std::string>> &datatable)
{
    try
    {
        if (datatable.empty())
        {
            // If nothing was previously defined
            if (house.empty())
            {
                throw std::length_error("No database was found.");
            }
        }
        else
        {
            //=====  Drop all Previous Data and Prepare for New  =====//
            _data_headers_.clear();
            for (const auto &x : datatable)
            {
                _data_headers_.push_back(x.first);
            }
            num_indiv = datatable.at(_data_headers_[0]).size();
            num_house = 0;
            house.clear();
            report_group_classes.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(_data_headers_.begin(), _data_headers_.end(), "house_id") == _data_headers_.end())
            {
                std::ostringstream oss;
                oss << "Impossible to associate a household structure: missing data column 'house_id'."
                    << " If it was intentional, ignore this message.";
                raise_warning(oss.str(), Unknown_Household_Data_Structure_Warning);
            }

            //===== Safe Load of Households =====//
            std::unordered_map<std::string, int> aux_hh_idx = {};
            for (int i = 0; i < num_indiv; ++i)
            {
                //----- Read individual data -----//
                std::string this_indiv_id = std::to_string(i);
                std::string this_house_id = std::to_string(num_house);
                int this_status = MISSING_STATUS;
                dictargs args = {};
                for (const std::string &key : _data_headers_)
                {
                    if (key == "id")
                    {
                        if (is_na(datatable.at("id")[i]))
                        {
                            throw std::invalid_argument("Missing id for individual in row n." + std::to_string(i+2) + ".");
                        }
                        else
                        {
                            this_indiv_id = datatable.at("id")[i];
                        }
                    }
                    else if (key == "house_id")
                    {
                        if (is_na(datatable.at("house_id")[i]))
                        {
                            throw std::invalid_argument("Missing house_id for individual in row n." + std::to_string(i+2) + ".");
                        }
                        else
                        {
                            this_house_id = datatable.at("house_id")[i];
                        }
                    }
                    else if (key == "status")
                    {
                        // NA status is a missing status
                        if (!is_na(datatable.at("status")[i]))
                        {
                            this_status = basics::stoul(datatable.at("status")[i], "individual status");
                        }
                    }
                    else if (is_na(datatable.at(key)[i]))
                    {
                        args.emplace(key, std::numeric_limits<double>::quiet_NaN());
                        std::ostringstream oss;
                        oss << "Found NA value in database. If intentional, ignore this warning.";
                        raise_warning(oss.str(), Unknown_Individual_Data_Warning);
                    }
                    else
                    {
                        args.emplace(key, basics::stod(datatable.at(key)[i], "data '" + key + "'"));
                    }
                }
                // create individual instance
                Individual this_indiv(this_indiv_id, this_status, args);

                //----- Assign and set rate per group -----//
                set_rate_per_group(this_indiv);

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (aux_hh_idx.find(this_house_id) == aux_hh_idx.end())
                {
                    // map data index to household number
                    aux_hh_idx.emplace(this_house_id, num_house);
                    // create a new empty household
                    house.emplace_back();
                    // assign the user data id
                    house[num_house].set_ID(this_house_id);
                    // set household number as its index
                    house[num_house].set_index(num_house);
                    // increment household number
                    ++num_house;
                }
                // assign the individual the respective household
                house[aux_hh_idx.at(this_house_id)].assign_individual(this_indiv);
            }

            //===== Set Contact Rates =====//
            for (Household &this_house : house)
            {
                set_rate_per_contact(this_house);
            }

            //===== Update Report for Data =====//
            update_report_data();

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const std::string &col : {"house_id", "id", "status"})
            {
                _data_headers_.erase(
                    std::remove(_data_headers_.begin(), _data_headers_.end(), col),
                    _data_headers_.end());
            }
        }
    }
    catch (const std::exception &e)
    {
        basics::raise_error(e.what());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Use the loaded data to fill-in the respective report. */
void MCMC::update_report_data()
{
    report_data = {{"n.cases", 0},
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

        report_data["n.cases"] += num_infec;
        report_data["n.miss.stat"] += this_house.num_missing_status;

        report_data["no.HH.case"] += (num_infec == 0);
        report_data["lst.1.case"] += (num_infec > 0);
        report_data["lst.2.case"] += (num_infec > 1);

        // max number of house cases
        if (report_data["max.n.HH.cases"] < num_infec)
        {
            report_data["n.HH.with.max.cases"] = 1;
            report_data["max.n.HH.cases"] = num_infec;
        }
        else if (report_data["max.n.HH.cases"] == num_infec)
        {
            report_data["n.HH.with.max.cases"] += 1;
        }

        // minimum house number of members
        if (report_data["min.HH.size"] > this_house_size)
        {
            report_data["n.HH.with.min.size"] = 1;
            report_data["min.HH.size"] = this_house_size;
        }
        else if (report_data["min.HH.size"] == this_house_size)
        {
            report_data["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (report_data["max.HH.size"] < this_house_size)
        {
            report_data["n.HH.with.max.size"] = 1;
            report_data["max.HH.size"] = this_house_size;
        }
        else if (report_data["max.HH.size"] == this_house_size)
        {
            report_data["n.HH.with.max.size"] += 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : set rate per group  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign group rate to an individual  */
void MCMC::assign_group_rate_to_individual(Individual &this_indiv)
{
    // this_indiv.assign_group_rate(EMPTY_RATE, EMPTY_RATE);
    UNUSED(this_indiv);
}

/* Update model report and set group rates associate to new individual. */
void MCMC::set_rate_per_group(Individual &this_indiv)
{
    try
    {
        assign_group_rate_to_individual(this_indiv);

        for (const auto &which : this_indiv.which_rate)
        {
            std::string rateID = which.first;
            if (rateID == CONTACT_RATE)
            {
                throw std::invalid_argument("The group rate name '" + CONTACT_RATE + "' is reserved. Choose a different one.");
            }

            std::string parID = which.second;

            //=====initialize model report=====//
            // if first individual, initialize rateID-report
            report_group_classes.emplace(rateID, std::map<std::string, int>());
            report_group_classes.at(rateID).emplace("[reference]", 0);

            //=====set rates and update report=====//
            if (is_na(parID))
            {
                // error : undefined rate
                std::invalid_argument("Rate '" + rateID + "' is not provided for individual " + this_indiv.get_ID() + ".");
            }
            else if (parID == REFERENCE)
            {
                // rate set to reference
                this_indiv.set_rate(rateID, get_reference());
                report_group_classes.at(rateID).at("[reference]") += 1;
            }
            else if (parID != REFERENCE)
            {
                // rate associated to parameter
                assert_param(parID);
                this_indiv.set_rate(rateID, get_param_pointer(parID));
                if (report_group_classes.at(rateID).find(parID) == report_group_classes.at(rateID).end())
                {
                    report_group_classes.at(rateID).emplace(parID, 1);
                }
                else
                {
                    report_group_classes.at(rateID).at(parID) += 1;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        basics::raise_error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
    this_house.assign_rate_per_contact(from, to, EMPTY_RATE);
}

/* Set contact rate according to user assignement. */
void MCMC::set_rate_per_contact(Household &this_house)
{
    try
    {
        // allocate contact patterns
        int hh_size = this_house.size();
        this_house.which_rate.resize(hh_size, std::vector<std::string>(hh_size, UNDEFINED));
        this_house.contact_rates.resize(hh_size, std::vector<const double *>(hh_size, &REFERENCE_VALUE));

        if (this_house.get_index() == 0)
        {
            // if first household, initialize contact rate report
            report_group_classes.emplace(CONTACT_RATE, std::map<std::string, int>());
            report_group_classes.at(CONTACT_RATE).emplace("[reference]", 0);
        }
        else
        {
            // all household must have a contact rate
            if (report_group_classes.find(CONTACT_RATE) == report_group_classes.end())
            {
                throw std::invalid_argument("Contact rate has not been provided for all households.");
            }
        }

        // loop over all couples
        for (int k = 0; k < hh_size; ++k) // from
        {
            for (int i = 0; i < hh_size; ++i) // to
            {
                // note : from 'k' --> to 'i'
                if (k != i)
                {
                    assign_contact_rate_to_household(this_house, k, i);
                    std::string parID = this_house.which_rate[k][i];

                    if (is_na(parID))
                    {
                        std::ostringstream oss;
                        oss << "Contact rate not provided for individuals " << this_house.indiv[k].get_ID()
                            << " and " << this_house.indiv[i].get_ID() << " within household " << this_house.get_ID() << ".";
                        throw std::invalid_argument(oss.str());
                    }
                    else if (parID == REFERENCE)
                    {
                        this_house.set_contact_rate(k, i, get_reference()); // reference rate
                        report_group_classes.at(CONTACT_RATE).at("[reference]") += 1;
                    }
                    else if (parID != EMPTY_RATE)
                    {
                        assert_param(parID);
                        this_house.set_contact_rate(k, i, get_param_pointer(parID)); // parameter rate
                        if (report_group_classes.at(CONTACT_RATE).find(parID) == report_group_classes.at(CONTACT_RATE).end())
                        {
                            report_group_classes.at(CONTACT_RATE).emplace(parID, 1);
                        }
                        else
                        {
                            report_group_classes.at(CONTACT_RATE).at(parID) += 1;
                        }
                    }

                    // count the total number of pairs
                    ++_n_contacts_;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        basics::raise_error(e.what());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : build augmented  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void MCMC::read_augm(const std::string &pathtofile)
{
    const std::vector<std::string> augm_file_headers = {
        "name", "status", "value", "valuetype",
        "model", "args", "argsval", "argstype", "rule",
        "rule_args", "rule_argsval", "rule_argstype", "rule_from", "rule_fromtype",
        "min", "mintype", "max", "maxtype", "init_rand"};

    read_csv_table(pathtofile, _augm_deftable_, augm_file_headers);
    build_augm();
    report_pathtofile["augm"] = pathtofile;
}

#ifdef USE_RCPP
/* Import data from R DataFrame and then run assignment. */
void MCMC::import_augm(Rcpp::DataFrame augm)
{
    rdataframe_to_map(augm, _augm_deftable_);
    build_augm();
    if (report_pathtofile.find("augm") == report_pathtofile.end())
    {
        report_pathtofile.erase("augm");
    }
}
#endif

/* Assign augmented variables to all individuals. */
void MCMC::build_augm()
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

                if (is_na(augmID))
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
                if (!is_na(_augm_deftable_.at("init_rand")[a]))
                {
                    is_init_rand = stobool(_augm_deftable_.at("init_rand")[a]);
                }

                // list of augmented variables to be inferred
                if (!is_na(_augm_deftable_.at("rule")[a]))
                {
                    augm_to_be_inferred.emplace(augmID);
                }

                // list of augmented variables with a model
                // note: no contribution of constant variables to Metropolis
                if (is_na(_augm_deftable_.at("model")[a]))
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
                    if (is_na(z))
                    {
                        if (mylist.size() > 1)
                        {
                            std::ostringstream oss;
                            oss << "Forbidden default and status definitions for augmented variable '"
                                << pair.first << "'.";
                            throw std::invalid_argument(oss.str());
                        }
                    }
                    else if (basics::stoi(z, "augmented status") < INFECTIOUS)
                    {
                        std::ostringstream oss;
                        oss << "Impossible to associate augmented '" << pair.first << "' for non infectious status.";
                        throw std::invalid_argument(oss.str());
                    }
                }
            }

            //===== Set internal state variable =====//
            _is_move_augm_ = !augm_to_be_inferred.empty();

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
        basics::raise_error(e.what());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : assign augmented variables  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::assign_augm(Individual &this_indiv)
{
    try
    {
        for (size_t a = 0; a < _augm_deftable_.at("name").size(); ++a)
        {
            // assign augmented variable of the corresponding status
            std::string augm_status = _augm_deftable_.at("status")[a];

            if ((is_na(augm_status)) || (this_indiv.get_status() == basics::stoi(augm_status, "augmented status")))
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
                if (stobool(_augm_deftable_.at("init_rand")[a]))
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
        basics::raise_error(e.what());
    }
}

/* Remove augmented definition of individual. */
void MCMC::clear_augm(Individual &this_indiv)
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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double MCMC::type_to_double(const std::string val, const std::string typ, const Individual &this_indiv)
{
    double output = INF_DOUBLE;

    try
    {
        if (is_na(val) && is_na(typ))
        {
            output = INF_DOUBLE; // default value
        }
        else if (!is_na(val) && is_na(typ))
        {
            output = basics::stod(val, "argument without type");
        }
        else if (is_na(val) && !is_na(typ))
        {
            raise_warning("Ignoring value type without a paired value.",
                          Non_Specified_Value_Warning);

            output = INF_DOUBLE; // FIXME default value
        }
        else
        {
            if (typ == TYPE_SCALAR)
            {
                output = basics::stod(val, "argument with scalar type");
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
                    raise_warning(oss.str(), Param_Infer_Fixed_Move_Warning);
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
        basics::raise_error(e.what());
    }

    return output;
}

/* Convert a type-std::string to a double-pointer */
const double *MCMC::type_to_doubleptr(const std::string val, const std::string typ, const Individual &this_indiv)
{
    const double *output = nullptr;

    try
    {
        if (is_na(val) && is_na(typ))
        {
            output = nullptr; // default value
        }
        else if (!is_na(val) && is_na(typ))
        {
            throw std::invalid_argument("Impossible to interpret default value of  '" + val + "'. Please provide a type.");
        }
        else if (is_na(val) && !is_na(typ))
        {
            raise_warning("Found a type without a paired value. Falling back to default.",
                          Non_Specified_Value_Warning);

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
                output = assign_scalar_pointer(basics::stod(val, "argument with scalar type"));
            }
            else
            {
                throw std::invalid_argument("Unsopported type '" + typ + "' found loading '" + val + "' for pointer value.");
            }
        }
    }
    catch (const std::exception &e)
    {
        basics::raise_error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
dictargs MCMC::typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv)
{
    std::vector<std::string> keys = parse_list_args(tkn_key);
    std::vector<std::string> vals = parse_list_args(tkn_val);
    std::vector<std::string> typs = parse_list_args(tkn_typ);

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
dictargs_ptr MCMC::typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv)
{
    std::vector<std::string> keys = parse_list_args(tkn_key);
    std::vector<std::string> vals = parse_list_args(tkn_val);
    std::vector<std::string> typs = parse_list_args(tkn_typ);

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