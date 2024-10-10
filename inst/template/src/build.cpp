/*
 * build.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

using namespace basics;

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Param_Infer_Fixed_Move_Warning = false;
bool Non_Specified_Type_Warning = false;
bool Non_Specified_Value_Warning = false;
bool Fix_Init_Mult_Chains_Warning = false; // fixed paramter values when multiple chains
bool Unkown_Setup_Entry_Warning = false;   // wrong entries in setup file
bool Unknown_Household_Data_Structure_Warning = false;

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
    print_decorator("Setup", add_file_info);
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
    print_decorator("Parameters", add_file_info);
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
    print_decorator("Database", add_file_info);
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
    print_decorator("Augmented", add_file_info);
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
            raise_warning(oss.str() + ".", Unkown_Setup_Entry_Warning);
        }

        //===== ASSIGN SETUP CONFIGURATIONS =====//

        //===== NUMBER OF THREADS =====//
        if (_setup_deflist_.find("num_threads") != _setup_deflist_.end())
        {
            num_threads = std::stoi(_setup_deflist_.at("num_threads"));
        }

        //===== PRNG SEED =====//
        if ((_setup_deflist_.find("seed") != _setup_deflist_.end()) && (!is_na(_setup_deflist_.at("seed"))))
        {
            setup_seed = std::stoull(_setup_deflist_.at("seed"));
        }

        //===== NUMBER OF CHAINS =====//
        if (_setup_deflist_.find("num_chains") != _setup_deflist_.end())
        {
            setup_num_chains = std::stoul(_setup_deflist_.at("num_chains"));
            if (setup_num_iter < 1)
            {
                raise_error("Number of chains must be greater than 0.");
            }
        }

        //===== NUMBER OF ITERATIONS =====//
        if (_setup_deflist_.find("num_iter") != _setup_deflist_.end())
        {
            setup_num_iter = std::stoul(_setup_deflist_.at("num_iter"));
            if (setup_num_iter < 2)
            {
                raise_error("Number of iterations must be greater than 1.");
            }
        }

        //===== NUMBER OF (FINITE) ADAPTIVE ITERATIONS =====//
        if (_setup_deflist_.find("num_adapt") != _setup_deflist_.end())
        {
            setup_num_adapt = std::stoul(_setup_deflist_.at("num_adapt"));
        }

        //===== THINNING =====//
        if (_setup_deflist_.find("thinning") != _setup_deflist_.end())
        {
            setup_thinning = std::stoul(_setup_deflist_.at("thinning"));
            if (setup_thinning >= setup_num_iter)
            {
                raise_error("Thinning must be smaller than number of iterations. ");
            }
            else if (setup_thinning == 0)
            {
                raise_error("Thinning must be greater than 0.");
            }
        }

        //===== BURN-IN =====//
        if (_setup_deflist_.find("burn_in") != _setup_deflist_.end())
        {
            setup_burn_in = std::stoul(_setup_deflist_.at("burn_in"));
            if (setup_burn_in >= setup_num_iter)
            {
                raise_error("Burn-in must be smaller than number of iterations.");
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
            set_inver_temper(std::stod(_setup_deflist_.at("inv_temper")));
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

/* Set inverse temperature and its relative functional methods. */
void MCMC::set_inver_temper(double inv_T)
{
    if ((inv_T < 0.) || (inv_T > 1.))
    {
        std::ostringstream oss;
        oss << "Inverse temperature " << inv_T << " out of range [0, 1]. Check your setup.";
        raise_error(oss.str());
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
    if (_param_deftable_.empty())
    {
        // error if nothing was previously defined
        if (param_to_be_static.empty() && param_to_be_inferred.empty())
        {
            raise_error("Undefined parameter definitions.");
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
                        oss << "Undefined value for parameter to infer '"
                            << parID << "' without random initialization.";
                        raise_error(oss.str());
                    }
                }
                else
                {
                    value = std::stod(_param_deftable_.at("value")[p]);
                }

                // load range
                dictargs range = {};
                if (!is_na(_param_deftable_.at("min")[p]))
                {
                    range.emplace("min", std::stod(_param_deftable_.at("min")[p]));
                }
                if (!is_na(_param_deftable_.at("max")[p]))
                {
                    range.emplace("max", std::stod(_param_deftable_.at("max")[p]));
                }

                // assign parameter to the set of variables
                Parameter new_param(parID, value,
                                    _param_deftable_.at("prior")[p],
                                    merge_list_args(_param_deftable_.at("args")[p],
                                                    _param_deftable_.at("argsval")[p]),
                                    _param_deftable_.at("rule")[p],
                                    merge_list_args(_param_deftable_.at("rargs")[p],
                                                    _param_deftable_.at("rargsval")[p]),
                                    range);
                assign_param(parID, new_param);
            }
            else
            {
                // load value
                if (is_na(_param_deftable_.at("value")[p]))
                {
                    raise_error("Undefined value for immutable parameter '" + parID + "'.");
                }
                else
                {
                    value = std::stod(_param_deftable_.at("value")[p]);
                }

                // assign immutable parameter to the set of variables
                Parameter new_param(parID, value);
                assign_param(parID, new_param);
            }
        }
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
    if (datatable.empty())
    {
        // If nothing was previously defined
        if (house.empty())
        {
            raise_error("No database was found.");
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

        //===== Assess the presence of indexes =====//
        bool is_index = (std::find(_data_headers_.begin(), _data_headers_.end(), "id") != _data_headers_.end());
        bool is_house_index = (std::find(_data_headers_.begin(), _data_headers_.end(), "house_id") != _data_headers_.end());
        if (!is_house_index)
        {
            std::ostringstream oss;
            oss << "Impossible to associate a household structure: missing data column 'house_index'."
                << " If it was intentional, ignore this message.";
            raise_warning(oss.str(), Unknown_Household_Data_Structure_Warning);
        }

        //===== Safe Load of Households =====//
        std::unordered_map<int, int> aux_hh_idx = {};
        for (int i = 0; i < num_indiv; ++i)
        {
            dictargs args = {};
            for (const std::string &key : _data_headers_)
            {
                if (is_na(datatable.at(key)[i]))
                {
                    args.emplace(key, std::numeric_limits<double>::quiet_NaN());
                }
                else
                {
                    args.emplace(key, std::stod(datatable.at(key)[i]));
                }
            }
            if (!is_index)
            {
                args.emplace("id", i);
            }

            Individual this_indiv(args);

            //----- Assign and set rate per group -----//
            set_rate_per_group(this_indiv);

            // user provided household index (to set structure or not)
            int hh_idx = (is_house_index) ? this_indiv.get_data("house_id") : num_house;

            // if household index doesn't exist, add new household
            // note : this is redundant if no structure
            if (aux_hh_idx.find(hh_idx) == aux_hh_idx.end())
            {
                // map data index to household number
                aux_hh_idx.emplace(hh_idx, num_house);
                // create a new empty household
                house.emplace_back();
                // assign the data index
                house[num_house].assign_data("house_id", hh_idx);
                // set household number as its id
                house[num_house].set_id(num_house);
                // increment household number
                ++num_house;
            }
            // assign the individual the respective household
            house[aux_hh_idx.at(hh_idx)].assign_individual(this_indiv);
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
    assign_group_rate_to_individual(this_indiv);

    for (const auto &which : this_indiv.which_rate)
    {
        std::string rateID = which.first;
        if (rateID == CONTACT_RATE)
        {
            std::ostringstream oss;
            oss << "The group rate name '" << CONTACT_RATE << "' is reserved."
                << "Please choose a different one.";
            raise_error(oss.str());
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
            std::ostringstream oss;
            oss << "Rate '" << rateID << "' is not provided for individual "
                << this_indiv.get_data("id") << ".";
            raise_error(oss.str());
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
    // allocate contact patterns
    int hh_size = this_house.size();
    this_house.which_rate.resize(hh_size, std::vector<std::string>(hh_size, UNDEFINED));
    this_house.contact_rates.resize(hh_size, std::vector<const double *>(hh_size, &REFERENCE_VALUE));

    if (this_house.get_id() == 0)
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
            raise_error("Contact rate has not been provided for all households.");
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
                    oss << "Contact rate not provided for individuals " << this_house.indiv[k].get_data("id")
                        << " and " << this_house.indiv[i].get_data("id") << " within household " << this_house.get_data("house_id") << ".";
                    raise_error(oss.str());
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
        "rargs", "rargsval", "rargstype", "rfrom", "rfromtype",
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
                raise_error(oss.str());
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
                    raise_error(oss.str());
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
                    raise_error("Impossible to initialize '" + augmID + "' without a model.");
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
                        raise_error(oss.str());
                    }
                }
                else if (std::stoi(z) < INFECTIOUS)
                {
                    std::ostringstream oss;
                    oss << "Impossible to associate augmented '" << pair.first << "' for non infectious status.";
                    raise_error(oss.str());
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : assign augmented variables  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::assign_augm(Individual &this_indiv)
{
    for (size_t a = 0; a < _augm_deftable_.at("name").size(); ++a)
    {
        // assign augmented variable of the corresponding status
        std::string augm_status = _augm_deftable_.at("status")[a];

        if ((is_na(augm_status)) || (this_indiv.get_status() == std::stoi(augm_status)))
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
                               typedict_to_sdargs(_augm_deftable_.at("rargs")[a],
                                                  _augm_deftable_.at("rargsval")[a],
                                                  _augm_deftable_.at("rargstype")[a], this_indiv),
                               range,
                               type_to_double(_augm_deftable_.at("rfrom")[a],
                                              _augm_deftable_.at("rfromtype")[a], this_indiv));
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

    if (is_na(val) && is_na(typ))
    {
        output = INF_DOUBLE; // default value
    }
    else if (!is_na(val) && is_na(typ))
    {
        try
        {
            output = std::stod(val);
        }
        catch (const std::invalid_argument &e)
        {
            std::ostringstream oss;
            oss << "Value '" + val + "' can't be cast to scalar value. Provide the respective type.";
            raise_warning(oss.str(), Non_Specified_Type_Warning);
        }
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
            output = std::stod(val);
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
                raise_error("Unknown command type '" + val + "'.");
            }
        }
        else
        {
            raise_error("Unknown variable type specifier '" + typ + "'.");
        }
    }

    return output;
}

/* Convert a type-std::string to a double-pointer */
const double *MCMC::type_to_doubleptr(const std::string val, const std::string typ, const Individual &this_indiv)
{
    const double *output = nullptr;

    if (is_na(val) && is_na(typ))
    {
        output = nullptr; // default value
    }
    else if (!is_na(val) && is_na(typ))
    {
        std::ostringstream oss;
        oss << "Impossible to interpret default value of  '" << val << "'. Please provide a type.";
        raise_error(oss.str());
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
            output = assign_scalar_pointer(std::stod(val));
        }
        else
        {
            std::ostringstream oss;
            oss << "Unsopported type '" + typ + "' found loading '" << val << "' for pointer value.";
            raise_error(oss.str());
        }
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