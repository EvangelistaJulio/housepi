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
const string TYPE_PARAM = "P";   // Parameter type declaration
const string TYPE_DATA = "D";    // Data type declaration
const string TYPE_SCALAR = "S";  // Scalar type declaration
const string TYPE_COMMAND = "C"; // Command type declaration

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

void MCMC::show_build_details()
{
    //----- SETUP -----//
    string add_file_info = UNDEFINED;
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
    show_parameters();

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
void MCMC::read_setup(const string &pathtofile)
{
    read_csv_list(pathtofile, _setup_deflist_);
    assign_setup();
    report_pathtofile["setup"] = pathtofile;
}

#ifdef USE_RCPP
void MCMC::import_setup(Rcpp::List setup)
{
    rlist_to_map(setup, _setup_deflist_);
    assign_setup();
    if (report_pathtofile.find("setup") == report_pathtofile.end())
    {
        report_pathtofile.erase("setup");
    }
}
#endif

/* Assign setup. */
void MCMC::assign_setup()
{
    // default number of threads
    unsigned num_threads = 1;

    if (!_setup_deflist_.empty())
    {
        //===== Check that all provided setups are correct =====//
        vector<string> file_entries = {"num_threads", "seed", "num_chains",
                                       "num_iter", "num_adapt", "thinning",
                                       "burn_in", "rnd_init_params",
                                       "num_rungs", "inv_temper"};
        vector<string> wrong_entries = {};
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
            setup_inv_T = std::stod(_setup_deflist_.at("inv_temper"));
            if ((setup_inv_T < 0.) || (setup_inv_T > 1.))
            {
                std::ostringstream oss;
                oss << "Inverse temperature " << setup_inv_T << " out of range [0, 1]. Check your setup.";
                raise_error(oss.str());
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
void MCMC::read_params(const string &pathtofile)
{
    read_csv_table(pathtofile, _param_deftable_, param_file_headers);
    assign_parameters();
    report_pathtofile["params"] = pathtofile;
}

#ifdef USE_RCPP
void MCMC::import_params(Rcpp::DataFrame params)
{
    rdataframe_to_map(params, _param_deftable_);
    assign_parameters();
    if (report_pathtofile.find("params") == report_pathtofile.end())
    {
        report_pathtofile.erase("params");
    }
}
#endif

/* Assign parameter definitions to dictionary. */
void MCMC::assign_parameters()
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
            // add entry to shared parameter dictionary
            if (stobool(_param_deftable_.at("infer")[p]))
            {
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
                assign_param(
                    _param_deftable_.at("name")[p],
                    std::stod(_param_deftable_.at("value")[p]),
                    _param_deftable_.at("prior")[p],
                    merge_list_args(_param_deftable_.at("args")[p], _param_deftable_.at("argsval")[p]),
                    _param_deftable_.at("rule")[p],
                    merge_list_args(_param_deftable_.at("rargs")[p], _param_deftable_.at("rargsval")[p]),
                    range);
            }
            else
            {
                assign_param(_param_deftable_.at("name")[p], stod(_param_deftable_.at("value")[p]));
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
void MCMC::read_data(const string &pathtofile)
{
    std::map<string, vector<string>> datatable = {};
    read_csv_table(pathtofile, datatable);
    assign_data(datatable);

    report_pathtofile["data"] = pathtofile;
}

#ifdef USE_RCPP
/* Import data from R DataFrame and then run assignment. */
void MCMC::import_data(Rcpp::DataFrame data)
{
    std::map<string, vector<string>> datatable = {};
    rdataframe_to_map(data, datatable);
    assign_data(datatable);
    if (report_pathtofile.find("data") == report_pathtofile.end())
    {
        report_pathtofile.erase("data");
    }
}
#endif

/* Assign data after pre-processing. */
void MCMC::assign_data(std::map<string, vector<string>> &datatable)
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
        bool is_index = (std::find(_data_headers_.begin(), _data_headers_.end(), "index") != _data_headers_.end());
        bool is_house_index = (std::find(_data_headers_.begin(), _data_headers_.end(), "house_index") != _data_headers_.end());
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
            for (const string &key : _data_headers_)
            {
                args.emplace(key, stod(datatable.at(key)[i]));
            }
            if (!is_index)
            {
                args.emplace("index", i);
            }

            Individual this_indiv(args);

            //----- Assign and set rate per group -----//
            set_rate_per_group(this_indiv);

            // user provided structure or no structure
            int hh_idx = (is_house_index) ? this_indiv.get_data("house_index") : num_house;

            // if household index doesn't exist, add new household
            // note : this is redundant if no structure
            if (aux_hh_idx.find(hh_idx) == aux_hh_idx.end())
            {
                house.emplace_back();                  // add a new empty household
                house[num_house].set_index(hh_idx);    // associate its index
                aux_hh_idx.emplace(hh_idx, num_house); // map index to household number
                ++num_house;                           // increment household number
            }
            house[aux_hh_idx.at(hh_idx)].assign_individual(this_indiv);
        }

        //===== Set Contact Rates =====//
        for (int h = 0; h < num_house; ++h)
        {
            set_rate_per_contact(h);
        }

        //===== Update Report for Data =====//
        update_report_data();

        //===== Clean up data entries for output files =====//
        // eventually remove standard headers to reorder them at the end
        for (const string &col : {"house_index", "index", "status"})
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
    for (int h = 0; h < num_house; ++h)
    {
        int num_infec = house[h].num_infectious;

        report_data["n.cases"] += num_infec;
        report_data["n.miss.stat"] += house[h].num_missing_status;

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
        if (report_data["min.HH.size"] > house[h].house_size)
        {
            report_data["n.HH.with.min.size"] = 1;
            report_data["min.HH.size"] = house[h].house_size;
        }
        else if (report_data["min.HH.size"] == house[h].house_size)
        {
            report_data["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (report_data["max.HH.size"] < house[h].house_size)
        {
            report_data["n.HH.with.max.size"] = 1;
            report_data["max.HH.size"] = house[h].house_size;
        }
        else if (report_data["max.HH.size"] == house[h].house_size)
        {
            report_data["n.HH.with.max.size"] += 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : set rate per group  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign group rate to an individual  */
void MCMC::assign_rate_per_group_to_individual(Individual &this_indiv)
{
    // this_indiv.assign_rate_per_group(EMPTY_RATE, EMPTY_RATE);
    UNUSED(this_indiv);
}

/* Update model report and set group rates associate to new individual. */
void MCMC::set_rate_per_group(Individual &this_indiv)
{
    assign_rate_per_group_to_individual(this_indiv);

    for (const auto &which : this_indiv.which_rate)
    {
        string rateID = which.first;
        if (rateID == CONTACT_RATE)
        {
            std::ostringstream oss;
            oss << "The group rate name '" << CONTACT_RATE << "' is reserved."
                << "Please choose a different one.";
            raise_error(oss.str());
        }

        string parID = which.second;

        //=====initialize model report=====//
        // if first individual, initialize rateID-report
        report_group_classes.emplace(rateID, std::map<string, int>());
        report_group_classes.at(rateID).emplace("[reference]", 0);

        //=====set rates and update report=====//
        if (is_na(parID))
        {
            // error : undefined rate
            std::ostringstream oss;
            oss << "Rate '" << rateID << "' is not provided for individual " << this_indiv.get_index() << ".";
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
void MCMC::assign_rate_per_contact_to_household(int h, int from, int to)
{
    house[h].assign_rate_per_contact(from, to, EMPTY_RATE);
}

/* Set contact rate according to user assignement. */
void MCMC::set_rate_per_contact(int h)
{
    // allocate contact patterns
    int s = house[h].house_size;
    house[h].which_rate.resize(s, vector<string>(s, UNDEFINED));
    house[h].contact_rates.resize(s, vector<const double *>(s, &REFERENCE_VALUE));

    if (h == 0)
    {
        // if first household, initialize contact rate report
        report_group_classes.emplace(CONTACT_RATE, std::map<string, int>());
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

    for (int k = 0; k < house[h].house_size; ++k)
    {
        for (int i = 0; i < house[h].house_size; ++i)
        {
            // note : from 'k' --> to 'i'
            if (k != i)
            {
                assign_rate_per_contact_to_household(h, k, i);
                string parID = house[h].which_rate[k][i];

                if (is_na(parID))
                {
                    std::ostringstream oss;
                    oss << "Contact rate not provided for individuals " << house[h].indiv[k].get_index()
                        << " and " << house[h].indiv[i].get_index() << " within household " << house[h].get_index() << ".";
                    raise_error(oss.str());
                }
                else if (parID == REFERENCE)
                {
                    house[h].set_rate(k, i, get_reference()); // reference rate
                    report_group_classes.at(CONTACT_RATE).at("[reference]") += 1;
                }
                else if (parID != EMPTY_RATE)
                {
                    assert_param(parID);
                    house[h].set_rate(k, i, get_param_pointer(parID)); // parameter rate
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
void MCMC::read_augm(const string &pathtofile)
{
    const vector<string> augm_file_headers = {
        "name", "value", "valuetype",
        "model", "args", "argsval", "argstype", "rule",
        "rargs", "rargsval", "rargstype", "rfrom", "rfromtype",
        "min", "mintype", "max", "maxtype", "init_rand"};

    read_csv_table(pathtofile, _augm_deftable_, augm_file_headers);
    assign_augmented();
    report_pathtofile["augm"] = pathtofile;
}

#ifdef USE_RCPP
/* Import data from R DataFrame and then run assignment. */
void MCMC::import_augm(Rcpp::DataFrame augm)
{
    rdataframe_to_map(augm, _augm_deftable_);
    assign_augmented();
    if (report_pathtofile.find("augm") == report_pathtofile.end())
    {
        report_pathtofile.erase("augm");
    }
}
#endif

/* Assign augmented variables to all individuals. */
void MCMC::assign_augmented()
{
    if (!_augm_deftable_.empty())
    {
        size_t table_length = _augm_deftable_.at("name").size();

        //===== Drop previous definitions =====//
        augm_full_set.clear();
        augm_to_be_inferred.clear();
        augm_with_model_set.clear();
        augm_rand_init_set.clear();

        //===== Augmented variables lists of usage =====//
        for (size_t a = 0; a < table_length; ++a)
        {
            string augmID = _augm_deftable_.at("name")[a];
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
            }
            else
            {
                std::ostringstream oss;
                oss << "Ambigous definitions for augmented variable '" << augmID << "'."
                    << " Check your augmented variable file.";
                raise_error(oss.str());
            }

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
                    raise_error("Impossible to randomly initialize '" + augmID + "' without a model.");
                }
            }
            else
            {
                augm_with_model_set.emplace(augmID);

                if (is_init_rand)
                {
                    augm_rand_init_set.emplace(augmID);
                }
            }
        }
        _is_data_augm_ = !augm_to_be_inferred.empty();

        //===== Eventually assign augmented variables to each individual =====//
        if (table_length > 0)
        {
            for (int h = 0; h < num_house; ++h)
            {
                for (int i = 0; i < house[h].house_size; ++i)
                {
                    for (size_t a = 0; a < table_length; ++a)
                    {
                        // load range
                        dictargs range = {};
                        double aux_min = type_to_double(_augm_deftable_.at("min")[a],
                                                        _augm_deftable_.at("mintype")[a], h, i);
                        if (aux_min != INF_DOUBLE)
                        {
                            range.emplace("min", aux_min);
                        }

                        double aux_max = type_to_double(_augm_deftable_.at("max")[a],
                                                        _augm_deftable_.at("maxtype")[a], h, i);
                        if (aux_max != INF_DOUBLE)
                        {
                            range.emplace("max", aux_max);
                        }

                        // assign augmented variable
                        house[h].indiv[i].assign_augm(
                            _augm_deftable_.at("name")[a],
                            type_to_double(_augm_deftable_.at("value")[a],
                                           _augm_deftable_.at("valuetype")[a], h, i),
                            _augm_deftable_.at("model")[a],
                            typedict_to_sdargs_ext(_augm_deftable_.at("args")[a],
                                                   _augm_deftable_.at("argsval")[a],
                                                   _augm_deftable_.at("argstype")[a], h, i),
                            _augm_deftable_.at("rule")[a],
                            typedict_to_sdargs(_augm_deftable_.at("rargs")[a],
                                               _augm_deftable_.at("rargsval")[a],
                                               _augm_deftable_.at("rargstype")[a], h, i),
                            range,
                            type_to_double(_augm_deftable_.at("rfrom")[a],
                                           _augm_deftable_.at("rfromtype")[a], h, i));
                    }
                }
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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed string to a double */
double MCMC::type_to_double(const string val, const string typ, int h, int i)
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
            output = get_data(h, i, val);
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

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(const string val, const string typ, int h, int i)
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
            output = house[h].indiv[i].get_data_pointer(val);
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

/* Convert a string with type-dictionary to string-to-double */
dictargs MCMC::typedict_to_sdargs(const string tkn_key, const string tkn_val, const string tkn_typ, int h, int i)
{
    vector<string> keys = parse_list_args(tkn_key);
    vector<string> vals = parse_list_args(tkn_val);
    vector<string> typs = parse_list_args(tkn_typ);

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
            output.emplace(keys[x], type_to_double(vals[x], typs[x], h, i));
        }
    }

    return output;
}

/* Convert a string with type-dictionary to string-to-double-pointer */
dictargs_ptr MCMC::typedict_to_sdargs_ext(const string tkn_key, const string tkn_val, const string tkn_typ, int h, int i)
{
    vector<string> keys = parse_list_args(tkn_key);
    vector<string> vals = parse_list_args(tkn_val);
    vector<string> typs = parse_list_args(tkn_typ);

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
            output.emplace(keys[x], type_to_doubleptr(vals[x], typs[x], h, i));
        }
    }

    return output;
}