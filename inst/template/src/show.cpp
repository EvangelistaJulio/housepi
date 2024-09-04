/*
 * show.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Inconmplete_Infer_Diagnostics_Warning = false; // wrong call to show diagnostics

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void MCMC::show_summary()
{
    using namespace screen;

    //===== PRINT =====//

    std::ostringstream oss;

    oss << "n. of individuals: " << num_indiv
        << "\n    * n. of cases: " << report_data["n.cases"]
        << " " << percent(report_data["n.cases"], num_indiv);
    if (report_data["n.miss.stat"] > 0)
    {
        oss << "\n    * n. of missing status : " << report_data["n.miss.stat"]
            << " " << percent(report_data["n.miss.stat"], num_indiv);
    }

    oss << "\nn. of households: " << num_househ;
    if (report_data["no.HH.case"] > 0)
    {
        oss << "\n    * with no cases: " << report_data["no.HH.case"]
            << " " << percent(report_data["no.HH.case"], num_househ);
    }

    oss << "\n    * with at least 1 case: " << report_data["lst.1.case"]
        << " " << percent(report_data["lst.1.case"], num_househ);
    if (report_data["lst.1.case"] > 0)
    {
        oss << "\n    * with at least 2 cases: " << report_data["lst.2.case"]
            << " " << percent(report_data["lst.2.case"], num_househ);
    }

    oss << "\nmax n. house cases: " << report_data["max.n.HH.cases"]
        << " " << percent(report_data["n.HH.with.max.cases"], num_househ)
        << "\nmin n. house members: " << report_data["min.HH.size"]
        << " " << percent(report_data["n.HH.with.min.size"], num_househ)
        << "\nmax n. house members: " << report_data["max.HH.size"]
        << " " << percent(report_data["n.HH.with.max.size"], num_househ);

    print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show group classes  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show model group classes. */
void MCMC::show_group_classes()
{
    using namespace screen;

    std::ostringstream oss;

    // format group classes
    for (const auto &entry : report_group_classes)
    {
        if ((entry.first != CONTACT_RATE) && (entry.first != EMPTY_RATE))
        {
            oss << entry.first << ":";
            for (const auto &x : entry.second)
            {
                oss << "\n    * " << x.first << ": " << x.second
                    << " " << percent(x.second, num_indiv);
            }
            oss << "\n";
        }
    }

    // format contact rates
    if (report_group_classes.at(CONTACT_RATE).find(EMPTY_RATE) != report_group_classes.at(CONTACT_RATE).end())
    {
        oss << CONTACT_RATE << ":";
        for (const auto &x : report_group_classes.at(CONTACT_RATE))
        {
            oss << "\n    * " << x.first << ": " << x.second
                << " " << percent(x.second, __n_contacts);
        }
        oss << "\n";
    }

    // print only if there's something to print
    if (oss.str().size() > 0)
    {
        print_stdout(oss.str());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : show setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the complete choices for the setup. */
void MCMC::show_setup()
{
    using namespace screen;

    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nnum_threads" << C_DELIM << get_num_threads()
        << "\nnum_chains" << C_DELIM << setup_num_chains
        << "\nnum_adapt" << C_DELIM << setup_num_adapt
        << "\nburn_in" << C_DELIM << setup_burn_in
        << "\nnum_iter" << C_DELIM << setup_num_iter
        << "\nthinning" << C_DELIM << setup_thinning
        << "\nrnd_init_params" << C_DELIM << truth(setup_rnd_init_params);

    // if an inverse temperature has been selected print it
    if (setup_inv_T != 1.)
    {
        oss << "\ninv_temper" << C_DELIM << setup_inv_T;
    }
    oss << "\n";

    print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::show_parameters()
{
    using namespace screen;

    std::ostringstream oss;
    oss << param_file_headers[0];
    for (size_t k = 1; k < param_file_headers.size(); k++)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    oss << "\n";

    for (const std::string &parID : param_to_be_inferred)
    {
        oss << show_param(parID) << "\n";
    }
    for (const std::string &parID : param_to_be_static)
    {
        oss << show_param(parID) << "\n";
    }

    print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void MCMC::show_results()
{
    using namespace screen;

    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";
    for (size_t k = 2; k < param_file_headers.size(); k++)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    oss << "\n";
    for (const std::string &parID : param_to_be_inferred)
    {
        oss << show_param_diagnostics(parID) << "\n";
    }

    print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show diagnostics  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show diagnostics as sampled during iterations. */
void MCMC::show_diagnostics()
{
    using namespace screen;

    if (!__is_infer_complete)
    {
        raise_warning("Ignored attempt to show diagnostics before inference.",
                      Inconmplete_Infer_Diagnostics_Warning);
    }
    else
    {
        std::ostringstream oss;
        oss << std::setprecision(2) << std::fixed
            << "loglikelihood_mean" << C_DELIM << dgn_avg << "\n"
            << "loglikelihood_sd" << C_DELIM << dgn_var() << "\n"
            << "log_pointwise_predictive_density" << C_DELIM << dgn_lppd << "\n"
            << "effective_parameters" << C_DELIM << dgn_eff_pars << "\n"
            << "WAIC" << C_DELIM << 2. * (dgn_eff_pars - dgn_lppd) << "\n";
        print_stdout(oss.str());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : init param file  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::inferstate_init(std::string output_filename)
{
    using namespace screen;

    if (output_filename == UNDEFINED)
    {
        raise_error("Undefined output file.");
    }
    else
    {
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
            __is_print_inferstate = true;
            __is_infer_complete = false;

            // print loglikelihood
            inferstate_ofstream << "loglikelihood";

            for (const std::string & parID : param_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << parID;
            }

            // which chain
            if (setup_num_chains > 1)
            {
                inferstate_ofstream << C_DELIM << "CHAIN__";
            }

            inferstate_ofstream << std::endl;
        }
        else
        {
            raise_error("Unable to open file: " + output_filename);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the model parameters at the given iteration time */
void MCMC::inferstate_print()
{
    using namespace screen;

    if (__is_print_inferstate)
    {
        if (inferstate_ofstream.is_open())
        {
            // print loglikelihood
            inferstate_ofstream << std::setprecision(2) << std::fixed << get_global_loglike();

            // print parameter values
            for (const std::string &parID : param_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << std::setprecision(5) << std::fixed << get_value(parID);
            }

            // print which chain
            if (setup_num_chains > 1)
            {
                inferstate_ofstream << C_DELIM << __chain;
            }

            inferstate_ofstream << std::endl;
        }
        else
        {
            raise_error("Output file closed unexpectedly.");
        }
    }
}

/* Close file with inferred parameters. */
void MCMC::inferstate_close()
{
    __is_print_inferstate = false;
    __is_infer_complete = true;
    inferstate_ofstream.close();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : print augm variables  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void MCMC::print_augm_variables(std::string save_augm_pathtofile)
{
    using namespace screen;

    if (save_augm_pathtofile != UNDEFINED)
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(save_augm_pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            //===== Print column headers =====//
            // default headers
            augm_ofstream << "index" << C_DELIM << "house_index" << C_DELIM << "status";

            // data headers
            for (const std::string &dataID : __data_headers)
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
            for (int h = 0; h < num_househ; h++)
            {
                for (int i = 0; i < househ[h].house_size; i++)
                {
                    augm_ofstream << househ[h].indiv[i].get_index()
                                  << C_DELIM << househ[h].get_index()
                                  << C_DELIM << househ[h].indiv[i].get_status();
                    for (const std::string &dataID : __data_headers)
                    {
                        augm_ofstream << C_DELIM << househ[h].indiv[i].get_data(dataID);
                    }
                    for (const std::string &augmID : augm_full_set)
                    {
                        augm_ofstream << C_DELIM << househ[h].indiv[i].get_augm(augmID);
                    }
                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            raise_error("Unable to open file: " + save_augm_pathtofile);
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
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_RCPP

/* Export current setup configurations to R List. */
Rcpp::List MCMC::export_setup()
{
    using namespace Rcpp;

    List output;
    output = List::create(
        _["seed"] = setup_seed,
        _["num_threads"] = get_num_threads(),
        _["num_chains"] = setup_num_chains,
        _["num_adapt"] = setup_num_adapt,
        _["burn_in"] = setup_burn_in,
        _["num_iter"] = setup_num_iter,
        _["thinning"] = setup_thinning,
        _["rnd_init_params"] = screen::truth(setup_rnd_init_params));

    if (setup_inv_T != 1.)
    {
        output["inv_temper"] = setup_inv_T;
    }

    return output;
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::export_params()
{
    using namespace Rcpp;
    DataFrame output;

    std::vector<std::string> vec_names = param_to_be_inferred;
    vec_names.insert(vec_names.end(), param_to_be_static.begin(), param_to_be_static.end());

    std::vector<double> vec_value = {};
    std::vector<bool> vec_infer = {};
    std::vector<std::string> vec_prior = {}, vec_args = {}, vec_argsval = {};
    std::vector<std::string> vec_rule = {}, vec_rargs = {}, vec_rargsval = {};
    std::vector<std::string> vec_min = {}, vec_max = {};

    for (const std::string &parID : vec_names)
    {
        vec_value.push_back(get_value(parID));
        vec_infer.push_back(get_infer(parID));
        vec_prior.push_back(show_param_distr_name(parID));
        vec_args.push_back(parse_list_args(show_param_distr_args(parID), screen::C_DELIM)[0]);
        vec_argsval.push_back(parse_list_args(show_param_distr_args(parID), screen::C_DELIM)[1]);
        vec_rule.push_back(show_param_prop_name(parID));
        vec_rargs.push_back(parse_list_args(show_param_prop_args(parID), screen::C_DELIM)[0]);
        vec_rargsval.push_back(parse_list_args(show_param_prop_args(parID), screen::C_DELIM)[1]);
        vec_min.push_back(parse_list_args(show_param_prop_range(parID), screen::C_DELIM)[0]);
        vec_max.push_back(parse_list_args(show_param_prop_range(parID), screen::C_DELIM)[1]);
    }

    output = DataFrame::create(_["name"] = vec_names,
                               _["value"] = vec_value,
                               _["infer"] = vec_infer,
                               _["prior"] = vec_prior,
                               _["args"] = vec_args,
                               _["argsval"] = vec_argsval,
                               _["rule"] = vec_rule,
                               _["rargs"] = vec_rargs,
                               _["rargsval"] = vec_rargsval,
                               _["min"] = vec_min,
                               _["max"] = vec_max);

    return output;
}

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::export_data()
{
    using namespace Rcpp;

    // initialize the dataframe as a list of named vectors
    std::map<std::string, std::vector<int>> df_int;
    std::map<std::string, std::vector<double>> df_dbl;

    DataFrame output;

    df_int.insert({"index", std::vector<int>(num_indiv, 0)});
    df_int.insert({"house_index", std::vector<int>(num_indiv, 0)});
    df_int.insert({"status", std::vector<int>(num_indiv, 0)});

    for (const std::string &dataID : __data_headers)
    {
        df_dbl.insert({dataID, std::vector<double>(num_indiv, 0.)});
    }
    for (const std::string &augmID : augm_full_set)
    {
        df_dbl.insert({augmID, std::vector<double>(num_indiv, 0.)});
    }

    size_t c = 0;
    for (int h = 0; h < num_househ; h++)
    {
        for (int i = 0; i < househ[h].house_size; i++)
        {
            df_int.at("index")[c] = househ[h].indiv[i].get_index();
            df_int.at("house_index")[c] = househ[h].get_index();
            df_int.at("status")[c] = househ[h].indiv[i].get_status();

            for (const std::string &dataID : __data_headers)
            {
                df_dbl.at(dataID)[c] = househ[h].indiv[i].get_data(dataID);
            }
            for (const std::string &augmID : augm_full_set)
            {
                df_dbl.at(augmID)[c] = househ[h].indiv[i].get_augm(augmID);
            }

            c++;
        }
    }

    output = DataFrame::create(_["index"] = df_int.at("index"),
                               _["house_index"] = df_int.at("house_index"),
                               _["status"] = df_int.at("status"));
    for (const std::string &dataID : __data_headers)
    {
        output.push_back(df_dbl.at(dataID), dataID);
    }
    for (const std::string &augmID : augm_full_set)
    {
        output.push_back(df_dbl.at(augmID), augmID);
    }

    return output;
}

#endif