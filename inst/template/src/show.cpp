/*
 * show.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void MCMC::show_summary() const
{
    //===== PRINT =====//

    std::ostringstream oss;

    oss << "n. of individuals: " << num_indiv
        << "\n    * n. of cases: " << report_data.at("n.cases")
        << " " << tools::percent(report_data.at("n.cases"), num_indiv);
    if (report_data.at("n.miss.stat") > 0)
    {
        oss << "\n    * n. of missing status : " << report_data.at("n.miss.stat")
            << " " << tools::percent(report_data.at("n.miss.stat"), num_indiv);
    }

    oss << "\nn. of households: " << num_house;
    if (report_data.at("no.HH.case") > 0)
    {
        oss << "\n    * with no cases: " << report_data.at("no.HH.case")
            << " " << tools::percent(report_data.at("no.HH.case"), num_house);
    }

    oss << "\n    * with at least 1 case: " << report_data.at("lst.1.case")
        << " " << tools::percent(report_data.at("lst.1.case"), num_house);
    if (report_data.at("lst.1.case") > 0)
    {
        oss << "\n    * with at least 2 cases: " << report_data.at("lst.2.case")
            << " " << tools::percent(report_data.at("lst.2.case"), num_house);
    }

    oss << "\nmax n. house cases: " << report_data.at("max.n.HH.cases")
        << " " << tools::percent(report_data.at("n.HH.with.max.cases"), num_house)
        << "\nmin n. house members: " << report_data.at("min.HH.size")
        << " " << tools::percent(report_data.at("n.HH.with.min.size"), num_house)
        << "\nmax n. house members: " << report_data.at("max.HH.size")
        << " " << tools::percent(report_data.at("n.HH.with.max.size"), num_house);

    tools::print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show group classes  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show model group classes. */
void MCMC::show_group_classes() const
{
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
                    << " " << tools::percent(x.second, num_indiv);
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
                << " " << tools::percent(x.second, _n_contacts_);
        }
        oss << "\n";
    }

    // print only if there's something to print
    if (oss.str().size() > 0)
    {
        tools::print_stdout(oss.str());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : show setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the complete choices for the setup. */
void MCMC::show_setup() const
{
    using tools::C_DELIM;

    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nnum_threads" << C_DELIM << get_num_threads()
        << "\nnum_chains" << C_DELIM << setup_num_chains
        << "\nnum_adapt" << C_DELIM << setup_num_adapt
        << "\nburn_in" << C_DELIM << setup_burn_in
        << "\nnum_iter" << C_DELIM << setup_num_iter
        << "\nthinning" << C_DELIM << setup_thinning
        << "\nrnd_init_params" << C_DELIM << tools::truth(setup_rnd_init_params);

    // if an inverse temperature has been selected print it
    if (setup_inv_T != 1.)
    {
        oss << "\ninv_temper" << C_DELIM << setup_inv_T;
    }
    oss << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::show_params() const
{
    std::ostringstream oss;
    oss << param_file_headers[0];
    for (size_t k = 1; k < param_file_headers.size(); ++k)
    {
        oss << tools::C_DELIM << param_file_headers[k];
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

    tools::print_stdout(oss.str());
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef USE_RCPP

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::String stoRString(const std::string &x)
{
    if (tools::is_na(x))
    {
        return NA_STRING;
    }
    else
    {
        return Rcpp::String(x);
    }
}

/* Export current setup configurations to R List. */
Rcpp::List MCMC::export_setup() const
{
    using namespace Rcpp;

    return List::create(
        _["seed"] = setup_seed,
        _["num_threads"] = get_num_threads(),
        _["num_chains"] = setup_num_chains,
        _["num_adapt"] = setup_num_adapt,
        _["burn_in"] = setup_burn_in,
        _["num_iter"] = setup_num_iter,
        _["thinning"] = setup_thinning,
        _["rnd_init_params"] = setup_rnd_init_params,
        _["inv_temper"] = setup_inv_T);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::export_params() const
{
    using namespace Rcpp;

    std::vector<std::string> V_names = param_to_be_inferred;
    V_names.insert(V_names.end(), param_to_be_static.begin(), param_to_be_static.end());

    NumericVector V_value = {};
    LogicalVector V_infer = {};
    CharacterVector V_prior = {}, V_args = {}, V_argsval = {};
    CharacterVector V_rule = {}, V_rule_args = {}, V_rule_argsval = {};
    NumericVector V_min = {}, V_max = {};

    for (const std::string &parID : V_names)
    {
        V_value.push_back(get_param(parID));
        V_infer.push_back(get_infer(parID));

        V_prior.push_back(stoRString(show_param_distr_name(parID)));
        V_args.push_back(stoRString(show_param_distr_arg_names(parID)));
        V_argsval.push_back(stoRString(show_param_distr_arg_values(parID)));

        V_rule.push_back(stoRString(show_param_prop_name(parID)));
        V_rule_args.push_back(stoRString(show_param_prop_arg_names(parID)));
        V_rule_argsval.push_back(stoRString(show_param_prop_arg_values(parID)));

        V_min.push_back(std::isnan(get_param_range_min(parID)) ? NA_REAL : get_param_range_min(parID));
        V_max.push_back(std::isnan(get_param_range_max(parID)) ? NA_REAL : get_param_range_max(parID));
    }

    return DataFrame::create(_["name"] = V_names, _["value"] = V_value, _["infer"] = V_infer,
                             _["prior"] = V_prior, _["args"] = V_args, _["argsval"] = V_argsval,
                             _["rule"] = V_rule, _["rule_args"] = V_rule_args, _["rule_argsval"] = V_rule_argsval,
                             _["min"] = V_min, _["max"] = V_max);
}

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::export_data() const
{
    using namespace Rcpp;

    // pre-allocate memories
    // initialize the dataframe as a list of named std::vectors
    CharacterVector df_index(num_indiv);
    CharacterVector df_house_index(num_indiv);
    IntegerVector df_status(num_indiv);
    std::map<std::string, NumericVector> df_dbl;

    for (const std::string &dataID : _data_headers_)
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

            for (const std::string &dataID : _data_headers_)
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
    for (const std::string &dataID : _data_headers_)
    {
        output.push_back(df_dbl.at(dataID), dataID);
    }
    for (const std::string &augmID : augm_full_set)
    {
        output.push_back(df_dbl.at(augmID), augmID);
    }

    return output;
}

/* Export MCMC traces with posterior samples. */
Rcpp::DataFrame MCMC::export_traces() const
{
    using namespace Rcpp;

    DataFrame output = DataFrame::create(
        _["loglikelihood"] = wrap(dgn_glob_loglike));

    for (const std::string &parID : param_to_be_inferred)
    {
        output.push_back(wrap(get_dgn_value(parID)), parID);
    }

    // retrieve chain identities
    if (setup_num_chains > 1)
    {
        ulong size = (setup_num_iter - setup_burn_in) / setup_thinning;
        IntegerVector df_chain(size * setup_num_chains);
        ulong k = 0;
        for (ulong c = 1; c <= setup_num_chains; ++c)
        {
            std::fill(df_chain.begin() + k, df_chain.begin() + k + size, c);
            k += size;
        }
        output.push_back(df_chain, "CHAIN__");
    }

    return output;
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix MCMC::export_loglike_matrix() const
{
    using namespace Rcpp;

    ulong sample_size = (setup_num_iter - setup_burn_in) / setup_thinning;
    sample_size *= setup_num_chains;

    // Create a NumericMatrix with N rows and K columns
    NumericMatrix output(sample_size, num_house);

    // Assign each NumericVector to a column of the matrix
    for (int h = 0; h < num_house; ++h)
    {
        NumericVector column = wrap(house[h].get_dgn_loglike());
        output(_, h) = column;
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
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void MCMC::show_results() const
{
    using tools::C_DELIM;

    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";
    for (size_t k = 2; k < param_file_headers.size(); ++k)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    oss << "\n";
    for (const std::string &parID : param_to_be_inferred)
    {
        oss << show_param_diagnostics(parID) << "\n";
    }

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::print_infer_state(const std::string &output_filename) const
{
    using tools::C_DELIM;

    if (!tools::is_na(output_filename))
    {
        std::ofstream inferstate_ofstream;
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
            ulong numsweeps = (setup_num_iter - setup_burn_in) / setup_thinning;

            //----- HEADERS -----//
            // header : loglikelihood
            inferstate_ofstream << "loglikelihood";
            // header : parameters
            for (const std::string &parID : param_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << parID;
            }
            // header : chain
            if (setup_num_chains > 1)
            {
                inferstate_ofstream << C_DELIM << "CHAIN__";
            }
            inferstate_ofstream << std::endl;

            //----- VALUES -----//

            // FIXME : // FLAG : // HERE :
            for (size_t c = 0; c < setup_num_chains; ++c)
            {
                for (size_t k = 0; k < numsweeps; k++)
                {
                    size_t i = k + c * numsweeps;
                    // print : loglikelihood
                    inferstate_ofstream << std::setprecision(2) << std::fixed << dgn_glob_loglike[i];
                    // print : parameter values
                    for (const std::string &parID : param_to_be_inferred)
                    {
                        inferstate_ofstream << C_DELIM << std::setprecision(5) << std::fixed << get_dgn_value(parID)[i];
                    }
                    // print : which chain
                    if (setup_num_chains > 1)
                    {
                        inferstate_ofstream << C_DELIM << (c + 1);
                    }
                    inferstate_ofstream << std::endl;
                }
            }

            inferstate_ofstream.close();
        }
        else
        {
            tools::error("Unable to open file: " + output_filename);
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
/////////////////////////////////  MCMC : print augm variables  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void MCMC::print_augm_variables(const std::string &save_augm_pathtofile) const
{
    using tools::C_DELIM;

    if (!tools::is_na(save_augm_pathtofile))
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(save_augm_pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            //===== Print column headers =====//
            // default headers
            augm_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

            // data headers
            for (const std::string &dataID : _data_headers_)
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
                    for (const std::string &dataID : _data_headers_)
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
            tools::error("Unable to open file: " + save_augm_pathtofile);
        }
    }
}

#endif