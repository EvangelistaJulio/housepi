/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/read.h"
#include "tools/string.h"

#include "mcmc.h"
#include "standalone.h"
#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Stdout Message Function  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    namespace str
    {
        std::string tolow(std::string s)
        {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c)
                           { return std::tolower(c); });
            return s;
        }
    }

    namespace msg
    {
        /// Throw an error message.
        [[noreturn]] void error(const std::string &txt)
        {
            std::cerr << "Error: " + txt + "\nExit..." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        /// Print an information message to stdout.
        void info(const std::string &txt) noexcept
        {
            std::cout << txt;
            std::cout.flush();
        }

        /// Display a warning message if not already displayed before.
        void warning(const std::string &txt)
        {
            if (is_new_warning(txt))
            {
                std::cerr << "Warning: " << txt << std::endl;
                add_new_warning(txt);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace read
    {
        /// Load files from the row with headers.
        bool read_row_by_headers_(std::ifstream &input,
                                  const std::vector<std::string> &headers,
                                  dict<std::string> &this_row)
        {
            std::string line, tkn;

            if (!std::getline(input, line))
                return false;

            std::istringstream iss(line);
            size_t n_cols = 0;

            while (std::getline(iss, tkn, DELIMITER))
                this_row.emplace(headers[n_cols++], strcleaner(tkn));

            // last column is empty
            for (size_t i = n_cols; i < headers.size(); i++)
                this_row.emplace(headers[i], "");

            return true;
        }

        /// Read column headers and compare to required headers.
        std::vector<std::string> read_column_headers_(std::ifstream &input)
        {
            // clean headers
            std::vector<std::string> headers = {};

            // read first line
            std::string line, tkn;
            std::getline(input, line);
            std::istringstream iss(line);

            // get headers from first line
            while (std::getline(iss, tkn, DELIMITER))
                headers.push_back(strcleaner(tkn));

            return headers;
        }

        /// Read csv table in pathtofile and store it parsed in csv_parsed_table.
        str_table csv_table(const std::string &pathtofile)
        {
            if (na::is(pathtofile))
                return {};

            std::ifstream input(pathtofile.c_str());
            if (!input.is_open())
                throw std::runtime_error(
                    "Unable to open file '" + pathtofile + "'.");

            // read and check for valid headers
            const std::vector<std::string> headers = read_column_headers_(input);

            str_table output = {};
            dict<std::string> this_row = {};
            while (read_row_by_headers_(input, headers, this_row))
            {
                for (const auto &[key, val] : this_row)
                    output[key].push_back(val);
                this_row.clear();
            }

            // ensure all columns have the same number of rows
            size_t n_rows = csv_n_rows(output);
            for (const auto &[col, vec] : output)
                if (vec.size() != n_rows)
                    throw std::runtime_error(
                        "Inconsistent number of rows found in column '" + col +
                        "' while reading file '" + pathtofile + "'.");

            return output;
        }
    }

#define STANDALONE_MAKE_READ_METHOD_(module, section_name)           \
    void Options::read_##module(MCMCBase &model)                     \
    {                                                                \
        try                                                          \
        {                                                            \
            model.module.build(read::csv_table(file.module));        \
        }                                                            \
        catch (const std::exception &e)                              \
        {                                                            \
            msg::error(std::string("While reading ") +               \
                       section_name + " definitions:\n" + e.what()); \
        }                                                            \
    }

    STANDALONE_MAKE_READ_METHOD_(setup, "setup")

    void Options::read_param(MCMCBase &model)
    {
        try
        {
            model.param.build(read::csv_table(file.param), model.translator);
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading parameters definitions:\n") + e.what());
        }
    }

    void Options::read_data_wide(MCMCBase &model)
    {
        try
        {
            model.data_wide().build(read::csv_table(file.data_wide));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading data-wide definitions:\n") + e.what());
        }
    }

    void Options::read_data_long(MCMCBase &model)
    {
        try
        {
            model.data_long().build(model.data_wide(), read::csv_table(file.data_long));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading data-long definitions:\n") + e.what());
        }
    }

    void Options::read_group_rate(MCMCBase &model)
    {
        try
        {
            model.group_rate.build(model.data_wide(), model.param, read::csv_table(file.group_rate));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading group rate definitions:\n") + e.what());
        }
    }

    void Options::read_contact_rate(MCMCBase &model)
    {
        try
        {
            model.contact_rate.build(model.data_wide(), model.param, read::csv_table(file.contact_rate));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading contact rate definitions:\n") + e.what());
        }
    }

    void Options::read_infec_profile(MCMCBase &model)
    {
        try
        {
            model.infec_profile.build(model.translator, read::csv_table(file.infec_profile));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading infectivity profile definitions:\n") + e.what());
        }
    }

    void Options::read_augm(MCMCBase &model)
    {
        try
        {
            model.augm.build(model.infec_profile, read::csv_table(file.augm));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While reading augmented definitions:\n") + e.what());
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the output file with inferred parameters with column names and initial values.
    void print_infer_state(MCMC &engine, const std::string &out_file)
    {
        engine.run_inference();

        if (na::is(out_file))
            return;

        const Setup &setup = engine.setup;
        const Posterior &post = engine.posterior;

        std::ofstream ofs;
        ofs.open(out_file);

        if (!ofs.is_open())
            throw std::runtime_error(
                "Unable to open file: " + out_file);

        // Headers
        // header log-likelihood
        ofs << "iteration" << read::DELIMITER << "log_lik";

        // header parameters
        for (const auto &name : post.param_names)
            ofs << read::DELIMITER << name;

        // header augmented variable of id trace
        for (const auto &name : post.augm_col_names)
            ofs << read::DELIMITER << name;

        // header chain
        ofs << read::DELIMITER << "_chain_" << std::endl;

        // Values
        for (size_t c = 0; c < setup.n_chains; ++c)
            for (size_t k = 0; k < setup.c_sample_size; k++)
            {
                size_t i = k + c * setup.c_sample_size;

                // print iteration
                unsigned long sweep = setup.n_burn + (k + 1) * setup.n_thin;
                ofs << sweep;

                // print transmission log-likelihood
                ofs << read::DELIMITER << std::setprecision(2)
                    << std::fixed << post.tot_transm_log_lik[0][i];

                // print parameter values
                for (size_t u = 0; u < post.param_values.size(); ++u)
                    ofs << read::DELIMITER << std::setprecision(5)
                        << std::fixed << post.param_values[u][i];

                // print augmented variable of id trace
                {
                    for (auto &member_matrix : post.augm_values)
                        for (auto &var_vec : member_matrix)
                            ofs << read::DELIMITER << std::setprecision(5)
                                << std::fixed << var_vec[i];
                }

                // print which chain
                ofs << read::DELIMITER << (c + 1) << std::endl;
            }

        ofs.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  MCMC::Impl : print simulation  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Print all augmented variables of each individual.
    std::ofstream open_ofstream_datawide(const std::string &f,
                                         const std::vector<std::string> &data_headers,
                                         const std::vector<std::string> &augm_headers)
    {
        std::ofstream ofs;
        ofs.open(f);

        if (!ofs.is_open())
            throw std::runtime_error("Unable to open file: '" + f + "'");

        ofs << "id"
            << read::DELIMITER << "house_id"
            << read::DELIMITER << "status";
        for (const std::string &datastr : data_headers)
            ofs << read::DELIMITER << datastr;
        for (const std::string &augmstr : augm_headers)
            ofs << read::DELIMITER << augmstr;
        ofs << std::endl;

        return ofs;
    }

    /// Print all augmented variables of each individual.
    void update_ofstream_datawide(std::ofstream &ofs,
                                  const Household &house,
                                  const DataWide &data_wide,
                                  const std::vector<size_t> &augm_export,
                                  const Augmented &augm)
    {
        if (!ofs.is_open())
            throw std::runtime_error("Output file closed unexpectedly.");

        for (const Individual &member : house)
        {
            ofs << member.id()
                << read::DELIMITER << house.house_id()
                << read::DELIMITER << member.status();
            for (size_t d = 0; d < data_wide.khold().size(); ++d)
                ofs << read::DELIMITER << data_wide.data(d, member);
            for (size_t a : augm_export)
            {
                const Variable &xaugm = augm.at(a, member);
                double x = xaugm.is_empty() ? na::DOUBLE : xaugm.value();
                ofs << read::DELIMITER << x;
            }
            ofs << std::endl;
        }
    }

    void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile)
    {
        DataWide &data_wide = engine.data_wide();
        const Augmented &augm = engine.augm;

        // initialization
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   data_wide.khold().all(),
                                                   augm.khold().as_string(augm.list_export));
        // model specific initialization
        engine.model->init_model_for_simulation();

        for (Household &house : data_wide)
            try
            {
                engine.model->simulate_single_house(house);
                update_ofstream_datawide(ofs, house,
                                         data_wide,
                                         augm.list_export, augm);
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" + e.what());
            }

        // close ofstream
        if (!ofs.is_open())
            throw std::runtime_error(
                "Output file '" + pathtofile + "' closed unexpectedly.");
        ofs.close();
    }

    Options::Options(int argc, char *argv[])
    {
        using str::tolow;
        using str::toup;

        int carg_i = 1;   // command line iterator
        std::string carg; // command line argument

        // check any command was supplied
        if (argc < 2)
            msg::error("No argument or command found. Check the usage through:\n" +
                       std::string(argv[0]) + " -h");

        // parse commands and options
        while (carg_i < argc)
        {
            carg = std::string(argv[carg_i]); // update command line argument

            // Command line argument asking for help
            // UPDATE : should I produce a manual for the standalone?
            if (toup(carg) == "-H" || tolow(carg) == "--help")
                msg::error(
                    "No Software Manual is available, refer to the Rcpp tutorials.\nApologies for the inconvenience.");

            else if (tolow(carg) == "--log-lik")
                is_log_lik = true;

            else if (tolow(carg) == "--infer")
                is_infer = true;

            else if (tolow(carg) == "--simulate")
                is_simulate = true;

            /*  --FUTURE-LOG-EVIDENCE--
            else if (tolow(carg) == "--log-evidence")
                num_rungs = str::toul(argv[++carg_i], "number of rungs");
            */

            else if (toup(carg) == "-A" || tolow(carg) == "--augmented")
                file.augm = argv[++carg_i];

            else if (toup(carg) == "-D" || tolow(carg) == "--dataset")
                file.data_wide = argv[++carg_i];

            else if (toup(carg) == "-L" || tolow(carg) == "--datalong")
                file.data_long = argv[++carg_i];

            else if (toup(carg) == "-G" || tolow(carg) == "--group-rates")
                file.group_rate = argv[++carg_i];

            else if (toup(carg) == "-C" || tolow(carg) == "--contact-rates")
                file.contact_rate = argv[++carg_i];

            else if (toup(carg) == "-I" || tolow(carg) == "--infec-profile")
                file.infec_profile = argv[++carg_i];

            else if (toup(carg) == "-O" || tolow(carg) == "--output")
                file.output = argv[++carg_i];

            else if (toup(carg) == "-P" || tolow(carg) == "--parameter")
                file.param = argv[++carg_i];

            else if (toup(carg) == "-S" || tolow(carg) == "--setup")
                file.setup = argv[++carg_i];

            else
                msg::error("Unrecognized option `" + carg + "`.");

            carg_i++; // move to the next command line argument
        }
        assert();
    }

    void Options::assert()
    {
        // check for single option
        int aux = is_log_lik + is_infer + is_simulate + is_log_evidence();

        if (aux == 0)
            is_log_lik = true; // default
        if (aux > 1)
            msg::error(
                "Incompatible commands: choose one between --log-lik, --infer, --simulate or --log-evidence.");

        if (!na::is(file.output))
            return;
        if (is_infer)
            msg::error("Missing output file for inference.");
        if (is_simulate)
            msg::error("Missing output file for simulation.");
    }

    /* --FUTURE-LOG-EVIDENCE--
    double log_evidence(MCMC &engine, std::string_view scheme, unsigned long n_rungs)
    {
    double output = 0;
    const auto pp = engine.power_posteriors(scheme, n_rungs);
    const std::vector<double> &inv_tem = pp[0];
    const std::vector<double> &avg_vec = pp[1];
    const std::vector<double> &var_vec = pp[2];

    for (unsigned long k = 0; k < engine.setup.n_rungs; ++k)
    {
        double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
        output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
        output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
    }

    return output;
    }
    */

    /// Main function for the standalone version of 'housepi'.
    int standalone(std::string_view model_name, int argc, char *argv[])
    {
        MCMC engine(model_name);

        ///////////////////////
        //////  OPTIONS  //////
        ///////////////////////

        Options opts(argc, argv);
        opts.read_data_wide(engine);     // 1) Load Data Wide
        opts.read_data_long(engine);     // ... and Data Long
        opts.read_setup(engine);         // 2) Load Setup
        opts.read_param(engine);         // 3) Load Parameters
        opts.read_infec_profile(engine); // 4) Load Infectivity Profile
        opts.read_group_rate(engine);    // 5) Load Group Rates
        opts.read_contact_rate(engine);  // ... and Contact Rates
        opts.read_augm(engine);          // 6) Load Augmented Variables

        //////////////////////////////
        ////////  EXECUTION  /////////
        //////////////////////////////

        try
        {
            if (opts.is_log_lik)
                std::cout << engine.evaluate_log_lik("") << std::endl;
            if (opts.is_simulate)
                print_simulated_data_wide(engine, opts.file.output);
            if (opts.is_infer)
                print_infer_state(engine, opts.file.output);

            /* --FUTURE-LOG-EVIDENCE--
            if (opts.is_log_evidence())
                std::cout << log_evidence(engine, "adaptive", num_rungs) << std::endl;
            */
        }
        catch (const std::exception &e)
        {
            engine.setup.crash(e);
        }

        return EXIT_SUCCESS;
    }
}
