/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <exception>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/string.h"
#include "standalone.h"
#include "mcmc.h"
#include "mcmc/setup.h"
#include "data/wide.h"
#include "data/long.h"
#include "parameter.h"
#include "risks.h"
#include "infec_profile.h"
#include "augmented.h"

//////////////////////////////  msg  //////////////////////////////

namespace housepi
{
    namespace str
    {
        constexpr char DELIMITER = ','; /// Delimiter of CSV columns.

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
        [[noreturn]] void stop(const std::string &txt)
        {
            std::cerr << "Error: " + txt + "\nExit..." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        /// Print message `txt` to stdout, with no line ending added.
        void info(const std::string &txt)
        {
            std::cout << txt;
            std::cout.flush();
        }

        /// Display a warning message.
        void warn(const std::string &txt)
        {
            if (WarningBuffer::hold(txt))
                return;
            std::cerr << "Warning: " << txt << std::endl;
        }
    }

    //////////////////////////////  read  //////////////////////////////

    namespace read
    {
        /// Read one row of `input` into `this_row`, keyed by `headers`.
        bool read_row_by_headers_(std::ifstream &input,
                                  const std::vector<std::string> &headers,
                                  dict<std::string> &this_row)
        {
            std::string line, tkn;

            if (!std::getline(input, line))
                return false;

            std::istringstream iss(line);
            size_t n_cols = 0;

            while (std::getline(iss, tkn, str::DELIMITER))
                this_row.emplace(headers[n_cols++], str::clean(tkn));

            // last column is empty
            for (size_t i = n_cols; i < headers.size(); i++)
                this_row.emplace(headers[i], "");

            return true;
        }

        /// Read the column headers from the first line of `input`.
        std::vector<std::string> read_column_headers_(std::ifstream &input)
        {
            // cleaned headers
            std::vector<std::string> headers = {};

            // read first line
            std::string line, tkn;
            std::getline(input, line);
            std::istringstream iss(line);

            // get headers from first line
            while (std::getline(iss, tkn, str::DELIMITER))
                headers.push_back(str::clean(tkn));

            return headers;
        }

        /// Read and parse the csv table at `pathtofile`.
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
            size_t n_rows = count_rows(output);
            for (const auto &[col, vec] : output)
                if (vec.size() != n_rows)
                    throw std::runtime_error(
                        "Inconsistent number of rows found in column '" + col +
                        "' while reading file '" + pathtofile + "'.");

            return output;
        }
    }

    void Options::read_setup(MCMCBase &model)
    {
        try
        {
            model.setup.build(read::csv_table(file.setup));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading setup definitions:\n") + e.what());
        }
    }


    void Options::read_param(MCMCBase &model)
    {
        try
        {
            model.param.build(read::csv_table(file.param), model.translator);
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading parameters definitions:\n") + e.what());
        }
    }

    void Options::read_data_wide(MCMCBase &model)
    {
        try
        {
            model.data.build_wide(read::csv_table(file.data_wide));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading data-wide definitions:\n") + e.what());
        }
    }

    void Options::read_data_long(MCMCBase &model)
    {
        try
        {
            model.data.build_long(read::csv_table(file.data_long));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading data-long definitions:\n") + e.what());
        }
    }

    void Options::read_group_risk(MCMCBase &model)
    {
        try
        {
            model.risks.build_group(model.data.wide(), model.param, read::csv_table(file.group_risk));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading group risk definitions:\n") + e.what());
        }
    }

    void Options::read_contact_risk(MCMCBase &model)
    {
        try
        {
            model.risks.build_contact(model.data.wide(), model.param, read::csv_table(file.contact_risk));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading contact risk definitions:\n") + e.what());
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
            msg::stop(std::string("While reading infectivity profile definitions:\n") + e.what());
        }
    }

    void Options::read_latent(MCMCBase &model)
    {
        try
        {
            model.augm.build(model.infec_profile, model.status_process, read::csv_table(file.latent));
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string("While reading latent definitions:\n") + e.what());
        }
    }

    //////////////////////////////  Inference trace  //////////////////////////////

    /// Run inference and write its trace to `out_file`.
    void print_infer_state(MCMC &engine, const std::string &out_file)
    {
        const InferenceRun run = run_inference({&engine}, 1, [](double) {});
        if (run.failure.error)
            std::rethrow_exception(run.failure.error);

        if (na::is(out_file))
            return;

        const Setup &setup = engine.setup;
        const Posterior &post = engine.posterior;

        std::ofstream ofs;
        ofs.open(out_file);

        if (!ofs.is_open())
            throw std::runtime_error(
                "Unable to open file: " + out_file);

        // header: log-likelihood
        ofs << "iteration" << str::DELIMITER << "log_lik";

        // header: parameters
        for (const auto &name : post.param_values.names)
            ofs << str::DELIMITER << name;

        // header: latents of the traced id
        for (const auto &name : post.latent_values.names)
            ofs << str::DELIMITER << name;

        // header: chain
        ofs << str::DELIMITER << "_chain_\n";

        // values: one row per draw
        for (size_t c = 0; c < setup.n_chains; ++c)
            for (size_t k = 0; k < setup.c_sample_size; k++)
            {
                size_t i = k + c * setup.c_sample_size;

                // print iteration
                std::uint64_t sweep = setup.n_burn + (k + 1) * setup.n_thin;
                ofs << sweep;

                // joint log-likelihood (transmission + observation), summed over households from pointwise draws
                double log_lik = 0.;
                for (size_t R : post.house_log_lik.source)
                    log_lik += post.house_log_lik.cols[R][i];
                ofs << str::DELIMITER << std::setprecision(2)
                    << std::fixed << log_lik;

                // print parameter values
                for (const auto &param_col : post.param_values.cols)
                    ofs << str::DELIMITER << std::setprecision(5)
                        << std::fixed << param_col[i];

                // print latents of the traced id
                for (const auto &latent_col : post.latent_values.cols)
                    ofs << str::DELIMITER << std::setprecision(5)
                        << std::fixed << latent_col[i];

                // print the chain
                ofs << str::DELIMITER << (c + 1) << '\n';
            }

        ofs.close();
    }

    //////////////////////////////  Simulation output  //////////////////////////////

    /// Open `f` and write its header: data columns, then latent columns.
    std::ofstream open_ofstream_datawide(const std::string &f,
                                         const std::vector<std::string> &data_headers,
                                         const std::vector<std::string> &latent_headers)
    {
        std::ofstream ofs;
        ofs.open(f);

        if (!ofs.is_open())
            throw std::runtime_error("Unable to open file: '" + f + "'");

        ofs << "rep"
            << str::DELIMITER << "id"
            << str::DELIMITER << "hh_id"
            << str::DELIMITER << "status";
        for (const std::string &datastr : data_headers)
            ofs << str::DELIMITER << datastr;
        for (const std::string &latentstr : latent_headers)
            ofs << str::DELIMITER << latentstr;
        ofs << '\n';

        return ofs;
    }

    /// Write the rows of `house` for replicate `rep`.
    void update_ofstream_datawide(std::ofstream &ofs,
                                  size_t rep,
                                  const Household &house,
                                  const Data &data,
                                  const std::vector<size_t> &latent_export,
                                  const Augmented &augm)
    {
        if (!ofs.is_open())
            throw std::runtime_error("Output file closed unexpectedly.");

        const data::Wide &data_wide = data.wide();
        for (const Individual &indiv : house)
        {
            ofs << rep
                << str::DELIMITER << data_wide.id(indiv)
                << str::DELIMITER << data_wide.hh_id(house)
                << str::DELIMITER;
            // StatusLvl is a byte: stream it as a number, not as a character
            if (indiv.is_missing_status())
                ofs << "NA";
            else
                ofs << static_cast<int>(indiv.status());
            for (size_t d = 0; d < data_wide.size(); ++d)
                ofs << str::DELIMITER << data.at(d, indiv).value();
            for (size_t a : latent_export)
            {
                const Variable &xlatent = augm.at(a, indiv);
                double x = xlatent.is_empty() ? na::DOUBLE : xlatent.value();
                ofs << str::DELIMITER << x;
            }
            ofs << '\n';
        }
    }

    void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile)
    {
        data::Wide &data_wide = engine.data.wide();
        const Augmented &augm = engine.augm;

        // no output file: simulate without serializing
        const bool is_writing = !na::is(pathtofile);

        engine.augm.init_house_for_simulation(); // build household-level observations
        engine.model->init_model_for_simulation(); // model-specific initialization
        engine.augm.init_members_for_simulation(engine.status_process); // bind member definitions once per run

        // header after initialization: a produced observation registers its data column there
        // the row loop reports every column the holder carries
        std::ofstream ofs;
        if (is_writing)
            ofs = open_ofstream_datawide(pathtofile,
                                         data_wide.all(),
                                         augm.consult(augm.scheme().list_latent_export));

        // the setup step count is the number of replicates
        for (size_t r = 0; r < engine.setup.n_steps; ++r)
        {
            // clear before recovering: clear_all() walks the infected members of the last replicate
            engine.augm.clear_all();
            data_wide.recover_all_status();

            for (Household &house : data_wide)
                try
                {
                    engine.model->simulate_single_house(house);
                    if (is_writing)
                        update_ofstream_datawide(ofs, r, house,
                                                 engine.data,
                                                 augm.scheme().list_latent_export, augm);
                }
                catch (const std::exception &e)
                {
                    throw std::runtime_error(
                        "While simulating " + data_wide.str(house) + ":\n" + e.what());
                }
        }

        // close ofstream
        if (!is_writing)
            return;
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

        // check a command is supplied
        if (argc < 2)
            msg::stop("No argument or command found. Check the usage through:\n" +
                      std::string(argv[0]) + " -h");

        // parse commands and options
        while (carg_i < argc)
        {
            carg = std::string(argv[carg_i]); // update command line argument

            // help argument
            // @todo Write a manual for the standalone.
            if (toup(carg) == "-H" || tolow(carg) == "--help")
                msg::stop(
                    "No Software Manual is available, refer to the Rcpp tutorials.\nApologies for the inconvenience.");

            else if (tolow(carg) == "--log-lik")
                is_log_lik = true;

            else if (tolow(carg) == "--infer")
                is_infer = true;

            else if (tolow(carg) == "--simulate")
                is_simulate = true;

            else if (toup(carg) == "-A" || tolow(carg) == "--latent")
                file.latent = argv[++carg_i];

            else if (toup(carg) == "-D" || tolow(carg) == "--dataset")
                file.data_wide = argv[++carg_i];

            else if (toup(carg) == "-L" || tolow(carg) == "--datalong")
                file.data_long = argv[++carg_i];

            else if (toup(carg) == "-G" || tolow(carg) == "--group-risks")
                file.group_risk = argv[++carg_i];

            else if (toup(carg) == "-C" || tolow(carg) == "--contact-risks")
                file.contact_risk = argv[++carg_i];

            else if (toup(carg) == "-I" || tolow(carg) == "--infec-profile")
                file.infec_profile = argv[++carg_i];

            else if (toup(carg) == "-O" || tolow(carg) == "--output")
                file.output = argv[++carg_i];

            else if (toup(carg) == "-P" || tolow(carg) == "--parameter")
                file.param = argv[++carg_i];

            else if (toup(carg) == "-S" || tolow(carg) == "--setup")
                file.setup = argv[++carg_i];

            else
                msg::stop("Unrecognized option `" + carg + "`.");

            carg_i++; // move to the next command line argument
        }
        assert();
    }

    void Options::assert()
    {
        // check for a single option
        int aux = is_log_lik + is_infer + is_simulate;

        if (aux == 0)
            is_log_lik = true; // default
        if (aux > 1)
            msg::stop(
                "Incompatible commands: choose one between --log-lik, --infer or --simulate.");

        // simulation without an output file runs but does not serialize
        if (!na::is(file.output) || is_simulate)
            return;
        if (is_infer)
            msg::stop("Missing output file for inference.");
    }

    /// Entry point of the standalone housepi build.
    int standalone(std::string_view model_name, int argc, char *argv[])
    {
        MCMC engine(model_name);

        //////  Options  //////

        Options opts(argc, argv);
        opts.read_data_wide(engine);     // 1) Load Data Wide
        opts.read_data_long(engine);     // ... and Data Long
        opts.read_setup(engine);         // 2) Load Setup
        opts.read_param(engine);         // 3) Load Parameters
        opts.read_infec_profile(engine); // 4) Load Infectivity Profile
        opts.read_group_risk(engine);    // 5) Load Group Risks
        opts.read_contact_risk(engine);  // ... and Contact Risks
        opts.read_latent(engine);        // 6) Load Latent Variables

        //////  Execution  //////

        try
        {
            if (opts.is_log_lik)
                std::cout << engine.evaluate_log_lik("") << std::endl;
            if (opts.is_simulate)
                print_simulated_data_wide(engine, opts.file.output);
            if (opts.is_infer)
                print_infer_state(engine, opts.file.output);
        }
        catch (const std::exception &e)
        {
            msg::stop(e.what());
        }

        return EXIT_SUCCESS;
    }
}
