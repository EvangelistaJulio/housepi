// main.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include <string>
#include <map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "mcmc.h"

// Header to the user defined model provided by the Make command
#include "model.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The Main Function. */
int main(int argc, char *argv[])
{
    using namespace tools;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool verbose = false; // verbose option
    bool progbar = true;  // progressbar option
    bool is_silent = false;

    bool is_infer = false;    // inference option
    bool is_simulate = false; // simulation option
    bool is_evaluate = false; // evaluate option
    bool is_marginal = false; // power posterior option
    ulong num_rungs = 0;

    // paths/to/files
    std::string params_file = UNDEFINED_STRING;
    std::string data_file = UNDEFINED_STRING;
    std::string augm_file = UNDEFINED_STRING;
    std::string setup_file = UNDEFINED_STRING;
    std::string output_file = UNDEFINED_STRING;
    std::string save_augm_file = UNDEFINED_STRING;

    //=====Test if some commands were supplied=====//
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        error(oss.str());
    }

    //=====Parse commands and options=====//
    while (carg_i < argc)
    {
        carg = std::string(argv[carg_i]); // update command line argument

        // Command line argument asking for help
        if (carg == std::string("-h") or carg == std::string("--help"))
        {
            // Show the manual using man and installed manpages
            // UPDATE : produce a manual
            std::ostringstream oss;
            oss << "No Software Manual is available, refer to the Rcpp tutorials.\n"
                << "Apologies for the inconvenience.";
            error(oss.str());
        }

        else if (carg == std::string("--marginal"))
        {
            is_marginal = true;
            num_rungs = tools::stoul(argv[++carg_i], "number of rungs");
        }
        
        else if ((carg == std::string("-v")) || carg == std::string("--verbose"))
            verbose = true;

        else if (carg == std::string("--silent"))
            is_silent = true;

        else if (carg == std::string("--evaluate"))
            is_evaluate = true;

        else if (carg == std::string("--infer"))
            is_infer = true;

        else if (carg == std::string("--simulate"))
            is_simulate = true;

        else if (carg == std::string("-s") || carg == std::string("--setup"))
            setup_file = argv[++carg_i];

        else if (carg == std::string("-p") || carg == std::string("--params"))
            params_file = argv[++carg_i];

        else if (carg == std::string("-d") || carg == std::string("--database"))
            data_file = argv[++carg_i];

        else if (carg == std::string("-a") || carg == std::string("--augmented"))
            augm_file = argv[++carg_i];

        else if (carg == std::string("-o") || carg == std::string("--output"))
            output_file = argv[++carg_i];

        else if (carg == std::string("--save-augm"))
            save_augm_file = argv[++carg_i];

        else
            error("Unrecognized option `" + carg + "`.");

        carg_i++; // Move to the next command line argument
    }

    //==========================//
    //===== INITIALIZATION =====//
    //==========================//

    if (is_silent)
    {
        verbose = false;
        progbar = false;
    }

    //----- Check for command -----//
    int aux = is_evaluate + is_infer + is_simulate + is_marginal;
    if (aux > 1)
        error(
            "Incompatible commands: choose one between --evaluate, --infer, --simulate or --marginal.");

    else if (aux == 0)
        is_evaluate = true; // default

    //----- Initialize the custom Model object -----//
    std::unique_ptr<Model> user_model = std::make_unique<Model>();

    // 1) Load Setup (optional)
    user_model->read_setup(setup_file);

    // 2) Load Parameters and Check (required)
    user_model->read_params(params_file);

    // 3) Load Data and Check (required)
    user_model->read_data(data_file);

    // 4) Load Augmented Variables (optional)
    user_model->read_augm(augm_file);

    if (verbose)
        user_model->show_full_build();

    //=====================//
    //===== EXECUTION =====//
    //=====================//

    if (is_infer)
    { //===== Inference =====//
        print_decorator("Inference", UNDEFINED_STRING, C_DECOR, verbose);
        user_model->infer(progbar);

        if (verbose)
        {
            print_decorator("Inference::results", UNDEFINED_STRING);
            user_model->show_results();
        }

        if (output_file != UNDEFINED_STRING)
        {
            print_decorator("Inference::infer state",
                            "# save: " + output_file,
                            C_DECOR, verbose);
            user_model->print_infer_state(output_file);
        }

        if (save_augm_file != UNDEFINED_STRING)
        {
            print_decorator("Inference::augmented dataframe",
                            "# save: " + save_augm_file,
                            C_DECOR, verbose);
            user_model->print_augm_variables(save_augm_file);
        }
    }

    else if (is_simulate)
    { //===== Simulation Study =====//
        if (output_file == UNDEFINED_STRING)
        {
            error("Missing output file for simulation.");
        }
        else
        {
            print_decorator("Simulation",
                            "# save: " + output_file,
                            C_DECOR, verbose);
            user_model->simulate(progbar);
            user_model->print_augm_variables(output_file);

            if (verbose)
            {
                print_decorator("Simulation::summary");
                user_model->show_summary();
            }
        }
    }

    else if (is_evaluate)
    { //===== Evaluate Loglikelihood =====//
        print_decorator("Loglikelihood", UNDEFINED_STRING, C_DECOR, verbose);
        print_stdout(std::to_string(user_model->evaluate()) + "\n");

        if (save_augm_file != UNDEFINED_STRING)
        {
            print_decorator("Loglikelihood::augmented dataframe",
                            "# save: " + save_augm_file,
                            C_DECOR, verbose);
            user_model->print_augm_variables(save_augm_file);
        }
    }

    else if (is_marginal)
    { //===== Power Posteriors =====//
        print_decorator("Marginal Loglikelihood", UNDEFINED_STRING, C_DECOR, verbose);
        print_stdout(std::to_string(user_model->evaluate_marginal(num_rungs)) + "\n");
    }

    return EXIT_SUCCESS;
}
