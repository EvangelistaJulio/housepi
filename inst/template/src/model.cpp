// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * UNDEFINED MODEL * *** *  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

/* Raise error for undefined model usage. */
void Model::init_model_for_inference()
{
	basics::raise_error("Attempt to compute the likelihood with undefined model (developer only).");
}

/* Raise error for undefined model usage. */
void Model::init_model_for_simulation()
{
	basics::raise_error("Attempt to perform a simulation with undefined model (developer only).");
}

