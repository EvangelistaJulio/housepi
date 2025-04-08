/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  * *** * DATA AUGMENTATION TO STUDY DELAYS * *** *  ///////////////////////
/////////////////////////////////  * Incubation Period example *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Global variables  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//----- Pointer to incubation period distribution -----//
std::shared_ptr<Distribution> incub_distr_ptr;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Model : init model for inference  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* [CHANGE-ME] Declare a common distribution. */
void Model::init_model_for_inference()
{
	// choose the incubation period distribution
	incub_distr_ptr = distribution_switchboard("gamma");

	// arguments are pointers to the respective parameters
	incub_distr_ptr->assign_pointer_arg("shape", param_pointer("shape"));
	incub_distr_ptr->assign_pointer_arg("scale", param_pointer("scale"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(Household &house) const
{
	double incub, output = 0.0;

	for (indiv i = 0; i < house.size(); ++i)
	{
		incub = house.augm(i, "onset") - house.augm(i, "exposure");

		if (incub_distr_ptr->valid(incub))
			output += incub_distr_ptr->logpdf(incub);
		else
			return -1e20;
	}

	return output;
}