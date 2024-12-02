/*
 * profiles.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__PROFILES__H
#define HOUSEPI__PROFILES__H

#include "function.h"

/* Base Class for Infectivity Profiles. */
class Base_Infectivity_Profile : public Base_Function
{
public:
    using Base_Function::Base_Function;

    virtual double pdf(double serial, const std::vector<double> &x = {}) const;
    virtual double cdf(double serial, const std::vector<double> &x = {}) const;
};

// Function to return a pointer to infectivity profile
std::unique_ptr<Base_Infectivity_Profile> infectivity_profile_switchboard(const std::string &name = UNDEFINED);


#endif // HOUSEPI__PROFILES__H
