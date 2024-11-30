/*
 * profiles.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__PROFILES__H
#define HOUSEPI__PROFILES__H

#include "varia.h"

class Base_Infectivity_Profile : public Base_Function
{
public:
    using Base_Function::Base_Function;

    virtual double pdf(const std::vector<double> &x) const;
    virtual double cdf(const std::vector<double> &x) const;
};

// Function to return a pointer to infectivity profile
std::unique_ptr<Base_Infectivity_Profile> infectivity_profile_ptr(const std::string &name = UNDEFINED);


#endif // HOUSEPI__PROFILES__H
