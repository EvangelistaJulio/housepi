
#include "function.h"
#include "function_registry.h"
#include "individual.h"
#include "household.h"
#include "model.h"

class MyMonomialFunction : public Function
{
public:
    MyMonomialFunction()                     // constructor
        : Function("my_monomial",            // function name used in registers
                   {{"a", 1.0}, {"b", 2.0}}, // names of arguments
                   0.0)
    {
    }

    void init_discretization() { init_discr_helper_(hsp::discr::LEFT); }

protected:
    // internal method called when the function is evaluated
    double raw_function_(double x) const override final
    {
        const double a = arg(0); // arguments order is by definition
        const double b = arg(1);
        return a * std::pow(x, b);
    }

    // internal method called when the function is evaluated
    double raw_primitive_(double x) const override final
    {
        const double a = arg(0), b_plus_1 = arg(1) + 1;
        return a * std::pow(x, b_plus_1) / (b_plus_1);
    }
};

// register user function so that it can be used in the analysis
REGISTER_USER_FUNCTION(MyMonomialFunction)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood within household.
double Model::compute_house_obs_log_lik(Household &house)
{
    double output = 0.0;
    for (const Individual &member : house.infected())
        output += std::log(member.aug("proba"));
    return output;
}

/// Augmented variable initialization.
void Model::simulate_single_house(Household &house)
{
    for (Individual &member : house.infected())
        init_missing_augm(member);
}
