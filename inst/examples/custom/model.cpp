
#include "housepi.h"
#include "mcmc.h"
#include "function.h"
#include "registry.h"

class Growth_ : public Function
{
public:
    Growth_() : Function({{"initial", 1.}, {"rate", 1.}}, 0.) {}
    HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Growth_)
    HOUSEPI_OVERRIDE_NAME_METHOD(growth)

    void init_discretization() override final { helper_init_range_min_(0.); }

protected:
    /// Internal method called when the function is evaluated.
    /// Arguments order is by definition.
    double raw_function_(double x) const override final
    {
        const double I0 = arg(0);
        const double r = arg(1);
        return I0 * std::exp(r * x);
    }
};

/// Register user function so that it can be used in the analysis.
HOUSEPI_REGISTER_FUNCTION(Growth_, Custom)

HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE()
HOUSEPI_DECLARE_REGISTER_KEYS(AugmKey, proba)
HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
HOUSEPI_DECLARE_EMPTY_KEYS(ParamKey)

class Model : public MCMC
{
public:
    HOUSEPI_MODEL_CONSTRUCTOR()

    /// Compute the log-likelihood within household.
    double compute_house_obs_log_lik(Household &house)
    {
        double output = 0.;
        for (size_t i : house.infected_hhindexes())
            output += std::log(house.member(i).aug(AugmKey::proba));
        return output;
    }

    /// Augmented variable initialization.
    void simulate_single_house(Household &house)
    {
        for (size_t i : house.infected_hhindexes())
            init_missing_augm(house.member(i));
    }
};

HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE()
HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE()
