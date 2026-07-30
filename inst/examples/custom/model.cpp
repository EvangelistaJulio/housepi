
#include "housepi.h"
#include "mcmc.h"
#include "function.h"
#include "registry.h"

HOUSEPI_OPEN_MODULE_NAMESPACE()

/// Custom function
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

    /// Internal method called when the primitive is evaluated.
    double raw_primitive_(double x) const override final
    {
        const double I0 = arg(0);
        const double r = arg(1);
        return I0 * std::exp(r * x) / r;
    }
};
/// Register user function so that it can be used in the analysis.
HOUSEPI_REGISTER_FUNCTION(Growth_, Custom)

/// Custom model.
class Model : public ModelBase
{
public:
    /// Model specific registries of mandatory definitions
    HOUSEPI_DECLARE_REGISTER_KEYS(AugmKey, proba)
    HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(GroupRiskKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHuintKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHuintvecKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
    HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha)

    /// Compute the log-likelihood within household.
    double compute_house_transm_log_lik(Household &house) override final
    {
        double output = -ParamKey::alpha;
        for (size_t i : house.infected_members())
            output += std::log(aug(AugmKey::proba, house.member(i)));
        return output;
    }

    /// Augmented variable initialization.
    void simulate_single_house(Household &house) override final
    {
        for (size_t i : house.infected_members())
            init_missing_augm(house.member(i));
    }
};

HOUSEPI_CLOSE_MODULE_NAMESPACE()
