
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
    HOUSEPI_DECLARE_REGISTER_KEYS(Augm, proba)
    HOUSEPI_DECLARE_EMPTY_KEYS(Data)
    HOUSEPI_DECLARE_EMPTY_KEYS(Long)
    HOUSEPI_DECLARE_EMPTY_KEYS(GroupRisk)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHuint)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHuintvec)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevec)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtx)
    HOUSEPI_DECLARE_REGISTER_KEYS(Param, alpha)

    /// Compute the log-likelihood within household.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        double output = -Param::alpha;
        for (const Individual &indiv : house.infected_indivs())
            output += std::log(aug(Augm::proba, indiv));
        return output;
    }

    /// Augmented variable initialization.
    void simulate_single_house(Household &house) override final
    {
		for (size_t i = 0; i < house.size(); ++i)
		{
			Individual &indiv = house.at(i);
			if (indiv.is_infected_status())
				init_missing_augm(indiv);
		}
    }
};

HOUSEPI_CLOSE_MODULE_NAMESPACE()
