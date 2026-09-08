
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

    void truncate() override final { helper_init_domain_min_(0.); }

    /// Values returned, as opposed to the domain set by truncate() above.
    /// Declaring it lets the engine prove a model requirement this function
    /// satisfies; leaving it out only means the engine assumes nothing.
    Codomain codomain() const noexcept override final
    {
        Codomain c;
        c.inf = 0.;
        c.below = Interval::Open;
        return c;
    }

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
    HOUSEPI_DECLARE_REGISTER_KEYS(Latent, proba)
    HOUSEPI_DECLARE_EMPTY_KEYS(Data)
    HOUSEPI_DECLARE_EMPTY_KEYS(Long)
    HOUSEPI_DECLARE_EMPTY_KEYS(GroupRisk)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHmember)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHmembervec)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevec)
    HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtx)
    HOUSEPI_DECLARE_REGISTER_KEYS(Param, alpha)

    /// Compute the log-likelihood within household.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        double output = -Param::alpha;
        for (const Individual &indiv : house.infected_indivs())
            output += std::log(aug(Latent::proba, indiv));
        return output;
    }

    /// Latent variable initialization.
    void simulate_single_house(Household &house) override final
    {
		for (MemberIdx i = 0; i < house.size(); ++i)
		{
			Individual &indiv = house.at(i);
			if (indiv.is_infected_status())
				init_missing_latent(indiv);
		}
    }
};

HOUSEPI_CLOSE_MODULE_NAMESPACE()
