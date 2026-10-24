
#include "function.h"
#include "housepi.h"
#include "mcmc.h"
#include "registry.h"

HOUSEPI_OPEN_MODULE_NAMESPACE()

/// Exponential growth `initial * exp(rate * x)`, a custom function.
class Growth_ final : public Cloneable<Growth_, Function>
{
public:
    static constexpr const char *NAME = "growth";
    Growth_() : Cloneable({{"initial", 1.}, {"rate", 1.}}, 0.) {}

    /// Set the domain to x >= 0.
    void truncate() override final { helper_init_domain_min_(0.); }

    /// Values returned.
    /// The engine can prove a model requirement this function satisfies;
    Codomain codomain() const noexcept override final
    {
        Codomain c;
        c.inf = 0.;
        c.below = Interval::Open;
        return c;
    }

protected:
    /// Evaluate the function at `x`.
    /// Arguments are read in declaration order.
    double raw_function_(double x) const override final
    {
        const double I0 = arg(0);
        const double r = arg(1);
        return I0 * std::exp(r * x);
    }

    /// Evaluate the primitive at `x`.
    double raw_primitive_(double x) const override final
    {
        const double I0 = arg(0);
        const double r = arg(1);
        return I0 * std::exp(r * x) / r;
    }
};
/// Register the function, available to definitions as `growth`.
HOUSEPI_REGISTER_FUNCTION(Growth_, Custom)

/// Custom model: one latent probability per infected member.
class Model : public ModelBase
{
public:
    /// Keys of the definitions the model requires.
    HOUSEPI_DECLARE_LATENT_KEYS(proba)
    HOUSEPI_DECLARE_PARAM_KEYS(alpha)

    /// Compute the transmission log-likelihood of `house`.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        double output = -par(Param::alpha);
        for (const Individual &indiv : house.infected_indivs())
            output += std::log(aug(Latent::proba, indiv));
        return output;
    }

    /// Initialize the latent variables of infected members.
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
