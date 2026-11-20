#include "housepi.h"

HOUSEPI_OPEN_MODULE_NAMESPACE()

/// Exponential growth `initial * exp(rate * x)`, a custom function.
class Growth final : public Cloneable<Growth, Function>
{
public:
    static constexpr const char *NAME = "growth";
    Growth() : Cloneable({{"initial", 1.}, {"rate", 1.}}, 0.) {}

    /// Set the domain to x >= 0.
    void truncate() override final { helper_init_domain_min_(0.); }

    /// Bound the values to (0, inf), which proves a positive model requirement.
    Codomain codomain() const noexcept override final
    {
        Codomain c;
        c.inf = 0.;
        c.below = Interval::Open;
        return c;
    }

protected:
    /// Require a positive `initial`.
    ArgRequirements arg_requirements_() const override final { return {{&arguments_(0), Codomain::positive()}}; }

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
/// Register the function and make it available by the name `growth`.
HOUSEPI_REGISTER_FUNCTION(Growth, Custom)

/// Community-only final size: each member is infected from outside with probability 1 - exp(-foi).
class Model final : public ModelBase
{
public:
    /// Keys of the definitions the model requires.
    HOUSEPI_DECLARE_PARAM_KEYS(foi)

    /// Require a positive force of infection.
    void init_model_requirements() override final { request_param_positive({Param::foi}); }

    /// Compute the transmission log-likelihood of `house`.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        const double log_esc = -par(Param::foi);
        const double log_inf = std::log(-std::expm1(log_esc));
        return house.n_infected_status() * log_inf + house.n_escaped_status() * log_esc;
    }

    /// Infect each member of `house` from outside.
    void simulate_single_house(Household &house) override final
    {
        const double foi = par(Param::foi);
        for (Individual &indiv : house)
            if (exponential() > foi)
                set_escaped(house, indiv);
            else
                set_infected(house, indiv);
    }
};

HOUSEPI_CLOSE_MODULE_NAMESPACE()
