#include <distribution.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// InfecProfile_Cortier25 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class InfecProfile_Cortier25 : public Distribution
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Cortier25()
		: Distribution("cortier25",
					   {"shape", "refer_vl", "base_vl", "max_vl",
						"prolif_period", "clear_period"},
					   {2.5, 0.5, -2, 8, 2.4, 3.1},
					   0., HSP_NA_DOUBLE, false, false) {}
	~InfecProfile_Cortier25() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		// FIXME: I don't want truncation renormalization since no distribution
		init_min_zero_reshaping(discrete, false);
	}

protected:
	double raw_pdf(double x) const override final
	{
		const double s = arg(0);				// 'shape'
		const double vref = arg(1);				// 'refer_vl'
		const double vl0 = arg(2) / vref;		// 'base_vl' / 'refer_vl'
		const double v = (arg(3) - vl0) / vref; // ('max_vl' - 'base_vl') / 'refer_vl'
		const double Tp = arg(4);				// 'prolif_period'
		const double Tc = arg(5);				// 'clear_period'

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v); // threshold VL > 0

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return pow(vl0 + v * dt, s);

		else if (dt < 1. + r * (1. - thres))
			return pow(vl0 + v * (1. + r - dt) / r, s);

		else
			return 0.;
	}

	double raw_cdf(double x) const override final
	{
		const double s = arg(0);				// 'shape'
		const double vref = arg(1);				// 'refer_vl'
		const double vl0 = arg(2) / vref;		// 'base_vl' / 'refer_vl'
		const double v = (arg(3) - vl0) / vref; // ('max_vl' - 'base_vl') / 'refer_vl'
		const double Tp = arg(4);				// 'prolif_period'
		const double Tc = arg(5);				// 'clear_period'

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v); // threshold VL > 0

		const double den = v * (s + 1.) / Tp;

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return (pow(vl0 + v * dt, s + 1.) - pow(vl0 + v * thres, s + 1.)) / den;

		else if (dt < 1. + r * (1. - thres))
		{
			const double F1 = pow(vl0 + v, s + 1.) - pow(vl0 + v * thres, s + 1.);
			return (F1 + r * (pow(vl0 + v, s + 1) - pow(vl0 + v * (1. + r - dt) / r, s + 1))) / den;
		}
		else
		{
			const double F1 = pow(vl0 + v, s + 1.) - pow(vl0 + v * thres, s + 1.);
			return F1 * (1 + r) / den;
		}
	}
};