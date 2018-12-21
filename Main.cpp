#include "Main.h"

double fx_minsum(const SParams& params, double x, double c)
{
	double p = params.pmax * params.pmax / 2 - PMIN * PMIN / 2;
	return (1 - x) * (params.b1 + params.b2 * p) + x * (params.a + p * getT(params, x));
}

// integral from pmin to p(x(c)) (b1+b2p)*x'(p)dp +
// integral from p(x(c)) to pmax (a+pT(x(c)))*x'(p)dp
double fx_minsum_integral(const SParams& params, double x, double c)
{
	double pmiddle = getP(params, x, c);
	double result = 0;
	const double STEP = 0.0001;
	for (double p = PMIN; p < pmiddle; p += STEP)
		result += STEP * (params.b1 + params.b2 * p) / (p * p);
	for (double p = pmiddle; p < params.pmax; p += STEP)
		result += STEP * (params.a + getT(params, x)) / (p * p);

	return result;
}

double fx_moreparkingmoney(const SParams& params, double x, double c)
{
	return -c * x;
}

double fx_onlyhalfcars(const SParams& params, double x, double c)
{
	return fabs(x - 0.5);
}

int main()
{
	SParams params;

	params.a = 120;
	params.gamma = 40;
	params.T0 = 10;
	params.b1 = 30;
	params.b2 = 95;
	params.pmax = 10;
	
	calculate("minsum", &fx_minsum, params);
	calculate("minsum_integral", &fx_minsum_integral, params);
	calculate("moreparkingmoney", &fx_moreparkingmoney, params);
	calculate("onlyhalfcars", &fx_onlyhalfcars, params);

	return 0;
}