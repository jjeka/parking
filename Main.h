#pragma once

#include <cstdio>
#include <cmath>
#include <cassert>
#include <cfloat>

#define PMIN 1

struct SParams
{
	// Ap = (a+c) + pT(x) = (a+c) + (p ^ -1) * (T0 + gamma * x^4);
	double a = 60;
	double gamma = 2;
	double T0 = 70;
	// Bp = b1 + p*b2 = b1 + (p ^ -1) * b2
	double b1 = 50;
	double b2 = 75;

	double pmax = 10;
};

bool checkParams(const SParams& f_params)
{
	// b1 < a
	if (!(f_params.b1 < f_params.a))
	{
		printf("Failed: b1 < a");
		return false;
	}
	// T(1) < b2
	if (!(f_params.T0 + f_params.gamma < f_params.b2))
	{
		printf("Failed: T(1) < b2");
		return false;
	}
	// a + T(0) > b1 + b2
	if (!(f_params.a + f_params.T0 > f_params.b1 + f_params.b2))
	{
		printf("Failed: a + T(0) > b1 + b2");
		return false;
	}
	// a + pmax*T(1) < b1 + pmax*b2
	if (!(f_params.a + f_params.pmax * (f_params.T0 + f_params.gamma) < f_params.b1 + f_params.pmax * f_params.b2))
	{
		printf("Failed: a + pmax*T(1) < b1 + pmax*b2");
		return false;
	}
	return true;
}

double getT(const SParams& f_params, double f_x)
{
	return f_params.T0 + f_params.gamma * pow(f_x, 4.0);
}

double getP(const SParams& f_params, double f_x, double f_c)
{
	return (f_params.a + f_c - f_params.b1) / (f_params.b2 - (f_params.T0 + f_params.gamma * pow(f_x, 4.0)));
}

double getX(const SParams& f_params, double c)
{
	// personal transport in first day
	double x0 = 0.5;

	const double EPS = 0.0000000001;
	double x = x0 + 2 * EPS;

	unsigned i = 0;
	for (; fabs(x - x0) > EPS; i++)
	{
		double oldx = x0;
		x0 = x;

		double p = getP(f_params, oldx, c);
		x = 1.0 / p;
	}
	return x;
}

typedef double fx_func_t(const SParams& f_params, double f_x, double f_c);

bool calculate(const char* f_name, fx_func_t* f_func, const SParams& f_params)
{
	if (!checkParams(f_params))
		return false;

	double cmin = 0.0;
	double cmax = f_params.b1 + f_params.b2 - f_params.T0;
	const unsigned int numSteps = 50;
	double cstep = (cmax - cmin) / numSteps;

	printf("%s\ncmin = %lf; cmax = %lf; numSteps = %u; step = %lf\n\n", 
		f_name, cmin, cmax, numSteps, cstep);
	printf("STARTING %s\n", f_name);

	char name[256];
	snprintf(name, 256, "%s.csv", f_name);
	FILE* file = fopen(name, "w");
	assert(file);
	fprintf(file, "%s; cmin = %lf; cmax = %lf; numSteps = %u; step = %lf\n\n",
		f_name, cmin, cmax, numSteps, cstep);
	fprintf(file, "C; X; F(X(C))\n");
	
	double min = 0;
	double minval = DBL_MAX;
	for (double c = cmin; c < cmax; c += cstep)
	{
		double x = getX(f_params, c);
		
		double y = f_func(f_params, x, c);
		if (y < minval)
		{
			min = c;
			minval = y;
		}
		
		printf("C = %-7.1lf| X=%-10.2lf| F(X(C)) = %-10.2lf\n", c, x, y);
		fprintf(file, "%lf; %lf; %lf\n", c, x, y);
	}

	// find closest
	cmin = (min - cstep < cmin) ? cmin : (min - cstep);
	cmax = (min + cstep > cmax) ? cmax : (min + cstep);
	
	const double EPS = 0.000001;
	double cleftval = f_func(f_params, getX(f_params, cmin), cmin);
	double crightval = f_func(f_params, getX(f_params, cmax), cmax);
	while (cmax - cmin >= EPS)
	{
		double c = (cmin + cmax) / 2.0;
		double ccenterval = f_func(f_params, getX(f_params, c), c);

		if (ccenterval > cleftval && ccenterval > crightval)
			printf("WARNING: center val > left val & right val\n");

		if (cleftval < crightval)
		{
			crightval = ccenterval;
			cmax = c;
		}
		else
		{
			cleftval = ccenterval;
			cmin = c;
		}
	}

	double x = getX(f_params, cmin);
	printf("\nC FOR MIN: %-7.1lf; X FOR MIN: %-10.2lf; MIN F(X(C))=%-10.2lf\n", cmin, x, cleftval);
	printf("===========================================\n\n");

	fprintf(file, "\nC FOR MIN; X FOR MIN; MIN F(X(C))\n%lf; %lf; %lf\n\n", cmin, x, cleftval);
	fclose(file);

	return true;
}