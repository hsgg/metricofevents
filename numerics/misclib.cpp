#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "global.h"
#include "myfloat.h"
#include "init.h"
#include "metric.h"
#include "misclib.h"

using namespace std;


static myfloat value_of_metric(const Metric metric, const myfloat* x,
		const int mu, const int nu)
{
	return (metric.*(metric.g[mu][nu]))(x);
}

// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric metric, const myfloat* x,
	const myfloat* y, const myfloat* z)
{
	myfloat l = 0.0;

	for (unsigned mu = 0; mu < metric.dim; mu++)
		for (unsigned nu = 0; nu < metric.dim; nu++)
			l += value_of_metric(metric, x, mu, nu) * y[mu] * z[nu];

	return l;
}

// used for spatial_projection_scalar()
static myfloat delta_distance(const Metric metric, const Particle& p,
		const myfloat* covar_u,
		const myfloat* x, const myfloat* y,
		const int mu, const int nu)
{
	return (covar_u[mu] * covar_u[nu] - value_of_metric(metric, p.x, mu, nu))
		* (y[mu] - x[mu]) * (y[nu] - x[nu]);
}

// Calculate spatial projection of x^mu * y^nu
myfloat spatial_projection_scalar(const Metric metric, const Particle& p,
		const myfloat* x, const myfloat* y)
{
	unsigned mu, nu;
	myfloat distancesq = 0.0; // "distance squared"
	myfloat const gamma = gammafactor(metric, p);

	myfloat covar_u[DIM];
	for (mu = 0; mu < metric.dim; mu++) {
		covar_u[mu] = 0.0;
		for (nu = 0; nu < metric.dim; nu++) {
			covar_u[mu] += value_of_metric(metric, p.x, mu, nu) * gamma * p.u[nu];
		}
	}

	// can optimize sum, since dx^mu * dx^nu is symmetric
	for (mu = 0; mu < metric.dim; mu++) {
		distancesq += delta_distance(metric, p, covar_u, x, y, mu, mu);
		for (nu = mu + 1; nu < metric.dim; nu++) {
			distancesq += 2.0 * delta_distance(metric, p, covar_u, x, y, mu, nu);
		}
	}

	return distancesq;
}


myfloat gammafactor(const Metric& metric, const Particle& particle)
{
	return sqrt(1.0 / scalar(metric, particle.x, particle.u, particle.u));
}


// Print some info
void info(const Metric& metric, const int taun, const myfloat tau, const myfloat dtau,
		const myfloat wrong, const Particle& p)
{
	cerr << taun
	     << ": tau = " << tau
	     << ", dtau = " << dtau
	     << ", wrong = " << wrong
	     << ", r = " << setprecision(12) << p.x[1]
	     << ", gamma = " << scientific << gammafactor(metric, p)
	     << endl
	     << setprecision(3);
}


// Calculate u
myfloat find_u(const Metric metric,
		const struct initializations* init,
		const myfloat* x,
		myfloat* u)
{
	myfloat min = init->umin;
	myfloat max = init->umax;
	myfloat un = u[init->change];
	myfloat old;

	myfloat fmin, fmax, fcurr;

	int i = 0;

	int error = 0;

	cout << scientific;

	do
	{
		old = un;

		u[init->change] = min;
		fmin = wrongness(metric, init, x, u);

		u[init->change] = max;
		fmax = wrongness(metric, init, x, u);

		u[init->change] = un;
		fcurr = wrongness(metric, init, x, u);


		if (fmin * fmax > 0.0)
		{
			cerr << "ERROR: fmin and fmax have same sign "
			     << "(i = " << i << "):" << endl
			     << "fmin = " << fmin << endl
			     << "fmax = " << fmax << endl;
			error = 1;
			break;
		}

		// reset boundaries
		if ( (absol(fcurr) < absol(fmax)) && (fcurr * fmax > 0.0) )
		{
			max = un;
		}
		else if ( (absol(fcurr) < absol(fmin)) && (fcurr * fmin > 0.0) )
		{
			min = un;
		}

		// reset un
		un = (min + max) / (myfloat) 2.0;

		if (i++ > 100000)
		{
			cerr << "WARNING: Calculation of u takes suspiciously "
			     << "long. Terminated early." << endl;
			break;
		}
	} while (un - old != 0.0);

	cout << "u" << init->change << "(" << i << ") = " << u[init->change]
	     << ", u^2(x) = " << scalar(metric, x, u, u)
	     << ", wrongness = " << fcurr << endl;

	if (error)
	{
		exit(-1);
	}

	return u[init->change];
}
