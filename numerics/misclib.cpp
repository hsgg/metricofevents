#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "global.h"
#include "myfloat.h"
#include "init.h"
#include "metric.h"
#include "misclib.h"

using namespace std;



// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric metric, const myfloat* x,
	const myfloat* y, const myfloat* z)
{
	myfloat l = 0.0;

	for (unsigned mu = 0; mu < metric.dim; mu++)
		for (unsigned nu = 0; nu < metric.dim; nu++)
			l += (metric.*(metric.g[mu][nu]))(x) * y[mu] * z[nu];

	return l;
}


// Print some info
void info(const int taun, const myfloat tau, const myfloat dtau,
		const myfloat wrong, const myfloat* x)
{
	cerr << taun
	     << ": tau = " << tau
	     << ", dtau = " << dtau
	     << ", wrong = " << wrong
	     << ", r = " << setprecision(12) << x[1] << setprecision(3)
	     << endl;
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
