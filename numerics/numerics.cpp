/*
 * HSGG Tü20070228
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

#include "global.h"
#include "myfloat.h"
#include "init.h"
#include "metric.h"
#include "particle.h"
#include "emfield.h"
#include "numerics.h"

using namespace std;



// set xpk = x + step * k
static inline void mk_cok(vector<myfloat>& xpk, const myfloat* x,
		const vector<myfloat>& xk_n, const myfloat step)
{
	for (unsigned mu = 0; mu < xpk.size(); mu++)
		xpk[mu] = x[mu] + step * xk_n[mu];
}


static inline myfloat christoffelsum(const Metric& metric, const unsigned& sigma,
	const vector<myfloat>& x, const vector<myfloat>& u)
{
	return (metric.*(metric.christoffelsum[sigma]))(x, u);
}

static inline myfloat acceleration3d(const Metric& metric, const unsigned& sigma,
	const vector<myfloat>& x, const vector<myfloat>& u)
{
	return (metric.*(metric.acceleration3d[sigma]))(x, u);
}



// electromagnetic acceleration
static inline myfloat emfieldforce(const EMField& emfield, const unsigned& sigma,
	const vector<myfloat>& x, const vector<myfloat>& u)
{
	myfloat force = 0.0;
	const Metric* metric = emfield.metric;

	for (unsigned nu = 0; nu < metric->dim; nu++)
	{
		for (unsigned mu = 0; mu < metric->dim; mu++)
		{
			force += (metric->*(metric->G[sigma][nu]))(x)
				* (emfield.*(emfield.F[nu][mu]))(x) * u[mu];
		}
	}

	return force;
}

static inline void mk_xk_uk(const Metric metric, const EMField emfield,
		const vector<myfloat> xpk, const vector<myfloat> upk,
		vector<vector<myfloat> >& xk, vector<vector<myfloat> >& uk,
		const int i, const myfloat dtau, const Particle particle)
{
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		xk[i][mu] = dtau * upk[mu];
		uk[i][mu] = dtau
			* (-acceleration3d(metric, mu, xpk, upk)
					- particle.q
					* emfieldforce(emfield, mu, xpk, upk));
	}
}



// Berechne x, u
// input: metric, x, u, dtau
// output:
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle)
{
	vector< vector<myfloat> > xk(4);
	vector< vector<myfloat> > uk(4);

	vector<myfloat> xpk(metric.dim);
	vector<myfloat> upk(metric.dim);

	for (unsigned j = 0; j < metric.dim; j++) {
		xk[j].resize(metric.dim);
		uk[j].resize(metric.dim);
		xk[0][j] = 0.0;
		uk[0][j] = 0.0;
	}

	//cerr << "in...";
	mk_cok(xpk, particle.x, xk[0], 0.0);
	mk_cok(upk, particle.u, uk[0], 0.0);
	//cerr << "with...";
	mk_xk_uk(metric, emfield, xpk, upk, xk, uk, 0, dtau, particle);
	//cerr << "in\n";

	mk_cok(xpk, particle.x, xk[0], (myfloat) 0.5);
	mk_cok(upk, particle.u, uk[0], (myfloat) 0.5);
	mk_xk_uk(metric, emfield, xpk, upk, xk, uk, 1, dtau, particle);


	mk_cok(xpk, particle.x, xk[1], (myfloat) 0.5);
	mk_cok(upk, particle.u, uk[1], (myfloat) 0.5);
	mk_xk_uk(metric, emfield, xpk, upk, xk, uk, 2, dtau, particle);

	mk_cok(xpk, particle.x, xk[2], 1.0);
	mk_cok(upk, particle.u, uk[2], 1.0);
	mk_xk_uk(metric, emfield, xpk, upk, xk, uk, 3, dtau, particle);

	// Berechne x, u
	for (unsigned mu = 0; mu < metric.dim; mu++) {
		particle.x[mu] = particle.x[mu]
			+ xk[0][mu]/6.0
			+ xk[1][mu]/3.0
			+ xk[2][mu]/3.0
			+ xk[3][mu]/6.0;
		particle.u[mu] = particle.u[mu]
			+ uk[0][mu]/6.0
			+ uk[1][mu]/3.0
			+ uk[2][mu]/3.0
			+ uk[3][mu]/6.0;
	}
}
