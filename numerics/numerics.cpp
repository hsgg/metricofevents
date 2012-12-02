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

static inline void mk_xk_uk(const Metric& metric, const EMField& emfield,
		const vector<myfloat>& xpk, const vector<myfloat>& upk,
		vector<vector<myfloat> >& xk, vector<vector<myfloat> >& uk,
		const int i, const myfloat dtau, const Particle& particle)
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


void init_rk_cache(struct rk_cache& rk_cache, unsigned dim)
{
	rk_cache.xk.resize(dim);
	rk_cache.uk.resize(dim);
	rk_cache.xpk.resize(dim);
	rk_cache.upk.resize(dim);
	for (unsigned j = 0; j < dim; j++) {
		rk_cache.xk[j].resize(dim);
		rk_cache.uk[j].resize(dim);
	}
}

// Berechne x, u
// input: metric, x, u, dtau
// output:
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle)
{
	struct rk_cache rk_cache;

	init_rk_cache(rk_cache, metric.dim);

	x_and_u(metric, emfield, dtau, particle, rk_cache);
}

void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle,
	struct rk_cache& rk_cache)
{
	for (unsigned j = 0; j < metric.dim; j++) {
		rk_cache.xk[0][j] = 0.0;
		rk_cache.uk[0][j] = 0.0;
	}

	//cerr << "in...";
	mk_cok(rk_cache.xpk, particle.x, rk_cache.xk[0], 0.0);
	mk_cok(rk_cache.upk, particle.u, rk_cache.uk[0], 0.0);
	//cerr << "with...";
	mk_xk_uk(metric, emfield, rk_cache.xpk, rk_cache.upk, rk_cache.xk, rk_cache.uk,
			0, dtau, particle);
	//cerr << "in\n";

	mk_cok(rk_cache.xpk, particle.x, rk_cache.xk[0], (myfloat) 0.5);
	mk_cok(rk_cache.upk, particle.u, rk_cache.uk[0], (myfloat) 0.5);
	mk_xk_uk(metric, emfield, rk_cache.xpk, rk_cache.upk, rk_cache.xk, rk_cache.uk,
			1, dtau, particle);


	mk_cok(rk_cache.xpk, particle.x, rk_cache.xk[1], (myfloat) 0.5);
	mk_cok(rk_cache.upk, particle.u, rk_cache.uk[1], (myfloat) 0.5);
	mk_xk_uk(metric, emfield, rk_cache.xpk, rk_cache.upk, rk_cache.xk, rk_cache.uk,
			2, dtau, particle);

	mk_cok(rk_cache.xpk, particle.x, rk_cache.xk[2], 1.0);
	mk_cok(rk_cache.upk, particle.u, rk_cache.uk[2], 1.0);
	mk_xk_uk(metric, emfield, rk_cache.xpk, rk_cache.upk, rk_cache.xk, rk_cache.uk,
			3, dtau, particle);

	// Berechne x, u
	for (unsigned mu = 0; mu < metric.dim; mu++) {
		particle.x[mu] = particle.x[mu]
			+ rk_cache.xk[0][mu]/6.0
			+ rk_cache.xk[1][mu]/3.0
			+ rk_cache.xk[2][mu]/3.0
			+ rk_cache.xk[3][mu]/6.0;
		particle.u[mu] = particle.u[mu]
			+ rk_cache.uk[0][mu]/6.0
			+ rk_cache.uk[1][mu]/3.0
			+ rk_cache.uk[2][mu]/3.0
			+ rk_cache.uk[3][mu]/6.0;
	}
}
