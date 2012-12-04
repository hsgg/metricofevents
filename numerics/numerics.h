#ifndef NUMERICS_H
#define NUMERICS_H

#include "global.h"
#include "myfloat.h"
#include "metric.h"
#include "particle.h"
#include "emfield.h"

using namespace std;


struct rk_cache {
	vector< vector<myfloat> > xk;
	vector< vector<myfloat> > uk;

	vector<myfloat> xpk;
	vector<myfloat> upk;
};

void init_rk_cache(struct rk_cache& rk_cache, unsigned dim);

// Berechne x, u
// input: metric, x, u, dtau
// output:
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle);
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, vector<Particle*>& particle,
	struct rk_cache& rk_cache);
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle,
	struct rk_cache& rk_cache);

#endif
