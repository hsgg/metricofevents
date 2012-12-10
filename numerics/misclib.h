#ifndef MISCLIB_H
#define MISCLIB_H

#include <vector>
#include "init.h"
#include "particle.h"

// absolute value of myfloat
inline myfloat absol(const myfloat number)
{
	return (number < 0.0) ? -number : number;
}

inline myfloat value_of_metric(const Metric metric, const myfloat* x,
		const int mu, const int nu)
{
	return (metric.*(metric.g[mu][nu]))(x);
}



// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric metric, const myfloat* x,
	const myfloat* y, const myfloat* z);

// Calculate spatial projection of x^mu * y^nu
myfloat spatial_projection_scalar(const Metric metric, const Particle& p,
		const myfloat* x, const myfloat* y);
myfloat spatial_projection_scalar(const Metric metric, const myfloat* xpos,
		const myfloat* covar_u, const myfloat* x, const myfloat* y);

myfloat gammafactor(const Metric& metric, const Particle& particle);

// Falschheit
inline myfloat wrongness(const Metric metric,
		const struct initializations* init,
		const myfloat* x, const myfloat* u)
{
	return scalar(metric, x, u, u) - init->teilchen_masse;
}


// Calculate u
myfloat find_u(const Metric metric,
		const struct initializations* init,
		const myfloat* x,
		myfloat* u);


// Print some info
void info(const Metric& metric, const int taun, const myfloat tau, const myfloat dtau,
		const myfloat wrong, const Particle& p);
void info(const Metric& metric, const int taun, const myfloat tau, const myfloat dtau,
		const myfloat wrong, const std::vector<Particle*>& p);


#endif
