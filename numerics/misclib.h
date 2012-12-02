#ifndef MISCLIB_H
#define MISCLIB_H

#include "init.h"
#include "particle.h"

// absolute value of myfloat
inline myfloat absol(const myfloat number)
{
	return (number < 0.0) ? -number : number;
}


// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric metric, const myfloat* x,
	const myfloat* y, const myfloat* z);

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


#endif
