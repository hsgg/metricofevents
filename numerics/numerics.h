#ifndef NUMERICS_H
#define NUMERICS_H

#include "global.h"
#include "myfloat.h"
#include "metric.h"
#include "particle.h"
#include "emfield.h"

using namespace std;


// Berechne x, u
// input: metric, x, u, dtau
// output:
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle);

#endif
