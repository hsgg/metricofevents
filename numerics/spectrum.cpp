#include <iostream>
#include <string>
#include <math.h>
#include "myfloat.h"
#include "metric.h"
#include "spectrum.h"
#include "misclib.h"


using namespace std;



// Konstruktor
Spectrum::Spectrum(const std::vector<myfloat> freqs, const Particle *part,
		const Metric *metric)
	: metric(metric), particle(part), freqmult(freqs)
{
	cnts.resize(freqmult.size());
	for(unsigned i = 0; i < cnts.size(); i++) {
		cnts[i] = 0.0;
	}
}


// Destruktor
Spectrum::~Spectrum(){}



// increase counts
void
Spectrum::inc_cnts()
{
	unsigned i = cnts.size();
	while (i--) {
		myfloat r = particle->x[1];
		// FIXME: We only consider a small belt
		if ((r < 4.0 * metric->m) && (r > 2.75 * metric->m)){
			// source movement
			myfloat src_u[DIM] = {
				1.0 / sqrt((metric->*metric->g[0][0])
						(particle->x)),
				0.0,
				0.0,
				0.0
			};
			myfloat E_over_k = scalar(*metric, particle->x,
					src_u, particle->u);

			// FIXME: This is only a very rudimentary gaussian
			myfloat energy = freqmult[i] * E_over_k;
			cnts[i] += exp(-105 * (energy - 5.9) * (energy - 5.9));
		}
	}
}
