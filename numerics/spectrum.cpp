#include <iostream>
#include <string>
#include <math.h>
#include "myfloat.h"
#include "metric.h"
#include "spectrum.h"


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
		// FIXME: This is only a very rudimentary gaussian in a band
		if ((r < 4.0 * metric->m) && (r > 2.75 * metric->m)){
			myfloat energy = freqmult[i] * particle->u[0];
			cnts[i] += exp(-1 * (energy - 5.9) * (energy - 5.9));
		}
	}
}
