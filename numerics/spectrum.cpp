#include <iostream>
#include <string>
#include <math.h>
#include "myfloat.h"
#include "metric.h"
#include "spectrum.h"
#include "misclib.h"


using namespace std;



// Konstruktor
Spectrum::Spectrum(const std::vector<myfloat> freqs, const Metric *metric)
	: metric(metric), freqmult(freqs)
{
	cnts.resize(freqmult.size());
	for(unsigned i = 0; i < cnts.size(); i++) {
		cnts[i] = 0.0;
	}
}


// Destruktor
Spectrum::~Spectrum(){}



// increase counts
int
Spectrum::inc_cnts(myfloat *x, myfloat *u)
{
	unsigned i = cnts.size();
	myfloat r = x[1];
	// FIXME: We only consider a small belt
	if ((r < 4.8 * metric->m) && (r > 3.3 * metric->m)){
		// source movement
		myfloat g00 = (metric->*metric->g[0][0])(x);
		myfloat g03 = (metric->*metric->g[0][3])(x);
		myfloat g33 = (metric->*metric->g[3][3])(x);
		myfloat u0 = 0.0, u3 = 0.0;
		// stationary src:
		u0 = 1.0 / sqrt(g00);
		// derived from schwarzschild, stable orbit:
		u0 = 1.0 / sqrt(1.0 - 3.0 * metric->m / r);
		u3 = (sqrt(g33 + u0 * u0 * (g03 * g03 - g00 * g33)) - g03 * u0) / g33;

		//u3 = 0.0;
		myfloat src_u[DIM] = { u0, 0.0, 0.0, u3 };
		myfloat E_over_k = scalar(*metric, x, src_u, u);
		cerr << "E_over_k = " << E_over_k << endl;
		/*cerr << "g00 = " << g00
			<< ", g03 = " << g03
			<< ", g33 = " << g33
			<< ", u0 = " << u0
			<< ", u3 = " << u3
			<< ", u^2 = " << scalar(*metric, x, src_u, src_u)
			<< endl;
			*/

		while (i--) {
			// FIXME: This is only a very rudimentary gaussian
			myfloat energy = freqmult[i] * E_over_k;
			size_t last = freqmult.size() - 1;
			myfloat dE = (freqmult[last] - freqmult[0]) / last;
			dE *= E_over_k;
			cnts[i] += dE * energy * exp(-10.0 * (energy - 5.9) * (energy - 5.9));
		}

		return 1;
	}
	return 0;
}
