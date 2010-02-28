#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>

#include "metric.h"
#include "myfloat.h"
#include "particle.h"


class Spectrum
{
    public:
	// Konstruktor
	Spectrum(const std::vector<myfloat> freqs, const Particle *part,
			const Metric *metric);
	// Destruktor
	~Spectrum();

	void inc_cnts();

	std::vector<myfloat> cnts;

    private:
	const Metric *metric;
	const Particle *particle;
	const std::vector<myfloat> freqmult; // E = u0 * freqmult
};

#endif
