#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>

#include "metric.h"
#include "myfloat.h"


class Spectrum
{
    public:
	// Konstruktor
	Spectrum(const std::vector<myfloat> freqs, const Metric *metric);
	// Destruktor
	~Spectrum();

	int inc_cnts(myfloat *x, myfloat *u);

	std::vector<myfloat> cnts;

    private:
	const Metric *metric;
	const std::vector<myfloat> freqmult; // E = u0 * freqmult
};

#endif
