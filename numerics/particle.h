#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>


#include "myfloat.h"
#include "metric.h"


class Particle
{
    public:
	// Konstruktor
	Particle(int unique, const myfloat& masse, const myfloat& ladung,
			const myfloat pos[DIM], const myfloat vel[DIM],
		       	const Metric& metric);
	// Destruktor
	~Particle();

	// Konstanten
	const myfloat m;
	const myfloat q;

	// Position, Vierergeschwindigkeit
	myfloat x[DIM];
	myfloat u[DIM];

	// for Runge-Kutta
	myfloat xk[4][DIM];
	myfloat uk[4][DIM];
	myfloat xpk[DIM];
	myfloat upk[DIM];

	void write_to_plotfile(const std::vector<Particle*>& particles);

    private:
	int const unique;
	std::ofstream plotfile;
	const Metric& metric;
};

#endif
