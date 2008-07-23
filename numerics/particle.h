/*
 * HSGG Tü2007.04.29
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include "myfloat.h"

class Particle
{
    public:
	// Konstruktor
	Particle(const myfloat& masse, const myfloat& ladung,
		const myfloat pos[DIM], const myfloat vel[DIM]);
	// Destruktor
	~Particle();

	// Konstanten
	const myfloat m;
	const myfloat q;

	// Position, Vierergeschwindigkeit
	myfloat x[DIM];
	myfloat u[DIM];

};

#endif
