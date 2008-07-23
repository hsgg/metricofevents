/*
 * HSGG Tü 2007.04.29
 */

#include "particle.h"


Particle::Particle(const myfloat& masse, const myfloat& ladung,
	const myfloat pos[DIM], const myfloat vel[DIM])
    : m(masse), q(ladung)
{
	for (int mu = 0; mu < DIM; mu++)
	{
		x[mu] = pos[mu];
		u[mu] = vel[mu];
	}
}


Particle::~Particle(){}
