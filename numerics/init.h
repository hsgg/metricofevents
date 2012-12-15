/*
 * HSGG Tü20070301
 */


#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H


#include "global.h"
#include "myfloat.h"

#define NUMPARTICLES 2


struct initializations {
	myfloat m;
	myfloat a;
	myfloat q;

	myfloat radius;

	myfloat teilchen_masse; // compatibility
	myfloat teilchen_ladung; // compatibility
	myfloat x[DIM]; // compatibility
	myfloat u[DIM]; // compatibility

	unsigned numparticles;
	myfloat mass[NUMPARTICLES];
	myfloat charge[NUMPARTICLES];
	myfloat xvec[NUMPARTICLES][DIM];
	myfloat uvec[NUMPARTICLES][DIM];

	int nrays;
	myfloat u3_inc;
	myfloat umin_inc;
	myfloat umax_inc;

	int change;
	myfloat umin;
	myfloat umax;

	myfloat dtau;
	myfloat tau_max;
	int max_n;
	myfloat max_wrongness;
	myfloat min_x1;
	myfloat max_x1;
};

struct initializations initials();
struct initializations initialize(char const*const filename);

#endif
