/*
 * HSGG Tü20070301
 */


#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H


#include "global.h"
#include "myfloat.h"



struct initializations {
	myfloat m;
	myfloat a;
	myfloat q;

	myfloat radius;

	myfloat teilchen_masse;
	myfloat teilchen_ladung;

	myfloat x[DIM];
	myfloat u[DIM];

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

#endif
