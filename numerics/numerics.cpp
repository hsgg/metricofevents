/*
 * HSGG Tü20070228
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

#include "global.h"
#include "myfloat.h"
#include "init.h"
#include "metric.h"
#include "emfield.h"
#include "numerics.h"

using namespace std;



// set xpk = x + step * k
static inline void mk_cok(const unsigned dim, myfloat* xpk, const myfloat* x,
		const myfloat* xk_n, const myfloat step)
{
	for (unsigned mu = 0; mu < dim; mu++)
		xpk[mu] = x[mu] + step * xk_n[mu];
}

static void mk_pk(vector<Particle*>& particle, const int k, const myfloat step)
{
	for (unsigned i = 0; i < particle.size(); i++) {
		// for x
		mk_cok(DIM, particle[i]->xpk, particle[i]->x,
				particle[i]->xk[k], step);
		// for u
		mk_cok(DIM, particle[i]->upk, particle[i]->u,
				particle[i]->uk[k], step);
	}
}


static inline myfloat christoffelsum(const Metric& metric, const unsigned& sigma,
	const myfloat* x, const myfloat* u)
{
	return (metric.*(metric.christoffelsum[sigma]))(x, u);
}

static inline myfloat acceleration3d(const Metric& metric, const unsigned& sigma,
	const myfloat* x, const myfloat* u)
{
	return (metric.*(metric.acceleration3d[sigma]))(x, u);
}

static myfloat local_gravity_acceleration3d(const Metric& metric, const unsigned mu,
		const vector<Particle*>& particle, const int i)
{
	return 0.0;
}


// electromagnetic acceleration
static inline myfloat emfieldforce(const EMField& emfield, const unsigned& sigma,
	const myfloat* x, const myfloat* u)
{
	myfloat force = 0.0;
	const Metric* metric = emfield.metric;

	for (unsigned nu = 0; nu < metric->dim; nu++)
	{
		for (unsigned mu = 0; mu < metric->dim; mu++)
		{
			force += (metric->*(metric->G[sigma][nu]))(x)
				* (emfield.*(emfield.F[nu][mu]))(x) * u[mu];
		}
	}

	return force;
}


static myfloat full_acceleration(const int mu, const Metric& metric, const EMField& emfield,
		vector<Particle*>& particle, unsigned i)
{
	Particle *const p = particle[i];
	return acceleration3d(metric, mu, p->xpk, p->upk)
		+ local_gravity_acceleration3d(metric, mu, particle, i)
		- p->q * emfieldforce(emfield, mu, p->xpk, p->upk);
}

static inline void mk_xk_uk(const Metric& metric, const EMField& emfield,
		const int k, const myfloat dtau, vector<Particle*> particle)
{
	for (unsigned i = 0; i < particle.size(); i++)
	{
		Particle* p = particle[i];
		for (unsigned mu = 0; mu < metric.dim; mu++)
		{
			p->xk[k][mu] = dtau * p->upk[mu];
			p->uk[k][mu] = dtau * full_acceleration(mu, metric, emfield, particle, i);
		}
	}
}



// Berechne x, u
// input: metric, x, u, dtau
// output:
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle)
{
	vector<Particle*> p(1);

	p[0] = &particle;

	x_and_u(metric, emfield, dtau, p);
}


// This is Runge-Kutta
void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, vector<Particle*>& particles)
{
	for (unsigned i = 0; i < particles.size(); i++) {
		for (unsigned j = 0; j < metric.dim; j++) {
			particles[i]->xk[0][j] = 0.0;
			particles[i]->uk[0][j] = 0.0;
		}
	}

	mk_pk(particles, 0, 0.0);
	mk_xk_uk(metric, emfield, 0, dtau, particles);

	mk_pk(particles, 0, 0.5);
	mk_xk_uk(metric, emfield, 1, dtau, particles);

	mk_pk(particles, 1, 0.5);
	mk_xk_uk(metric, emfield, 2, dtau, particles);

	mk_pk(particles, 2, 1.0);
	mk_xk_uk(metric, emfield, 3, dtau, particles);

	// Calculate x, u
	for (unsigned i = 0; i < particles.size(); i++) {
		Particle *const p = particles[i];
		for (unsigned mu = 0; mu < metric.dim; mu++) {
			p->x[mu] = p->x[mu]
				+ p->xk[0][mu]/6.0
				+ p->xk[1][mu]/3.0
				+ p->xk[2][mu]/3.0
				+ p->xk[3][mu]/6.0;
			p->u[mu] = p->u[mu]
				+ p->uk[0][mu]/6.0
				+ p->uk[1][mu]/3.0
				+ p->uk[2][mu]/3.0
				+ p->uk[3][mu]/6.0;
		}
	}
}
