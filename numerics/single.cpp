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
#include "particle.h"
#include "emfield.h"
#include "numerics.h"

using namespace std;


// absolute value of myfloat
inline myfloat absol(const myfloat number)
{
	return (number < 0.0) ? -number : number;
}


// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric metric, const myfloat* x,
	const myfloat* y, const myfloat* z)
{
	myfloat l = 0.0;

	for (unsigned mu = 0; mu < metric.dim; mu++)
		for (unsigned nu = 0; nu < metric.dim; nu++)
			l += (metric.*(metric.g[mu][nu]))(x) * y[mu] * z[nu];

	return l;
}


// Falschheit
inline myfloat wrongness(const Metric metric,
		const struct initializations* init,
		const myfloat* x, const myfloat* u)
{
	return scalar(metric, x, u, u) - init->teilchen_masse;
}


// Calculate u
myfloat find_u(const Metric metric,
		const struct initializations* init,
		const myfloat* x,
		myfloat* u)
{
	myfloat min = init->umin;
	myfloat max = init->umax;
	myfloat un = u[init->change];
	myfloat old;

	myfloat fmin, fmax, fcurr;

	int i = 0;

	int error = 0;

	cout << scientific;

	do
	{
		old = un;

		u[init->change] = min;
		fmin = wrongness(metric, init, x, u);

		u[init->change] = max;
		fmax = wrongness(metric, init, x, u);

		u[init->change] = un;
		fcurr = wrongness(metric, init, x, u);


		if (fmin * fmax > 0.0)
		{
			cerr << "ERROR: fmin and fmax have same sign "
			     << "(i = " << i << "):" << endl
			     << "fmin = " << fmin << endl
			     << "fmax = " << fmax << endl;
			error = 1;
			break;
		}

		// reset boundaries
		if ( (absol(fcurr) < absol(fmax)) && (fcurr * fmax > 0.0) )
		{
			max = un;
		}
		else if ( (absol(fcurr) < absol(fmin)) && (fcurr * fmin > 0.0) )
		{
			min = un;
		}

		// reset un
		un = (min + max) / (myfloat) 2.0;

		if (i++ > 100000)
		{
			cerr << "WARNING: Calculation of u takes suspiciously "
			     << "long. Terminated early." << endl;
			break;
		}
	} while (un - old != 0.0);

	cout << "u" << init->change << "(" << i << ") = " << u[init->change]
	     << ", u^2(x) = " << scalar(metric, x, u, u)
	     << ", wrongness = " << fcurr << endl;

	if (error)
	{
		exit(-1);
	}

	return u[init->change];
}



// Print some info
inline void info(const int taun, const myfloat tau, const myfloat dtau,
		const myfloat wrong, const myfloat* x)
{
	cerr << taun
	     << ": tau = " << tau
	     << ", dtau = " << dtau
	     << ", wrong = " << wrong
	     << ", r = " << setprecision(12) << x[1] << setprecision(3)
	     << endl;
}



// ********* MAIN ************
int main()
{
	struct initializations init = initials();

	// spacetime
	Metric metric(init.m, init.a, init.q);

	// emfield
	EMField emfield(metric);

	// particle
	init.u[init.change] = find_u(metric, &init, init.x, init.u);
	Particle particle(init.teilchen_masse, init.teilchen_ladung,
		init.x, init.u);

	cout << endl;
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		cout << "x" << mu << " = " << particle.x[mu] << endl;
	}
	cout << endl;
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		cout << "u" << mu << " = " << particle.u[mu] << endl;
	}
	cout << endl;


	// Test antisymmetriy of fieldstrengthtensor
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		for (unsigned nu = 0; nu < metric.dim; nu++)
		{
			vector<myfloat> teilchenx(metric.dim);
			for (unsigned i = 0; i < metric.dim; i++)
				teilchenx[i] = particle.x[i];
			cout << mu << nu << ": "
			     << (emfield.*(emfield.F[mu][nu]))(teilchenx)
			     << " ---- "
			     << (emfield.*(emfield.F[nu][mu]))(teilchenx)
			     << endl;
		}
	}



	// Iterate
	myfloat dtau = init.dtau;
	myfloat wrong = wrongness(metric, &init, particle.x, particle.u);
	myfloat wrong_old = wrong;
	ofstream plotfile("plot.dat");
	ofstream dtaufile("dtau.dat");
	ofstream wrongfile("wrong.dat");
	plotfile << scientific;
	dtaufile << scientific;
	wrongfile << scientific;
	cout << endl << metric.name << ":" << endl
	     << "m = " << metric.m << endl
	     << "a = " << metric.a << endl
	     << "q = " << metric.q << endl
	     << endl;

	// Print important information to plotfile
	plotfile << "# Metric: " << metric.name << endl;
	plotfile << "#	m = " << metric.m << endl
		<< "#	a = " << metric.a << endl
		<< "#	q = " << metric.q << endl;
	plotfile << "# particle:" << endl
		<< "#	m = " << init.teilchen_masse << endl
		<< "#	q = " << init.teilchen_ladung << endl;
	plotfile << "# length_4velocity = "
		<< scalar(metric, particle.x, particle.u, particle.u) << endl
		<< "# dtau = " << dtau << endl
		<< "# tau_max = " << init.tau_max << endl
		<< "# max_wrongness = " << init.max_wrongness << endl
		<< "# max_x1 = " << init.max_x1 << endl;
	plotfile << "# x3 x1 x2 x0	wrong	u3 u1 u2 u0" << endl;


	// Iterate
	int taun = -1;
	myfloat tau = 0.0;
	while ((absol(wrong) <= init.max_wrongness)
		&& (tau <= init.tau_max + dtau))
	{
		dtau = init.dtau * exp(-1e15 * absol(wrong - wrong_old) / init.max_wrongness);

		if (++taun % 1000 == 0) {
			info(taun, tau, dtau, wrong, particle.x);
			if (taun >= 10000) {
				cerr << "WARNING: Aborting prematurely: taking too long" << endl;
				break;
			}
		}

		// print to file
		plotfile<< particle.x[3] << ' '
			<< particle.x[1] << ' '
			<< particle.x[2] << ' '
			<< particle.x[0] << '\t'
			<< wrong	 << '\t'
			<< particle.u[3] << ' '
		        << particle.u[1] << ' '
		        << particle.u[2] << ' '
		        << particle.u[0] << '\t'
			<< tau << '\t'
			<< dtau
			<< endl;

		dtaufile << tau << "\t" << dtau << endl;
		wrongfile << tau << "\t" << wrong << endl;

		tau += dtau;


		// Berechne x, u
		x_and_u(metric, emfield, dtau, particle);

		wrong_old = wrong;
		wrong = wrongness(metric, &init, particle.x, particle.u);

		if ((particle.x[1] > init.max_x1)
				|| (particle.x[1] < init.min_x1))
			break;
	}
	plotfile << endl;
	plotfile << "# length_4velocity = "
		 << scalar(metric, particle.x, particle.u, particle.u)
		 << endl;
	plotfile.close();
	dtaufile.close();
	wrongfile.close();


	info(taun, tau, dtau, wrong, particle.x);

	cout << endl;

	return 0;
}
