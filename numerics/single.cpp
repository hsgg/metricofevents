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
#include "misclib.h"
#include "spectrum.h"



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

	// spectrum
	myfloat emin = 3.5;
	myfloat emax = 7.5;
	myfloat start_u0 = particle.u[0];
	vector<myfloat> freqmult(100);
	for (unsigned i = 0; i < freqmult.size(); i++){
		// E = u0 * freqmult
		myfloat energy = i * (emax - emin) / freqmult.size() + emin;
		freqmult[i] = energy / start_u0;
	}
	Spectrum spec(freqmult, &metric);

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
	ofstream specfile("spec.dat");
	plotfile << scientific;
	dtaufile << scientific;
	wrongfile << scientific;
	specfile << scientific;
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
		//dtau = init.dtau * exp(-1e15 * absol(wrong - wrong_old) / init.max_wrongness);
		//dtau = init.dtau * pow(absol(particle.x[1] - 2 * init.m), sqrt(3.0));
			//* (1.0 - wrong / init.max_wrongness);

		//dtau = dtau * (1.0 - wrong / init.max_wrongness);

		if (++taun % 1000 == 0) {
			info(taun, tau, dtau, wrong, particle.x);
			if (taun >= init.max_n) {
				cerr << "WARNING: Aborting prematurely: taking too long" << endl;
				break;
			}
		}

		// print to file
		if (taun % 100) {
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
		}

		tau += dtau;


		// Berechne x, u
		x_and_u(metric, emfield, dtau, particle);

		// Update spectrum
		spec.inc_cnts(particle.x, particle.u);

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

	for (unsigned i = 0; i < spec.cnts.size(); i++){
		specfile << freqmult[i] * start_u0 << '\t'
			<< spec.cnts[i] << endl;
	}

	plotfile.close();
	dtaufile.close();
	wrongfile.close();
	specfile.close();


	info(taun, tau, dtau, wrong, particle.x);

	cout << endl;

	return 0;
}
