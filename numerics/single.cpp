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



// ********* MAIN ************
int main(int argc, char *argv[])
{
	if (argc != 2) {
		cerr << "Error: Must have a filename on the command line." << endl;
		return 2;
	}

	struct initializations init = initialize(argv[1]);

	// spacetime
	Metric metric(init.m, init.a, init.q);

	// emfield
	EMField emfield(metric);

	// particle
	if (init.u[0] != 1) {
		cerr << "Error: u[0] = " << init.u[0] << " must be exactly 1!" << endl;
		return 1;
	}
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
	ofstream plotfile("plot.dat");
	ofstream dtaufile("dtau.dat");
	plotfile << scientific;
	dtaufile << scientific;
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
	plotfile << "# dtau = " << dtau << endl
		<< "# tau_max = " << init.tau_max << endl
		<< "# max_x1 = " << init.max_x1 << endl;
	plotfile << "# x3 x1 x2 x0	u3 u1 u2 u0" << endl;


	// Iterate
	int taun = -1;
	myfloat tau = 0.0;
	while (tau <= init.tau_max + dtau)
	{
		if (++taun % 1000 == 0) {
			info(taun, tau, dtau, 0, particle.x);
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
				<< particle.u[3] << ' '
				<< particle.u[1] << ' '
				<< particle.u[2] << ' '
				<< particle.u[0] << '\t'
				<< tau << '\t'
				<< dtau
				<< endl;

			dtaufile << tau << "\t" << dtau << endl;
		}

		tau += dtau;


		// Berechne x, u
		x_and_u(metric, emfield, dtau, particle);

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


	info(taun, tau, dtau, 0, particle.x);

	cout << endl;

	return 0;
}
