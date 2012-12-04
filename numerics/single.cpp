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

	// rk_cache
	struct rk_cache rk_cache;
	init_rk_cache(rk_cache, metric.dim);

	// particle
	if (init.uvec[0][0] != 1 || init.uvec[1][0] != 1) {
		cerr << "Error: u[0][0] = " << init.uvec[0][0] << " must be exactly 1!" << endl;
		cerr << "Error: u[1][0] = " << init.uvec[1][0] << " must be exactly 1!" << endl;
		return 1;
	}
	unsigned numparticles = 2;
	vector<Particle*> particle;
	for (unsigned i = 0; i < numparticles; i++) {
		particle.push_back(new Particle(i, init.mass[i], init.charge[i],
				init.xvec[i], init.uvec[i], metric));
		myfloat gammacheck = scalar(metric, particle[i]->x, particle[i]->u, particle[i]->u);
		if (gammacheck <= 0.0) {
			cerr << "Error: gammacheck = " << gammacheck << " must be positive!" << endl;
			return 1;
		}

		cout << endl;
		for (unsigned mu = 0; mu < metric.dim; mu++)
		{
			cout << "x" << mu << " = " << particle[i]->x[mu] << endl;
		}
		cout << endl;
		for (unsigned mu = 0; mu < metric.dim; mu++)
		{
			cout << "u" << mu << " = " << particle[i]->u[mu] << endl;
		}
		cout << "gamma = " << gammafactor(metric, *particle[i]) << endl;
		cout << endl;
	}


	// Test antisymmetriy of fieldstrengthtensor
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		for (unsigned nu = 0; nu < metric.dim; nu++)
		{
			vector<myfloat> teilchenx(metric.dim);
			for (unsigned i = 0; i < metric.dim; i++)
				teilchenx[i] = particle[0]->x[i];
			cout << mu << nu << ": "
			     << (emfield.*(emfield.F[mu][nu]))(teilchenx)
			     << " ---- "
			     << (emfield.*(emfield.F[nu][mu]))(teilchenx)
			     << endl;
		}
	}



	// Iterate
	myfloat dtau = init.dtau;
	ofstream dtaufile("dtau.dat");
	dtaufile << scientific;
	cout << endl << metric.name << ":" << endl
	     << "m = " << metric.m << endl
	     << "a = " << metric.a << endl
	     << "q = " << metric.q << endl
	     << endl;


	// Iterate
	int taun = -1;
	myfloat tau = 0.0;
	while (tau <= init.tau_max + dtau)
	{
		if (++taun % 100 == 0) {
			info(metric, taun, tau, dtau, 0, *particle[0]);
			if (taun >= init.max_n) {
				cerr << "WARNING: Aborting prematurely: taking too long" << endl;
				break;
			}
		}

		// print to file
		if (taun % 100) {
			for (unsigned i = 0; i < particle.size(); i++)
				particle[i]->write_to_plotfile();
			dtaufile << tau << "\t" << dtau << endl;
		}

		tau += dtau;


		// Berechne x, u
		x_and_u(metric, emfield, dtau, particle, rk_cache);

		if ((particle[0]->x[1] > init.max_x1)
				|| (particle[0]->x[1] < init.min_x1))
			break;
	}
	dtaufile.close();


	info(metric, taun, tau, dtau, 0, *particle[0]);

	cout << endl;

	return 0;
}
