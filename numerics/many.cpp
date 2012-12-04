#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <string>

#include "global.h"
#include "myfloat.h"
#include "init.h"
#include "metric.h"
#include "particle.h"
#include "emfield.h"
#include "numerics.h"
#include "misclib.h"
#include "spectrum.h"


static char const *filename = "calclen.dat";
static void truncate_file()
{
	ofstream tpfile(filename, ios::trunc);
	tpfile << "# nray	taun" << endl;
	tpfile.close();
}
static void append_file(int taun)
{
	static int i = -250;
	i++;
	ofstream tpfile(filename, ios::app);
	tpfile << scientific;
	tpfile << i << '\t' << taun << endl;
	tpfile.close();
}

int absorb;
int const absorb_max = 1;

static void go_ray(struct initializations &init, Spectrum &spec,
		Metric &metric, Particle &particle, EMField &emfield,
		string &plotfilename)
{
	myfloat dtau = init.dtau;

	ofstream plotfile(plotfilename.c_str());
	ofstream dtaufile("globaldtau.dat");
	ofstream wrongfile("globalwrong.dat");
	plotfile << scientific;
	dtaufile << scientific;
	wrongfile << scientific;

	// Print important information to plotfile
	plotfile << "# Metric: " << metric.name << endl;
	plotfile << "#	m = " << metric.m << endl
		<< "#	a = " << metric.a << endl
		<< "#	q = " << metric.q << endl;
	plotfile << "# particle:" << endl
		<< "#	m = " << init.teilchen_masse << endl
		<< "#	q = " << init.teilchen_ladung << endl;
	plotfile << "# length_4velocity = "
		<< scalar(metric, particle.x, particle.u, particle.u)
		<< endl
		<< "# dtau = " << dtau << endl
		<< "# tau_max = " << init.tau_max << endl
		<< "# max_wrongness = " << init.max_wrongness << endl
		<< "# max_x1 = " << init.max_x1 << endl;
	plotfile << "# x3 x1 x2 x0	wrong	u3 u1 u2 u0" << endl;

	myfloat wrong = wrongness(metric, &init, particle.x,
			particle.u);
	myfloat wrong_old = wrong;

	// Iterate
	int absorb_last = 0;
	int taun = -1;
	myfloat tau = 0.0;
	while ((absol(wrong) <= init.max_wrongness)
		&& (tau <= init.tau_max + dtau))
	{
		dtau = init.dtau * pow(absol(particle.x[1] - 2 * init.m), sqrt(3.0));

		if (++taun % 256 == 0) {
			info(metric, taun, tau, dtau, wrong, particle);
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

		// Update spectrum
		absorb += spec.inc_cnts(particle.x, particle.u);

		wrong_old = wrong;
		wrong = wrongness(metric, &init, particle.x, particle.u);

		if (absorb) {
			if (absorb >= absorb_max)
				break;
			if (absorb == absorb_last)
				// error, currently
				break;
			absorb_last = absorb;
		}
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

	info(metric, taun, tau, dtau, wrong, particle);

	append_file(taun);
}


// ********* MAIN ************
int main()
{
	struct initializations init = initials();

	truncate_file();

	// spacetime
	Metric metric(init.m, init.a, init.q);
	cout << endl << metric.name << ":" << endl
		<< "m = " << metric.m << endl
		<< "a = " << metric.a << endl
		<< "q = " << metric.q << endl
		<< endl;

	// emfield
	EMField emfield(metric);

	// particle
	init.u[init.change] = find_u(metric, &init, init.x, init.u);
	Particle initparticle(0, init.teilchen_masse, init.teilchen_ladung,
		init.x, init.u, metric);

	// spectrum
	myfloat emin = 2.0;
	myfloat emax = 13.0;
	myfloat start_u0 = initparticle.u[0];
	vector<myfloat> freqmult(10000);
	for (unsigned i = 0; i < freqmult.size(); i++){
		// E = u0 * freqmult
		myfloat energy = i * (emax - emin) / freqmult.size() + emin;
		freqmult[i] = energy / start_u0;
	}
	Spectrum spec(freqmult, &metric);

	ofstream specfile("globalspec.dat");
	specfile << scientific;


	// per-ray data
	int nrays;
	for (nrays = init.nrays - 1; nrays >= 0; nrays--) {
		cout << "===========================================" << endl
			<< "Ray index: " << nrays << endl
			<< endl;
		// particle
		Particle particle(0, init.teilchen_masse, init.teilchen_ladung,
				init.x, init.u, metric);
		init.change = 1;
		init.umin = init.umin_inc;
		init.umax = init.umax_inc;
		particle.u[0] = start_u0;
		particle.u[3] = init.u[3] + nrays * init.u3_inc;
		particle.u[init.change] = find_u(metric, &init, particle.x,
				particle.u);

		cout << endl;
		for (unsigned mu = 0; mu < metric.dim; mu++) {
			cout << "x" << mu << " = " << particle.x[mu] << '\t';
			cout << "u" << mu << " = " << particle.u[mu] << endl;
		}
		cout << endl;

		stringstream nrays_str;
		nrays_str << nrays;
		string plotfilename;
		plotfilename = "globalplot.dat/plot";
		plotfilename += nrays_str.str();
		plotfilename += ".dat";

		// Iterate
		absorb = 0;
		go_ray(init, spec, metric, particle, emfield, plotfilename);

		cout << endl;
		for (unsigned mu = 0; mu < metric.dim; mu++) {
			cout << "x" << mu << " = " << particle.x[mu] << '\t';
			cout << "u" << mu << " = " << particle.u[mu] << endl;
		}
		cout << endl;

		if (absorb && (absorb != absorb_max)){
			cerr << "Error detected" << endl;
			exit(1);
		}

		cout << endl;
	}


	for (unsigned i = 0; i < spec.cnts.size(); i++){
		specfile << freqmult[i] * start_u0 << '\t'
			<< spec.cnts[i] << endl;
	}
	specfile.close();


	return 0;
}
