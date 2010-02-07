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

using namespace std;


// absolute value of myfloat
inline myfloat abs(const myfloat number)
{
	return (number < 0) ? -number : number;
}


// set xpk = x + step * k
inline void mk_cok(vector<myfloat>& xpk, const myfloat* x,
		const vector<myfloat>& xk_n, const myfloat step)
{
	for (unsigned mu = 0; mu < xpk.size(); mu++)
	{
		xpk[mu] = x[mu] + step * xk_n[mu];
	}
}


// contract twice
myfloat christoffelsum(const Metric& metric, const unsigned& sigma,
	const vector<myfloat>& x, const vector<myfloat>& u)
{
	/*
	myfloat result = 0.0;
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		for (unsigned nu = 0; nu < metric.dim; nu++)
		{
			result +=
				(metric.*(metric.christoffel[sigma][mu][nu]))(x)
				* u[mu] * v[nu];
		}
	}*/

	return (metric.*(metric.christoffelsum[sigma]))(x, u);
}



// electromagnetic acceleration
myfloat emfieldforce(const EMField& emfield, const unsigned& sigma,
	const vector<myfloat>& x, const vector<myfloat>& u)
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


// calculate scalarproduct of y and z at x.
myfloat scalar(const Metric& metric, const myfloat* x,
	const myfloat* y, const myfloat* z)
{
	myfloat l = 0.0;

	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		for (unsigned nu = 0; nu < metric.dim; nu++)
		{
			l += (metric.*(metric.g[mu][nu]))(x) * y[mu] * z[nu];
		}
	}

	return l;
}


// Falschheit
inline myfloat wrongness(const Metric& metric,
		const struct initializations* init,
		const myfloat* x, const myfloat* u)
{
	return scalar(metric, x, u, u) - init->teilchen_masse;
}


// Calculate u
myfloat find_u(const Metric& metric,
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
		if ( (abs(fcurr) < abs(fmax)) && (fcurr * fmax > 0.0) )
		{
			max = un;
		}
		else if ( (abs(fcurr) < abs(fmin)) && (fcurr * fmin > 0.0) )
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


// Berechne x, u
// input: metric, x, u, dtau
inline void x_and_u(const Metric& metric, const EMField& emfield,
	const myfloat dtau, Particle& particle)
{
	vector< vector<myfloat> > xk(4);
	vector< vector<myfloat> > uk(4);

	vector<myfloat> xpk(metric.dim);
	vector<myfloat> upk(metric.dim);

	// Berechne xk, uk
	for (unsigned i = 0; i < 4; i++)
	{
		xk[i].resize(metric.dim);
		uk[i].resize(metric.dim);

		// set cok = co + k, depending on how far we are
		switch (i)
		{
			case 0:
				for (unsigned j = 0; j < metric.dim; j++)
				{
					xk[0][j] = 0.0;
					uk[0][j] = 0.0;
				}
				mk_cok(xpk, particle.x, xk[0], 0.0);
				mk_cok(upk, particle.u, uk[0], 0.0);
				break;
			case 1:
				mk_cok(xpk, particle.x, xk[0], (myfloat) 0.5);
				mk_cok(upk, particle.u, uk[0], (myfloat) 0.5);
				break;
			case 2:
				mk_cok(xpk, particle.x, xk[1], (myfloat) 0.5);
				mk_cok(upk, particle.u, uk[1], (myfloat) 0.5);
				break;
			case 3:
				mk_cok(xpk, particle.x, xk[2], 1.0);
				mk_cok(upk, particle.u, uk[2], 1.0);
			break;
		}

		// Berechne xk, uk
		for (unsigned mu = 0; mu < metric.dim; mu++)
		{
			xk[i][mu] = dtau * upk[mu];
			uk[i][mu] = dtau
				* (-christoffelsum(metric, mu, xpk, upk)
				- particle.q
				* emfieldforce(emfield, mu, xpk, upk));
		}
	}

	// Berechne x, u
	for (unsigned mu = 0; mu < metric.dim; mu++)
	{
		particle.x[mu] = particle.x[mu]
			+ xk[0][mu]/6.0
			+ xk[1][mu]/3.0
			+ xk[2][mu]/3.0
			+ xk[3][mu]/6.0;
		//cout << co[c.x[mu]] << endl;
		particle.u[mu] = particle.u[mu]
			+ uk[0][mu]/6.0
			+ uk[1][mu]/3.0
			+ uk[2][mu]/3.0
			+ uk[3][mu]/6.0;
		//cout << co[c.u[mu]] << endl;
	}
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

	/*
	// Print important information to plotfile
	plotfile << "# Metric: " << c.name_of_metric << endl << endl;
	for (unsigned i = 0; i < metric.dim; i++)
	{
		plotfile << "# u" << i << "_pkt = " << u_pkt[i] << endl;
	}
	plotfile << endl;
	plotfile << "# " << c.m << " = " << co[c.m] << endl
		<< "# " << c.a << " = " << co[c.a] << endl
		<< "# teilchen_masse = " << init.teilchen_masse << endl
		<< endl;
	for (unsigned i = 0; i < metric.dim; i++)
	{
		plotfile << "# " << c.x[i] << " = " << co[c.x[i]] << endl;
	}
	plotfile << endl;
	for (unsigned i = 0; i < metric.dim; i++)
	{
		plotfile << "# " << c.u[i] << " = " << co[c.u[i]] << endl;
	}
	plotfile << endl;
	plotfile << "# length_4velocity = " << length_4velocity(c).subs(co)
		<< endl << endl
		<< "# dtau = " << dtau << endl
		<< "# tau_max = " << taun_max * dtau << endl
		<< endl;
	*/

	//cout << "taun_max = " << taun_max << endl;
	//cout << "tau_max = " << init.tau_max << endl << endl;
	//cerr << scientific;

	// Iterate
	myfloat tau_old = 1.0;
	int taun = -1;
	myfloat tau = 0.0;
	//while (tau != tau_old)
	while
	(
		//(dtau >= pow(init.max_wrongness, 2))
		//&&
		(abs(wrong) <= init.max_wrongness)
		&&
		(tau <= init.tau_max + dtau)
	)
	{
		dtau = init.dtau * exp(-abs(wrong) / init.max_wrongness);

		if (++taun % 1000 == 0)
		{
			info(taun, tau, dtau, wrong, particle.x);
		}

		// print to file
		//plotfile<< x[1] * cos(x[3]) * sin(x[2]) << "\t"
		//	<< x[1] * sin(x[3]) * sin(x[2]) << "\t"
		//	<< x[1] 	    * cos(x[2]) << "\t";
		plotfile<< particle.x[3] << ' '
			<< particle.x[1] << ' '
			<< particle.x[2] << ' '
			<< particle.x[0] << '\t' ;
		plotfile<< wrong << '\t';
		plotfile<< particle.u[3] << ' '
		        << particle.u[1] << ' '
		        << particle.u[2] << ' '
		        << particle.u[0]
			<< endl;

		dtaufile << tau << "\t" << dtau << endl;
		wrongfile << tau << "\t" << wrong << endl;

		tau_old = tau;
		tau += dtau;


		// Berechne x, u
		x_and_u(metric, emfield, dtau, particle);

		wrong = wrongness(metric, &init, particle.x, particle.u);

		if (particle.x[1] > 20)
		{
			break;
		}
	}
	plotfile << endl;
	plotfile << "# length_4velocity = "
		 << scalar(metric, particle.x, particle.u, particle.u)
		 << endl;
	plotfile.close();
	dtaufile.close();
	wrongfile.close();


	//cout << "u0 = " << u[0] << endl;
	//cout << "u3 = " << u[3] << endl;


	info(taun, tau, dtau, wrong, particle.x);

	cout << endl;

	return 0;
}
