/*
 * HSGG Tü 2007.04.29
 */

#include <sstream>
#include <cmath>

#include "misclib.h"
#include "particle.h"
#include "numerics.h"


using namespace std;


Particle::Particle(int unique, const myfloat& masse, const myfloat& ladung,
	const myfloat pos[DIM], const myfloat vel[DIM], const Metric& metric)
    : m(masse), q(ladung), unique(unique), metric(metric)
{
	for (int mu = 0; mu < DIM; mu++)
	{
		x[mu] = pos[mu];
		u[mu] = vel[mu];
	}

	stringstream s;
	s << "plot" << unique << ".dat";
	plotfile.open(s.str().c_str());
	plotfile << scientific << setprecision(12);
	/*
	plotfile << "# Metric: " << metric.name << endl;
	plotfile << "#	m = " << metric.m << endl
		<< "#	a = " << metric.a << endl
		<< "#	q = " << metric.q << endl;
		*/
	plotfile << "# particle " << unique << ":" << endl
		<< "#	m = " << m << endl
		<< "#	q = " << q << endl;
	/*
	plotfile << "# dtau = " << dtau << endl
		<< "# tau_max = " << init.tau_max << endl
		<< "# max_x1 = " << init.max_x1 << endl;
	plotfile << "# x3 x1 x2 x0	gamma u3 u1 u2 u0" << endl;
	*/
}



Particle::~Particle()
{
	plotfile << endl;
	plotfile.close();
}


void Particle::write_to_plotfile(const vector<Particle*>& particles)
{
	plotfile<< x[3] << ' '
		<< x[1] << ' '
		<< x[2] << ' '
		<< x[0] << '\t'
		<< gammafactor(metric, *this) << ' '
		<< u[3] << ' '
		<< u[1] << ' '
		<< u[2] << ' '
		<< u[0] << '\t';
	for (unsigned i = 0; i < particles.size(); i++) {
		plotfile<< sqrtl(spatial_projection_scalar(metric, *this, x, particles[i]->x))
			<< ' ';
	}
	plotfile<< endl;
}
