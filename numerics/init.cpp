/*
 * HSGG Tü20070301
 */

#include <math.h>

#include "myfloat.h"
#include "init.h"

struct initializations initials()
{
	struct initializations i;

	i.radius = 10.0;

	//myfloat c = 299792458;

	//i.m = (6.6726e-11) * (1.99e+30) / (c * c);
	i.m = 1.0;
	i.a = 0.1;
	i.q = 0.15;


	i.teilchen_masse = 1.0;
	i.teilchen_ladung = 0.0;

	i.x[0] = 0.0;
	//i.x[1] = 1.496e+11;
	i.x[1] = 3.0 * i.m;
	i.x[2] = (myfloat) (M_PI / 2.0);
	i.x[3] = 0.0;

	//myfloat E = 1.0070;
	i.change = 0; // welches u ueberschrieben wird
	i.umin = 0;
	i.umax = 10;
	i.u[0] = 0.0; //sqrt( E / (1.0 - 2.0 * i.m / i.x[1] - i.m / i.x[1]) );
	i.u[1] = -1.1; //-0.7;
	i.u[2] = 0.0;
	i.u[3] = 0.0; // -sqrt(i.m / (i.x[1] * i.x[1] * i.x[1])) * i.u[0] * i.a;
	//i.u[3] = L / (i.x[1] * i.x[1])
	//	* (myfloat) ( pow(sin(i.x[2]),-0.999999999992) );

	//i.dtau = c * 86400 * 160;
	//i.tau_max = 1000 * 31.536e+6 * c;
	i.dtau = 0.001;
	i.tau_max = 65;
	i.max_wrongness = 1e-3;

	return i;
}
