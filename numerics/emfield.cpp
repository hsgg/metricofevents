/*
 * HSGG Tü, 2007.05.05
 */

#include <iostream>
#include <string>
#include <math.h>
#include "myfloat.h"
#include "metric.h"
#include "emfield.h"


using namespace std;


inline const myfloat EMField::rho4(const vector<myfloat>& x) const
{
	return pow((x[1] * x[1] + pow(metric->a * cos(x[2]), 2)), 2);
}


// F_00
const myfloat EMField::F_00(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_01
const myfloat EMField::F_01(const vector<myfloat>& x) const
{
	return -F_10(x);
}

// F_02
const myfloat EMField::F_02(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_03
const myfloat EMField::F_03(const vector<myfloat>& x) const
{
	return 0.0;
}


// F_10
const myfloat EMField::F_10(const vector<myfloat>& x) const
{
	return 2.0 * (metric->q / rho4(x))
		* (x[1] * x[1] - pow(metric->a * cos(x[2]), 2));
}

// F_11
const myfloat EMField::F_11(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_12
const myfloat EMField::F_12(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_13
const myfloat EMField::F_13(const vector<myfloat>& x) const
{
	return 2.0 * (metric->q / rho4(x))
		* (x[1] * x[1] - pow(metric->a * cos(x[2]), 2))
		* (-metric->a * pow(sin(x[2]), 2));
}


// F_20
const myfloat EMField::F_20(const vector<myfloat>& x) const
{
	return 2.0 * (2.0 * metric->q / rho4(x))
		* (metric->a * x[1] * cos(x[2]) * sin(x[2]))
		* (-metric->a);
}

// F_21
const myfloat EMField::F_21(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_22
const myfloat EMField::F_22(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_23
const myfloat EMField::F_23(const vector<myfloat>& x) const
{
	return 2.0 * (2.0 * metric->q / rho4(x))
		* (metric->a * x[1] * cos(x[2]) * sin(x[2]))
		* (x[1] * x[1] + metric->a * metric->a);
}


// F_30
const myfloat EMField::F_30(const vector<myfloat>& x) const
{
	return 0.0;
}

// F_31
const myfloat EMField::F_31(const vector<myfloat>& x) const
{
	return -F_13(x);
}

// F_32
const myfloat EMField::F_32(const vector<myfloat>& x) const
{
	return -F_23(x);
}

// F_33
const myfloat EMField::F_33(const vector<myfloat>& x) const
{
	return 0.0;
}




// Konstruktor
EMField::EMField(const Metric& gmunu)
    : metric(&gmunu)
{
	if (metric->name != std::string("Kerr-Newman"))
	{
		std::cerr << "Possibly wrong metric!!" << std::endl;
	}

	F[0][0] = &EMField::F_00;
	F[0][1] = &EMField::F_01;
	F[0][2] = &EMField::F_02;
	F[0][3] = &EMField::F_03;
	F[1][0] = &EMField::F_10;
	F[1][1] = &EMField::F_11;
	F[1][2] = &EMField::F_12;
	F[1][3] = &EMField::F_13;
	F[2][0] = &EMField::F_20;
	F[2][1] = &EMField::F_21;
	F[2][2] = &EMField::F_22;
	F[2][3] = &EMField::F_23;
	F[3][0] = &EMField::F_30;
	F[3][1] = &EMField::F_31;
	F[3][2] = &EMField::F_32;
	F[3][3] = &EMField::F_33;
}


// Destruktor
EMField::~EMField(){}
