/*
 * HSGG Tü20070301
 */

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "myfloat.h"
#include "init.h"

using namespace std;


struct initializations initials()
{
	struct initializations i;

	i.radius = 10.0;

	//myfloat c = 299792458;

	//i.m = (6.6726e-11) * (1.99e+30) / (c * c);
	i.m = 1.0;
	i.a = 0.9;
	i.q = 0.0; //0.15;


	i.teilchen_masse = 1.0;
	i.teilchen_ladung = 0.0;

	i.x[0] = 0.0;
	//i.x[1] = 1.496e+11;
	i.x[1] = 9.0 * i.m;
	i.x[2] = (myfloat) (M_PI_2);
	i.x[3] = 0.0;

	//myfloat E = 1.0070;
	i.change = 0; // welches u ueberschrieben wird
	i.umin = 0;
	i.umax = 10;
	i.u[0] = 1.0; //sqrt( E / (1.0 - 2.0 * i.m / i.x[1] - i.m / i.x[1]) );
	i.u[1] = -1.1; //-0.7;
	i.u[2] = 0.0;
	//i.u[3] = 0.13; //0.21; // -sqrt(i.m / (i.x[1] * i.x[1] * i.x[1])) * i.u[0] * i.a;
	i.u[3] = 0.06;
	//i.u[3] = L / (i.x[1] * i.x[1])
	//	* (myfloat) ( pow(sin(i.x[2]),-0.999999999992) );

	i.nrays = 100;
	i.u3_inc = -(2.0 * i.u[3]) / (i.nrays - 1);
	i.umin_inc = -10.0;
	i.umax_inc = 0.0;

	//i.dtau = c * 86400 * 160;
	//i.tau_max = 1000 * 31.536e+6 * c;
	i.dtau = 0.001;
	i.tau_max = 65;
	i.max_n = 150000;
	i.max_wrongness = 1e-4;
	i.min_x1 = 0.0;
	i.max_x1 = 50.1;

	return i;
}

// Eliminate beginning Whitespace
static void eliminate_beginning_whitespace(string& s)
{
	unsigned i = 0;
	while ( (s[i] == ' ') || (s[i] == '	') )
		i++;

	s.erase(0, i);
}

// Eliminate ending whitespace
static void eliminate_ending_whitespace(string& s)
{
	int i = s.length() - 1;
	while ( (s[i] == ' ') || (s[i] == '	') )
		i--;

	s.erase(i + 1);
}


// Eliminate all comments (everything on a line starting from '#').
static void eliminate_comments(string& s)
{
	unsigned i = 0;
	while ( (s[i] != '#') && (i < s.length()) )
		i++;

	s.erase(i);
}


// Returns the number of occurrences of the character 'ch' in 's'.
static unsigned contains(const char& ch, const string& s)
{
	unsigned num = 0;
	unsigned i = 0; 
	while (i < s.length())
	{
		if (s[i++] == ch)
			num++;
	}

	return num;
}


// split_at_equalsign
static void split_at_equalsign(const string& s,
		string& before, string& after)
{
	unsigned i = s.find_first_of('=');

	before.assign(s, 0, i);
	after.assign(s, i + 1, s.length() - i - 1);

	eliminate_beginning_whitespace(after);
	eliminate_ending_whitespace(before);
}



struct initializations initialize(char *filename)
{
	struct initializations i;
	bool version_ok = false;

	ifstream fin(filename);
	string s;
	string name ("");

	unsigned lineno = 0;

	while (getline(fin, s))
	{
		lineno++;

		eliminate_beginning_whitespace(s);
		eliminate_comments(s);
		eliminate_ending_whitespace(s); // Must come after eliminate_comments
		if (s.empty()) continue;

		cerr << lineno << ": " << s << "." << endl;

		// Name? (Starts with '"'.)
		if (s[0] == '"') {
			name += s.substr(1, s.length() - 2);
			continue;
		}

		if (contains('=', s) != 1) {
			cerr << lineno << ": Syntax error: Must contain exactly one '='."
				<< endl;
			exit(1);
		}


		string before;
		string after;
		split_at_equalsign(s, before, after);

		if (before == string("version")) {
			if (after != string("0.1"))
				cerr << "Warning: We only support version 0.1! Good luck!" << endl;
			else
				version_ok = true;

		} else if (before == string("radius")) {
			i.radius = atof(after.c_str());

		} else if (before == string("m")) {
			i.m = atof(after.c_str());
		} else if (before == string("a")) {
			i.a = atof(after.c_str());
		} else if (before == string("q")) {
			i.q = atof(after.c_str());

		} else if (before == string("change")) {
			i.change = atof(after.c_str());
		} else if (before == string("umin")) {
			i.umin = atof(after.c_str());
		} else if (before == string("umax")) {
			i.umax = atof(after.c_str());

		} else if (before == string("mass[0]")) {
			i.mass[0] = atof(after.c_str());
		} else if (before == string("charge[0]")) {
			i.charge[0] = atof(after.c_str());
		} else if (before == string("x[0][0]")) {
			i.xvec[0][0] = atof(after.c_str());
		} else if (before == string("x[0][1]")) {
			i.xvec[0][1] = atof(after.c_str());
		} else if (before == string("x[0][2]")) {
			i.xvec[0][2] = atof(after.c_str());
		} else if (before == string("x[0][3]")) {
			i.xvec[0][3] = atof(after.c_str());
		} else if (before == string("u[0][0]")) {
			i.uvec[0][0] = atof(after.c_str());
		} else if (before == string("u[0][1]")) {
			i.uvec[0][1] = atof(after.c_str());
		} else if (before == string("u[0][2]")) {
			i.uvec[0][2] = atof(after.c_str());
		} else if (before == string("u[0][3]")) {
			i.uvec[0][3] = atof(after.c_str());

		} else if (before == string("mass[1]")) {
			i.mass[1] = atof(after.c_str());
		} else if (before == string("charge[1]")) {
			i.charge[1] = atof(after.c_str());
		} else if (before == string("x[1][0]")) {
			i.xvec[1][0] = atof(after.c_str());
		} else if (before == string("x[1][1]")) {
			i.xvec[1][1] = atof(after.c_str());
		} else if (before == string("x[1][2]")) {
			i.xvec[1][2] = atof(after.c_str());
		} else if (before == string("x[1][3]")) {
			i.xvec[1][3] = atof(after.c_str());
		} else if (before == string("u[1][0]")) {
			i.uvec[1][0] = atof(after.c_str());
		} else if (before == string("u[1][1]")) {
			i.uvec[1][1] = atof(after.c_str());
		} else if (before == string("u[1][2]")) {
			i.uvec[1][2] = atof(after.c_str());
		} else if (before == string("u[1][3]")) {
			i.uvec[1][3] = atof(after.c_str());

		} else if (before == string("nrays")) {
			i.nrays = atof(after.c_str());
		} else if (before == string("u3_inc")) {
			i.u3_inc = atof(after.c_str());
		} else if (before == string("umin_inc")) {
			i.umin_inc = atof(after.c_str());
		} else if (before == string("umax_inc")) {
			i.umax_inc = atof(after.c_str());

		} else if (before == string("dtau")) {
			i.dtau = atof(after.c_str());
		} else if (before == string("tau_max")) {
			i.tau_max = atof(after.c_str());
		} else if (before == string("max_n")) {
			i.max_n = atof(after.c_str());
		} else if (before == string("max_wrongness")) {
			i.max_wrongness = atof(after.c_str());
		} else if (before == string("min_x1")) {
			i.min_x1 = atof(after.c_str());
		} else if (before == string("max_x1")) {
			i.max_x1 = atof(after.c_str());
		}
	}

	if (!version_ok)
		cerr << "Warning: We only support version 0.1! Good luck!" << endl;

	return i;
}
