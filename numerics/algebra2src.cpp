/*
 * HSGG Tü20070228
 * HSGG Tü20070707
 * HSGG Tü20070810
 */

#include <iostream>
#include <fstream>
#include <sstream>
//#include <math.h>
#include <ginac/ginac.h>

#include "../calculus/coordinate_system.h"


using namespace std;
using namespace GiNaC;


// ********* MAIN - algebra ************
int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "Usage: " << argv[0] << " <infile.metric>" << endl;
		cerr << "Output is written to \"metric.cpp\" and \"metric.h\"."
			<< endl;
		exit(2);
	}

	// Calculate equations
	Coordinate_System c(argv[1]);


	exmap co;
	vector<symbol> u(c.dim);
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		stringstream s;
		s << "x[" << mu << "]";
		co[c.coord[mu]] = symbol(s.str().c_str());
	}
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		stringstream s;
		s << "u[" << mu << "]";
		u[mu] = symbol(s.str().c_str());
	}


	// Write metric.cpp and metric.h
	ofstream out("metric.cpp");
	ofstream outh("metric.h");

	// Set output format
	out << csrc_double << endl;
	outh << csrc_double << endl;

	out << "/* HSGG: Created by " << argv[0] << ". */" << endl << endl;
	outh << "/* HSGG: Created by " << argv[0] << ". */" << endl << endl;

	outh << "#ifndef METRIC_H" << endl
	     << "#define METRIC_H" << endl << endl;

	out << "#include <math.h>" << endl;
	out << "#include <iostream>" << endl;
	out << "#include \"metric.h\"" << endl << endl << endl;

	outh << "#include <string>" << endl;
	outh << "#include <vector>" << endl;
	outh << "#include \"global.h\"" << endl;
	outh << "#include \"myfloat.h\"" << endl;
	outh << endl;


	// Define class Metric
	outh << "class Metric" << endl
	     << "{" << endl
	     << "    public:" << endl
	     << "\t// Konstruktor" << endl
	     << "\tMetric(const myfloat& mass"
		<< ", const myfloat& ang"
		<< ", const myfloat& charge);" << endl
	     << "\t// Destruktor" << endl
	     << "\t~Metric();" << endl
	     << endl
	     << "\tconst std::string name;" << endl
	     << "\tconst unsigned dim;" << endl
	     << "\tconst myfloat m;" << endl
	     << "\tconst myfloat a;" << endl
	     << "\tconst myfloat q;" << endl
	     << "\tconst myfloat (Metric::*g[" << c.dim << "][" << c.dim << "])"
		<< "(const myfloat x[" << c.dim << "]) const;"
		<< endl
	     << "\tconst myfloat (Metric::*G[" << c.dim << "][" << c.dim << "])"
		<< "(const myfloat* x) const;"
		<< endl
	     << "\tconst myfloat (Metric::*christoffelsum[" << c.dim << "])"
		<< "(const myfloat* x, const myfloat* u) const;"
		<< endl
	     << "\tconst myfloat (Metric::*acceleration3d[" << c.dim << "])"
		<< "(const myfloat* x, const myfloat* u) const;"
		<< endl
	     << endl
	     << "    private:" << endl;
	


	// g_mu_nu functions
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		for (unsigned nu = 0; nu < c.dim; nu++)
		{
			outh << "\tconst myfloat g" << mu << nu
			     << "(const myfloat x[" << c.dim << "]) const"
			     << ";" << endl;
			out << "const myfloat Metric::g" << mu << nu
			    << "(const myfloat x[" << c.dim << "]) const" << endl
			    << "{" << endl
			    << "\treturn " << c.g(mu,nu).subs(co) << ";"
				<< endl
			    << "}" << endl << endl;
		}
	}
	outh << endl;

	// g^mu^nu functions (G^mu^nu)
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		for (unsigned nu = 0; nu < c.dim; nu++)
		{
			outh << "\tconst myfloat G" << mu << nu
			     << "(const myfloat* x) const"
			     << ";" << endl;
			out << "const myfloat Metric::G" << mu << nu
			    << "(const myfloat* x) const" << endl
			    << "{" << endl
			    << "\treturn " << c.G(mu,nu).subs(co) << ";"
				<< endl
			    << "}" << endl << endl;
		}
	}
	outh << endl;

	/*
	// christoffel functions
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		for (unsigned nu = 0; nu < c.dim; nu++)
		{
			for (unsigned sigma = 0; sigma < c.dim; sigma++)
			{
				outh << "\tconst myfloat christoffel"
				     << mu << "_" << nu << sigma
				     << "(const myfloat x[" << c.dim << "]) const"
				     << ";" << endl;
				out << "const myfloat Metric::christoffel"
				    << mu << "_" << nu << sigma
				    << "(const myfloat x[" << c.dim << "]) const" << endl
				    << "{" << endl
				    << "\treturn "
					<< c.christoffel[mu][nu][sigma].subs(co)
					<< ";" << endl
				    << "}" << endl << endl;
			}
		}
	}
	*/


	// equation of motion (sigma, x[c.dim], u[c.dim])
	cerr << "Calculating christoffelsum... ";
	vector<ex> csum(c.dim);
	for (unsigned sigma = 0; sigma < c.dim; sigma++)
	{
		csum[sigma] = 0;
		for (unsigned mu = 0; mu < c.dim; mu++)
		{
			for (unsigned nu = 0; nu < c.dim; nu++)
			{
				// with respect to tau:
				csum[sigma] = csum[sigma]
					+ c.christoffel(sigma,mu,nu)
					* u[mu] * u[nu];
			}
		}
		cerr << sigma;
		csum[sigma] = csum[sigma].normal();
		cerr << ":";
		outh << "\tconst myfloat chsum"
		     << sigma
		     << "(const myfloat* x, const myfloat* u) const"
		     << ";" << endl;
		out << "const myfloat Metric::chsum"
		    << sigma
		    << "(const myfloat* x, const myfloat* u) const"
		    << endl
		    << "{" << endl
		    << "\treturn "
			<< csum[sigma].subs(co)
			<< ";" << endl
		    << "}" << endl << endl;
	}
	cerr << "done" << endl;


	// 3-acceleration
	exmap m;
	m[u[0]] = 1;
	cerr << "Calculating 3-acceleration... ";
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		ex ddotx = 0;
		for (unsigned sigma = 0; sigma < c.dim; sigma++)
		{
			for (unsigned rho = 0; rho < c.dim; rho++)
			{
				ddotx += (c.christoffel(0, sigma, rho) * u[mu]
						- c.christoffel(mu, sigma, rho))
					* u[sigma] * u[rho];
			}
		}
		cerr << mu;
		ddotx = ddotx.subs(m).normal();
		cerr << ":";
		outh << "\tconst myfloat acc3d"
		     << mu
		     << "(const myfloat* x, const myfloat* u) const"
		     << ";" << endl;
		out << "const myfloat Metric::acc3d"
		    << mu
		    << "(const myfloat* x, const myfloat* u) const"
		    << endl
		    << "{" << endl
		    << "\treturn "
			<< ddotx.subs(co)
			<< ";" << endl
		    << "}" << endl << endl;
	}
	cerr << "done" << endl;


	/* Define constructor */
	out << "Metric::Metric(const myfloat& mass"
	    << ", const myfloat& ang"
	    << ", const myfloat& charge)" << endl
	    << "    : ";
	// Konstants
	out << "name(\"" << c.name << "\")" << endl
	    << "\t, dim(" << c.dim << ")" << endl
	    << "\t, m(mass)" << endl
	    << "\t, a(ang)" << endl
	    << "\t, q(charge)" << endl
	    ;
	out << "{" << endl;

	// g_mu_nu
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
	    	for (unsigned nu = 0; nu < c.dim; nu++)
		{
			out << "\tg[" << mu << "][" << nu << "]"
			    << " = &Metric::g" << mu << nu << ";" << endl;
		}
	}
	out << endl;
	// g^mu^nu (G^mu^nu)
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
	    	for (unsigned nu = 0; nu < c.dim; nu++)
		{
			out << "\tG[" << mu << "][" << nu << "]"
			    << " = &Metric::G" << mu << nu << ";" << endl;
		}
	}
	out << endl;
	/* christoffel^mu_nu_sigma
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		for (unsigned nu = 0; nu < c.dim; nu++)
		{
			for (unsigned sigma = 0; sigma < c.dim; sigma++)
			{
				out << "\tchristoffel["
				    << mu << "][" << nu << "][" << sigma << "]"
				    << " = &Metric::christoffel"
				    << mu << "_" << nu << sigma << ";" << endl;
			}
		}
	}
	out << endl; */
	// christoffelsum, acceleration3d
	for (unsigned mu = 0; mu < c.dim; mu++)
	{
		out << "\tchristoffelsum["
		    << mu << "]"
		    << " = &Metric::chsum"
		    << mu << ";" << endl;
		out << "\tacceleration3d["
		    << mu << "]"
		    << " = &Metric::acc3d"
		    << mu << ";" << endl;
	}
	out << endl;
	out << "}" << endl << endl; // end of constructor


	/* Define destructor */
	out << "Metric::~Metric(){}" << endl << endl;


	outh << "};" << endl << endl
	     << "#endif" << endl;


	out.close();
	outh.close();

	cout << "Created metric.h und metric.cpp." << endl;



	return 0;
}

