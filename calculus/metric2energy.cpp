/*
 * HSGG Tü, 2007.06.25: Started from kerr5/HSGG Tü20070228
 */

#include <iostream>
#include <fstream>
#include <sstream>
//#include <math.h>
#include <ginac/ginac.h>

#include "mytime.h"
#include "coordinate_system.h"
#include "functions.h"



using namespace std;
using namespace GiNaC;

// metric()
struct basic_metric metric();

#define FOR(mu) for (unsigned (mu) = 0; (mu) < c.dim; (mu)++)
#define FOR2(m,n) FOR(m) FOR(n)
#define FOR3(m,n,s) FOR(m) FOR(n) FOR(s)


// line-end-finish
#define LEF(mu, nu) out << "."; \
	if ( ((mu) != c.dim - 1) || ((nu) != c.dim - 1) ) \
	{ out << " \\\\*"; } \
	out << endl;
#define LEF1(mu) out << "."; \
	if ((mu) != c.dim - 1) \
	{ out << " \\\\*"; } \
	out << endl;


#define ENVBEGIN out << "\\begin{align*}" << endl;
#define EQUALSIGN " & = "
#define ENVEND out << "\\end{align*}" << endl << endl;




// ********* MAIN - algebra ************
int main(int argc, char* argv[])
{
	const clock_t ta = clock();

	if (argc != 3)
	{
		cerr << "Usage: " << argv[0]
			<< " <infile.metric> <outfile.tex>" << endl;
		exit(2);
	}

	Coordinate_System c(argv[1]);
	string outfilename(argv[2]);


	ofstream out(outfilename.c_str());

	// Set output format
	out << latex << endl;

	out << "%/* HSGG: Created by " << argv[0] << ". */" << endl << endl;


	// LaTeX-Kopf
	out << "\\documentclass[fleqn]{article}" << endl
		<< "\\usepackage[utf8x]{inputenc}" << endl
		<< "\\usepackage{ucs}" << endl
		<< "\\usepackage{amsmath}" << endl
		<< "\\usepackage{amsfonts}" << endl
		<< "\\usepackage{amssymb}" << endl
		<< endl
		<< "\\usepackage[a4paper, left=1.5cm, right=1.5cm, top=1.5cm, bottom=4cm]{geometry}" << endl
		<< endl
		<< "\\allowdisplaybreaks[1]" << endl
		<< endl
		<< "\\begin{document}" << endl
		<< endl
		<< endl;



	// Name
	out << "{ \\Large" << endl
		<< c.name << ":"<< endl
		<< "}" << endl
		<< "\\bigskip" << endl
		<< "\\bigskip" << endl
		<< endl;

	// x^m
	out << "\\framebox{$x^\\mu$}" << endl;
	ENVBEGIN
	FOR(m)
	{
		out << "x^{" << m << "}"
			<< EQUALSIGN
			<< c.coord[m];
		LEF1(m);
	}
	ENVEND
	out << endl;

	// g_mu_nu functions
	out << "\\framebox{$g_{\\mu \\nu}$}" << endl;
	ENVBEGIN
	FOR2(mu,nu)
	{
		out << "g_{" << mu << nu << "}"
		    << EQUALSIGN
		    << c.g(mu,nu);
		LEF(mu, nu)

		cerr << "g" << mu << nu << " = " << c.g(mu,nu) << endl;
	}
	ENVEND
	cerr << "Defined metric." << endl;


	// surd
	out << "\\framebox{$\\surd = \\sqrt{-\\det(g_{\\mu \\nu})}$}" << endl;
	ENVBEGIN
		out << "\\surd" << EQUALSIGN
			<< sqrt(-c.det_g()).expand().normal()
			<< ".";
	ENVEND



	// g^mu^nu functions (G^mu^nu)
	out << "\\framebox{$g^{\\mu \\nu}$}" << endl;
	out << "\\begin{align*}" << endl;
	FOR2(mu,nu)
	{
		out << "g^{" << mu << nu << "}"
		    << EQUALSIGN
		    << c.G(mu,nu);
		LEF(mu, nu)
	}
	out << "\\end{align*}" << endl;
	out << endl;
	cerr << "Inverted metric." << endl;

	// christoffel functions Gamma^mu_nu_sigma
	cerr << "Christoffels:";
	out << "\\framebox{$\\Gamma^{\\sigma}_{\\mu \\nu}$}" << endl;
	FOR(mu)
	{
		out << "\\begin{align*}" << endl;
		FOR2(nu,sigma)
		{
			cerr << mu << nu << sigma << ":";
			out << "\\Gamma^"
			    << mu << "_{" << nu << sigma << "}"
			    << EQUALSIGN
			    << c.christoffel(mu,nu,sigma);
			LEF(nu, sigma)
		}
		out << "\\end{align*}" << endl;
	}
	out << endl;
	cerr << "done" << endl;

	/*
	// chritoffel_komma
	out << "\\framebox{$\\Gamma^\\alpha_{\\mu \\alpha , \\nu}$}" << endl;
	out << "\\begin{align*}" << endl;
	FOR(mu)
	{
		FOR(nu)
		{
			out << "\\Gamma^\\alpha_{"
			    << mu << "\\alpha , " << nu << "}"
			    << EQUALSIGN
			    << c.christoffel_komma[mu][nu];
			LEF(mu, nu)
		}
	}
	out << "\\end{align*}" << endl;
	out << endl;

	// christoffel sum komma
	out << "\\framebox{$\\Gamma^\\alpha_{\\mu \\nu , \\alpha}$}" << endl;
	ENVBEGIN
	FOR(mu)
	{
		FOR(nu)
		{
			out << "\\Gamma^\\alpha_{"
			    << mu << nu << " , \\alpha}"
			    << EQUALSIGN
			    << c.christoffel_sum_komma[mu][nu];
			LEF(mu, nu)
		}
	}
	ENVEND
	*/



	// Riemann R_mu_nu
	cerr << "Ricci:";
	out << "\\framebox{$R_{\\mu \\nu}$}" << endl;
	ENVBEGIN
	FOR2(mu,nu)
	{
		cerr << mu << nu << ":";
		out << "R_{" << mu << nu << "}"
		    << EQUALSIGN
		    << c.riemann(mu,nu);
		LEF(mu, nu)
	}
	ENVEND
	cerr << "done" << endl;

	// Riemann R^mu_nu
	cerr << "Ricci:";
	out << "\\framebox{$R^\\mu_{\\;\\nu}$}" << endl;
	ENVBEGIN
	FOR(mu)
	{
		FOR(nu)
		{
			cerr << mu << nu << ":";
			out << "R^" << mu << "_{\\;" << nu << "}"
			    << EQUALSIGN
			    << c.ricci(mu,nu);
			LEF(mu, nu)
		}
	}
	ENVEND
	cerr << "done" << endl;


	// R
	out << "\\framebox{$R$}" << endl;
	ENVBEGIN
	out << "R" << EQUALSIGN << c.gauss();
	out << "." << endl;
	ENVEND

	// Einstein G^mu_nu
	cerr << "Einstein:";
	out << "\\framebox{$G^\\mu_{\\;\\nu}$}" << endl;
	ENVBEGIN
	FOR(mu)
	{
		FOR(nu)
		{
			cerr << mu << nu << ":";
			out << "G^" << mu << "_{\\;" << nu << "}"
			    << EQUALSIGN
			    << c.einstein(mu,nu);
			LEF(mu, nu)
		}
	}
	ENVEND
	cerr << "done" << endl;

	// G
	out << "\\framebox{$G$}" << endl;
	ENVBEGIN
	out << "G" << EQUALSIGN << c.einstein_scalar();
	out << "." << endl;
	ENVEND
	cerr << "Calculated Einstein scalar." << endl;

	// G^mu_nu:mu
	out << "\\framebox{$G^\\mu_{\\; \\nu : \\mu} = 0$}" << endl;
	ENVBEGIN
	FOR(nu)
	{
		out << "G^\\mu_{\\;" << nu << ":\\mu}"
		    << EQUALSIGN
		    << c.div_einstein(nu);
		LEF1(nu)
	}
	ENVEND
	cerr << "Calculated Einstein divergence." << endl;



	// Harmony
	cerr << "Harmony:";
	out << "\\framebox{$g^{\\mu \\nu}\\,\\Gamma^\\lambda_{\\mu\\nu} = 0$?}"
	    << endl;
	ENVBEGIN
	FOR(m)
	{
		cerr << m << ":";
		out << "g^{\\mu \\nu}\\,\\Gamma^" << m << "_{\\mu \\nu}"
		    << EQUALSIGN
		    << c.harmony(m);
		LEF1(m)
	}
	ENVEND
	cerr << "done" << endl;


	// LaTeX-Fuß
	out << endl
		<< "\\end{document}" << endl
		<< endl;

	out.close();
	cout << "Created \"" << outfilename << "\"." << endl;

	cerr << "====================="
	     << "   Time taken: " << timediff(ta, clock()) << " s   "
	     << "=====================" << endl;

	return 0;
}

