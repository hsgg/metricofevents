/*
 * HSGG Tü20070629
 * HSGG Tü20070706
 */

#ifndef COORDINATE_SYSTEM_H
#define COORDINATE_SYSTEM_H

#include <ginac/ginac.h>
#include <string>
#include <vector>



// Information about space-time
class Coordinate_System
{

	// basic_equation
	typedef struct {
		GiNaC::ex e;
		bool calculated;
	} eq0;
	typedef std::vector<eq0> eq1;
	typedef std::vector<eq1> eq2;
	typedef std::vector<eq2> eq3;

public:
	// Constructor
	Coordinate_System(const char * filename);

	// Destructor
	~Coordinate_System();


	std::string name; // name
	unsigned dim;

	std::vector<GiNaC::symbol> coord; // x^mu

	GiNaC::ex g(const int m, const int n);
	GiNaC::ex G(const int m, const int n);

	GiNaC::ex det_g();

	GiNaC::ex christoffel(const int s, const int m, const int n);

	GiNaC::ex riemann(const int m, const int n);
	GiNaC::ex ricci(const int m, const int n);
	GiNaC::ex gauss();

	GiNaC::ex einstein(const int m, const int n);
	GiNaC::ex einstein_scalar();
	GiNaC::ex div_einstein(const int m);

	GiNaC::ex harmony(const int m);

private:
	eq2 g_mn; // g_mu_nu
	eq2 G_mn; // g^mu^nu

	eq3 g_mn_s; // g_mn,s
	GiNaC::ex g_diff(const int m, const int n, const int s);

	eq3 christoffel_s_mn; // christoffel^sigma_mu_nu

	eq1 chris_a_na; // chris^a_na
	GiNaC::ex christoffel_a_ma(const int m);

	eq2 chris_a_na_m; // chris^a_na,m
	GiNaC::ex christoffel_a_ma_n(const int m, const int n);
	eq2 chris_a_mn_a; // chris^a_mu_nu,a
	GiNaC::ex christoffel_a_mn_a(const int m, const int n);
	eq2 chris_a_mn_chris_b_ab; // chris^a_mn chris^b_ab
	GiNaC::ex christoffel_a_mn_christoffel_b_ab(const int m, const int n);
	eq2 chris_a_mb_chris_b_na; // chris^a_mb chris^b_na
	GiNaC::ex christoffel_a_mb_christoffel_b_na(const int m, const int n);

	eq2 riemann_mn; // R_mu_nu
	eq2 ricci_m_n; // R^mu_nu
	eq0 gauss_a_a; // R

	eq2 einstein_m_n; // G^mu_nu
	eq0 einstein_trace; // G = G^mu_mu
	eq1 div_einstein_m; // G^mu_nu:mu

	eq1 harmony_m; // g^mn * L^a_mn



	// compute inverse metric G from g
	void calc_inverse_metric();


	// Metric reading help.
	void eliminate_beginning_whitespace(std::string& s);
	void eliminate_ending_whitespace(std::string& s);
	void eliminate_comments(std::string& s);
	bool long_comment_here(const std::string& s);
	unsigned contains(const char& ch, const std::string& s);
	void split_at_equalsign(const std::string& s,
			std::string& before, std::string& after);
	void inlst(const unsigned& lineno, const std::string& sym,
			const GiNaC::lst& symbollist);
	void split_atoi(const std::string& s, unsigned& a, unsigned& b);
	void resize_g_mni(const unsigned& i, const unsigned& dim);

	// init
	void init(const unsigned& dimension);



};


#endif
