/*
 * HSGG Tü, 2007.06.25: Started from kerr5/HSGG Tü20070228
 * HSGG Tü, 2007.07.05: Started from einstein_tensor algebra2src.cpp
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>

#include <ginac/ginac.h>

#include "coordinate_system.h"


using namespace std;
using namespace GiNaC;



#define FOR(mu) for (unsigned (mu) = 0; (mu) < dim; (mu)++)

#define FOR2(m,n) FOR(m) FOR(n)

#define FOR3(m,n,s) FOR(m) FOR(n) FOR(s)


/* For the metric */
#define NORMG .normal()

/* For riemann and up */
#define NORMR .normal().expand()

/* For everything else */
#define NORM .normal()





// const int argument macro
#define ARG0 
#define ARG1(mu) const int mu
#define ARG2(mu,nu) const int mu, const int nu
#define ARG3(mu,nu,sigma) const int mu, const int nu, const int sigma

// calc_function macro
#define EXFUNC_HEAD(name, arguments) \
	GiNaC::ex Coordinate_System::name(arguments)



// metric
EXFUNC_HEAD(g, ARG2(m,n))
{
	if (!g_mn[m][n].calculated)
	{
		g_mn[m][n].e = (g_mn[m][n].e)NORMG;

		// symmetric?
		if ( (g_mn[n][m].calculated)
			&& (g_mn[m][n].e != g_mn[n][m].e) )
		{
			cerr << "Metric not symmetric:" << endl
				<< "\tg" << m << n
				<< " = " << g_mn[m][n].e
				<< endl
				<< "\tg" << n << m
				<< " = " << g_mn[n][m].e
				<< endl;
			exit (1);
		}
		else
		{
			g_mn[m][n].calculated = true;
		}
	}

	return g_mn[m][n].e;
}



// Inverse metric
EXFUNC_HEAD(G, ARG2(m,n))
{
	static bool calculated = false;

	if (!calculated)
	{
		calc_inverse_metric(); // TODO: test for symmetry
		calculated = true;
	}

	return G_mn[m][n].e;
}


// g_mn,s
EXFUNC_HEAD(g_diff, ARG3(m,n,s))
{
	if (!g_mn_s[m][n][s].calculated)
	{
		if (g_mn_s[n][m][s].calculated)
		{
			g_mn_s[m][n][s].e = g_mn_s[n][m][s].e;
		}
		else
		{
			g_mn_s[m][n][s].e = g(m,n).diff(coord[s])NORM;
		}

		g_mn_s[m][n][s].calculated = true;
	}

	return g_mn_s[m][n][s].e;
}


// christoffel
EXFUNC_HEAD(christoffel, ARG3(s,m,n))
{
	if (!christoffel_s_mn[s][m][n].calculated)
	{
		if (christoffel_s_mn[s][n][m].calculated)
		{
			christoffel_s_mn[s][m][n].e
				= christoffel_s_mn[s][n][m].e;
		}
		else
		{
			ex toffel = 0;
			FOR(lambda)
			{
				toffel = toffel + G(s,lambda)
					* (
			 		g_diff(lambda,n,m)
					+ g_diff(m,lambda,n)
					- g_diff(m,n,lambda)
					);
			}
			christoffel_s_mn[s][m][n].e = (toffel/2)NORM;
		}
		christoffel_s_mn[s][m][n].calculated = true;
	}

	return christoffel_s_mn[s][m][n].e;
}


// Christoffel^a_na
EXFUNC_HEAD(christoffel_a_ma, ARG1(n))
{
	if (!chris_a_na[n].calculated)
	{
		// TODO: Use 1/2 * g^lm * g_lm,n instead.
		ex chris = 0;
		FOR(a)
		{
			chris = chris + christoffel(a, n,a);
		}
		chris_a_na[n].e = (chris)NORM;
		chris_a_na[n].calculated = true;
	}

	return chris_a_na[n].e;
}


// chris^a_ma_n
EXFUNC_HEAD(christoffel_a_ma_n, ARG2(m,n))
{
	if (!chris_a_na_m[m][n].calculated)
	{
		if (chris_a_na_m[n][m].calculated)
		{
			chris_a_na_m[m][n].e = chris_a_na_m[n][m].e;
		}
		else
		{
			chris_a_na_m[m][n].e
				= christoffel_a_ma(m).diff(coord[n])NORM;
		}
		chris_a_na_m[m][n].calculated = true;
	}

	return chris_a_na_m[m][n].e;
}


// chris^a_mn_a
EXFUNC_HEAD(christoffel_a_mn_a, ARG2(m,n))
{
	if (!chris_a_mn_a[m][n].calculated)
	{
		if (chris_a_mn_a[n][m].calculated)
		{
			chris_a_mn_a[m][n].e = chris_a_mn_a[n][m].e;
		}
		else
		{
			ex toffel = 0;
			FOR(a)
			{
				toffel = toffel
					+ christoffel(a, m,n).diff(coord[a]);
			}
			chris_a_mn_a[m][n].e = (toffel)NORM;
		}
		chris_a_mn_a[m][n].calculated = true;
	}
	return chris_a_mn_a[m][n].e;
}


// chris^a_mn * chris^b_ab
EXFUNC_HEAD(christoffel_a_mn_christoffel_b_ab, ARG2(m,n))
{
	if(!chris_a_mn_chris_b_ab[m][n].calculated)
	{
		if (chris_a_mn_chris_b_ab[n][m].calculated)
		{
			chris_a_mn_chris_b_ab[m][n].e
				= chris_a_mn_chris_b_ab[n][m].e;
		}
		else
		{
			ex tof = 0;
			FOR(a)
			{
				tof = tof + christoffel(a,m,n)
					* christoffel_a_ma(a);
			}
			chris_a_mn_chris_b_ab[m][n].e = (tof)NORM;
		}
		chris_a_mn_chris_b_ab[m][n].calculated = true;
	}
	return chris_a_mn_chris_b_ab[m][n].e;
}


// chris^a_mb * chris^b_na
EXFUNC_HEAD(christoffel_a_mb_christoffel_b_na, ARG2(m,n))
{
	if(!chris_a_mb_chris_b_na[m][n].calculated)
	{
		if (chris_a_mb_chris_b_na[n][m].calculated)
		{
			chris_a_mb_chris_b_na[m][n].e
				= chris_a_mb_chris_b_na[n][m].e;
		}
		else
		{
			ex tof = 0;
			FOR2(a,b)
			{
				tof = tof + christoffel(a,m,b)
					* christoffel(b,n,a);
			}
			chris_a_mb_chris_b_na[m][n].e = (tof)NORM;
		}
		chris_a_mb_chris_b_na[m][n].calculated = true;
	}
	return chris_a_mb_chris_b_na[m][n].e;
}


// riemann_mn
EXFUNC_HEAD(riemann, ARG2(m,n))
{
	if(!riemann_mn[m][n].calculated)
	{
		ex rie;
		if (riemann_mn[n][m].calculated)
		{
			rie = riemann_mn[n][m].e;
		}
		else
		{
			rie = christoffel_a_ma_n(m,n)
				- christoffel_a_mn_a(m,n)
				- christoffel_a_mn_christoffel_b_ab(m,n)
				+ christoffel_a_mb_christoffel_b_na(m,n);
			rie = (rie)NORMR;
		}
		riemann_mn[m][n].e = rie;
		riemann_mn[m][n].calculated = true;
	}
	return riemann_mn[m][n].e;
}


// Ricci: R^m_n
EXFUNC_HEAD(ricci, ARG2(m,n))
{
	if (!ricci_m_n[m][n].calculated)
	{
		ex ric = 0;
		FOR(a)
		{
			ric = ric + G(m,a) * riemann(a,n);
		}
		ricci_m_n[m][n].e = (ric)NORMR;
		ricci_m_n[m][n].calculated = true;
	}
	return ricci_m_n[m][n].e;
}


// Gauss: R
EXFUNC_HEAD(gauss, ARG0)
{
	if (!gauss_a_a.calculated)
	{
		ex gau = 0;
		FOR(a)
		{
			gau = gau + ricci(a,a);
		}
		gauss_a_a.e = (gau)NORMR;
		gauss_a_a.calculated = true;
	}
	return gauss_a_a.e;
}


// Einstein: G^m_n
EXFUNC_HEAD(einstein, ARG2(m,n))
{
	if (!einstein_m_n[m][n].calculated)
	{
		ex ein = 0;
		if (m == n)
		{
			ein = ricci(m,n) - gauss()/2;
			ein = (ein)NORMR;
		}
		else
		{
			ein = ricci(m,n);
		}
		einstein_m_n[m][n].e = ein;
		einstein_m_n[m][n].calculated = true;
	}
	return einstein_m_n[m][n].e;
}


// Einstein: G
EXFUNC_HEAD(einstein_scalar, ARG0)
{
	if (!einstein_trace.calculated)
	{
		ex ein = 0;
		FOR(a)
		{
			ein = ein + einstein(a,a);
		}
		einstein_trace.e = (ein)NORMR;
		einstein_trace.calculated = true;
	}
	return einstein_trace.e;
}


// Einstein: G^a_m;a
EXFUNC_HEAD(div_einstein, ARG1(m))
{
	if (!div_einstein_m[m].calculated)
	{
		ex div = 0;
		FOR(a)
		{
			div = div + einstein(a,m).diff(coord[a]);
			FOR(b)
			{
				div = div + christoffel(a,a,b) * einstein(b,m);
				div = div - christoffel(b,m,a) * einstein(a,b);
			}
		}
		div_einstein_m[m].e = (div)NORMR;
		div_einstein_m[m].calculated = true;
	}
	return div_einstein_m[m].e;
}




// check harmony g^mn * L^a_mn
EXFUNC_HEAD(harmony, ARG1(m))
{
	if (!harmony_m[m].calculated)
	{
		ex har = 0;
		FOR2(a,b)
		{
			har = har + g(a,b) * christoffel(m,a,b);
		}
		harmony_m[m].e = (har)NORMR;
		harmony_m[m].calculated = true;
	}
	return harmony_m[m].e;
}



ex Coordinate_System::det_g()
{
	matrix gmatrix(dim,dim);
	FOR2(m,n)
	{
		gmatrix(m,n) = g(m,n);
	}
	return gmatrix.determinant();
}




// compute inverse metric G from g
void Coordinate_System::calc_inverse_metric()
{
	// Dummy symbols
	symbol **gmn = new symbol* [dim];
	FOR(a)
	{
		gmn[a] = new symbol [dim];
	}

	// Create substitution map
	exmap substitute;
	FOR2(i,j)
	{
		substitute[gmn[i][j]] = g(i,j);
	}


	// Construct dummy matrix
	matrix gmatrix(dim, dim);
	FOR2(i,j)
	{
		gmatrix(i, j) = gmn[i][j];
	}


	// Invert dummy matrix
	matrix Gmatrix = gmatrix.inverse();


	// Substitute and copy to output
	FOR2(i,j)
	{
		if (j < i)
			G_mn[i][j].e = G_mn[j][i].e;
		else
			G_mn[i][j].e
				= Gmatrix(i, j).subs(substitute)NORMG;
	}

	// unalloc dummy symbols
	FOR(a)
	{
		//delete gmn[a];
	}
	delete gmn;
}





// init the elements
void Coordinate_System::init(const unsigned& dimension)
{
	eq0 e0;
	e0.e = 0;
	e0.calculated = false;
	eq1 e1(dimension, e0);
	eq2 e2(dimension, e1);
	eq3 e3(dimension, e2);

	G_mn = e2;

	g_mn_s = e3;

	christoffel_s_mn = e3;

	chris_a_na = e1;

	chris_a_na_m = e2;
	chris_a_mn_a = e2;
	chris_a_mn_chris_b_ab = e2;
	chris_a_mb_chris_b_na = e2;

	riemann_mn = e2;
	ricci_m_n = e2;
	gauss_a_a = e0;

	einstein_m_n = e2;
	einstein_trace = e0;
	div_einstein_m = e1;

	harmony_m = e1;
}


// Eliminate beginning Whitespace
void Coordinate_System::eliminate_beginning_whitespace(string& s)
{
	unsigned i = 0;
	while ( (s[i] == ' ') || (s[i] == '	') )
		i++;

	s.erase(0, i);
}

// Eliminate ending whitespace
void Coordinate_System::eliminate_ending_whitespace(string& s)
{
	int i = s.length() - 1;
	while ( (s[i] == ' ') || (s[i] == '	') )
		i--;

	s.erase(i + 1);
}


// Eliminate all comments (everything on a line starting from '#').
void Coordinate_System::eliminate_comments(string& s)
{
	unsigned i = 0;
	while ( (s[i] != '#') && (i < s.length()) )
		i++;

	s.erase(i);
}


// test for a long comment
bool Coordinate_System::long_comment_here(const string& s)
{
	string long_comment_str("##########");
	if (s.compare(0, long_comment_str.size(), long_comment_str))
		return false;
	return true;
}


// Returns the number of occurrences of the character 'ch' in 's'.
unsigned Coordinate_System::contains(const char& ch, const string& s)
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
void Coordinate_System::split_at_equalsign(const string& s,
		string& before, string& after)
{
	unsigned i = s.find_first_of('=');

	before.assign(s, 0, i);
	after.assign(s, i + 1, s.length() - i - 1);

	eliminate_beginning_whitespace(after);
	eliminate_ending_whitespace(before);
}



// Already in list?
void Coordinate_System::inlst(const unsigned& lineno, const string& sym,
		const lst& symbollist)
{
	for (lst::const_iterator i = symbollist.begin();
			i != symbollist.end();
			i++)
	{
		stringstream ss;
		ss << *i;
		if (ss.str() == sym)
		{
			cerr << lineno << ": ERROR: Symbol '" << sym
				<< "' previously defined." << endl;
			exit(1);
		}
	}
}


// Split a string of numbers in the middle in two, convert each to an integer.
void Coordinate_System::split_atoi(const string& s, unsigned& a, unsigned& b)
{
	if (s.size() % 2)
	{
		cerr << "Cannot split \"" << s << "\" in two. "
			<< "Odd number of elements." << endl;
		exit(1);
	}
	string aa(s, 0, s.size() / 2);
	string bb(s, s.size() / 2, s.size() / 2);
	a = atoi(aa.c_str());
	b = atoi(bb.c_str());
}


// Resize g_mn[i] to the specified size at least.
void Coordinate_System::resize_g_mni(const unsigned& i, const unsigned& dim)
{
	if (g_mn[i].size() < dim)
	{
		eq0 e0;
		e0.e = 0;
		e0.calculated = false;
		g_mn[i].resize(dim, e0);
	}
}




// Constructor, the second
Coordinate_System::Coordinate_System(const char * filename)
{
	ifstream fin(filename);
	string s;

	lst symbollist;
	exmap symbolmap;

	unsigned lineno = 0;
	bool long_comment = false;

	while (getline(fin, s))
	{
		lineno++;

		eliminate_beginning_whitespace(s);
		if (long_comment_here(s)) long_comment = !long_comment;
		if (long_comment) continue;
		eliminate_comments(s);
		eliminate_ending_whitespace(s); // Must after eliminate_comments
		if (s.empty()) continue;
		if (contains('=', s) > 1)
		{
			cerr << lineno << ": Syntax error: Too many '='."
				<< endl;
			exit(1);
		}


		cerr << lineno << ": " << s << "." << endl;

		// Name? (Starts with '"'.)
		if (s[0] == '"')
		{
			name += s.substr(1, s.length() - 2);
			continue;
		}

		// Is it a declaration? (Contains no '='.)
		if (!contains('=', s))
		{
			// Already declared?
			inlst(lineno, s, symbollist);

			// Is it a constant?
			if (!contains('(', s))
			{
				symbollist.append(symbol(s));
			}
			else if (contains(',', s) == 1)
			{
				// Is it a function of 1 variable?
				// TODO
			}
			else
			{
				// TODO
				// ..of 2 variables?
				// ..of 3 variables?
				// ...
			}
		}
		else // contains '='
		{
			string before;
			string after;
			split_at_equalsign(s, before, after);

			// Already defined?
			inlst(lineno, before, symbollist);

			symbol symname(before);
			symbollist.append(symname);

			// Is it a coordinate (starts with 'x')?
			if ( (before[0] == 'x') && (isdigit(before[1])) )
			{
				before.erase(0, 1);
				// TODO: stringstream?
				unsigned mu = atoi(before.c_str());
				if (coord.size() < mu + 1)
					coord.resize(mu + 1);
				// else error

				coord[mu] = symbol(after);
				symbollist.append(coord[mu]);
			}


			// All equations, regardless of race and ethnicity,
			// go in symbolmap.
			try {
				symbolmap[symname] = ex(after, symbollist);
			} catch (exception &p) {
				cerr << p.what() << endl;
				exit(1);
			}


			// Is it the metric (starts with 'g')?
			if ( (before[0] == 'g')
					&& isdigit(before[1])
					&& isdigit(before[2]) )
			{
				before.erase(0, 1);
				unsigned m, n;
				split_atoi(before, m, n);

				if (g_mn.size() < m + 1)
					g_mn.resize(m + 1);
				resize_g_mni(m, n + 1);

				g_mn[m][n].e = symbolmap[symname]
					.subs(symbolmap)
					.subs(symbolmap)
					.subs(symbolmap)
					.subs(symbolmap)
					.subs(symbolmap)
					.subs(symbolmap)
					.subs(symbolmap);
			}
		} // contains '='
	}

	// check dimensions
	dim = coord.size();
	if (dim != g_mn.size())
	{
		cerr << "Dimensions do not match. "
			<< "dim = " << dim << ", "
			<< "g_mn.size() = " << g_mn.size() << endl;
		exit(1);
	}
	for (unsigned i = 0; i < dim; i++)
	{
		if (dim < g_mn[i].size())
		{
			cerr << "Dimensions do not match. "
				<< "dim = " << dim << ", "
				<< "g_mn[" << i << "].size() = "
				<< g_mn[i].size() << endl;
			exit(1);
		}
		else if (g_mn[i].size() < dim)
		{
			resize_g_mni(i, dim);
		}
	}

	init(dim);

	cerr << "dim = " << dim << endl;
}


// Destruktor
Coordinate_System::~Coordinate_System()
{
}
