/*
 * HSGG Tü2007062x
 * HSGG Tü20070810
 */

#include <ginac/ginac.h>
#include "functions.h"

using namespace GiNaC;


// declare and register a function of 1 variable
#define MAKE_FUNC1(name, latex, funcprime, \
		name0, latex0, func0prime, \
		name00, latex00) \
	\
	DECLARE_FUNCTION_1P(name00) \
	REGISTER_FUNCTION(name00, latex_name(latex00)) \
	\
	\
	static ex func0prime(const ex & x, unsigned diff_param) \
	{ \
		return name00(x); \
	} \
	DECLARE_FUNCTION_1P(name0) \
	REGISTER_FUNCTION(name0, latex_name(latex0). \
			derivative_func(func0prime)) \
	\
	\
	static ex funcprime(const ex & x, unsigned diff_param) \
	{ \
		return name0(x); \
	} \
	DECLARE_FUNCTION_1P(name) \
	REGISTER_FUNCTION(name, latex_name(latex). \
			derivative_func(funcprime))


// declare and register a function of 2 variables
#define MAKE_FUNC2(name, latex, funcprime, \
		name1, latex1, func1prime, \
		name2, latex2, func2prime, \
		name11, latex11, \
		name12, latex12, \
		name22, latex22) \
	\
	DECLARE_FUNCTION_2P(name11) \
	REGISTER_FUNCTION(name11, latex_name(latex11)) \
	\
	DECLARE_FUNCTION_2P(name12) \
	REGISTER_FUNCTION(name12, latex_name(latex12)) \
	\
	DECLARE_FUNCTION_2P(name22) \
	REGISTER_FUNCTION(name22, latex_name(latex22)) \
	\
	\
	DECLARE_FUNCTION_2P(name1) \
	static ex func1prime(const ex & x, const ex & y, unsigned diff_param) \
	{ \
		ex retval = 0; \
		if (diff_param == 0) \
			retval = name11(x,y); \
		if (diff_param == 1) \
			retval = name12(x,y); \
		return retval; \
	} \
	REGISTER_FUNCTION(name1, latex_name(latex1). \
			derivative_func(func1prime)) \
	\
	DECLARE_FUNCTION_2P(name2) \
	static ex func2prime(const ex & x, const ex & y, unsigned diff_param) \
	{ \
		ex retval = 0; \
		if (diff_param == 0) \
			retval = name12(x,y);   /* for symmetry */ \
		if (diff_param == 1) \
			retval = name22(x,y); \
		return retval; \
	} \
	REGISTER_FUNCTION(name2, latex_name(latex2). \
			derivative_func(func2prime)) \
	\
	\
	DECLARE_FUNCTION_2P(name) \
	static ex funcprime(const ex & x, const ex & y, unsigned diff_param) \
	{ \
		ex retval = 0; \
		if (diff_param == 0) \
			retval = name1(x,y); \
		if (diff_param == 1) \
			retval = name2(x,y); \
		return retval; \
	} \
	REGISTER_FUNCTION(name, latex_name(latex). \
			derivative_func(funcprime))



// define a function with partial derivatives of 3 variables
#define FUNC3_PRIME(name, latex, nameprime, \
		name1, name2, name3) \
	DECLARE_FUNCTION_3P(name) \
	static ex nameprime(const ex & x, const ex & y, const ex & z, \
			unsigned diff_param) \
	{ \
		ex retval = 0; \
		if (diff_param == 0) \
			retval = name1(x,y,z); \
		if (diff_param == 1) \
			retval = name2(x,y,z); \
		if (diff_param == 2) \
			retval = name3(x,y,z); \
		return retval; \
	} \
	REGISTER_FUNCTION(name, latex_name(latex). \
			derivative_func(nameprime))

// declare and register a function of 3 variables
#define MAKE_FUNC3(name, latex, funcprime, \
		name1, latex1, func1prime, \
		name2, latex2, func2prime, \
		name3, latex3, func3prime, \
		name11, latex11, \
		name12, latex12, \
		name13, latex13, \
		name22, latex22, \
		name23, latex23, \
		name33, latex33) \
	\
	DECLARE_FUNCTION_3P(name11) \
	REGISTER_FUNCTION(name11, latex_name(latex11)) \
	\
	DECLARE_FUNCTION_3P(name12) \
	REGISTER_FUNCTION(name12, latex_name(latex12)) \
	\
	DECLARE_FUNCTION_3P(name13) \
	REGISTER_FUNCTION(name13, latex_name(latex13)) \
	\
	DECLARE_FUNCTION_3P(name22) \
	REGISTER_FUNCTION(name22, latex_name(latex22)) \
	\
	DECLARE_FUNCTION_3P(name23) \
	REGISTER_FUNCTION(name23, latex_name(latex23)) \
	\
	DECLARE_FUNCTION_3P(name33) \
	REGISTER_FUNCTION(name33, latex_name(latex33)) \
	\
	\
	FUNC3_PRIME(name1, latex1, func1prime, \
			name11, name12, name13) \
	\
	FUNC3_PRIME(name2, latex2, func2prime, \
			name12, name22, name23) \
	\
	FUNC3_PRIME(name3, latex3, func3prime, \
			name13, name23, name33) \
	\
	\
	FUNC3_PRIME(name, latex, funcprime, \
			name1, name2, name3)


/* 4 Variables */
// define a function of 4 variables
#define FUNC4(name, latex) \
	DECLARE_FUNCTION_4P(name) \
	REGISTER_FUNCTION(name, latex_name(latex))

// define a function of 4 variables with partial derivatives
#define FUNC4_PRIME(name, latex, nameprime, \
		name0prime, name1prime, name2prime, name3prime) \
	DECLARE_FUNCTION_4P(name) \
	static ex nameprime(const ex & x, const ex & y, const ex & z, \
			const ex & w, unsigned diff_param) \
	{ \
		ex retval = 0; \
		if (diff_param == 0) \
			retval = name0prime(x,y,z,w); \
		if (diff_param == 1) \
			retval = name1prime(x,y,z,w); \
		if (diff_param == 2) \
			retval = name2prime(x,y,z,w); \
		if (diff_param == 3) \
			retval = name3prime(x,y,z,w); \
		return retval; \
	} \
	REGISTER_FUNCTION(name, latex_name(latex). \
			derivative_func(nameprime))

// declare and register a function of 4 variables
#define MAKE_FUNC4(name, latex, funcprime, \
		name0, latex0, func0prime, \
		name1, latex1, func1prime, \
		name2, latex2, func2prime, \
		name3, latex3, func3prime, \
		name00, latex00, \
		name01, latex01, \
		name02, latex02, \
		name03, latex03, \
		name11, latex11, \
		name12, latex12, \
		name13, latex13, \
		name22, latex22, \
		name23, latex23, \
		name33, latex33) \
	\
	FUNC4(name00, latex00) \
	FUNC4(name01, latex01) \
	FUNC4(name02, latex02) \
	FUNC4(name03, latex03) \
	\
	FUNC4(name11, latex11) \
	FUNC4(name12, latex12) \
	FUNC4(name13, latex13) \
	\
	FUNC4(name22, latex22) \
	FUNC4(name23, latex23) \
	\
	FUNC4(name33, latex33) \
	\
	\
	FUNC4_PRIME(name0, latex0, func0prime, \
			name00, name01, name02, name03) \
	\
	FUNC4_PRIME(name1, latex1, func1prime, \
			name01, name11, name12, name13) \
	\
	FUNC4_PRIME(name2, latex2, func2prime, \
			name02, name12, name22, name23) \
	\
	FUNC4_PRIME(name3, latex3, func3prime, \
			name03, name13, name23, name33) \
	\
	\
	FUNC4_PRIME(name, latex, funcprime, \
			name0, name1, name2, name3)



// a(t)
MAKE_FUNC1(a, "a", adot_func,
	ad_func, "\\dot{a}", addot_func,
	add_func, "\\ddot{a}")

// c(t) = sol(t) (speed-of-light)
MAKE_FUNC1(c, "c", cdot,
	cd, "\\dot{c}", cddot,
	cdd, "\\ddot{c}")


// R(t)
MAKE_FUNC1(R, "R", Rdot,
        Rd, "\\dot{R}", Rddot,
        Rdd, "\\ddot{R}")

// H(t)
MAKE_FUNC1(Hu, "H", Hdot,
        Hd, "\\dot{H}", Hddot,
        Hdd, "\\ddot{H}")

// A(r)
MAKE_FUNC1(A, "A", Ap,
	A0, "A'", A0p,
	A00, "A''")

// B(r)
MAKE_FUNC1(B, "B", Bp,
	B0, "B'", B0p,
	B00, "B''")

// C(r)
MAKE_FUNC1(C, "C", Cp,
	C0, "C'", C0p,
	C00, "C''")

// D(r)
MAKE_FUNC1(D, "D", Dp,
	D0, "D'", D0p,
	D00, "D''")

// F(r)
MAKE_FUNC1(F, "F", Fp,
	F0, "F'", F0p,
	F00, "F''")

// G(r)
MAKE_FUNC1(Gm, "G", Gp,
	G0, "G'", G0p,
	G00, "G''")

// m(r)
//MAKE_FUNC1(mass, "m",
//	mp, "m'", mprime,
//	mpp, "m''", mpprime)

// m(t,r)
/*MAKE_FUNC2(mass, "m", mprime,
	m1p, "\\dot{m}", m1prime,
	m2p, "m'", m2prime,
	m11pp, "\\ddot{m}",
	m12pp, "\\dot{m}'",
	m22pp, "m''")
*/


// r(R)
MAKE_FUNC1(	r, "r", rp,
		r0, "r'", r0p,
		r00, "r''"	)

// rho(r)
MAKE_FUNC1(	rho, "\\rho", rhop,
		rho0, "\\rho'", rho0p,
		rho00, "\\rho''"	)

// epsilon(r)
MAKE_FUNC1(	epsilon, "\\epsilon", epsilonp,
		epsilon0, "\\epsilon'", epsilon0p,
		epsilon00, "\\epsilon''"	)

// rho(t,r)
MAKE_FUNC2(rhotr, "\\rho", rhotr_n,
	rhotr0, "\\dot{\\rho}", rhotr0prime,
	rhotr1, "\\rho'", rhotr1prime,
	rhotr11, "\\ddot{\\rho}",
	rhotr12, "\\dot{\\rho \\mathbf{NOTICE}}'",
	rhotr22, "\\rho''")

// aa(t,x)
MAKE_FUNC2(aa, "a", aa_n,
	aa0, "\\dot{a}", aa0prime,
	aa1, "a'", aa1prime,
	aa11, "\\ddot{a}",
	aa12, "\\dot{a}'",
	aa22, "a''")

// bb(t,x)
MAKE_FUNC2(bb, "b", bb_n,
	bb0, "\\dot{b}", bb0prime,
	bb1, "b'", bb1prime,
	bb11, "\\ddot{b}",
	bb12, "\\dot{b \\mathbf{NOTICE}}'",
	bb22, "b''")

// cc(t,x)
MAKE_FUNC2(cc, "c", cc_n,
	cc0, "\\dot{c}", cc0prime,
	cc1, "c'", cc1prime,
	cc11, "\\ddot{c}",
	cc12, "\\dot{c \\mathbf{NOTICE}}'",
	cc22, "c''")

// dd(t,x)
MAKE_FUNC2(dd, "d", dd_n,
	dd0, "\\dot{d}", dd0prime,
	dd1, "d'", dd1prime,
	dd11, "\\ddot{d}",
	dd12, "\\dot{d \\mathbf{NOTICE}}'",
	dd22, "d''")

// E(rho, z)
MAKE_FUNC2(E, "E", E_n,
	E0, "\\dot{E}", E0prime,
	E1, "E'", E1prime,
	E11, "\\ddot{E}",
	E12, "\\dot{E'}",
	E22, "E''")

// masse(x,y,z)
MAKE_FUNC3(masse, "m", massep,
	masse1, "m'", masse1p,
	masse2, "\\tilde{m}", masse2p,
	masse3, "\\hat{m}", masse3p,
	masse11, "m''",
	masse12, "\\tilde{m}'",
	masse13, "\\hat{m}'",
	masse22, "\\tilde{\\tilde{m}}",
	masse23, "\\tilde{\\hat{m}}",
	masse33, "\\hat{\\hat{m}}")

// massig(t,r,theta,phi)
MAKE_FUNC4(massig, "m", mp,
	m0, "\\dot{m}", m0p,
	m1, "m'", m1p,
	m2, "\\tilde{m}", m2p,
	m3, "\\hat{m}", m3p,
	m00, "\\ddot{m}",
	m01, "\\dot{m}'",
	m02, "\\tilde{\\dot{m}}",
	m03, "\\hat{\\dot{m}}",
	m11, "m''",
	m12, "\\tilde{m}'",
	m13, "\\hat{m}'",
	m22, "\\tilde{\\tilde{m}}",
	m23, "\\tilde{\\hat{m}}",
	m33, "\\hat{\\hat{m}}")

// P(t,r,theta,phi)
/*MAKE_FUNC4(Potential, "\\Phi", Pp, \
	P0, "\\dot{\\Phi}", P0p, \
	P1, "\\Phi'", P1p, \
	P2, "\\tilde{\\Phi}", P2p, \
	P3, "\\hat{\\Phi}", P3p, \
	P00, "\\ddot{\\Phi}", \
	P01, "\\dot{\\Phi}'", \
	P02, "\\tilde{\\dot{\\Phi}}", \
	P03, "\\hat{\\dot{\\Phi}}", \
	P11, "\\Phi''", \
	P12, "\\tilde{\\Phi}'", \
	P13, "\\hat{\\Phi}'", \
	P22, "\\tilde{\\tilde{\\Phi}}", \
	P23, "\\tilde{\\hat{\\Phi}}", \
	P33, "\\hat{\\hat{\\Phi}}")
*/

// P(r,theta,phi)
MAKE_FUNC3(Potential, "\\Phi", Pp,
	P1, "\\Phi'", P1p,
	P2, "\\tilde{\\Phi}", P2p,
	P3, "\\hat{\\Phi}", P3p,
	P11, "\\Phi''",
	P12, "\\tilde{\\Phi}'",
	P13, "\\hat{\\Phi}'",
	P22, "\\tilde{\\tilde{\\Phi}}",
	P23, "\\tilde{\\hat{\\Phi}}",
	P33, "\\hat{\\hat{\\Phi}}")


// g_mn(t,x,y,z)
#define MAKE_G_MN(munu) \
	MAKE_FUNC4(g_##munu, "g_{"#munu"}", g_##munu##p, \
			g_##munu##_0, "g_{"#munu",0}", g_##munu##_0p, \
			g_##munu##_1, "g_{"#munu",1}", g_##munu##_1p, \
			g_##munu##_2, "g_{"#munu",2}", g_##munu##_2p, \
			g_##munu##_3, "g_{"#munu",3}", g_##munu##_3p, \
			g_##munu##_00, "g_{"#munu",00}", \
			g_##munu##_01, "g_{"#munu",01}", \
			g_##munu##_02, "g_{"#munu",02}", \
			g_##munu##_03, "g_{"#munu",03}", \
			g_##munu##_11, "g_{"#munu",11}", \
			g_##munu##_12, "g_{"#munu",12}", \
			g_##munu##_13, "g_{"#munu",13}", \
			g_##munu##_22, "g_{"#munu",22}", \
			g_##munu##_23, "g_{"#munu",23}", \
			g_##munu##_33, "g_{"#munu",33}")
MAKE_G_MN(00)
MAKE_G_MN(01)
MAKE_G_MN(02)
MAKE_G_MN(03)
MAKE_G_MN(11)
MAKE_G_MN(12)
MAKE_G_MN(13)
MAKE_G_MN(22)
MAKE_G_MN(23)
MAKE_G_MN(33)
