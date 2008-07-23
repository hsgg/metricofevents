/*
 * HSGG Tü20070319
 */


#ifndef MYFLOAT_H
#define MYFLOAT_H

#include "global.h"


#ifdef USE_CLN
	#include <cln/cln.h>
#endif




#ifdef USE_CLN
	typedef cln::cl_R myfloat;
#else
	typedef long double myfloat;
#endif



#endif
