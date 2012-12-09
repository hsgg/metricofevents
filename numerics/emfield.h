/*
 * HSGG Tü, 2007.05.05
 */

#ifndef EMFIELD_H
#define EMFIELD_H

#include <vector>

#include "metric.h"
#include "myfloat.h"

class EMField
{
    public:
	// Konstruktor
	EMField(const Metric& metric);
	// Destruktor
	~EMField();

	// Konstanten, Metric
	const Metric* metric;

	// EM Fieldstrengthtensor
	const myfloat (EMField::*F[DIM][DIM])(const myfloat* x) const;

    private:
	inline const myfloat rho4(const myfloat* x) const;

	const myfloat F_00(const myfloat* x) const;
	const myfloat F_01(const myfloat* x) const;
	const myfloat F_02(const myfloat* x) const;
	const myfloat F_03(const myfloat* x) const;
	const myfloat F_10(const myfloat* x) const;
	const myfloat F_11(const myfloat* x) const;
	const myfloat F_12(const myfloat* x) const;
	const myfloat F_13(const myfloat* x) const;
	const myfloat F_20(const myfloat* x) const;
	const myfloat F_21(const myfloat* x) const;
	const myfloat F_22(const myfloat* x) const;
	const myfloat F_23(const myfloat* x) const;
	const myfloat F_30(const myfloat* x) const;
	const myfloat F_31(const myfloat* x) const;
	const myfloat F_32(const myfloat* x) const;
	const myfloat F_33(const myfloat* x) const;
};

#endif
