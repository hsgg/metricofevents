#!/bin/sh

cut -f 1-2 -d' ' <plot0.dat >pl0.dat
grep ^3.14159 pl0.dat >pl.dat
grep ^9.42477 pl0.dat >pl9.dat

