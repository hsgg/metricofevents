/*
 * TÃ¼, 2006/01/29 HSGG, GNULin
 * Tü, 2006/01/29 HSGG, WinXP
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

double timediff(const clock_t ta, const clock_t tb)
{
	return (tb - ta) / (double) CLOCKS_PER_SEC;
}

clock_t msg_starting()
{
	fprintf(stderr, "Starting...");
	return clock();
}

void msg_finished(const clock_t ta)
{
	fprintf(stderr, "done. Time taken: %6.2f\n", timediff(ta, clock()));
}
