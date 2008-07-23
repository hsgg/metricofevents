/*
 * HSGG Tü20070228
 */

#include <math.h>
#include <stdio.h>

#include "global.h"
#include "init.h"
#include "myfloat.h"

using namespace std;


// ********* MAIN ************
int main()
{
	#ifdef USE_CLN
		const double m = double_approx(initials().m);
		const double r = double_approx(initials().radius);
	#else
		const double m = initials().m;
		const double r = initials().radius;
	#endif

	FILE* circ_file = fopen("sphere.dat", "w");
	FILE* circ2_file = fopen("sphere2.dat", "w");

	fprintf(circ_file, "# m = %f\n", m);
	fprintf(circ2_file, "# r = %f\n", r);

	// Radius ist 2 * m
	for (int i = 0; i < 36; i++)
	{
		for (int j = 0; j < 18; j++)
		{
			/*
			fprintf(circ_file, "%e\t%e\t%e\n",
				2 * m * cos(i * M_PI / 18) * sin(j * M_PI / 18),
				2 * m * sin(i * M_PI / 18) * sin(j * M_PI / 18),
				2 * m			   * cos(j * M_PI / 18)
			);
			*/
			fprintf(circ_file, "%e\t%e\t%e\n",
				i * M_PI / 18,
				2 * m * j / 17.0,
				j * M_PI / 18);
			fprintf(circ2_file, "%e\t%e\t%e\n",
				r * cos(i * M_PI / 18)	* sin(j * M_PI / 18),
				r * sin(i * M_PI / 18)	* sin(j * M_PI / 18),
				r			* cos(j * M_PI / 18)
			);
		}
	}
	fclose(circ_file);
	fclose(circ2_file);

	return 0;
}
