/*
 * HSGG Tü20070318
 */


	// Define Kerr-metric:

	c.name_of_metric = "Kerr";

	// Hilfsexpressions
	ex delta = pow(r,2) + pow(a,2) - 2 * m * r;
	ex rho2 = pow(r,2) + pow(a * cos(theta),2);
	ex sigma2 = pow(pow(r,2) + pow(a,2),2) - pow(a * sin(theta),2) * delta;
	ex omega = 2 * m * r * a / sigma2;

	// metric
	c.g[0][0] = rho2 * delta / sigma2
		- pow(omega * sin(theta),2) * sigma2 / rho2 ;
	c.g[0][1] = 0;
	c.g[0][2] = 0;
	c.g[0][3] = omega * sigma2 / rho2 * pow(sin(theta),2);
	c.g[1][0] = c.g[0][1];
	c.g[1][1] = -rho2 / delta;
	c.g[1][2] = 0;
	c.g[1][3] = 0;
	c.g[2][0] = c.g[0][2];
	c.g[2][1] = c.g[1][2];
	c.g[2][2] = -rho2;
	c.g[2][3] = 0;
	c.g[3][0] = c.g[0][3];
	c.g[3][1] = c.g[1][3];
	c.g[3][2] = c.g[2][3];
	c.g[3][3] = -sigma2 / rho2 * pow(sin(theta),2);

