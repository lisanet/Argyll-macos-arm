
/*
 *  Roots3And4.c
 *
 *  Utility functions to find cubic and quartic roots,
 *  coefficients are passed like this:
 *
 *	  c[0] + c[1]*x + c[2]*x^2 + c[3]*x^3 + c[4]*x^4 = 0
 *
 *  The functions return the number of non-complex roots and
 *  put the values into the s array.
 *
 *  Author:		 Jochen Schwarze (schwarze@isa.de)
 *               (From Graphics Gems I)
 */

int SolveQuadric(double c[3], double s[2]);

int SolveCubic(double c[4], double s[3]);

int SolveQuartic(double c[5], double s[4]);


