
/* Test the CIE delta E 2000 code */

#include <stdio.h>
#include <math.h>
#include "icc.h"


/* Reference data */

#define NTESTS 34

/* From the Sharma, Wu and Dalal "Implementation Notes" etc. paper: */
struct {
	double Lab1[3];
	double Lab2[3];
	double de;
} ref[NTESTS] = {
	{ 50.0000, 2.6772, -79.7751, 50.0000, 0.0000, -82.7485, 2.0425 },
	{ 50.0000, 3.1571, -77.2803, 50.0000, 0.0000, -82.7485, 2.8615 },
	{ 50.0000, 2.8361, -74.0200, 50.0000, 0.0000, -82.7485, 3.4412 },
	{ 50.0000, -1.3802, -84.2814, 50.0000, 0.0000, -82.7485, 1.0000 },
	{ 50.0000, -1.1848, -84.8006, 50.0000, 0.0000, -82.7485, 1.0000 },
	{ 50.0000, -0.9009, -85.5211, 50.0000, 0.0000, -82.7485, 1.0000 },
	{ 50.0000, 0.0000, 0.0000, 50.0000, -1.0000, 2.0000, 2.3669 },
	{ 50.0000, -1.0000, 2.0000, 50.0000, 0.0000, 0.0000, 2.3669 },
	{ 50.0000, 2.4900, -0.0010, 50.0000, -2.4900, 0.0009, 7.1792 },
	{ 50.0000, 2.4900, -0.0010, 50.0000, -2.4900, 0.0010, 7.1792 },
	{ 50.0000, 2.4900, -0.0010, 50.0000, -2.4900, 0.0011, 7.2195 },
	{ 50.0000, 2.4900, -0.0010, 50.0000, -2.4900, 0.0012, 7.2195 },
	{ 50.0000, -0.0010, 2.4900, 50.0000, 0.0009, -2.4900, 4.8045 },
	{ 50.0000, -0.0010, 2.4900, 50.0000, 0.0010, -2.4900, 4.8045 },
	{ 50.0000, -0.0010, 2.4900, 50.0000, 0.0011, -2.4900, 4.7461 },
	{ 50.0000, 2.5000, 0.0000, 50.0000, 0.0000, -2.5000, 4.3065 },
	{ 50.0000, 2.5000, 0.0000, 73.0000, 25.0000, -18.0000, 27.1492 },
	{ 50.0000, 2.5000, 0.0000, 61.0000, -5.0000, 29.0000, 22.8977 },
	{ 50.0000, 2.5000, 0.0000, 56.0000, -27.0000, -3.0000, 31.9030 },
	{ 50.0000, 2.5000, 0.0000, 58.0000, 24.0000, 15.0000, 19.4535 },
	{ 50.0000, 2.5000, 0.0000, 50.0000, 3.1736, 0.5854, 1.0000 },
	{ 50.0000, 2.5000, 0.0000, 50.0000, 3.2972, 0.0000, 1.0000 },
	{ 50.0000, 2.5000, 0.0000, 50.0000, 1.8634, 0.5757, 1.0000 },
	{ 50.0000, 2.5000, 0.0000, 50.0000, 3.2592, 0.3350, 1.0000 },
	{ 60.2574, -34.0099, 36.2677, 60.4626, -34.1751, 39.4387, 1.2644 },
	{ 63.0109, -31.0961, -5.8663, 62.8187, -29.7946, -4.0864, 1.2630 },
	{ 61.2901, 3.7196, -5.3901, 61.4292, 2.2480, -4.9620, 1.8731 },
	{ 35.0831, -44.1164, 3.7933, 35.0232, -40.0716, 1.5901, 1.8645 },
	{ 22.7233, 20.0904, -46.6940, 23.0331, 14.9730, -42.5619, 2.0373 },
	{ 36.4612, 47.8580, 18.3852, 36.2715, 50.5065, 21.2231, 1.4146 },
	{ 90.8027, -2.0831, 1.4410, 91.1528, -1.6435, 0.0447, 1.4441 },
	{ 90.9257, -0.5406, -0.9208, 88.6381, -0.8985, -0.7239, 1.5381 },
	{ 6.7747, -0.2908, -2.4247, 5.8714, -0.0985, -2.2286, 0.6377 },
	{ 2.0776, 0.0795, -1.1350, 0.9033, -0.0636, -0.5514, 0.9082 }
};

double icmCIE2K(double *Lab1, double *Lab2);

int main(void) {
	int rv = 0;
	int i;

	printf("Starting Test:\n");

#ifdef NEVER
	for (i = 0; i < NTESTS; i++) {
		printf("Lab1 = %f %f %f\n", ref[i].Lab1[0], ref[i].Lab1[1], ref[i].Lab1[2]);
		printf("Lab2 = %f %f %f\n", ref[i].Lab2[0], ref[i].Lab2[1], ref[i].Lab2[2]);
		printf("de = %f\n",ref[i].de);
	}
#endif

	/* Test it all out */
	for (i = 0; i < NTESTS; i++) {
		double de;

		de = icmCIE2K(ref[i].Lab1, ref[i].Lab2);
		if (fabs(de - ref[i].de) > 0.0001) {
			printf("Error at index %d:\n",i);
			printf("Lab1 = %f %f %f\n", ref[i].Lab1[0], ref[i].Lab1[1], ref[i].Lab1[2]);
			printf("Lab2 = %f %f %f\n", ref[i].Lab2[0], ref[i].Lab2[1], ref[i].Lab2[2]);
			printf("DeltaE is %f, should be %f\n\n",de,ref[i].de);
			rv = 1;
		}
		de = icmCIE2K(ref[i].Lab2, ref[i].Lab1);
		if (fabs(de - ref[i].de) > 0.0001) {
			printf("Error at index %d:\n",i);
			printf("Lab1 = %f %f %f\n", ref[i].Lab2[0], ref[i].Lab2[1], ref[i].Lab2[2]);
			printf("Lab2 = %f %f %f\n", ref[i].Lab1[0], ref[i].Lab1[1], ref[i].Lab1[2]);
			printf("DeltaE is %f, should be %f\n\n",de,ref[i].de);
			rv = 1;
		}
	}

	printf("Test Finished\n");
	return rv;
}

#ifdef NEVER		/* Test implementation in icc.c */

/* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
/* Supplementary Test Data, and Mathematical Observations", by */
/* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
/* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */

/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
double icmCIE2Ksq(double *Lab0, double *Lab1) {
	double C1, C2;
	double h1, h2;
	double dL, dC, dH;
	double dsq;

	/* The trucated value of PI is needed to ensure that the */
	/* test cases pass, as one of them lies on the edge of */
	/* a mathematical discontinuity. The precision is still */
	/* enough for any practical use. */
#define RAD2DEG(xx) (180.0/3.14159265358979 * (xx))
#define DEG2RAD(xx) (3.14159265358979/180.0 * (xx))

	/* Compute Cromanance and Hue angles */
	{
		double C1ab, C2ab;
		double Cab, Cab7, G;
		double a1, a2;

		C1ab = sqrt(Lab0[1] * Lab0[1] + Lab0[2] * Lab0[2]);
		C2ab = sqrt(Lab1[1] * Lab1[1] + Lab1[2] * Lab1[2]);
		Cab = 0.5 * (C1ab + C2ab);
		Cab7 = pow(Cab,7.0);
		G = 0.5 * (1.0 - sqrt(Cab7/(Cab7 + 6103515625.0)));
		a1 = (1.0 + G) * Lab0[1];
		a2 = (1.0 + G) * Lab1[1];
		C1 = sqrt(a1 * a1 + Lab0[2] * Lab0[2]);
		C2 = sqrt(a2 * a2 + Lab1[2] * Lab1[2]);

		if (C1 < 1e-9)
			h1 = 0.0;
		else {
			h1 = RAD2DEG(atan2(Lab0[2], a1));
			if (h1 < 0.0)
				h1 += 360.0;
		}

		if (C2 < 1e-9)
			h2 = 0.0;
		else {
			h2 = RAD2DEG(atan2(Lab1[2], a2));
			if (h2 < 0.0)
				h2 += 360.0;
		}
	}

	/* Compute delta L, C and H */
	{
		double dh;

		dL = Lab1[0] - Lab0[0];
		dC = C2 - C1;
		if (C1 < 1e-9 || C2 < 1e-9) {
			dh = 0.0;
		} else {
			dh = h2 - h1;
			if (dh > 180.0)
				dh -= 360.0;
			else if (dh < -180.0)
				dh += 360.0;
		}

		dH = 2.0 * sqrt(C1 * C2) * sin(DEG2RAD(0.5 * dh));
	}

	{
		double L, C, h, T;
		double hh, ddeg;
		double C7, RC, L50sq, SL, SC, SH, RT;
		double dLsq, dCsq, dHsq, RCH;

		L = 0.5 * (Lab0[0]  + Lab1[0]);
		C = 0.5 * (C1 + C2);
		if (C1 < 1e-9 || C2 < 1e-9) {
			h = h1 + h2;
		} else {
			h = h1 + h2;
			if (fabs(h1 - h2) > 180.0) {
				if (h < 360.0)
					h += 360.0;
				else if (h >= 360.0)
					h -= 360.0;
			}
			h *= 0.5;
		}
		T = 1.0 - 0.17 * cos(DEG2RAD(h-30.0)) + 0.24 * cos(DEG2RAD(2.0 * h))
		  + 0.32 * cos(DEG2RAD(3.0 * h + 6.0)) - 0.2 * cos(DEG2RAD(4.0 * h - 63.0));
		hh = (h - 275.0)/25.0;
		ddeg = 30.0 * exp(-hh * hh);
		C7 = pow(C,7.0);
		RC = 2.0 * sqrt(C7/(C7 + 6103515625.0));
		L50sq = (L - 50.0) * (L - 50.0);
		SL = 1.0 + (0.015 * L50sq)/sqrt(20.0 + L50sq);
		SC = 1.0 + 0.045 * C;
		SH = 1.0 + 0.015 * C * T;
		RT = -sin(DEG2RAD(2 * ddeg)) * RC;

		dLsq = dL/SL;
		dCsq = dC/SC;
		dHsq = dH/SH;

		RCH = RT * dCsq * dHsq;

		dLsq *= dLsq;
		dCsq *= dCsq;
		dHsq *= dHsq;

		dsq = dLsq + dCsq + dHsq + RCH;
	}

	return dsq;

#undef RAD2DEG
#undef DEG2RAD
}

/* Return the CIE2DE000 Delta E color difference measure for two Lab values */
double icmCIE2K(double *Lab0, double *Lab1) {
	return sqrt(icmCIE2Ksq(Lab0, Lab1));
}

#endif /* NEVER */
