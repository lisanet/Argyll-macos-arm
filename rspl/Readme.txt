!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rspl now supports different resolution grids in each dimension.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This is the second generation Regular Spline library.

It contains scattered data point to regular grid interpolation,
as well as spline smoothing, and the reverse interpolation
code. This version is more modular, and uses better solution
algorithms than the earlier REGSPL, and generally replaces it.

The reverse interpolation algorithms support features needed
for devices like CMYK printers, such as total ink limiting,
black locus selection, gamut boundary detection, vector
and nearest gamut clipping.

It has been written with operation with 6 color printing
devices in mind (ie. 3 extra degrees of freedom, and hence
a 3 dimensional inking locus), although this usage is likely
to be unwealdy.


Misc test files:

c1.c		Test 1D curves. First test is to check tracking at multiple resolutions 
c1df.c		Test 1D curve with weak default function.
c1i.c		Test 1D curve with incremental points
sm1.c		Discover 1D smoothness factor vs resolution tracking factor
sm2.c		Discover 2D smoothness factor vs resolution tracking factor
sm3.c		Discover 3D smoothness factor vs resolution tracking factor
t2d.c		Test 2D fitting. Test againt two resolutions.
t2ddf.c		Test 2D fitting with weak default function.
t3d.c		Test 3D fitting. Test againt two resolutions.
t3ddf.c		Test 3D fitting with weak default function.
tnd.c		Simple test of 4D
trnd.c		Simple test of 4D reverse lookup.
smtnd.c		Sythetic function, nD multi-parameter fitting test function.
smtmpp.c	MPP function, nD multi-parameter fitting test function.

