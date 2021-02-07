This directory contains the gamut boundary
finding, and gamut mapping code.
It also contains utilities for generating
3D visualisations of gamut boundaries.

	iccgamut:
		Take an icc profile, and derive a gamut surface from one
		of the forward or reverse transforms, using the desired
		intent. Output the result as a CGATS format .gam file,
		and optionaly as a VRML model file.

	tiffgamut:
		Take an tiff file and an icc profile, and derive a gamut
		surface from the forward absolute profile and the tiff
		image.  Output the result as a CGATS format .gam file,
		and optionaly as a VRML model file.

	viewgam:

		This program reads one or more CGATS format triangular gamut
		surface descriptions (.gam), and combines them into a VRML file,
		so that the gamuts can be visually compared.
		Allowance is made for representing the gamuts as solid or
		tranparent surfaces, wire frames, and coloring the gamut
		with natural or uniform colors.
