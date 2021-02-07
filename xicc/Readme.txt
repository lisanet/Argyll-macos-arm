This directory holds the "expansion" icc libraries.

These suplement the base icc library with enhanced
profile functionality, such as smoothed interpolation,
reverse interpolation, table creation from scattered
data etc.

Most of this functionality is based on the rspl and icc libraries.
This is where ink limiting and black generation policies
for CMYK devices is implemented.

xicclu.exe	is the analog of the icclib utility icclu, but
			expands the capability to reverse lookup the
			Lut tables.

