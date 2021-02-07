# default Makefile to invoke Jam

all::
	jam -q -fJambase -j 3

clean::
	jam -fJambase clean

install::
	jam -q -fJambase -j 3 install
