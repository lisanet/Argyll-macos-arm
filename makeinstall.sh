#!/bin/sh
echo "Script to invoke Jam from the top and install"
echo "binaries in the ./bin directory and samples into the ref ./directory"

if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.sh
fi
jam -q -fJambase -j${NUMBER_OF_PROCESSORS:-2} install
