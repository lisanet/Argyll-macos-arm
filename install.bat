echo Simple batch file to invoke Jam install from the top
jam -q -fJambase -j%NUMBER_OF_PROCESSORS% install
rem If you have trouble with the parallel build, try the
rem version with only one thread.
rem jam -q -fJambase install
