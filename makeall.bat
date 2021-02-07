echo Simple batch file to invoke Jam from the top
jam -q -fJambase -j%NUMBER_OF_PROCESSORS% %1
rem If you have trouble with the parallel build, try the
rem version with only one thread.
rem jam -q -fJambase %1
