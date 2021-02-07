# Argyll CMS for macOS with Apple Silicon cpu   #
This is a friendly fork of [ArgyllCMS][https://www.argyllcms.com] with patches to run on the new Apple Silicon cpu with the M1 chip set.

The original distribution cannot detect the display and therefor sadly crashes on the new M1 chip. This repo includes some small patches to
get the display detection working again and to smoothly compile as arm64 binaries. 

## Installation ##

### Binary installation, e.g. to use with DisplayCAL ###
The easiest way is to download the binary package and unzip the archiv to a directory of your choice. 

To use this with [displayCAL](https://displaycal.net), open displayCAL and choose from the menu 'Locate ArgyllCMS executables ...' and point to this directory.

If you want to use argyll on the command line, make sure to put the bin directory into your PATH environment by doing

export PATH=/our-custom-path/Argyll_V2.1.2/bin/:$PATH


### Install from source code ###

If you want to compile it from source code, you need the JAM build tool. You can get it form one of my other repos here .....

After building and installing JAM, make sure you have set the following environment variables:

export HOSTTYPE=arm64
export MACHTYPE=arm64

Then clone the repo, go to the source directory and type

./makeall.sh
./makeinstall.sh

The binaries will be located in the bin subdirectory. To make it easier to install argyll, you can build a binary package with 

./makepackagebin.sh





