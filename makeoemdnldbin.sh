#!/bin/sh
echo "Script to invoke Jam and then package spectro/oemdnld binary."

#   Typical environment variables:
#   (NOTE some systems don't export these ENV vars. by default !!!)
#
#   Platform                        $OSTYPE      $MACHTYPE                $HOSTTYPE
#
#   Win2K [CMD.EXE]                 (none)       (none)                   (none)        
#
#   Cygwin Win2K [bash]             cygwin       i686-pc-cygwin           i686
#
#   OS X PPC 10.3 [zsh]             darwin7.0    powerpc                  (none)
#
#   OS X i386 10.4 [bash]           darwin8.0    i386-apple-darwin8.0     i386
#
#   OS X i386 10.5 [bash]           darwin9.0    i386-apple-darwin9.0     i386
#
#   OS X i386 10.6 [bash]           darwin10.0   x86_64-apple-darwin10.0  x86_64
#
#   OS X i386 10.7 [bash]           darwin11     x86_64-apple-darwin11    x86_64
#
#   Linux RH 4.0 [bash]             linux-gnu    i686-redhat-linux-gnu    i686
#
#   Linux Fedora 7.1 [bash]         linux-gnu    i386-redhat-linux-gnu    i386
#   Linux Ubuntu  ??7               linux-gnu    i486-pc-linux-gnu        i686
#
#   Linux Fedora 7.1 64 bit [bash]  linux-gnu    x86_64-redhat-linux-gnu  x86_64
#   Ubuntu 12.10 64 bit [bash]      linux-gnu    x86_64-pc-linux-gnu      x86_64
#
#   FreeBSD 9.1 64 bit [bash]       freebsd9.1   amd64-portbld-freebsd9.1 amd64
#

# We don't mark the binaries with the version, so that we don'y
# haveto keep the web page in sync with ArgyllCMS.

echo "About to make OEMdnld binary distribution"

TOPDIR=OEMdist

# Make sure that some environment variable are visible to Jam:
export OSTYPE MACHTYPE HOSTTYPE

# Make sure it's built and installed
if ! jam -q -fJambase -j${NUMBER_OF_PROCESSORS:-2} -sBUILTIN_TIFF=true -sBUILTIN_JPEG=true "<spectro>oemdnld" ; then
	echo "Build failed!"
	exit 1
fi 

# Maybe we could get Jam to do the following ?

if [ X$OS = "XWindows_NT" ] ; then
	echo "We're on MSWindows!"
	# Hack cross comile
	if [ X$COMPILER = "XMINGW64" ] ; then
		echo "We're cross compiling to MSWin 64 bit !"
		PACKAGE=oemdnld_win64_exe.zip
		EXE=.exe
		unset USETAR
	else
		# ~~ should detect native 64 bit here ~~
		echo "We're on MSWin 32 bit !"
		PACKAGE=oemdnld_win32_exe.zip
		EXE=.exe
		unset USETAR
	fi
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
	echo "We're on OSX 10.3 PPC!"
	PACKAGE=oemdnld_osx10.3_ppc_bin.tgz
	USETAR=true
	EXE=
else if [ X$OSTYPE = "Xdarwin8.0" ] ; then
	if [ X$MACHTYPE = "Xi386-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 i386!"
		PACKAGE=oemdnld_osx10.4_i86_bin.tgz
	else if [ X$MACHTYPE = "Xpowerpc-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 PPC!"
		PACKAGE=oemdnld_osx10.4_ppc_bin.tgz
	fi
	fi
	EXE=
	USETAR=true
else if [ X$OSTYPE = "Xdarwin10.0" \
       -o X$OSTYPE = "Xdarwin11" ] ; then
	if [ X$HOSTTYPE = "Xx86_64" ] ; then
		echo "We're on OSX 10.6 x86_64!"
		PACKAGE=oemdnld_osx10.6_x86_64_bin.tgz
	fi
	EXE=
	USETAR=true
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
	if [[ "$MACHTYPE" = x86_64-*-linux-gnu ]] ; then
		echo "We're on Linux x86_64!"
		PACKAGE=oemdnld_linux_x86_64_bin.tgz
	else if [[ "$MACHTYPE" = *86-*-linux-gnu ]] ; then
		echo "We're on Linux x86!"
		PACKAGE=oemdnld_linux_x86_bin.tgz
	fi
	fi
	EXE=
	USETAR=true
fi
fi
fi
fi
fi

if [ X$PACKAGE = "X" ] ; then
	echo "Unknown host - build failed!"
	exit 1
fi 

echo "Making OEMdnld binary distribution $PACKAGE"

rm -rf $TOPDIR
mkdir $TOPDIR

# Collect the names of all the files that we're going to package

allfiles="spectro/oemdnld${EXE}"

# Copy all the files to the package top directory
for i in ${allfiles}; do
	path=${i%/*}		# extract path without filename
	file=${i##*/}		# extract filename
	if [ $path = $i ] ; then
		path=
	fi
	cp $i $TOPDIR/${file}
done

# Create the package
rm -f $PACKAGE
if [ X$USETAR = "Xtrue" ] ; then
	cd $TOPDIR
	tar -czvf ../$PACKAGE *
	# tar -xzf to extract
	# tar -tzf to list
	cd ..
else
	cd $TOPDIR
	zip -9 -r ../$PACKAGE *
	# unzip to extract
	# unzip -l to list
	cd ..
fi
rm -rf $TOPDIR
echo "Done GNU Argyll binary distribution $PACKAGE"

exit 0

