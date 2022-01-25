#!/bin/sh
echo "Script to invoke Jam and then package the binary release."

# Must use this rather than "jam -q" to ensure builtin libraries are used.

# Set the environment string VERSION from the #define, ie 1.0.0
VERSION=`grep ARGYLL_VERSION_STR h/aconfig.h | head -1 | sed 's/# define ARGYLL_VERSION_STR //' | sed 's/"//g'`

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
#   OS X i386 10.11 [bash]          darwin16     x86_64-apple-darwin16    x86_64
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

echo "About to make Argyll binary distribution $VERSION"

TOPDIR=Argyll_V$VERSION

if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.sh
fi

# Make sure that some environment variable are visible to Jam:
export OSTYPE MACHTYPE HOSTTYPE
unset USETARPREFIX

# .sp come from profile, .cht from scanin and .ti3 from spectro
rm -f bin/*.exe bin/*.dll
rm -f ref/*.sp ref/*.cht ref/*.ti2

# Make sure it's built and installed
if ! jam -q -fJambase -j${NUMBER_OF_PROCESSORS:-2} -sBUILTIN_TIFF=true -sBUILTIN_JPEG=true -sBUILTIN_PNG=true -sBUILTIN_Z=true -sBUILTIN_SSL=true install ; then
	echo "Build failed!"
	exit 1
fi 

# Maybe we could get Jam to do the following ?

if [ X$OS = "XWindows_NT" ] ; then
	echo "We're on MSWindows!"
	# Hack cross comile
	if [ X$COMPILER = "XMINGW64" ] ; then
		echo "We're cross compiling to MSWin 64 bit !"
		PACKAGE=Argyll_V${VERSION}_win64_exe.zip
		USBDIRS="usb"
		USBBINFILES="binfiles.msw"
		unset USETAR
	else
		# ~~ should detect native 64 bit here ~~
		echo "We're on MSWin 32 bit !"
		PACKAGE=Argyll_V${VERSION}_win32_exe.zip
		USBDIRS="usb"
		USBBINFILES="binfiles.msw"
		unset USETAR
	fi
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
	echo "We're on OSX 10.3 PPC!"
	PACKAGE=Argyll_V${VERSION}_osx10.3_ppc_bin.tgz
	USBDIRS="usb"
	USBBINFILES="binfiles.osx"
	USETAR=true
else if [ X$OSTYPE = "Xdarwin8.0" ] ; then
	if [ X$MACHTYPE = "Xi386-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 i386!"
		PACKAGE=Argyll_V${VERSION}_osx10.4_i86_bin.tgz
	else if [ X$MACHTYPE = "Xpowerpc-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 PPC!"
		PACKAGE=Argyll_V${VERSION}_osx10.4_ppc_bin.tgz
	fi
	fi
	USBDIRS="usb"
	USBBINFILES="binfiles.osx"
	USETAR=true
else if [ X$OSTYPE = "Xdarwin10.0" \
       -o X$OSTYPE = "Xdarwin11"   \
       -o X$OSTYPE = "Xdarwin12"   \
       -o X$OSTYPE = "Xdarwin13"   \
       -o X$OSTYPE = "Xdarwin14"   \
       -o X$OSTYPE = "Xdarwin15"   \
       -o X$OSTYPE = "Xdarwin16"   \
       -o X$OSTYPE = "Xdarwin17"   \
       -o X$OSTYPE = "Xdarwin18"   \
       -o X$OSTYPE = "Xdarwin19"   \
       -o X$OSTYPE = "Xdarwin20" ] ; then
	if [ X$HOSTTYPE = "Xx86_64" ] ; then
		echo "We're on OSX 10.X x86_64!"
		PACKAGE=Argyll_V${VERSION}_osx10.6_x86_64_bin.tgz
	else if [ X$HOSTTYPE = "Xarm64" ] ; then
		echo "We're on OSX 11.X arm64!"
		PACKAGE=Argyll_V${VERSION}_macOS11_arm64_bin.tgz
	fi
	fi
	USBDIRS="usb"
	USBBINFILES="binfiles.osx"
	USETAR=true
	USETARPREFIX=true
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
	if [[ "$MACHTYPE" = x86_64-*-linux-gnu ]] ; then
		echo "We're on Linux x86_64!"
		PACKAGE=Argyll_V${VERSION}_linux_x86_64_bin.tgz
	else if [[ "$MACHTYPE" = *86-*-linux-gnu ]] ; then
		echo "We're on Linux x86!"
		PACKAGE=Argyll_V${VERSION}_linux_x86_bin.tgz
	fi
	fi
	USBDIRS="usb"
	USBBINFILES="binfiles.lx"
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

echo "Making GNU Argyll binary distribution $PACKAGE for Version $VERSION"

rm -rf $TOPDIR
mkdir $TOPDIR

# Collect the names of all the files that we're going to package
unset topfiles; for i in `cat binfiles`; do topfiles="$topfiles ${i}"; done
unset docfiles; for i in `cat doc/afiles`; do docfiles="$docfiles doc/${i}"; done
unset usbfiles;
for j in ${USBDIRS}; do
	if [ ${j} ]; then
		for i in `cat ${j}/${USBBINFILES}`; do usbfiles="$usbfiles ${j}/${i}"; done
	fi
done

allfiles="${topfiles} bin/* ref/* ${docfiles} ${usbfiles}"

# Copy all the files to the package top directory
for i in ${allfiles}; do
	path=${i%/*}		# extract path without filename
	file=${i##*/}		# extract filename
	if [ $path = $i ] ; then
		path=
	fi
	if [ X$path != "X" ] ; then
		mkdir -p $TOPDIR/${path}
	fi
	cp $i $TOPDIR/$i
done

# Create the package
rm -f $PACKAGE
if [ X$USETAR = "Xtrue" ] ; then
	if [ X$USETARPREFIX = "Xtrue" ] ; then
		# Don't save ._* files...
		COPYFILE_DISABLE=1 tar -czvf $PACKAGE $TOPDIR
	else
		tar -czvf $PACKAGE $TOPDIR
	fi
	# tar -xzf to extract
	# tar -tzf to list
	# to update a file:
	#  gzip -d archive.tgz
	#  tar -uf application.tar file
	#  gzip application.tar
	#  mv application.tar.gz application.tgz
	# or "tgzupdate.sh application fullpath/file1 fullpath/file2 fullpath/file3"
	# Should we use "COPYFILE_DISABLE=1 tar .." on OS X ??
else
	zip -9 -r $PACKAGE $TOPDIR
	# unzip to extract
	# unzip -l to list
	# zip archive.zip path/file to update
fi
rm -rf $TOPDIR
echo "Done GNU Argyll binary distribution $PACKAGE"

exit 0

