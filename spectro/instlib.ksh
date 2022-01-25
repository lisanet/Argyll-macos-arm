#!/bin/sh

# Package up the source and other files for the standalone GPLv2 Instrument library.

# Copyright 2007 - 2013 Graeme W. Gill
# This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 2 or later :-
# see the License2.txt file for licencing details.

echo "Making standalone GPLv2 instrument archive instlib.zip "

H_FILES="
	../h/sort.h
	"

NUMLIB_FILES="
	../numlib/numsup.h
	../numlib/numsup.c
	"

CGATS_FILES="
	../cgats/pars.h
	../cgats/pars.c
	../cgats/parsstd.c
	../cgats/cgats.h
	../cgats/cgats.c
	../cgats/cgatsstd.c
	"

XICC_FILES="
	../xicc/xspect.h
	../xicc/xspect.c
	../xicc/ccss.h
	../xicc/ccss.c
	../xicc/ccmx.h
	../xicc/ccmx.c
	"

RSPL_FILES="
	../rspl/rspl1.h
	../rspl/rspl1.c
	"

USB_FILES="
	../usb/driver/driver_api.h
	"

SPECTRO_FILES="
	License2.txt
	spotread.c
	Makefile.OSX
	Makefile.UNIX
	Makefile.WNT
	pollem.h
	pollem.c
	conv.h
	conv.c
	sa_conv.h
	sa_conv.c
	aglob.c
	aglob.h
	hidio.h
	hidio.c
	icoms.h
	dev.h
	inst.h
	inst.c
	insttypes.c
	insttypes.h
	insttypeinst.h
	instappsup.c
	instappsup.h
	disptechs.h
	disptechs.c
	dtp20.c
	dtp20.h
	dtp22.c
	dtp22.h
	dtp41.c
	dtp41.h
	dtp51.c
	dtp51.h
	dtp92.c
	dtp92.h
	ss.h
	ss.c
	ss_imp.h
	ss_imp.c
	i1disp.c
	i1disp.h
	i1d3.h
	i1d3.c
	i1pro.h
	i1pro.c
	i1pro_imp.h
	i1pro_imp.c
	i1pro3.h
	i1pro3.c
	i1pro3_imp.h
	i1pro3_imp.c
	munki.h
	munki.c
	munki_imp.h
	munki_imp.c
	hcfr.c
	hcfr.h
	huey.c
	huey.h
	colorhug.c
	colorhug.h
	spyd2.c
	spyd2.h
	spydX.c
	spydX.h
	specbos.h
	specbos.c
	kleink10.h
	kleink10.c
	ex1.c
    ex1.h
    smcube.h
	smcube.c
    cubecal.h
	oemarch.c
	oemarch.h
	oeminst.c
	vinflate.c
	inflate.c
	LzmaDec.c
	LzmaDec.h
	LzmaTypes.h
	icoms.c
	icoms_nt.c
	icoms_ux.c
	iusb.h
	usbio.h
	usbio.c
	usbio_nt.c
	usbio_ox.c
	usbio_lx.c
	rspec.h
	rspec.c
	xdg_bds.c
	xdg_bds.h
	base64.h
	base64.c
	xrga.h
	xrga.c
	"

FILES=" $H_FILES $CGATS_FILES $NUMLIB_FILES $RSPL_FILES $XICC_FILES $USB_FILES $SPECTRO_FILES "

rm -f instlib.zip
rm -rf _zipdir
rm -f _ziplist
mkdir _zipdir
mkdir _zipdir/instlib


# Archive the Argyll files needed
for j in $FILES
do
	if [ ! -e ${j} ] ; then
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!! Can't find file ${j} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		NOTFOUND="$NOTFOUND ${j}"
	else 
		cp ${j} _zipdir/instlib/${j##*/}
		echo instlib/${j##*/} >> _ziplist
	fi
done

# Plus renamed files
cp IntsLib_Readme.txt _zipdir/instlib/Readme.txt
echo instlib/Readme.txt >> _ziplist
cp Jamfile.SA _zipdir/instlib/Jamfile
cp Makefile.SA _zipdir/instlib/Makefile
echo instlib/Jamfile >> _ziplist
echo instlib/Makefile >> _ziplist
cp ../h/aconfig.h _zipdir/instlib/sa_config.h
echo instlib/sa_config.h >> _ziplist

# Create usb archive

for j in `cat ../usb/afiles | grep -E -v 'afiles|binfiles.msw|binfiles.osx|binfiles.lx|Jamfile|ArgyllCMS.inf.t|ArgyllCMS.inf.d'`
do
	echo "File ${j}"

	# Create any needed temporary directories

	tt=usb/${j}
	path=${tt%/*}		# extract path without filename
		
	echo "path ${path}"

	if [ ! -e _zipdir/instlib/${path} ] ; then                     # if not been created
		echo "Creating directory _zipdir/instlib/${path}"
		mkdir -p _zipdir/instlib/${path}
	fi

	tt=../${tt}

	if [ ! -e ${tt} ] ; then
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!! Can't find file ${tt} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		NOTFOUND="$NOTFOUND ${tt}"
	else 
		cp ${tt} _zipdir/instlib/usb/${j}
		echo instlib/usb/${j} >> _ziplist
	fi
done

cd _zipdir
zip -9 -m ../instlib.zip `cat ../_ziplist`
cd ..
rm -rf _zipdir
rm -f _ziplist

if [ "X$NOTFOUND" != "X" ] ; then
	echo "!!!!!! Didn't find $NOTFOUND !!!!!!"
fi

echo "Created instlib.zip"
