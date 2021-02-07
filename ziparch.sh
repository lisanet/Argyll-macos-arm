#!/bin/sh

# Make complete source distribution archive.

echo "Making Complete Argyll source archive argyll.zip... "

rm -f argyll.zip
rm -rf _zipdir
mkdir _zipdir
NOTFOUND=

# Split on lines, not spaces
OIFS="$IFS"
IFS='
'
for i in `cat adirs bdirs`
do
	echo
	echo "#### Doing Directory $i ####"
	if [ ! -e ${i}/afiles ] ; then
		if [ ! -e ${i}/bfiles ] ; then
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Can't find ${i}/afiles or ${i}/bfiles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			NOTFOUND="$NOTFOUND ${i}/afiles ${i}/bfiles"
		fi
	fi

	if [ -e ${i}/afiles ] ; then
		rm -f _ziplist

		for j in `cat $i/afiles`
		do
			# Create any needed temporary directories
			tt=${i}/${j}
			path=${tt%/*}		# extract path without filename
			
			if ! expr _zipdir/${path} : '\b\.\b' > /dev/null  ; then   # if not "."
				if [ ! -e _zipdir/${path} ] ; then                     # if not been created
					mkdir -p _zipdir/${path}
				fi
			fi

			if [ ! -e "${i}/${j}" ] ; then
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Can't find file ${i}/${j} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				NOTFOUND=$NOTFOUND ${i}/${j}
			else 
				dos2unix ${i}/${j}
				cp ${i}/${j} _zipdir/${i}/${j}
				echo ${i}/${j} >> _ziplist
			fi
		done

		cd _zipdir
		zip -9 -m ../argyll.zip `cat ../_ziplist`
		cd ..
		#if ! expr ${i} : '\b\.\b' > /dev/null ; then
		if ! expr ${i} : '\.' > /dev/null ; then
			rm -r _zipdir/${i}
		fi
	fi

	# same as above, but for "bfiles", if it exists
	if [ -e ${i}/bfiles ] ; then
		rm -f _ziplist

		for j in `cat $i/bfiles`
		do

			# Create any needed temporary directories
			tt=${i}/${j}
			path=${tt%/*}		# extract path without filename
			
			if ! expr _zipdir/${path} : '\b\.\b' > /dev/null  ; then   # if not "."
				if [ ! -e _zipdir/${path} ] ; then                     # if not been created
					mkdir -p _zipdir/${path}
				fi
			fi

			if [ ! -e ${i}/${j} ] ; then
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Can't find file ${i}/${j} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				NOTFOUND=$NOTFOUND ${i}/${j}
			else 
				dos2unix ${i}/${j}
				cp ${i}/${j} _zipdir/${i}/${j}
				echo ${i}/${j} >> _ziplist
			fi
		done

		cd _zipdir
		zip -9 -m ../argyll.zip `cat ../_ziplist`
		cd ..
		#if ! expr ${i} : '\b\.\b' > /dev/null ; then
		if ! expr ${i} : '\.' > /dev/null ; then
			rm -r _zipdir/${i}
		fi
	fi
done
rm -r _zipdir
rm _ziplist
if [ "X$NOTFOUND" != "X" ] ; then
    echo "!!!!!! Didn't find $NOTFOUND !!!!!!"
fi
echo "Finished Complete Argyll source archive argyll.zip... "

IFS="$OIFS"
