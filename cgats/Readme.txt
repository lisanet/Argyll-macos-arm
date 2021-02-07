
CGATS file I/O library, V2.01 README file
-----------------------------------------

Package contents:
-----------------
cgatslib.zip   ZIP archive of the following files
Readme.txt     This file.
License4.txt   Important! - Permissions for use of this package.
cgats.c        CGATS Library source code.
cgatsstd.c     I/O and malloc source code.
cgats.h        CGATS Library include file. Note machine dependent defines.
pars.c         Parser source code.
parsstd.c      I/O and malloc source code.
pars.h         Parser include file. Note machine dependent defines.
Jampfile       JAM style "makefile" see <http://www.perforce.com/jam/jam.html>
Makefile       Makefile. Modify this to include one of the following rule sets.
Makefile.WNT   Makefile defines for Microsoft C++ on Windows NT
Makefile.IBMNT Makefile defines for IBM C++ on Windows NT
Makefile.UNIX  Makefile defines for generic UNIX system
Makefile.OSX   Makefile defines for Apple MAC OSX

Changes:
--------

Changes since V2.00

	Removed all exit()s from code - now return error values
	from all functions.

	Added abstract objects for File I/O and memory allocation
	to improve system compatibility.

	Separated the implimentations of the abstract I/O and
	memory objects that ise stdio and malloc into separate
	files so that a library can be compiled without reference
	to these system calls.


---------------------------------------------------------------
                DOCUMENTATION

	This library has been implemented from the CGATS.5 Data Exchange Format
	specification, in Annex J, of the ANSI CGATS.5-1993 standard.
	See <http://www.npes.org/standards/description.htm>
	
    This module attempts to make reading CGATS.5 and IT8.7 files
    easy and convenient. Since the standard is a less than clear
    on some points it is hard to know how much compatibility is
    to be expected.

    The module supports non-standard keywords and fields automatically.
    It does not support reading of comments.
    It supports reading and writing multiple tables within one file.
	Abreviated tables may be written.
	Non-standard fields written by this module will be recognized
	correcty when read by this module (data types are written
	unambiguously - reals all have a decimal point, no non-quoted strings
	are written that could be interpreted as a real or integer), but
	there is no certainty that non-standard fields written by other
	software will be recognized corectly (e.g. - how to you tell
	whether 1234 is an integer, real, or non-quoted string ?)

    To creat a file each element needs to be built up in turn:

    Create an empty cgats structure:
        new_cgats().

    A non-standard memory allocator can be specified for use
    by the cgats object by passing an object that inherits
    from the cgatsAlloc class defined in parse.h, to the
    new_cgats_al(cgatsAlloc *al) constructor.


    Use the module methods to:


	Add a user defined file identifier to augment the standard identifiers:
        add_other(cgats *p, char *osym)
	This can be used to read or write table that are compatible with CGATS.5
	syntax but have different file identifiers.
	Use a zero length string (ie. just "") for wildcard.
	Normaly returns 0, returns -ve if there was a system error,
	& sets p->errc & p->err appropriately.


    create an empty table entry:
        add_table(cgats *p, table_type tt, int oi)
    tt is the table type: it8_7_1, it8_7_2, it8_7_3, it8_7_4, cgats_5, cgats_X,
	tt_other or tt_none.
	cgat_X represents any "CGATS.XXX" table type.
	If tt_other is used for creating a file, then oi must be set to the index of the
	user defined table identifiers set by calls to p->add_other();
	The table type read can be identified by looking at the index
	p->t[table_number].tt, which will be set to the table type.
	If the type is cgats_X type, the type actually found will be in p->cgats_type.
	If tt == tt_other, then the "other index" value is p->t[table_numberx].oi.
	For a wildcard tt_other, the actual type found will be at p->others[oi].

    Tables are added in turn, and are indexed from 0 in other functions.
	Returns 0 normally.
	It returns -2 if there was a system error, & sets p->errc & p->err appropriately.


	Suppress the writing of the file identifier string, standard keyword definitions
	and/or field definitions for a second or subsequent table:
		set_table_flags(cgats *p, int table, int sup_id, int sup_kwords, int sup_fields);
	This only makes sense if the subsequent table has the same identifier type, keywords
	and field definitions as the preceding table. Returns 0 normally.
	It returns -1 if there was an error, & sets p->errc & p->err appropriately.


    Add keywords and their values:
        int add_kword(cgats *p, int table, char *ksym, char *kdata, char *kcom);
	The return value is the index of the new keyword, or -1, errc & err on error.
    Any non-standard keywords will automatically be declared in the file.
	The comment is optional and NULL should be passed if no comment is to be used
	after the keyword/value pair.

    Standard keywords are:
        ORIGINATOR          System, organization or individual that created data
        DESCRIPTOR          Purpose or contents of the data file
        CREATED             Date of creation of the file
        MANUFACTURER        Manufacturer of physical target
        PROD_DATE           Year and month of physical target production yyyy:mm
        SERIAL              Unique physical target serial number
        MATERIAL            Material physical target was produced on
        INSTRUMENTATION     Manufacturer and model number of measuring instrument
        MEASUREMENT_SOURCE  Illumination used for spectral measurements
        PRINT_CONDITIONS    Characteristics of printed sheet being reported

    Standard keywords will be created automatically if necessary for a legal file format.

    The following kewords are supplied automatically by the module,
    and cannot be used for other things:    
        NUMBER_OF_FIELDS
        BEGIN_DATA_FORMAT
        END_DATA_FORMAT
        NUMBER_OF_SETS
        BEGIN_DATA
        END_DATA
        KEYWORD


    Add fields:
        int add_field(cgats *p, int table, char *fsym, data_type ftype);

	The return value is the index of the new field, or return
	-1, errc & err on error, -2, errc & err on system error.

    ftype defines the data type from: r_t, i_t, cs_t, nqcs_t.
        r_t is the real (double) type,
        i_t is the integer (int) type,
        cs_t is the character string (char*) type,
        nqcs_t is the same as cs_t except that it will be non-quoted if possible.
    Note that the type must agree with the standard type if the field is
    from the set default data format identifiers:
    
        SAMPLE_ID    nqcs_t   Identifies sample which data represents
        STRING       cs_t     Identifies label, or other non-machine readable value.
        CMYK_C       r_t      Cyan percentage of CMYK
        CMYK_M       r_t      Magenta percentage of CMYK
        CMYK_Y       r_t      Yellow percentage of CMYK
        CMYK_K       r_t      Black percentage of CMYK
        D_RED        r_t      Red filter reflection density
        D_GREEN      r_t      Green filter reflection density
        D_BLUE       r_t      Blue filter reflection density
        D_VIS        r_t      Visual filter reflection density
        RGB_R        r_t      Red component of RGB data
        RGB_G        r_t      Green component of RGB data
        RGB_B        r_t      Blue component of RGB data
        SPECTRAL_NM  r_t      Wavelength of measurement in nanometers
        SPECTRAL_PCT r_t      Precentage reflectance/transmittance
        XYZ_X        r_t      X component of tristimulus data
        XYZ_Y        r_t      Y component of tristimulus data
        XYZ_Z        r_t      Z component of tristimulus data
        XYY_X        r_t      x component of chromaticity data
        XYY_Y        r_t      y component of chromaticity data
        XYY_CAPY     r_t      Y component of chromaticity data
        LAB_L        r_t      L* component of Lab data
        LAB_A        r_t      a* component of Lab data
        LAB_B        r_t      b* component of Lab data
        LAB_C        r_t      C*ab component of Lab data
        LAB_H        r_t      hab component of Lab data
        LAB_DE       r_t      CIA delta E
        STDEV_X      r_t      Standard deviation of X (tristimulous data)
        STDEV_Y      r_t      Standard deviation of Y (tristimulous data)
        STDEV_Z      r_t      Standard deviation of Z (tristimulous data)
        STDEV_L      r_t      Standard deviation of L*
        STDEV_A      r_t      Standard deviation of a*
        STDEV_B      r_t      Standard deviation of b*
        STDEV_DE     r_t      Standard deviation of CIE delta E

    Add a set of data:
        add_set(cgats *p, int table, ...)
    The data should be supplied as a varargs list in the appropropriate
    data format [char*, double or int].
	Returns 0 normally, -1 errc & err if parameter error,
	-2 errc & err if system error.

    Add a set of data from union array:
		add_setarr(cgats *p, int table, cgats_set_elem *args);
    The data should be supplied as an array of cgats_set_elem unions.
	Returns 0 normally, -1 errc & err if parameter error,
	-2 errc & err if system error.

    Write the data out to a file.
        write_name(cgats *p, char *fname);
    The method will return non-zero on an error, with an error
    description in the err location of the structure.


    A non-standard destination of data can be read by passing an
    object that inherits from the cgatsFile class defined in
    parse.h to the write(cgats *p, cgatsFile *fp) method.

    To read in a data file, the cgats structure should be created
    as usual. The read method can then be called to read in the file:
        read_name(cgats *p, char *fname)
    Returns 0 normally, and -ve on an error, with an error
    description in the err location of the structure,
	and errc set with the return code.

    A non-standard source of data can be read by passing an
    object that inherits from the cgatsFile class defined in
    parse.h to the read(cgats *p, cgatsFile *fp) method.

	The reader will deal automaticaly with carry over of keywords
	and/or field definitions from one table to another, making each
	table appear independent once read.

	The data is accessed by refering to the following read-only
	structure entries:

    The number of tables will be in p->ntables
    
    The number of keywords will be in p->t[table_number].nkwords
    The number of fields will be in p->t[table_number].nfields
    The number of sets will be in p->t[table_number].nsets

	Tables, keywords, fields and sets index from 0.

    The keywords will be in
        p->t[table_number].ksym[keyword_index]

    The keywords character string value be in
        p->t[table_number].kdata[keyword_index]

    The field format identifiers of each field will be in
        p->t[table_number].fsym[field_index]

    The data type of each field will be in
        p->t[table_number].ftype[field_index]

    A void pointer to the data of each field of each set will be in
        p->t[table_number].fdata[set_index][field_index]
	Cast the void pointet according to its type to retrieve the data.
	Alternatively the
	p->get_setarr(struct _cgats *p, int table, int set_index, cgats_set_elem *ary)
	method can be used to fill in a suitable sized cgats_set_elem array with
	the value of all the fields at a particular index. Any character string
	type will be a pointer to the data in p->t[table_number].fdata[set_index][field_index]. 

    To find the index to a particular keyword, use:
        find_kword(cgats *p, int table, char *ksym)
    -1 will be returned if no match is found.
    -2 will be returned, p->errc & p->err will be set if table is out of range.

    To find the index to a particular field, use:
        find_field(cgats *p, int table, char *fsym);
    -1 will be returned if no match is found.
    -2 will be returned, p->errc & p->err will be set if table is out of range.

	Rather than checking the error return codes from every method,
    the first error is "sticky", and recorded in the object.
	This can be checked after a series of operations by calling
	the error(cgats *p, char **mes) method, which will return the
	error code, and (optionaly) the error message.

    Once operations are finished, the object can be deleted by calling
    the delete method:
        del(cgats *p)

Graeme Gill.

---------------------------------------------------------------
