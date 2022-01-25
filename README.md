# Argyll CMS for macOS with Apple Silicon cpu   #
This is a friendly fork of [ArgyllCMS](https://www.argyllcms.com) with patches to run on the new Apple Silicon cpu with the M1 chip set.

The original distribution cannot detect the display and therefor sadly crashes on the new M1 chip. This repo includes some small patches to
get the display detection working again and to smoothly compile as arm64 binaries. 

## Installation ##

### Binary installation, e.g. to use with DisplayCAL ###
The easiest way is to download the binary package and unzip the archiv to a directory of your choice. 

The binary package is not signed with a developer certificate so you need to remove the quarantine bit to allow the system to run them. 
Just open Termonal.app and go to the directory where you've unzipped the archiv and type the following commands. The directory is named with the version number e.g. Argyll_V2.3.0 for Argyll version 2.3.0. You need to replace the version number in the following commands with the one you use.

```
cd Argyll_V2.3.0
xattr -dr com.apple.quarantine *
```

To use this with [displayCAL](https://displaycal.net), open displayCAL and choose from the menu 'Locate ArgyllCMS executables ...' and point to this directory.

If you want to use argyll on the command line, make sure to put the bin directory into your PATH environment by doing
```
export PATH=/our-custom-path/Argyll_V2.3.0/bin/:$PATH
```

### Install from source code ###

If you want to compile it from source code, you need the JAM build tool. You can get it from one of my other repos here .....

After building and installing JAM, make sure you have set the following environment variables:

```
export HOSTTYPE=arm64
export MACHTYPE=arm64
```

The easiest way to build Argyll is to clone the repo, then go to the source directory and type
```
./makepackagebin.sh
```

This will build a binary package which you can install as described above. If you want, for some reason, to build without packaging, then after cloning, type 
```
./makeall.sh
./makeinstall.sh
```

This way the binaries will be located in the bin subdirectory. 

