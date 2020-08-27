# FEANBEAM: Finite Element Analysis of Beams

## Description

FEANBEAM is a fortran application that can be used to analyse beams and slender structures.

Current features:

* Geometrically nonlinear beam theory (Simo and Vu-Quoc)

* Linear elastic material model

* Initally straight elements

* Quasi-static analysis

* Newton-Raphson and Arc-Length solution algorithm

## Usage

Application functions are written in modules which need to be compiled only the first time. For now the input is also a fortran file, which needs to be compiled every time it is changed. The goal is to completely separate the input from application. This will be done some time in the future.

For now, there are two possibilities:

1. Input using a script, where you specify everything. This is convinient for simple first attempts.

2. Input using a script and support files, where you specify parameters in a separate file. This is more suitable for analysis that is run multiple times with slighlty tweaked parameters.

Both examples can be seen in *Examples*.

## Requirements

FEANBEAM provides you with source code, which needs to be compiled in order to run. In the folder `Build` you should find some scripts to automate this process written for different combinations of operating systems and compilers. If you have another set up you can probably mesh something together or direct the question to me as it might be useful to the others as well.

## Installation

### Installing on Windows 10

1. Install [GNU Fortran Compiler](https://gcc.gnu.org/wiki/GFortranBinaries) or [Intel Fortran Compiler](https://software.intel.com/content/www/us/en/develop/tools/compilers/fortran-compilers.html).

 * MinGW Installation Options: Under all packages select *mingw32-make-bin* and *mingw32-gcc-fortran-bin* options. Then click *Installation* in the upper-right corner and select *Apply changes*.

 * If you installed MinGW, add it to the `PATH` variables. First go to Windows Settings (Control Panel). In the search box start typing "environment" and select *Edit environment variables for your account*. In the first box select *PATH* and click *EDIT*. Click *New* and then *Browse* and navigate to *MinGW* installation folder and select bin. An example would be `C:\MinGW\bin`. Finally click *OK* and *OK* again.

2. Instal LAPACK

 * Install [Visual Studio](https://visualstudio.microsoft.com/downloads/) (Community edition is free)

 * Install CMake by obtaining the latest binary distribution from: https://cmake.org/download/

 * Add CMake to the `PATH` variables, same as with MinGW. This time, the *bin* folder is probably under `C:\Program Files\CMake\bin`

 * Get LAPACK as source code (.zip) from their [Github page](https://github.com/Reference-LAPACK/lapack/releases)

 * Unzip the archive to `lapack-version` and move the folder to a convinient location, example:  `C:\LAPACK\lapack-version`

 * Create a *build* folder to install LAPACK. A good choice would be `C:\LAPACK\build`.

 * Open CMake application and choose *Source* (unzipped folder `lapack-version`) and *Build* folders (previosuly created).

 * Click *Configure*.

 * Select the installed *MinGW Makefiles*. Click next.

 * Select the Fortran and C compilers in MinGW *bin* folder. Probably at `C:\MinGW\bin\gfortran.exe` and `C:\MinGW\bin\gcc.exe`.

 * Click continue. Some files may be falsely recognised as threats by windows deffender in which case you must click on the notification and under *Actions* click *Allow*. A lot of red options appear. Additionaly check *BUILD_SHARED_LIBS* and *CMAKE_GNUtoMS* and change *CMAKE_INSTALL_PREFIX* to `C:/LAPACK/LAPACK`.

 * Press again configure until everything becomes white.

 * Click *Generate*. This will generate build and after it is finished, you may close CMake.

 * Now open command promt and navigate to the *build* folder with command `cd C:\LAPACK\build`.

* Execute command `mingw32-make`.

* Execute command `mingw32-make test` to test if everything is installed correctly.

### Installing on MacOS

1. Install GNU Fortran compiler.

``` zsh
brew install gcc
```

2. Instal LAPACK

``` zsh
brew install lapack
```



## Contribution

If you like this project and think you can contribute, please do not hesitate to contact me.


