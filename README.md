Description
===========

flattorus is a program that generates three corrugations on a standard
torus embedding. See "Flat tori in three dimensional space and convex integration". PNAS 109:7218-7223, 2012 by Borrelli, V.; Jabrane, S.; Lazarus, F. & Thibert, B. for more details. 
The program can ouput the resulting corrugated torus as a 3D mesh in VRML, VTK, POVRAY or OFF formats. The VRML format can be used for 3D printing.

The archive contains a directory named "hevea" and
- a subdirectory `src/` with the main file flattorus.cpp 
- a subdirectory `bin/` with the Makefile and 
- 8 subdirectories
`argstream/` `cylindre/`  `integration/`  `iso/`  `local/`  `matrices/`  `tore/`  `utils/`
containing header files.

Installation on UNIX like architectures 
==========================
(The same source files should work on other platforms with an appropriate compilation command)

Prerequisites:
* g++ >= 2.4.0 (check with g++ --version)
* blitz++ >= 0.9 (http://sourceforge.net/projects/blitz/. There is a Ubuntu package)
The standard blitz installation directory is `/usr/local`. You can use any other directory if you prefer. 
* openmp if you want to benefit from parallel programming on multi-core CPUs

1) Unpack the source archive `flat_torus.tgz` in the directory of your choice. 
2) cd into this directory.
3) cd bin.
4) to compile flattorus, run:
make blitzprefix=[your prefered Blitz++ prefix directory. Default is /usr/local.]
you can also disable openMP by setting 
`OPENMPMODE=no`
or use gcc debug compilation by setting
`DEBUGMODE=yes`

Examples:
* If blitz-09 is installed in the standard directory `/usr/local` and you want to compile flattorus with openMP, you should just type:
```bash
$ make
```

* If you installed the 0.10 version of blitz++ in the hevea directory and want to compile in debug mode without openMP, you should type
```bash
$ make blitzprefix=$PWD/../blitz-0.10 OPENMPMODE=no DEBUGMODE=yes
```



RUN flattorus
=============
To run flattorus from the hevea directory just type:
```bash
$ ./bin/flattorus 
```

This will apply 3 corrugations with 12, 80 and 500 oscillations to the standard torus of revolution with radii r=0.2 and R=0.5. The computation is performed on a square grid with 10,000\*10,000 nodes. The program outputs after each corrugation
- a subsampled embedding of size 1,000\*1,000 in VTK format for easy handling
- a local piece of size 500\*500 in VTK format 
- a VTK file of size 10,000\*10,000 
- a subsampled VRML file of size 3500\*3500 that can be used for high quality 3D printing

If you want to use a 3D printer for the third corrugated torus it seems accurate to  
use the default grid size (10,000\*10,1000) for the computations and an  
undersampled grid of size 3,500\*3,500 for the VRML ouput (larger grids may 
not be accepted by 3D printers). 

Several parameters can be modified by the user. Type
```bash
$ ./bin/flattorus --help
```
to list all the possible arguments. 


