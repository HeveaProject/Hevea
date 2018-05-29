/***************************************************************************
                           NonStiffIntegratorT.h
                             -------------------
    written by:                  : Blake Ashby
    last modified                : Nov 15, 2002
    email                        : bmashby@stanford.edu

This code computes the numerical solution of a system of first order ordinary
differential equations y'=f(x,y). It uses an explicit Runge-Kutta method of
order (4)5 due to Dormand & Prince with step size control and dense output.

This code is written in C++ and is a modification of the code written originally
in Fortran (Version of April 28, 1994) by:

	E. Hairer & G. Wanner
	Universite de Geneve, dept. de Mathematiques
	CH-1211 GENEVE 4, SWITZERLAND
	E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

and adapted for C by:
	J.Colinge (COLINGE@DIVSUN.UNIGE.CH).

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).


USER PROVIDED FUNCTION:
----------------------

Function	function defining the differential equation. It must have the following
			prototype:

			void Function(NB_TYPE x, NB_TYPE *y, NB_TYPE *f);

			where x is the value of the independent variable x at which you want y
			integrated. y is the array of the current values of the state
			vector. The array f will be filled with the function result y' = f

INPUT PARAMETERS
----------------

	n			Dimension of the system

	x			Initial value of dependent variable (usually time)

	*y			Initial y values (NB_TYPE y[n]).

	xend		Final x value (xend-x may be positive or negative).

	*rtoler		Relative and absolute error tolerances. They can be both scalars or
	*atoler		vectors of length n (in the scalar case pass the addresses of
	 			variables where you have placed the tolerance values). If set as
	 			NULL on input, they are set to 1.0e-7 and itoler is set to 0.

	itoler		Switch for atoler and rtoler :
				itoler = 0:	Both atoler and rtoler are scalars, the code keeps
							roughly the local error of y[i] below rtoler*fabs(y[i])+atoler.
				itoler = 1: Both rtoler and atoler are vectors, the code keeps the
							local error of y[i] below rtoler[i]*fabs(y[i])+atoler[i].

	iout		Switch for calling solout :
				iout = 0: 	no call,
				iout = 1: 	solout only used for output,
				iout = 2: 	dense output is performed in solout (in this case nrdens
							must be greater than 0).

	icont		***To utilize this feature SolutionOutput routine must be rewritten**
				An array containing the indexes of components for which dense
				output is required. If no dense output is required, pass NULL.

Sophisticated setting of parameters
-----------------------------------

Several parameters have a default value (if set to 0) but, to better adapt the
code to your problem, you can specify particular initial values.

	uround		The rounding unit, default 1.0e-16.

	safe		Safety factor in the step size prediction, default 0.9.

	facl		Parameters for step size selection; the new step size is chosen
	facr		subject to the restriction  facl <= hnew/hold <= facr.
				Default values are facl = 0.2 and facr = 10.0.

	beta		The "beta" for stabilized step size control (see section IV.2 of
				Hairer and Wanner's book). Larger values for beta ( <= 0.1 ) make
				the step size control more stable. This program needs a larger
				beta than Higham & Hall. Negative initial value provoke beta = 0.0.
				Default beta = 0.04.

	hmax		Maximal step size, default xend - x.

	h			Initial step size, default is a guess computed by the function hinit.

	nmax		Maximal number of allowed steps, default 100000.

	meth		Switch for the choice of the method coefficients; at the moment the
				only possibility and default value are 1.

	nstiff		Test for stiffness is activated when the current step number is a
				multiple of nstiff. A negative value means no test and the default
				is 1000.

	nrdens		***To utilize the dense output feature for particular components of y
				***SolutionOutput routine must be rewritten
				***For now, all components of y are output
				Number of components for which dense outpout is required, default 0.
				For 0 < nrdens < n, the components have to be specified in icont[0],
				icont[1], ... icont[nrdens-1]. Note that if nrdens=0 or nrdens=n, no
				icont is needed, specify NULL.

 ***************************************************************************/
//    THIS CODE HAS BEEN SLIGHTLY MODIFIED BY THE HEVEA TEAM TO FIT OUR DATA
//    Email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr

#ifndef _NON_STIFF_INTEGRATOR_T_H_
#define _NON_STIFF_INTEGRATOR_T_H_

#include "IntegratorT.h"
#include<vector>

class NonStiffIntegratorT : public IntegratorT
{

public:

  // Constructors

  NonStiffIntegratorT( HAI_function & f,
		       const int nin, NB_TYPE xin, NB_TYPE xendin, 
		       NB_TYPE dxin, int ntotin, int nrdens, int itolerin = 0, 
		       NB_TYPE *rtolerin = 0, NB_TYPE *atolerin = 0,
		       const int ioutin = 2, NB_TYPE hin = 0.0, 
		       NB_TYPE hmaxin = 0.0, int nmaxin = 0.0, 
		       NB_TYPE uroundin = 0.0, NB_TYPE safein = 0.0, 
		       NB_TYPE faclin = 0.0, NB_TYPE facrin = 0.0, 
		       NB_TYPE betain = 0.0, int nstiffin = 0, 
		       unsigned *icontin = NULL);
  
  // Still need to implement copy constructor and default constructor

  ~NonStiffIntegratorT();

  //void Integrate();
  void Integrate(NB_TYPE yin[], 
		 std::vector< std::vector<NB_TYPE> > &flot);//ADDED BY BORIS

  // The Function to Integrate
  HAI_function & Function;
//    void Function( NB_TYPE x, NB_TYPE *y, NB_TYPE *z) {}


private:

  // core integrater for dopri5
  virtual int CoreIntegrator(NB_TYPE yin[], 
			     std::vector< std::vector<NB_TYPE> > &flot);

  // determine initial step size, h
  NB_TYPE hinit(NB_TYPE yin[]);

  // calculates intermediate values for y
  virtual NB_TYPE ContinuousOutput(unsigned i);

  // member variables

  // for stabilized step-size control
  NB_TYPE beta;
  // switch for the choice of the coefficients
  int meth;
  // test for stiffness
  int nstiff;
  // number of components for which dense output is required
	
  // To utilize the dense output feature for particular components of y
  // SolutionOutput routine must be rewritten
  // For now, all components of y are output
  int nrdens;
  // indices for components for which dense output is required
  unsigned *icont;
  // array used for dense output
  int *indir;
	
  // working arrays used for dense output
  NB_TYPE *rcont1;
  NB_TYPE *rcont2;
  NB_TYPE *rcont3;
  NB_TYPE *rcont4;
  NB_TYPE *rcont5;
	
  // vectors used in integration steps
  NB_TYPE *f;
  NB_TYPE *yy1;
  NB_TYPE *k1;
  NB_TYPE *k2;
  NB_TYPE *k3;
  NB_TYPE *k4;
  NB_TYPE *k5;
  NB_TYPE *k6;
  NB_TYPE *ysti;
	
};

#endif // _NON_STIFF_INTEGRATOR_T_H_
