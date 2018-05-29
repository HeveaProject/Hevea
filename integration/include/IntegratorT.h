/***************************************************************************
                               IntegratorT.h
                             -------------------
    written by           : Blake Ashby
    last updated         : Nov 15, 2002
    email                : bmashby@stanford.edu

This is the parent class for StiffIntegratorT (implicit integrator based on
RADAU5) and NonStiffIntegratorT (explicit integrator based on DOPRI5). The
code has been modified from the code written originally in Fortran by

         E. Hairer and G. Wanner
         Universite de Geneve, Dept. de Mathematiques
         Ch-1211 Geneve 24, Switzerland
         E-mail:  ernst.hairer@math.unige.ch
                  gerhard.wanner@math.unige.ch

See the header files StiffIntegratorT.h and NonStiffIntegratorT.h for
more details.

 ***************************************************************************/

#ifndef _INTEGRATOR_T_H_
#define _INTEGRATOR_T_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <limits.h>
#include <decsol.h>
#include <cmath>
#include <cstring>
#include <vector>

using namespace std;

class HAI_function
{
public: // pure virtual operator. Must be implemented in a derived class.
  virtual void operator() (NB_TYPE x, NB_TYPE *y, NB_TYPE *f) = 0;
};

class IntegratorT
{
public:

  // Constructor	
  IntegratorT(const int nin, NB_TYPE xin, const NB_TYPE xendin, 
	      NB_TYPE dxin, int ntotin = 1, int itolerin = 0, 
	      NB_TYPE *rtolerin = 0, NB_TYPE *atolerin = 0, 
	      const int ioutin = 1, NB_TYPE hin = 0.0, NB_TYPE hmaxin = 0.0, 
	      int nmaxin = 0, NB_TYPE uroundin = 0.0, NB_TYPE safein = 0.0, 
	      NB_TYPE faclin = 0.0, NB_TYPE facrin = 0.0);
  
  // Still need to implement copy constructor and default constructor
  
  // Destructor
  virtual ~IntegratorT();
  
  virtual void Integrate(NB_TYPE yin[], 
			 std::vector< std::vector<NB_TYPE> > &flot) = 0;
  
  // Function that controls the output of the results.
  //Modify this routine according to your needs
  int SolutionOutput(std::vector< std::vector<NB_TYPE> > &flot);
  
  // get number of function evaluations
  int NumFunction() const { return nfcn; }
  // get number of attempted steps
  int NumStep() const { return nstep; }
  // get number of accepted steps
  int NumAccept() const { return naccpt; }
  // get number of rejected steps
  int NumReject() const { return nrejct; }
  
protected:
  
  // CoreIntegrator
  virtual int CoreIntegrator(NB_TYPE yin[], // vector for y initial values
			     std::vector< std::vector<NB_TYPE> > &flot) = 0;
  
  virtual NB_TYPE ContinuousOutput(unsigned i) = 0;
  
  // Member variables
  
  // dimension of system
  const int n;
  // independent variable (usually time)
  NB_TYPE x;
  // final value for independent variable
  const NB_TYPE xend;
  // time step for intermediate output
  NB_TYPE dx;
  // total number of desired output values. Should be 1 + (xend - x)/dx
  int ntot;
  // switch for rtol and atol; if itol = 0, rtol and atol are scalars
  // if itol = 1, rtol and atol are vectors
  int itoler;
  // relative error tolerance
  NB_TYPE *rtoler;
  // absolute error tolerance
  NB_TYPE *atoler;
  // variables that indicate whether rtoler and atoler are NULL on input
  // or not--needed for proper memory management in destructor
  bool rtolerNULL;
  bool atolerNULL;
  // routine for dense output at every time step is called if iout = 1
  const int iout;
  // integration step length
  NB_TYPE h;
  
  // Derived variables
  
  // maximal step size
  NB_TYPE hmax;
  // maximal number of steps
  int nmax;
  // smallest number satisfying 1.0 + uround > 1.0
  NB_TYPE uround;
  // safety factor in step size prediction
  NB_TYPE safe;
  // facl, facr--parameters for step size selection
  NB_TYPE facl;
  NB_TYPE facr;
  
  // Counting variables
  
  // number of function evaluations (not counting those in numerical
  // Jacobian calculations)
  int nfcn;
  // number of attempted steps
  int nstep;
  // number of accepted steps
  int naccpt;
  // number of rejected steps
  int nrejct;
  
  // stores past value of x
  NB_TYPE xold;
  // stores past value of h
  NB_TYPE hold;
  // x at discrete points specified by dx interval
  NB_TYPE xd;
};

#endif // _INTEGRATOR_T_H_
