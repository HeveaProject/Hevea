/***************************************************************************
                               IntegratorT.cpp
                             -------------------
    written by           : Blake Ashby
    last updated         : Nov 15, 2002
    email                : bmashby@stanford.edu
 ***************************************************************************/


#include "IntegratorT.h"

// constructors

IntegratorT::IntegratorT(const int nin, NB_TYPE xin, 
			 const NB_TYPE xendin, NB_TYPE dxin, int ntotin, 
			 int itolerin, NB_TYPE *rtolerin, NB_TYPE *atolerin, 
			 const int ioutin, NB_TYPE hin, NB_TYPE hmaxin, 
			 int nmaxin, NB_TYPE uroundin, NB_TYPE safein,
			 NB_TYPE faclin, NB_TYPE facrin) :
  n(nin), x(xin), xend(xendin), dx(dxin), ntot(ntotin), 
  itoler(itolerin), rtoler(rtolerin), atoler(atolerin), rtolerNULL(false), 
  atolerNULL(false), iout(ioutin), h(hin), hmax(hmaxin), nmax(nmaxin), 
  uround(uroundin), safe(safein), facl(faclin), facr(facrin), nfcn(0), 
  nstep(0), naccpt(0), nrejct(0), xold(xin), hold(hin), xd(xin)
{
  // n, the dimension of the system
  if (n == UINT_MAX) {
    cout << "System too big, max. n = " << UINT_MAX - 1 << endl;
    throw -1;
  }
  
  // rtoler, the relative tolerance of the integration
  if (!rtoler) {
    itoler = 0;
    rtoler = new NB_TYPE;
    *rtoler = 1.0e-7;
    rtolerNULL = true;
  }
  
  // atoler, the absolute tolerance of the integration
  if (!atoler) {
    itoler = 0;
    atoler = new NB_TYPE;
    *atoler = 1.0e-7;
    atolerNULL = true;
  }
  
  // -------- maximal step size
  if (hmax == 0.0) hmax = xend - x;
  
  // -------- nmax--maximal number of steps
  if (nmax == 0) nmax = 100000;
  if (nmax <= 0) {
    cout << " wrong input, nmax = " << nmax << endl;
    throw -1;
  }
  
  // -------- uround--smallest number satisfying 1.0 + uround > 1.0
  if (uround == 0.0) uround = 1.0e-16;
  if ((uround <= 1.0e-19) || (uround >= 1.0)) {
    cout << " coefficients have 20 digits, uround = " << uround << endl;
    throw -1;
  }
  
  // --------- safe--safety factor in step size prediction
  if (safe == 0.0) safe = 0.9;
  if ((safe <= 0.001) || (safe >= 1.0)) {
    cout << " curious input for safety factor, safe = " << safe << endl;
    throw -1;
  }
  
}  // Constructor

// Destructor
IntegratorT::~IntegratorT()
{
  if (rtolerNULL) delete rtoler;
  if (atolerNULL) delete atoler;
}

// Function that controls the output of the results.
// Modify this routine according to your needs

int IntegratorT::SolutionOutput(std::vector< std::vector<NB_TYPE> > &flot)
{
  cout << setiosflags(ios::showpoint);// | ios::fixed);
  
  if (naccpt == 0) xd = xold;
  
  NB_TYPE posneg = sign(1.0, dx);           // Ajouté par Francis.
  while ((posneg > 0 && xd < x) || (posneg < 0 && xd > x)) {
    if ( (posneg > 0 && xold <= xd) ||  (posneg < 0 && xold >= xd)) {   
      std::vector<NB_TYPE> flot_t(n);
      for ( int i=0; i < n; i++ ) flot_t[i] = ContinuousOutput(i);
      flot.push_back(flot_t);
      xd += dx;
    }
  }
  
  return 0;
  
}  // SolutionOutput
