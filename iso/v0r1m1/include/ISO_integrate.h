//************************************************************************
//                                                   ISO_integrate.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************
// In this file, we provide an interface for Blake Ashby's code 
// which is in the directory HEVEA_C++/integration

#ifndef _ISO_INTEGRATE_H_
#define _ISO_INTEGRATE_H_

//-----------------------------------------------------------//
// class HAI_function
// 
// Used to define Function in the Blake Ashby's code
//-----------------------------------------------------------//

class ISO_zeta: public HAI_function 
{
  typedef NB_TYPE Num_type;

  TOR_embedding<Num_type> const  & F_;
  UTI_2vector<Num_type> V_, U_;
public:
  ISO_zeta(TOR_embedding<Num_type> const & embedding, 
	    UTI_2vector<Num_type> d, 
	    UTI_2vector<Num_type> d_perp): 
    F_(embedding), V_(d), U_(d_perp) {}

  //FUNCTION needed by B. Ashby's code. X'=f(t,X) = zeta(X,t) with X=y[0]
  void operator() (Num_type t, Num_type *y, Num_type *f)
  {
    // We need to know zeta at the point (y[0],t) in the basis
    // (V[i],U[i]).
    UTI_2vector<Num_type> param = y[0]*V_ + t*U_;
    UTI_3vector<Num_type> V_F = F_.deriv(V_,param);
    // zeta = - (U(F_) * V(F_)) / (V(F_) * V(F_))
    f[0] = - (F_.deriv(U_,param) * V_F) / (V_F*V_F);;    
    return;
  }
};

class ISO_h_s: public HAI_function 
{
  typedef NB_TYPE Num_type;

  ISO_embedding<Num_type> const & iso_;
  TOR_embedding<Num_type> const & F_;
  int const stage_;
  int const dir_;
  CYL_function_S1<Num_type> const & phi_X_;
  UTI_2vector<Num_type> const & V_;
  UTI_2vector<Num_type> const & U_;
  Num_type s_; // h_s(t) = h(s,t)
  int n_oscil_;

public:
  ISO_h_s(Num_type s,
	  ISO_embedding<Num_type> const & iso, 
	  TOR_embedding<Num_type> const & F,
	  int const k, 
	  int const i,
	  CYL_function_S1<Num_type> const & phi_X,
	  UTI_2vector<Num_type> const & d,
	  UTI_2vector<Num_type> const & d_perp,
	  int n_oscil) : 
    s_(s), iso_(iso), F_(F), stage_(k), dir_(i), 
    phi_X_(phi_X), V_(d), U_(d_perp), n_oscil_(n_oscil) {}

  //FUNCTION  needed by B. Ashby's code. X' = f(t,X) = h(t)
  void operator() (Num_type t, Num_type *x, Num_type *f)
  {
    UTI_2vector<Num_type> phi_s_t = phi_X_(s_,t) * V_ + t * U_;
    
    UTI_2vector<Num_type> W = iso_.W(phi_s_t, V_, U_, F_);
    UTI_3vector<Num_type> z = F_.deriv(W, phi_s_t);
    UTI_3vector<Num_type> n = F_.normal(phi_s_t);

    Num_type r = iso_.radius(phi_s_t, stage_, dir_, V_, U_, F_);

    Num_type alpha = J0_inv( z.norm() / r );
    Num_type beta = cos( 2*M_PI*(n_oscil_*t - floor(n_oscil_*t)));
    
    z.normalize();
    UTI_3vector<Num_type> h = r*( cos(alpha*beta)*z + sin(alpha*beta)*n);
    
    for (int k=0; k<3; k++) 
      f[k] = h[k];    
    return;
  }
};

//------------------------------------------------------------------------
//  ISO_embedding< Num_type > flow()
//
// Computation of the fisrt coordinate of the flow in the basis (V[i],U[i]))
//-----------------------------------------------------------------------
template <typename Num_type >  CYL_function_S1<Num_type> 
ISO_embedding<Num_type>::
flow(TOR_embedding<Num_type> const & f,
     UTI_2vector<Num_type> const & V,
     UTI_2vector<Num_type> const & U) const
{
  // Function to integrate in B. Ashby's code
  ISO_zeta zeta_function(f, V, U);

  // FOR DEBUG
  if (debug_level_ > 2) {
    for (int i = 0; i < 10; i++) {
      Num_type y_tmp[1] = {i*i/100.};
      Num_type z_tmp[1];
      zeta_function(i/10., y_tmp, z_tmp);
      cout << " zeta(" << y_tmp[0] << ", " << i/10. << ") = " << z_tmp[0] 
	   << endl;
    }
  }
  // END FOR DEBUG

  //PARAMETERS OF B ASHBY'S CODE
  // dimension of problem
  const int dim = 1;

  // initial value for x (correspond au temps t=0)
  Num_type xbeg = 0.0;
  // final value for x (a quel temps t on arrete)
  const Num_type xend = 1.;
  // interval of x for printing output
  Num_type dy = 1 / ((Num_type) n_y_);
  //ERROR PARAMETERS
  int nrdens = dim;
  // rtoler and atoler are scalars
  int itoler = 0;
  // relative tolerance
  Num_type rtoler = 1.0e-10;
  // absolute tolerance
  Num_type atoler = 1.0e-10;


  // There are other hidden parameters (see header files of Hairer code) 
  // for these parameters.

  // Creation of flot_X
  CYL_function_S1< Num_type > flot_X(n_x_, n_y_); //cautious, it means that...
                              // ..there are n_x_ times n_y_+1 possible values.
    
  int i;
#pragma omp parallel for shared(flot_X, zeta_function) private(i)
  for (i=0; i < n_x_; i++) {  
    Num_type s = (Num_type)i / ((Num_type) n_x_);
    // Initialisation at the point (s,0)
    // initial values for y
    Num_type y[1] = {s};

    // We build the nonstiffT integrator
    NonStiffIntegratorT nonstiffT(zeta_function, dim, xbeg, xend, dy, 
				  n_y_+1, nrdens, itoler, &rtoler, &atoler );
    
    // We integrate. The result : flot[i] contains y[0].
    std::vector< std::vector<Num_type> > flot;
    flot.reserve(n_y_ + 1);
    nonstiffT.Integrate(y, flot);
    assert ( flot.size() == n_y_ + 1 );
    // We write the result in flot_X
    for (int j=0; j < n_y_ + 1; j++) { // Cautious, +1 for the cylinder.
      flot_X.mat(i,j)=flot[j][0];     
    }
  }
  return flot_X;
}

//------------------------------------------------------------------------
//  ISO_embedding< Num_type > integrate_func()
//
// Computation of the integral curve of X' = f, with initial condition init_cond
// The result is stored in flot. The curve is discretised with Ny + 1 
// points.
//-----------------------------------------------------------------------
template <typename Num_type >  void 
ISO_embedding<Num_type>::integrate_func( 
			   HAI_function & f,
			   int const dim, // dimension of problem
			   Num_type *init_cond,    
			   std::vector< std::vector<Num_type> > & flot,
			   int Ny ) const
{
  //BLAKE'S PARAMETERS

  // initial value for x (corresponding to time t=0)
  Num_type xbeg = 0.0;
  // final value for x at what time do we stop)
  const Num_type xend = 1;
  // interval of x for printing output
  Num_type dy = 1 / ((Num_type) Ny);

    //ERROR PARAMETERS
  int nrdens = dim;
  // rtoler and atoler are scalars
  int itoler = 0;
  // relative tolerance
  Num_type rtoler = 1.0e-10;
  // absolute tolerance
  Num_type atoler = 1.0e-10;
  // There are other hidden parameters (see header files of Hairer code) 
  // for these parameters.

  // We build nonstiffT
  NonStiffIntegratorT nonstiffT(f, dim, xbeg, xend, dy, Ny+1, 
				nrdens, itoler, &rtoler, &atoler );
    
  nonstiffT.Integrate(init_cond, flot);
  assert ( flot.size() == Ny + 1 );
}

//------------------------------------------------------------------------
//  ISO_embedding< Num_type > integrate_h_hai()
//
// Integration along flow lines with Ashby's method 
//-----------------------------------------------------------------------
template <typename Num_type>  CYL_embedding<Num_type> 
ISO_embedding<Num_type>::
integrate_h_hai(CYL_function_S1< Num_type > const & phi_X,
		TOR_embedding<Num_type> const & f,
		UTI_2vector<Num_type> const & V,
		UTI_2vector<Num_type> const & U,
		int const stage,
		int const dir,
		int const N_oscil) const
{
  int nx(n_x_), ny(n_y_);

  CYL_embedding<Num_type> ftilde(nx, ny);
  int i;

#pragma omp parallel for shared(ftilde) private(i)
  for(i=0; i < nx; i++) {
    // Value of s  for the initial condition
    Num_type s = (Num_type)i/(Num_type) nx;
    //Coordonnees (x,y) du point du flot au debut.
    UTI_2vector<Num_type> phi_i0 = phi_X(s,0) * V;
    UTI_3point<Num_type> f_phi_i0 =  f(phi_i0);


    ISO_h_s h_func(s, *this, f, stage, dir, phi_X, V, U, N_oscil);
    std::vector< std::vector<Num_type> > flot;
    flot.reserve(ny + 1);
    Num_type init_cond[3] = {0., 0., 0.};
    integrate_func( h_func, 3, init_cond, flot, ny );

    for (int j = 0; j <= ny; j++)
      for (int k = 0; k < 3; k++)
	ftilde.mat(i,j)[k] = f_phi_i0[k] + flot[j][k];
  }

  return ftilde;
}

//------------------------------------------------------------------------
//  ISO_embedding< Num_type > integrate_h_trap()
//
// Integration along flow lines with trapeze's method 
//-----------------------------------------------------------------------
template <typename Num_type>  CYL_embedding<Num_type> 
ISO_embedding<Num_type>::
integrate_h_trap(CYL_function_S1< Num_type > const & phi_X,	  
		 TOR_embedding<Num_type> const & f,
		 UTI_2vector<Num_type> const & V,
		 UTI_2vector<Num_type> const & U,
		 int const stage,
		 int const dir,
		 int const N_oscil) const
{
  CYL_embedding<Num_type> ftilde(n_x_, n_y_);

  for(int i=0; i < n_x_; i++) {
 
    Num_type s = (Num_type)i/(Num_type) n_x_;
    ISO_h_s h_func(s, *this, f, stage, dir, phi_X, V, U, N_oscil);
    Num_type y_tmp[3] = {0,0,0};
    Num_type h_tmp[3];

    //Coordinate of the flow at the beginning.
    UTI_2vector<Num_type> phi_i0 = phi_X.mat(i,0) * V; 
    ftilde.mat(i,0) = f(phi_i0);

    h_func(0, y_tmp, h_tmp);
    UTI_3vector<Num_type> h_vieux(h_tmp[0],h_tmp[1], h_tmp[2]);
    for (int j = 1; j < n_y_ + 1; j++) {
      // Value of t
      Num_type t = (Num_type)j/((Num_type) n_y_);
      h_func(t, y_tmp, h_tmp); // We use the same function than Ashby's code
      UTI_3vector<Num_type> h(h_tmp[0],h_tmp[1], h_tmp[2]);// Value dof h
      ftilde.mat(i,j) = ftilde.mat(i,j-1) + (h+h_vieux) / (2.*(Num_type) n_y_);
      h_vieux=h;
    }
  }
  return ftilde;
}

#endif // _ISO_INTEGRATE_H_
