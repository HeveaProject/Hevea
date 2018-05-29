//************************************************************************
//                                                    CYL_cylinder.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : Aug 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains the implementation of the class members of 
// CYL_function_S1 or embedding
// =======================================================

#ifndef _CYL_IMPLEMENTATION_H_
#define _CYL_IMPLEMENTATION_H_

#include <MAT_matrices.h>
#include <CYL_cylinder.h>
#include <cassert>

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::normal()
//
// Returns the normal to the embedding.
// i.e. n = df/dx ^ df/dy / || df/dx ^ df/dy ||
//-------------------------------------------------------------
template < typename Num_type > UTI_3vector<Num_type>
CYL_embedding <Num_type>::normal( Num_type x, Num_type y ) const
{
  UTI_3vector<Num_type> vx = deriv_x(x,y);
  UTI_3vector<Num_type> vy = deriv_y(x,y);

  UTI_3vector<Num_type> n = cross_prod( vx, vy );
  return n.normalize();
} 

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::pull_back()
//
// Returns the pullback of the embedding at the point (x,y).
//-------------------------------------------------------------
template < typename Num_type > UTI_sym_2form<Num_type>
CYL_embedding <Num_type>::pull_back( Num_type x, Num_type y ) const
{
  UTI_3vector<Num_type> vx = deriv_x(x,y);
  UTI_3vector<Num_type> vy = deriv_y(x,y);

  return UTI_sym_2form<Num_type>( vx * vx, vx * vy, vy * vy );
}

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::gluing()
//
// Here, we glue the cylinder so that the embedding of the two 
// boundaries of the cylinder coincide. The space parameter is
// deformed so that the point of parameter (s,t) corresponds 
// to the point of parameter  (phi_X(s,t),t)  expressed in the 
// basis (V,U) 
//-------------------------------------------------------------
template < typename Num_type > void CYL_embedding <Num_type>::
gluing( TOR_embedding<Num_type> const & f,
	 CYL_function_S1<Num_type> const & phi_X,
	 UTI_2vector<Num_type> const & d_2,
	 UTI_2vector<Num_type> const & d_2_perp )
{
  int nx(n_x()), ny(n_y());
  assert( phi_X.n_x() == nx && phi_X.n_y() == ny );

  int i;
#pragma omp parallel for private(i)
  for (i=0; i< nx; i++) {
    UTI_2vector<Num_type> param = phi_X.mat(i,ny)*d_2 + d_2_perp;
    UTI_3vector<Num_type> delta_f_i = f( param ) - m_(i, ny);
    for (int j=0;j < ny + 1; j++) {
      m_(i,j) += interpolant((Num_type)j/(Num_type)ny)* delta_f_i;
    }
  }
}

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::cyl_to_torus()
//
// Here, we copy an embedding of a cylinder to an embedding of a torus.
// We assume that the two boundaries of the cylinder coincide.
// We also suppose that the parameter space of the cylinder is deformed 
// so that the point of parameter (s,t) corresponds to the point of 
// parameter (phi_X(s,t),t)  expressed in the basis (V,U). The basis
// (d_2=V,d_3) is used so as to sweep the point samples. The parameter 
// lambda is used for the change of basis.
//-------------------------------------------------------------
template < typename Num_type > TOR_embedding<Num_type> & 
CYL_embedding <Num_type>::cyl_to_torus(CYL_function_S1<Num_type> const & phi_X,
				       Num_type const lambda,
				       UTI_2vector<Num_type> const & d_2,
				       UTI_2vector<Num_type> const & d_3,
				       TOR_embedding<Num_type> & f) const
{
  int nx(n_x()), ny(n_y());
  assert( nx == ny );

  UTI_2vector<int> d_2_int( (int) d_2.x(), (int) d_2.y() ); // The int version. 
  UTI_2vector<int> d_3_int( (int) d_3.x(), (int) d_3.y() ); // The int version.

  for (int l=0; l < ny; l++) { // l=0 is theoretically useless since f is not changed.
    Num_type x_delta = lambda*l; // Note: nx/ny = 1;
    for(int k=0; k < nx; k++) {
      Num_type x_k = nx * phi_X.mat(k,l) +  x_delta;  // n_x * xin the basis 
                                                      // (d_2, d_3).
      Num_type x_k_plus_1 = x_delta; 
      if (k != nx - 1) x_k_plus_1 += nx * phi_X.mat(k+1, l); 
      else x_k_plus_1 += nx * (phi_X.mat(0, l) + 1.); 

      Num_type x_k_moins_1 = x_delta; 
      if (k != 0) x_k_moins_1 += nx * phi_X.mat(k - 1, l);
      else x_k_moins_1 += nx * (phi_X.mat(nx - 1, l) - 1.); 

      Num_type x_k_plus_2 = x_delta; 
      if (k < nx - 2) x_k_plus_2 += nx * phi_X.mat(k+2, l); 
      else if (k == nx - 1) x_k_plus_2 += nx * (phi_X.mat(1, l) + 1.); 
      else if (k == nx - 2) x_k_plus_2 += nx * (phi_X.mat(0, l) + 1.); 

      // Approximation of derivatives in x_k and x_k+1 by mean of left and right derivatives.
      UTI_3vector<Num_type> v_k = mat((k+1) % nx,l) - mat(k,l) +
	(x_k_plus_1 - x_k) / (x_k - x_k_moins_1) * 
	(mat(k,l) - mat((k+nx-1) % nx,l));
      v_k /= 2;

      UTI_3vector<Num_type> v_k_plus_1 = mat((k+1) % nx,l) - mat(k,l) +
	(x_k_plus_1 - x_k) / (x_k_plus_2 - x_k_plus_1) * 
	(mat((k+2) % nx,l) - mat((k+1) % nx,l));
      v_k_plus_1 /= 2;

      int q_beg = ceil( x_k );
      int q_end = ceil( x_k_plus_1 );
      for (int q = q_beg; q < q_end; q++) {
	Num_type t = (q - x_k) / (x_k_plus_1 - x_k);

	UTI_2vector<int> param = q * d_2_int + l * d_3_int;
	int i = param.x() % nx; // can be < 0.
	if (i < 0) i += nx; assert( i >= 0 && i < nx );
	int j = param.y() % ny;
	if (j < 0) j += ny; assert( j >= 0 && j < ny );
	f.mat(i,j) = c_spline_interp( t, 
				      mat(k,l), 
				      mat((k+1) % nx,l), 
				      v_k, 
				      v_k_plus_1 );
      }
    }
  }
  return f;
}

#endif // _CYL_IMPLEMENTATION_H_
