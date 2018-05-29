//************************************************************************
//                                                  CYL_interp_S1.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains interpolation functions for bidimentional grids.
//
// =======================================================

#ifndef _CYL_INTERP_S1_H_
#define _CYL_INTERP_S1_H_

//--------------------------------------------------------------------
//          cyl_bilinear_interp_S1()
//
// generic bilinear interpolation formula for a point (x,y) of [0,1]^2
// in the square [i,i+1]*[j,j+1] in the m grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
// We assume that the image of the line segment beteen two successive 
// samples is a shortest path  (!= géodésique) between the image of the
// two samples.
//--------------------------------------------------------------------

template < typename S1_type, typename Num_type > 
S1_type
cyl_bilinear_interp_S1(Num_type x, Num_type y, // coordinates of the point in the
		       int i, int j,           // square of lower left vertex (i,j)
		       MAT_matrix_2<S1_type> const & m )    
{
  int nx = m.nb_row();
  
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;

  S1_type f_ll = m( i, j ), 
             f_rl = m( i_plus_1, j ), 
             f_lu = m( i, j_plus_1 ), 
             f_ru = m( i_plus_1, j_plus_1 );

  Num_type delta_l = f_rl - f_ll;
  Num_type delta_u = f_ru - f_lu;
  if (i_plus_1 == 0) { // We assume the only discontinuity is 
    delta_l += 1;      // in between the samples nx - 1 et 0.
    delta_u += 1;
  }
  assert( delta_l < .5 && delta_l > -.5 && delta_u < .5 && delta_u > -.5 );

  S1_type temp_l = f_ll + x * delta_l;
  S1_type temp_u = f_lu + x * delta_u;

  return temp_l + y * ( temp_u - temp_l );
}

//--------------------------------------------------------
//          cyl_bicubic_interp_S1()
//
// generic bicubic interpolation formula for a point (x,y) of [0,1]^2
// in the square [i,i+1]*[j,j+1] in the m grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
// The formulae is based on a bicubic interpolation with C-splines.
//--------------------------------------------------------

template < typename S1_type, typename Num_type > 
S1_type
cyl_bicubic_interp_S1(Num_type x, Num_type y, // coordinates of the point in the
		      int i, int j,  // square of lower left vertex (i,j)
		      MAT_matrix_2<S1_type> const & m ) 
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1; // Caution : the cylinder's size is nx*(ny+1) !!
  assert( m.nb_col() > 5 && j < ny );
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type v_i_shift = ( i_plus_2 - i_moins_2 - 4 + 
			 8*(i_moins_1 - i_plus_1 + 2) )/ nx;
  Num_type v_i_plus_1_shift = ( i_plus_3 - i_moins_1 - 4 + 
				8*(i - i_plus_2 + 2) )/ nx;
  Num_type diff_i_shift = (i + 1 - i_plus_1)/nx;

  Num_type h_01_x = x*x*(3 - 2*x);  // Hermite polynomials
  Num_type h_10_x = x*(x*(x - 2) + 1); 
  Num_type h_11_x = x*x*(x - 1); 
  
  S1_type f_j_moins_2[6];
  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 derivative approximation
    Num_type v_i = v_i_shift + 
                   ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
              8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= 12.;

    Num_type v_i_plus_1 = v_i_plus_1_shift + 
           ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
      8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= 12.;

    Num_type diff_i = diff_i_shift + 
      m(i_plus_1, jk_moins_2) - m(i, jk_moins_2);

    f_j_moins_2[k] =  m(i, jk_moins_2) + 
      h_01_x * diff_i + h_10_x * v_i + h_11_x * v_i_plus_1;
  }

  Num_type v_j = finite_diff_derivative( f_j_moins_2, j_aux - j, 0 );
  Num_type v_j_plus_1 = finite_diff_derivative( 
					    f_j_moins_2, j_aux - j, 1 );
  S1_type f_j = f_j_moins_2[2 + j - j_aux];
  S1_type f_j_plus_1 = f_j_moins_2[3 + j - j_aux];

  return c_spline_interp(y, f_j, f_j_plus_1, v_j, v_j_plus_1);
}

//--------------------------------------------------------
//          cyl_generic_interp_S1()
//
// generic interpolation formula for a point (x,y) of [0,1]^2
//replaced in grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//--------------------------------------------------------

template < typename S1_type, typename Num_type > 
S1_type
cyl_generic_interp_S1( Num_type x, 
		       Num_type y, 
		       MAT_matrix_2<S1_type> const & m,
		       int degree )  // the interpolation degree
{
  // x and y are supposed to be in [0,1[. x is actually of the
  // form t - floor(t) for some t when cyl_interp is called. It seems
  // that on some processors we may have t - floor(t) >= 1. We perform
  // some preprocessing to ensure x in [0,1[.

  if (x < 0) x += 1;
  if (x >= 1) x -= 1;
  assert( x >= 0 && x < 1);
  assert( y >= 0 && y <= 1);

  int nx = m.nb_row();
  int ny = m.nb_col() - 1; // Caution : the cylinder's size is nx*(ny+1) !!
  Num_type x_scaled = x * nx;
  Num_type y_scaled = y * ny;
  
  int i = floor( x_scaled );
  int j = floor( y_scaled );

  Num_type u = x_scaled - i;
  Num_type v = y_scaled - j;

  if (j == ny) { 
    assert( v == 0. );
    v = 1;
    j = ny -1;
  }

  if (degree == 1)
    return cyl_bilinear_interp_S1<S1_type, Num_type> (u, v, i, j, m);
  else 
    return cyl_bicubic_interp_S1<S1_type, Num_type> (u, v, i, j, m);
}

//------------------------------------------------------------
//          cyl_bilinear_deriv_x_S1()
//
// d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bilinear interpolation. 
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Num_type > 
Num_type
cyl_bilinear_deriv_x_S1(Num_type x, 
			Num_type y,   // coordinates of the point in the
			int i, int j, // square of lower left vertex (i,j)
			MAT_matrix_2<S1_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1;// Caution : the cylinder's size is nx*(ny+1) !!
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;
  assert( j < ny ); 

  // So as to take into account the modulo 1, we assume that the only 
  // discontinuity is a gap of +1 between the samples  nx - 1 et 0.
  int delta = i + 1 - i_plus_1;

  // Multiply by nx to take reparameterization into account. Add delta
  // to take discontinuity into account.
  return delta + (Num_type)nx * ( (1.- y)*(m(i_plus_1,j) - m(i,j)) +
			  y * (m(i_plus_1,j_plus_1) - m(i,j_plus_1)) );  
}

//------------------------------------------------------------
//          cyl_bilinear_deriv_y_S1()
//
// d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bilinear interpolation. 
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Num_type > 
Num_type
cyl_bilinear_deriv_y_S1(Num_type x, 
			Num_type y,    // coordinates of the point in the
			int i, int j, // square of lower left vertex (i,j)
			MAT_matrix_2<S1_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1; // Caution : the cylinder's size is nx*(ny+1) !!
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;
  assert( j < ny ); 

  // Multiply by ny to take reparameterization into account
  return (Num_type)ny * ( (1.- x)*(m(i,j_plus_1) - m(i,j)) +
			  x * (m(i_plus_1,j_plus_1) - m(i_plus_1,j)) );  
}

//------------------------------------------------------------
//          cyl_bicubic_deriv_x_S1()
//
// d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bicubic interpolation. 
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bicubic_deriv_x_S1( Num_type x, 
			Num_type y,  // coordinates of the point in the
			int i, int j,  // square of lower left vertex (i,j)
			MAT_matrix_2<S1_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1;  // Caution : the cylinder's size is nx*(ny+1) !!
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  // So as to take into account the modulo 1, we assume that the only 
  // discontinuity is a gap of +1 between the samples  nx - 1 et 0.
  Vector_type v_i_shift = ( i_plus_2 - i_moins_2 - 4 + 
			 8*(i_moins_1 - i_plus_1 + 2) )/ nx;
  Vector_type v_i_plus_1_shift = ( i_plus_3 - i_moins_1 - 4 + 
				8*(i - i_plus_2 + 2) )/ nx;
  Vector_type diff_i_shift = (i + 1 - i_plus_1)/nx;

  Num_type dh_01_dx = 6*x*(1-x);  // derivatives of Hermite polynomials
  Num_type dh_10_dx = (x-1)*(3*x - 1); 
  Num_type dh_11_dx = x*(3*x - 2); 
  
  Vector_type df_dx_j_moins_2[6];
  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 derivative approximation
    Vector_type v_i = v_i_shift + 
                      ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
                 8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= (Num_type)12.;

    Vector_type v_i_plus_1 = v_i_plus_1_shift + 
                      ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
                 8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= (Num_type)12.;

    Vector_type diff_i = diff_i_shift + 
      m(i_plus_1, jk_moins_2) - m(i, jk_moins_2);
    df_dx_j_moins_2[k] = dh_01_dx * diff_i + 
                         dh_10_dx * v_i + dh_11_dx * v_i_plus_1;
  }

  Vector_type v_j = finite_diff_derivative( df_dx_j_moins_2, j_aux - j, 0 );
  Vector_type v_j_plus_1 = finite_diff_derivative( 
					    df_dx_j_moins_2, j_aux - j, 1 );
  Vector_type df_dx_j = df_dx_j_moins_2[2 + j - j_aux];
  Vector_type df_dx_j_plus_1 = df_dx_j_moins_2[3 + j - j_aux];

  // Multiply by nx to take reparameterization into account
  return (Num_type)nx * c_spline_interp( y, 
					 df_dx_j, 
					 df_dx_j_plus_1, 
					 v_j, 
					 v_j_plus_1);
}

//------------------------------------------------------------
//          cyl_bicubic_deriv_y_S1()
//
// d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bicubic interpolation. 
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bicubic_deriv_y_S1( Num_type x, 
			Num_type y, // coordinates of the point in the
			int i, int j,  // square of lower left vertex (i,j)
			MAT_matrix_2<S1_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1; // Caution : the cylinder's size is nx*(ny+1) !!
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  // So as to take into account the modulo 1, we assume that the only 
  // discontinuity is a gap of +1 between the samples  nx - 1 et 0.
  Vector_type v_i_shift = ( i_plus_2 - i_moins_2 - 4 + 
			 8*(i_moins_1 - i_plus_1 + 2) )/ nx;
  Vector_type v_i_plus_1_shift = ( i_plus_3 - i_moins_1 - 4 + 
				8*(i - i_plus_2 + 2) )/ nx;
  Vector_type diff_i_shift = (i + 1 - i_plus_1)/nx;

  Num_type h_01_x = x*x*(3 - 2*x);  // Hermite polynomials
  Num_type h_10_x = x*(x*(x - 2) + 1); 
  Num_type h_11_x = x*x*(x - 1); 

  S1_type f_j_moins_2[6];
  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 derivative approximation
    Vector_type v_i = v_i_shift + 
                      ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
                 8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= 12.;

    Vector_type v_i_plus_1 = v_i_plus_1_shift + 
              ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
         8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= 12.;

    Vector_type diff_i = diff_i_shift + 
      m(i_plus_1, jk_moins_2) - m(i, jk_moins_2);
    f_j_moins_2[k] =  m(i, jk_moins_2) + 
      h_01_x * diff_i + h_10_x * v_i + h_11_x * v_i_plus_1;
  }

  Vector_type v_j = finite_diff_derivative( f_j_moins_2, j_aux - j, 0 );
  Vector_type v_j_plus_1 = finite_diff_derivative( 
					    f_j_moins_2, j_aux - j, 1 );
  S1_type f_j = f_j_moins_2[2 + j - j_aux];
  S1_type f_j_plus_1 = f_j_moins_2[3 + j - j_aux];

  Num_type dh_01_dy = 6*y*(1-y);  // derivatives of Hermite polynomials
  Num_type dh_10_dy = (y-1)*(3*y - 1); 
  Num_type dh_11_dy = y*(3*y - 2); 
  
  // Multiply by ny to take reparameterization in to account
  return (Num_type)ny * (dh_01_dy*( f_j_plus_1 - f_j ) + 
			 dh_10_dy*v_j + dh_11_dy*v_j_plus_1);
}

//------------------------------------------------------------
//          cyl_generic_deriv_x_S1()
//
// generic d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Num_type > 
Num_type
cyl_generic_deriv_x_S1( Num_type x, 
			Num_type y, 
			MAT_matrix_2<S1_type> const & m,
			int degree )
{
  // x and y are supposed to be in [0,1[. x is actually set to some 
  // t - floor(t) before calling. It seems that on some processors we may have 
  // t - floor(t) >= 1. We perform some preprocessing to ensure x,y in [0,1[.

  if (x < 0) x += 1;
  if (x >= 1) x -= 1;

  assert( x >= 0 && x <= 1 && y >= 0 && y <= 1 );

  int nx = m.nb_row();
  int ny = m.nb_col() - 1;
  Num_type x_scaled = x * nx;
  Num_type y_scaled = y * ny;
  
  int i = floor( x_scaled );
  int j = floor( y_scaled );

  Num_type u = x_scaled - i;
  Num_type v = y_scaled - j;

  if (j == ny) { 
    assert( v == 0. );
    v = 1;
    j = ny -1;
  }
  if (degree == 1) 
    return cyl_bilinear_deriv_x_S1<S1_type, Num_type>( u, v, i, j, m);
  else 
    return cyl_bicubic_deriv_x_S1<S1_type, Num_type, Num_type >( 
							       u, v, i, j, m );
}

//------------------------------------------------------------
//          cyl_generic_deriv_y_S1()
//
// generic d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. The image is in S1=R/Z. 
//------------------------------------------------------------

template < typename S1_type, typename Num_type > 
Num_type
cyl_generic_deriv_y_S1( Num_type x, 
			Num_type y, 
			MAT_matrix_2<S1_type> const & m,
			int degree )
{
  // x and y are supposed to be in [0,1[. x is actually set to some 
  // t - floor(t). It seems that on some processors we may have 
  // t - floor(t) >= 1. We perform some preprocessing to ensure x,y in [0,1[.

  if (x < 0) x += 1;
  if (x >= 1) x -= 1;

  assert( x >= 0 && x <= 1 && y >= 0 && y <= 1 );

  int nx = m.nb_row();
  int ny = m.nb_col() - 1;
  Num_type x_scaled = x * nx;
  Num_type y_scaled = y * ny;
  
  int i = floor( x_scaled );
  int j = floor( y_scaled );

  Num_type u = x_scaled - i;
  Num_type v = y_scaled - j;

  if (j == ny) { 
    assert( v == 0. );
    v = 1;
    j = ny -1;
  }
  if (degree == 1) 
    return cyl_bilinear_deriv_y_S1<S1_type, Num_type> ( u, v, i, j, m);
  else 
    return cyl_bicubic_deriv_y_S1<S1_type, Num_type>( u, v, i, j, m );
}

#endif // _CYL_INTERP_S1_H_
