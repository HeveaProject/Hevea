//************************************************************************
//                                                    CYL_interp.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains interpolation functions for bidimentional grids.
//
// =======================================================

#ifndef _CYL_INTERP_H_
#define _CYL_INTERP_H_


//-------------------------------------------------------------------
//          finite_diff_derivative()
//
// Order 4 finite difference derivative formulas for point k (=0,1)
// including shifts for boundary evaluation. Formulas taken from
// reference.wolfram.com/mathematica/tutorial/NDSolvePDE.html Prend en
// compte le remplissage de f[] avec les décalages d'indices
// (cf. cyl_bicubic_interp() ).
// //----------------------------------------------------------------
template < typename Vector_type > 
Vector_type
finite_diff_derivative(Vector_type *f_k_moins_2, int shift, int k )
{
  assert ( k == 0 || k == 1 );

  Vector_type v;

  if ( shift == 0 ) {  
    v = (f_k_moins_2[k] - f_k_moins_2[k+4]) +
      8. * (f_k_moins_2[k+3] - f_k_moins_2[k+1]);
  }
  else if ( shift == 2 && k == 0 ) { 
    v = -25. *f_k_moins_2[0] + 48.*f_k_moins_2[1] -
      36.*f_k_moins_2[2] + 16.*f_k_moins_2[3] -3.*f_k_moins_2[4];
  }
  else if ( (shift == 2 && k == 1) || (shift == 1 && k == 0) ) { 
    v = -3. *f_k_moins_2[0] - 10.*f_k_moins_2[1] + 
      18.*f_k_moins_2[2] - 6.*f_k_moins_2[3] + f_k_moins_2[4];
  }
  else if ( shift == 1 && k == 1 || (shift == -1 && k == 0) ) { 
    v = (f_k_moins_2[1-k] - f_k_moins_2[5-k]) +
      8. * (f_k_moins_2[4-k] - f_k_moins_2[2-k]);
  }
  else if ( (shift == -1 && k == 1) || (shift == -2 && k == 0) ) { 
    v = -1.*f_k_moins_2[1] + 6.*f_k_moins_2[2] - 
      18.*f_k_moins_2[3] + 10.*f_k_moins_2[4] + 3.*f_k_moins_2[5];
  }
  else if ( shift == -2 && k == 1 ) { 
    v = 3.*f_k_moins_2[1] - 16.*f_k_moins_2[2] + 
      36.*f_k_moins_2[3] - 48.*f_k_moins_2[4] + 25.*f_k_moins_2[5];
  }
  return v / 12.;
}

//--------------------------------------------------------
//          cyl_bilinear_interp()
//
// generic bilinear interpolation formula for a point (x,y) of [0,1]^2
// in the square [i,i+1]*[j,j+1] in the m grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. 
//--------------------------------------------------------

template < typename Point_type, typename Num_type > 
Point_type
cyl_bilinear_interp(Num_type x, Num_type y,    // coordinates of the point in the
		    int i, int j,              // square of lower left vertex (i,j)
		    MAT_matrix_2<Point_type> const & m )     
{
  int nx = m.nb_row();
  
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;

  Point_type f_ll = m( i, j ), 
             f_rl = m( i_plus_1, j ), 
             f_lu = m( i, j_plus_1 ), 
             f_ru = m( i_plus_1, j_plus_1 );

  Point_type temp_l = f_ll + x * ( f_rl - f_ll );
  Point_type temp_u = f_lu + x * ( f_ru - f_lu );

  return temp_l + y * ( temp_u - temp_l );
}

//--------------------------------------------------------
//          cyl_bicubic_interp()
//
// generic bicubic interpolation formula for a point (x,y) of [0,1]^2
// in the square [i,i+1]*[j,j+1] in the m grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. 
// The formulae is based on a bicubic interpolation with C-splines.
//--------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Point_type
cyl_bicubic_interp(Num_type x, Num_type y,// coordinates of the point in the
		   int i, int j,   // square of lower left vertex (i,j)
		   MAT_matrix_2<Point_type> const & m ) 
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1;
  assert( m.nb_col() > 5 && j < ny );
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type h_01_x = x*x*(3 - 2*x);  // Hermite polynomials
  Num_type h_10_x = x*(x*(x - 2) + 1); 
  Num_type h_11_x = x*x*(x - 1); 
  
  Vector_type f_j_moins_2[6]; // Point -> Vector to allow
                              // affine combinations.
  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 for derivative approximation
    Vector_type v_i = ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
                 8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= 12.;

    Vector_type v_i_plus_1 = 
      ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
      8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= 12.;

    // instead of calling c_spline_interp we use the computations on
    // Hermite polynoms.
    f_j_moins_2[k] =  ( m(i, jk_moins_2) - Point_type(0,0,0) ) + 
      h_01_x * ( m(i_plus_1, jk_moins_2) - m(i, jk_moins_2) ) +
      h_10_x * v_i + h_11_x * v_i_plus_1;
  }

  Vector_type v_j = finite_diff_derivative( f_j_moins_2, j_aux - j, 0 );
  Vector_type v_j_plus_1 = finite_diff_derivative( 
					    f_j_moins_2, j_aux - j, 1 );
  Vector_type f_j = f_j_moins_2[2 + j - j_aux];
  Vector_type f_j_plus_1 = f_j_moins_2[3 + j - j_aux];

  return Point_type(0,0,0) + c_spline_interp( y, 
					      f_j, 
					      f_j_plus_1, 
					      v_j, 
					      v_j_plus_1);
}

//--------------------------------------------------------
//          cyl_interp()
//
// generic interpolation formula for a point (x,y) of [0,1]^2
//replaced in grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. 
//--------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Point_type
cyl_generic_interp( Num_type x, 
		    Num_type y, 
		    MAT_matrix_2<Point_type> const & m,
		    int degree ) // The interpolation degree
{
  // x and y are supposed to be in [0,1[. They are actually of the
  // form t - floor(t) for some t when cyl_interp is called. It seems
  // that on some processors we may have t - floor(t) >= 1. We perform
  // some preprocessing to ensure x,y in [0,1[.

  if (x < 0) x += 1;
  if (x >= 1) x -= 1;
  assert( x >= 0 && x < 1);
  assert( y >= 0 && y <= 1);

  int nx = m.nb_row();
  int ny = m.nb_col() - 1;  // Caution : the cylinder's size is nx*(ny+1) !!
  Num_type x_scaled = x * nx;
  Num_type y_scaled = y * ny;
  
  int i = floor( x_scaled );
  int j = floor( y_scaled );

  Num_type u = x_scaled - i;
  Num_type v = y_scaled - j;

  if (j == ny) { 
    assert( v == 0. );
    v = 1;
    j = ny - 1;
  }

  if (degree == 1)
    return cyl_bilinear_interp<Point_type, Num_type> (u, v, i, j, m);
  else
    return cyl_bicubic_interp<Point_type, Vector_type, Num_type> 
      (u, v, i, j, m);
}


//------------------------------------------------------------
//          cyl_bilinear_deriv_x()
//
// d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bilinear interpolation. 
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bilinear_deriv_x(Num_type x, 
		     Num_type y,   // coordinates of the point in the
		     int i, int j, // square of lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1;
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;
  assert( j < ny );

  // Multiply by nx to take reparameterization into account
  return (Num_type)nx * ( (1.- y)*(m(i_plus_1,j) - m(i,j)) +
			  y * (m(i_plus_1,j_plus_1) - m(i,j_plus_1)) );  
}


//------------------------------------------------------------
//          cyl_bicubic_deriv_x()
//
// d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in grid of size nb_col * nb_row. The formula is a 
// simple derivation of the bicubic interpolation. 
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bicubic_deriv_x( Num_type x, 
		     Num_type y, // coordinates of the point in the
		     int i, int j,   // square of lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1; // Caution : the cylinder's size is nx*(ny+1) !!
  assert( m.nb_col() > 5 );
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type dh_01_dx = 6*x*(1-x);  // derivatives of Hermite polynomials
  Num_type dh_10_dx = (x-1)*(3*x - 1); 
  Num_type dh_11_dx = x*(3*x - 2); 
  
  Vector_type df_dx_j_moins_2[6];

  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 derivative approximation
    Vector_type v_i = ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
         (Num_type)8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= (Num_type)12.;

    Vector_type v_i_plus_1 = 
           ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
      (Num_type)8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= (Num_type)12.;

    df_dx_j_moins_2[k] =  
      dh_01_dx * ( m(i_plus_1, jk_moins_2) - m(i, jk_moins_2) ) +
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
//          cyl_generic_deriv_x()
//
// generic d/dx derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_generic_deriv_x( Num_type x, 
		     Num_type y, 
		     MAT_matrix_2<Point_type> const & m,
		     int degree = CYL_deriv_order )
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
    return cyl_bilinear_deriv_x<Point_type, Vector_type, Num_type> 
                                                       ( u, v, i, j, m );
  else 
    return cyl_bicubic_deriv_x<Point_type, Vector_type, Num_type>
                                                      ( u, v, i, j, m );
}

//------------------------------------------------------------
//          cyl_bilinear_deriv_y()
//
// generic d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bilinear_deriv_y(Num_type x, 
		     Num_type y,   
		     int i, int j, 
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1;
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = j + 1;
  assert( j < ny );

  // Multiply by nx to take reparameterization into account
  return (Num_type)ny * ( (1.- x)*(m(i,j_plus_1) - m(i,j)) +
			  x * (m(i_plus_1,j_plus_1) - m(i_plus_1,j)) );  
}

//------------------------------------------------------------
//          cyl_bicubic_deriv_y()
//
// generic d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row. This formula is a 
// simple derivation of the bicubic interpolation.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_bicubic_deriv_y( Num_type x, 
		     Num_type y, 
		     int i, int j,  
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col() - 1; 
  assert( m.nb_col() > 5 );
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type h_01_x = x*x*(3 - 2*x);  // Hermite polynomials
  Num_type h_10_x = x*(x*(x - 2) + 1); 
  Num_type h_11_x = x*x*(x - 1); 
  
  Vector_type f_j_moins_2[6]; // Vector plutôt que Point pour les combinaisons
                              // affines.
  int j_aux = min(max(j, 2), ny - 3);

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = j_aux + k - 2;

    // Uses order 4 approximation for derivative
    Vector_type v_i = ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
                 8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= 12.;

    Vector_type v_i_plus_1 = 
           ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
      8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= 12.;

    // plutôt que d'appeler c_spline_interp on réutilise les calculs
    // des polynômes de Hermite.
    f_j_moins_2[k] =  (m(i, jk_moins_2) - Point_type(0,0,0)) + 
      h_01_x * ( m(i_plus_1, jk_moins_2) - m(i, jk_moins_2) ) +
      h_10_x * v_i + h_11_x * v_i_plus_1;
  }

  Vector_type v_j = finite_diff_derivative( f_j_moins_2, j_aux - j, 0 );
  Vector_type v_j_plus_1 = finite_diff_derivative( 
					    f_j_moins_2, j_aux - j, 1 );
  Vector_type f_j = f_j_moins_2[2 + j - j_aux];
  Vector_type f_j_plus_1 = f_j_moins_2[3 + j - j_aux];

  Num_type dh_01_dy = 6*y*(1-y);  // derivatives of Hermite polynomials
  Num_type dh_10_dy = (y-1)*(3*y - 1); 
  Num_type dh_11_dy = y*(3*y - 2); 

  // Multiply by ny to take reparameterization in to account
  return (Num_type)ny * (dh_01_dy*( f_j_plus_1 - f_j ) + 
			 dh_10_dy*v_j + dh_11_dy*v_j_plus_1);
}

//------------------------------------------------------------
//          cyl_generic_deriv_y()
//
// generic d/dy derivation formula for a point (x,y) of [0,1]^2
// replaced in a grid of size nb_col * nb_row.
// The grid is supposed to be periodic in x. 
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
cyl_generic_deriv_y( Num_type x, 
		     Num_type y, 
		     MAT_matrix_2<Point_type> const & m,
		     int degree = CYL_deriv_order )
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
    return cyl_bilinear_deriv_y< Point_type, Vector_type, Num_type> ( 
							     u, v, i, j, m);
  else 
    return cyl_bicubic_deriv_y<Point_type, Vector_type, Num_type>( 
							     u, v, i, j, m );
}

#endif // _CYL_INTERP_H_
