//************************************************************************
//                                                    TOR_interp.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : Dec 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains interpolation functions for bidimentional and
// biperiodic grids (on a torus).
//
// =======================================================

#ifndef _TOR_INTERP_H_
#define _TOR_INTERP_H_

//----------------------------------------------------------------------//
//  Fonction c_spline_interp()
//
// Interpolation with C-spline basis.
//----------------------------------------------------------------------//

template < typename Point_type, typename Vector_type, typename Num_type > 
Point_type c_spline_interp( Num_type const t, 
			    Point_type const & p0, 
			    Point_type const & p1,
			    Vector_type const & v0,
			    Vector_type const & v1 ) 
{
  Num_type h_01 = t*t*(3 - 2*t);  // Hermite polynomials
  Num_type h_10 = t*(t*(t - 2) + 1); 
  Num_type h_11 = t*t*(t - 1); 
 

  return p0 + h_01 * (p1 -p0) + h_10 * v0 + h_11 * v1;
}

//--------------------------------------------------------
//          tor_bilinear_interp()
//
// generic bilinear interpolation formula for a point 
// (x,y) of [0,1]^2 in the square  [i,i+1]*[j,j+1]
// in the m grid of size  nb_col * nb_row.
// The grid is supposed to be doubly periodic.
//--------------------------------------------------------

template < typename Point_type, typename Num_type > 
Point_type
tor_bilinear_interp( Num_type x, Num_type y, // coordinates of the point in the  
		     int i, int j,   // square with lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m ) 
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = (j + 1) % ny;

  Point_type f_ll = m( i, j ), 
             f_rl = m( i_plus_1, j ), 
             f_lu = m( i, j_plus_1 ), 
             f_ru = m( i_plus_1, j_plus_1 );

  Point_type temp_l = f_ll + x * ( f_rl - f_ll );
  Point_type temp_u = f_lu + x * ( f_ru - f_lu );

  return temp_l + y * ( temp_u - temp_l );
}

//--------------------------------------------------------
//          tor_bicubic_interp()
//
// generic bicubic interpolation formula for a point 
// (x,y) of [0,1]^2 in the square  [i,i+1]*[j,j+1]
// in the m grid of size  nb_col * nb_row.
// The grid is supposed to be doubly periodic. The formulae is 
// using a bicubic interpolation based on C-splines.
//--------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Point_type
tor_bicubic_interp( Num_type x, Num_type y, // coordinates of the point in the   
		    int i, int j,    // square with lower left vertex (i,j)
		    MAT_matrix_2<Point_type> const & m ) 
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type h_01_x = x*x*(3 - 2*x);  // Hermite polynomials
  Num_type h_10_x = x*(1-x)*(1-x); 
  Num_type h_11_x = x*x*(x - 1); 
  
  Point_type f_j_moins_2[6];

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = (j + k + ny - 2) % ny;

    // Uses order 4 derivative approximation
    Vector_type v_i = ( m(i_moins_2, jk_moins_2) - m(i_plus_2, jk_moins_2) ) +
       (Num_type)8. * ( m(i_plus_1, jk_moins_2) - m(i_moins_1, jk_moins_2) );
    v_i /= (Num_type)12.;

    Vector_type v_i_plus_1 = 
           ( m(i_moins_1, jk_moins_2) - m(i_plus_3, jk_moins_2) ) +
      (Num_type)8. * ( m(i_plus_2, jk_moins_2) - m(i, jk_moins_2) );
    v_i_plus_1 /= (Num_type)12.;

    f_j_moins_2[k] =  m(i, jk_moins_2) + 
      h_01_x * ( m(i_plus_1, jk_moins_2) - m(i, jk_moins_2) ) +
      h_10_x * v_i + h_11_x * v_i_plus_1;
  }

  // Uses order 4 derivative approximation
  Vector_type v_j = (f_j_moins_2[0] - f_j_moins_2[4]) +
     (Num_type)8. * (f_j_moins_2[3] - f_j_moins_2[1]);
  v_j /= (Num_type)12.;

  Vector_type v_j_plus_1 = (f_j_moins_2[1] - f_j_moins_2[5]) +
            (Num_type)8. * (f_j_moins_2[4] - f_j_moins_2[2]);
  v_j_plus_1 /= (Num_type)12.;

  return c_spline_interp(y, f_j_moins_2[2], f_j_moins_2[3], v_j, v_j_plus_1);
}

//--------------------------------------------------------
//          tor_generic_interp()
//
// generic interpolation formula for a point (x,y)
// of [0,1]^2 replaced a grid of size nb_col * nb_row.
// The grid is supposed to be doubly periodic. 
//--------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Point_type
tor_generic_interp( Num_type x, 
		    Num_type y, 
		    MAT_matrix_2<Point_type> const & m,
		    int degree = TOR_interp_order ) // The interpolation degree
{
  // x and y are supposed to be in [0,1[. They are actually set to some 
  // t - floor(t). It seems that on some processors we may have 
  // t - floor(t) >= 1. We perform some preprocessing to ensure x,y in [0,1[.

  if (x < 0) x += 1;
  if (x >= 1) x -= 1;
  if (y < 0) y += 1;
  if (y >= 1) y -= 1;

  assert( x >= 0 && x <= 1 && y >= 0 && y <= 1 );

  int nx = m.nb_row();
  int ny = m.nb_col();
  Num_type x_scaled = x * nx;
  Num_type y_scaled = y * ny;
  
  int i = floor( x_scaled );
  int j = floor( y_scaled );

  if (degree == 1) return tor_bilinear_interp<Point_type, Num_type> (
					  x_scaled - i, y_scaled - j, i, j, m);
  else return tor_bicubic_interp<Point_type, Vector_type, Num_type>(
			                 x_scaled - i, y_scaled - j, i, j, m );
}

//------------------------------------------------------------
//          tor_bilinear_deriv_x()
//
// derivation formulae d/dx for a point (x,y) in [0,1]^2
// replaced in a a grid of size nb_col * nb_row.
// This formulae is a simple derivation of the bilinear interpolation.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_bilinear_deriv_x(Num_type x, 
		     Num_type y, 
		     int i, int j,   // square with lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = (j + 1) % ny;

  // Multiply by nx to take reparameterization into account
  return (Num_type)nx * ( (1.- y)*(m(i_plus_1,j) - m(i,j)) +
			  y * (m(i_plus_1,j_plus_1) - m(i,j_plus_1)) );  
}

//------------------------------------------------------------
//          tor_bilinear_deriv_y()
//
// derivation formulae d/dy for a point (x,y) in [0,1]^2
// replaced in a a grid of size nb_col * nb_row.
// This formulae is a simple derivation of the bilinear interpolation.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_bilinear_deriv_y(Num_type x, 
		     Num_type y, 
		     int i, int j,   // square with lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  int i_plus_1 = (i + 1) % nx;
  int j_plus_1 = (j + 1) % ny;

  // Multiply by nx to take reparameterization into account
  return (Num_type)ny * ( (1.- x)*(m(i,j_plus_1) - m(i,j)) +
			  x * (m(i_plus_1,j_plus_1) - m(i_plus_1,j)) );  
}

//------------------------------------------------------------
//          tor_bicubic_deriv_x()
//
// derivation formulae d/dx for a point (x,y) in [0,1]^2
// replaced in a a grid of size nb_col * nb_row.
// This formulae is a simple derivation of the bicubic interpolation.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_bicubic_deriv_x( Num_type x, 
		     Num_type y, 
		     int i, int j,   // square with lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  
  int i_moins_2 = (i + nx - 2) % nx;
  int i_moins_1 = (i_moins_2 + 1) % nx;
  int i_plus_1 = (i + 1) % nx;
  int i_plus_2 = (i_plus_1 + 1) % nx;
  int i_plus_3 = (i_plus_2 + 1) % nx;

  Num_type dh_01_dx = 6*x*(1-x);  // derivatives of Hermite polynomials
  Num_type dh_10_dx = (x-1)*(3*x - 1); 
  Num_type dh_11_dx = x*(3*x - 2); 
  
  Vector_type df_dx_j_moins_2[6];

  for (int k=0; k < 6; k++) {
    int jk_moins_2 = (j + k + ny - 2) % ny;

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

  // Uses order 4 derivative approximation
  Vector_type v_j = (df_dx_j_moins_2[0] - df_dx_j_moins_2[4]) +
     (Num_type)8. * (df_dx_j_moins_2[3] - df_dx_j_moins_2[1]);
  v_j /= (Num_type)12.;

  Vector_type v_j_plus_1 = (df_dx_j_moins_2[1] - df_dx_j_moins_2[5]) +
            (Num_type)8. * (df_dx_j_moins_2[4] - df_dx_j_moins_2[2]);
  v_j_plus_1 /= (Num_type)12.;

  // Multiply by nx to take reparameterization into account
  return (Num_type)nx * c_spline_interp( y, 
					 df_dx_j_moins_2[2], 
					 df_dx_j_moins_2[3], 
					 v_j, 
					 v_j_plus_1);
}

//------------------------------------------------------------
//          tor_bicubic_deriv_y()
//
// derivation formulae d/dy for a point (x,y) in [0,1]^2
// replaced in a a grid of size nb_col * nb_row.
// This formulae is a simple derivation of the bicubic interpolation.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_bicubic_deriv_y( Num_type x, 
		     Num_type y, 
		     int i, int j,   // square with lower left vertex (i,j)
		     MAT_matrix_2<Point_type> const & m )
{
  int nx = m.nb_row();
  int ny = m.nb_col();
  
  int j_moins_2 = (j + ny - 2) % ny;
  int j_moins_1 = (j_moins_2 + 1) % ny;
  int j_plus_1 = (j + 1) % ny;
  int j_plus_2 = (j_plus_1 + 1) % ny;
  int j_plus_3 = (j_plus_2 + 1) % ny;

  Num_type dh_01_dy = 6*y*(1-y);  // derivatives of Hermite polynomials
  Num_type dh_10_dy = (y-1)*(3*y - 1); 
  Num_type dh_11_dy = y*(3*y - 2); 
  
  Vector_type df_dy_i_moins_2[6];

  for (int k=0; k < 6; k++) {
    int ik_moins_2 = (i + k + nx - 2) % nx;

    // Uses order 4 derivative approximation
    Vector_type v_j = ( m(ik_moins_2, j_moins_2) - m(ik_moins_2, j_plus_2) ) +
         (Num_type)8. * ( m(ik_moins_2, j_plus_1) - m(ik_moins_2, j_moins_1) );
    v_j /= (Num_type)12.;

    Vector_type v_j_plus_1 = 
           ( m(ik_moins_2, j_moins_1) - m(ik_moins_2, j_plus_3) ) +
      (Num_type)8. * ( m(ik_moins_2, j_plus_2) - m(ik_moins_2, j) );
    v_j_plus_1 /= (Num_type)12.;

    df_dy_i_moins_2[k] =  
      dh_01_dy * ( m(ik_moins_2, j_plus_1) - m(ik_moins_2, j) ) +
      dh_10_dy * v_j + dh_11_dy * v_j_plus_1;
  }

  // Uses order 4 derivative approximation
  Vector_type v_i = (df_dy_i_moins_2[0] - df_dy_i_moins_2[4]) +
               (Num_type)8. * (df_dy_i_moins_2[3] - df_dy_i_moins_2[1]);
  v_i /= (Num_type)12.;

  Vector_type v_i_plus_1 = (df_dy_i_moins_2[1] - df_dy_i_moins_2[5]) +
                  (Num_type)8. * (df_dy_i_moins_2[4] - df_dy_i_moins_2[2]);
  v_i_plus_1 /= (Num_type)12.;

  // Multiply by ny to take reparameterization into account
  return (Num_type)ny * c_spline_interp( x, 
					 df_dy_i_moins_2[2], 
					 df_dy_i_moins_2[3], 
					 v_i, 
					 v_i_plus_1 );
}

//------------------------------------------------------------
//          tor_generic_deriv_x()
//
// generic  d/dx derivation formulae for a point (x,y) of 
// [0,1]^2 replaced in a grid of size nb_col * nb_row.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_generic_deriv_x( Num_type u, 
		     Num_type v, 
		     MAT_matrix_2<Point_type> const & m,
		     int degree )
{
  // x and y are supposed to be in [0,1[. They are actually set to some 
  // t - floor(t). It seems that on some processors we may have 
  // t - floor(t) >= 1. We perform some preprocessing to ensure x,y in [0,1[.

  if (u < 0) u += 1;
  if (u >= 1) u -= 1;
  if (v < 0) v += 1;
  if (v >= 1) v -= 1;

  assert( u >= 0 && u <= 1 && v >= 0 && v <= 1 );

  int nx = m.nb_row();
  int ny = m.nb_col();
  Num_type u_scaled = u * nx;
  Num_type v_scaled = v * ny;
  
  int i = floor( u_scaled );
  int j = floor( v_scaled );

  if (degree == 1) 
    return tor_bilinear_deriv_x<Point_type, Vector_type, Num_type> (
					  u_scaled - i, v_scaled - j, i, j, m);
  else 
    return tor_bicubic_deriv_x<Point_type, Vector_type, Num_type>(
			                 u_scaled - i, v_scaled - j, i, j, m );
}

//------------------------------------------------------------
//          tor_generic_deriv_y()
//
// generic  d/dy derivation formulae for a point (x,y) of 
// [0,1]^2 replaced in a grid of size nb_col * nb_row.
//------------------------------------------------------------

template < typename Point_type, typename Vector_type, typename Num_type > 
Vector_type
tor_generic_deriv_y( Num_type u, 
		     Num_type v, 
		     MAT_matrix_2<Point_type> const & m,
		     int degree )
{
  // x and y are supposed to be in [0,1[. They are actually set to some 
  // t - floor(t). It seems that on some processors we may have 
  // t - floor(t) >= 1. We perform some preprocessing to ensure x,y in [0,1[.

  if (u < 0) u += 1;
  if (u >= 1) u -= 1;
  if (v < 0) v += 1;
  if (v >= 1) v -= 1;

  assert( u >= 0 && u <= 1 && v >= 0 && v <= 1 );

  int nx = m.nb_row();
  int ny = m.nb_col();
  Num_type u_scaled = u * nx;
  Num_type v_scaled = v * ny;
  
  int i = floor( u_scaled );
  int j = floor( v_scaled );

  if (degree == 1) return 
		     tor_bilinear_deriv_y<Point_type, Vector_type, Num_type> (
					  u_scaled - i, v_scaled - j, i, j, m);
  else return tor_bicubic_deriv_y<Point_type, Vector_type, Num_type>(
			                 u_scaled - i, v_scaled - j, i, j, m );
}

#endif // _TOR_INTERP_H_
