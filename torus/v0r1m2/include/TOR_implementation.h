//************************************************************************
//                                               TOR_implementation.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : July 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains the implementation of the class members of 
// various TOR_ functions, vector fields, 2forms,...
//
// =======================================================

#ifndef _TOR_IMPLEMENTATION_H_
#define _TOR_IMPLEMENTATION_H_

#include <TOR_torus.h>
#include <cassert>


//--------------------------------------------------------
//          TOR_function<Num_type>::operator/=
//
// Divides by another function.
//--------------------------------------------------------

template < typename Num_type > TOR_function<Num_type> &
TOR_function <Num_type>::
operator/=( TOR_function<Num_type> const & quotient )
{
  assert( ( quotient.n_x() == this->n_x() ) && 
	  ( quotient.n_y() == this->n_y() ) );

  this->m_ /= quotient.m_;

  return *this;
}

//--------------------------------------------------------
//          TOR_function<Num_type>::operator/
//
// Divise une fonction par une autre.
//--------------------------------------------------------

template < typename Num_type > TOR_function<Num_type>
TOR_function <Num_type>::
operator/ ( TOR_function<Num_type> const & quotient ) const
{
  assert( ( quotient.n_x() == this->n_x() ) && 
	  ( quotient.n_y() == this->n_y() ) );

  return m_ / quotient.m_;
}

//--------------------------------------------------------
//          TOR_function<Num_type>::operator-
//
// Change the signe of a function.
//--------------------------------------------------------

template < typename Num_type > TOR_function<Num_type> &
TOR_function <Num_type>::
operator-()
{
  this->m_ = - this->m_;
  return *this;
}

//--------------------------------------------------------
//          TOR_function<Num_type>::operator+=
//
// Add another function.
//--------------------------------------------------------
template < typename Num_type > TOR_function<Num_type> &
TOR_function <Num_type>::
operator+=( TOR_function<Num_type> const & other )
{
  assert( ( other.n_x() == n_x() ) && ( other.n_y() == n_y() ) );

  m_ += other.m_;
  return *this;
}


//-------------------------------------------------------------
//          TOR_2vector_field<Num_type>::apply()
//
// Apply the vector field (X = *this) to an embedding f.
// Returns the field (TOR_3vector_field) df.X.
//-------------------------------------------------------------

template < typename Num_type > TOR_3vector_field <Num_type> &
TOR_2vector_field <Num_type>::
apply( TOR_embedding<Num_type> const & f, 
       TOR_3vector_field<Num_type> & out_field ) const
{
  int nx( this->n_x() ), ny( this->n_y() );

  assert( ( f.n_x() == nx ) && ( out_field.n_x() == nx ) &&
	  ( f.n_y() == ny ) && ( out_field.n_y() == ny ) );
  // relies on a simple approximation of the jacobian of f
  // where df/dx is obtained by finite differences of order 2 or 4

  Num_type nx_prime(nx), ny_prime(ny);
  if (TOR_deriv_order == 2) {   // The deriv_x,y formula should be normalized
    nx_prime /= 2.;
    ny_prime /= 2.;
  }
  else {
    assert( TOR_deriv_order == 4);
    nx_prime /= 12.;
    ny_prime /= 12.;
  }
  
  TOR_3vector_field<Num_type> tmp_field(nx, ny); // for the calculation of derivatives

  if ( ! this->is_constant_ ) {
    out_field = f.deriv_x(out_field);
    tmp_field = f.deriv_y(tmp_field);
	
    for( int i = 0; i < nx; i++ )
      for( int j = 0; j < ny; j++ ) {
	out_field.mat(i,j) *= (this->mat(i,j)).x() * nx_prime; 
	out_field.mat(i,j) += ((this->mat(i,j)).y() * ny_prime) * 
	                      tmp_field.mat(i,j); 
      }
  }
  else { // this is a constant field
    out_field = UTI_3vector<Num_type>(0,0,0); // initialisation to 0;

    if (const_vec_[0] != 0) {
      out_field = f.deriv_x(out_field);
	out_field *= const_vec_.x() * nx_prime;
    }
    if (const_vec_[1] != 0) {
      tmp_field = f.deriv_y(tmp_field);
      tmp_field *= const_vec_.y() * ny_prime;
      out_field += tmp_field;
    }
  }

  return out_field;
}

//-------------------------------------------------------------
//          TOR_sym_2form_field<Num_type>:: decompose_alpha()
//
// returns the alpha coefficient for the decomposition in  
// primitive metrics.
//-------------------------------------------------------------
template < typename Num_type >   TOR_function<Num_type> & 
TOR_sym_2form_field<Num_type>:: 
decompose_alpha( UTI_3vector<Num_type> const & alpha_filter, 
		 TOR_function<Num_type> & alpha ) const
{
  int nx( n_x() ), ny( n_y() );
  assert( (alpha.n_x() == nx) && (alpha.n_y() == ny) );

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      alpha.mat(i,j) = mat(i,j) * alpha_filter;

  return alpha;
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>:: TOR_embedding()
//
// Constructor of TOR_embedding.
// 
//-------------------------------------------------------------

template < typename Num_type > void
TOR_embedding <Num_type>::
init_standard_( Num_type r1, Num_type r2 )
{
  assert( r1 + r2 < 1 && r1 > 0 && r2 > 0 ); // The embedding is strictly short !

  int nx( n_x() ), ny( n_y() );

  // We make a two passes computations (lines-columns) then
  // tcolumns-lines) to avoid to calculate the same sine and
  // cosine several times.
  
  for (int j = 0; j < ny; j++) {

    Num_type cos_2piy = cos(j*2*M_PI/ny);
    Num_type sin_2piy = sin(j*2*M_PI/ny);

    for (int i = 0; i < nx; i++) {

      m_(i,j).x() = cos_2piy;
      m_(i,j).y() = sin_2piy;
    }
  }
  for (int i = 0; i < nx; i++) {

    Num_type r2_plus_r1_cos_2pix = (r2 + r1*cos(i*2.*M_PI/nx)) * M_1_PI/2.;
    Num_type r1_sin_2pix = r1 * sin(i*2.*M_PI/nx) * M_1_PI / 2.;

    for (int j = 0; j < ny; j++) {

      m_(i,j).x() *= r2_plus_r1_cos_2pix;
      m_(i,j).y() *= r2_plus_r1_cos_2pix;
      m_(i,j).z() = r1_sin_2pix;
    }
  }
}
 
//-------------------------------------------------------------
//          TOR_embedding<Num_type>:: init_standard_bis_()
//
// Constructor BIS of TOR_embedding. Exchange of x and y 
// in the standard embedding.
// 
//-------------------------------------------------------------

template < typename Num_type > void
TOR_embedding <Num_type>::
init_standard_bis_( Num_type r1, Num_type r2 )
{
  assert( r1 + r2 < 1 && r1 > 0 && r2 > 0 ); // The embedding is strictly !

  int nx( n_x() ), ny( n_y() );

  // We make a two passes computations (lines-columns) then 
  // (columns-lines) to avoid to calculate the same sine and
  // cosine several times.
  
  for (int i = 0; i < nx; i++) {

    Num_type cos_2pix = cos(i*2*M_PI/nx);
    Num_type sin_2pix = sin(i*2*M_PI/nx);

    for (int j = 0; j < ny; j++) {

      m_(i,j).x() = cos_2pix;
      m_(i,j).y() = sin_2pix;
    }
  }
  for (int j = 0; j < ny; j++) {

    Num_type r2_plus_r1_cos_2piy = (r2 + r1*cos(j*2.*M_PI/ny)) * M_1_PI/2.;
    Num_type r1_sin_2piy = r1 * sin(j*2.*M_PI/ny) * M_1_PI / 2.;

    for (int i = 0; i < nx; i++) {

      m_(i,j).x() *= r2_plus_r1_cos_2piy;
      m_(i,j).y() *= r2_plus_r1_cos_2piy;
      m_(i,j).z() = r1_sin_2piy;
    }
  }
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::normal_field()
//
// Returns the normals of the embedding (TOR_3vector_field).
// i.e. n = df/dx ^ df/dy / || df/dx ^ df/dy ||
//-------------------------------------------------------------
template < typename Num_type > TOR_3vector_field<Num_type> &
TOR_embedding <Num_type>::
normal_field( TOR_3vector_field<Num_type> & normals ) const
{
  int nx( this->n_x() ), ny( this->n_y() );

  assert( ( normals.n_x() == nx ) && ( normals.n_y() == ny ) );

  TOR_3vector_field<Num_type> df_x(nx, ny), df_y(nx, ny);
  df_x = this->deriv_x( df_x );
  df_y = this->deriv_y( df_y );

  normals = df_x.cross_prod_mat( df_y );

  normals /= normals.norm();
  return normals;
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::normal()
//
// Returns the embedding's normal at the point (x,y)
// i.e. n = df/dx ^ df/dy / || df/dx ^ df/dy ||
//-------------------------------------------------------------
template < typename Num_type > UTI_3vector<Num_type>
TOR_embedding <Num_type>::normal( Num_type x, Num_type y ) const
{
  UTI_3vector<Num_type> vx = deriv_x(x,y);
  UTI_3vector<Num_type> vy = deriv_y(x,y);

  UTI_3vector<Num_type> n = cross_prod( vx, vy );
  return n.normalize();
} 

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::pull_back()
//
// Returns the pullback of the embedding.
//-------------------------------------------------------------
template < typename Num_type > TOR_sym_2form_field<Num_type> &
TOR_embedding <Num_type>::
pull_back( TOR_sym_2form_field<Num_type> & g ) const
{
  int nx( this->n_x() ), ny( this->n_y() );
  
  assert( ( g.n_x() == nx ) && ( g.n_y() == ny ) );

  TOR_3vector_field<Num_type> df_x(nx, ny), df_y(nx, ny);
  df_x = this->deriv_x( df_x );
  df_y = this->deriv_y( df_y );

  Num_type nx_prime(nx), ny_prime(ny);
  if (TOR_deriv_order == 2) {   // The deriv_x,y formula should be normalized
    nx_prime /= 2.;
    ny_prime /= 2.;
  }
  else {
    assert(TOR_deriv_order == 4);
    nx_prime /= 12.;
    ny_prime /= 12.;
  }
  

  for( int i = 0; i < nx; i++ ) {
    for( int j = 0; j < ny; j++ ) {
      UTI_3vector<Num_type> vx, vy;

      vx = nx_prime * df_x.mat(i,j);
      vy = ny_prime * df_y.mat(i,j);

      g.E(i,j) = vx * vx;
      g.F(i,j) = vx * vy;
      g.G(i,j) = vy * vy;
    }
  }
  return g;
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::pull_back()
//
// Returns the pullback of the embedding at (x,y).
//-------------------------------------------------------------
template < typename Num_type > UTI_sym_2form<Num_type>
TOR_embedding <Num_type>::pull_back( Num_type x, Num_type y ) const
{
  UTI_3vector<Num_type> vx = deriv_x(x,y);
  UTI_3vector<Num_type> vy = deriv_y(x,y);

  return UTI_sym_2form<Num_type>( vx * vx, vx * vy, vy * vy );
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::deriv_x()
//
// Derives in the x direction. Uses a finite differences formula. 
// The sampling step is supposed to be equal to 1/2.
//
//  If order  = 2 : df/dx(i) = f(i+1)-f(i-1)
//  If order  = 4 : df/dx(i) = 8*(f(i+1)-f(i-1)) + (f(i-2) - f(i+2))
//-------------------------------------------------------------
template < typename Num_type > TOR_3vector_field<Num_type> &
TOR_embedding<Num_type>::
deriv_x( TOR_3vector_field<Num_type> & df_x ) const
{
  int nx(n_x());

  if ( TOR_deriv_order == 2 ) {
    for (int i = 0; i < nx; i++) {
      int i_pred(i-1), i_succ(i+1);
      if (i==0) i_pred = nx-1;
      else if (i==nx-1) i_succ = 0;

      df_x.m_(i, Range::all()) = substract( m_(i_succ, Range::all()),
					    m_(i_pred, Range::all()) );
    }
  }
  else if ( TOR_deriv_order == 4 ) {
    for (int i = 0; i < nx; i++) {
      int i_moins_2 = (nx + i - 2) % nx;
      int i_moins_1 = (i_moins_2 + 1) % nx;
      int i_plus_1 = (i + 1) % nx;
      int i_plus_2 = (i + 2) % nx;

      df_x.m_(i, Range::all()) = 
	(Num_type)8. * substract( m_(i_plus_1, Range::all()), 
				  m_(i_moins_1, Range::all())) + 
	substract(m_(i_moins_2, Range::all()),
		  m_(i_plus_2, Range::all()));
    }
  }
  return df_x;
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::deriv_y()
//
// Derives in the y direction. Uses a finite differences formula. 
// The sampling step is supposed to be equal to 1/2.
//
//  If ordre  = 2 : df/dy(i) = f(i+1)-f(i-1)
//  If ordre  = 4 : df/dy(i) = 8*(f(i+1)-f(i-1)) - (f(i+2) - f(i-2)
//-------------------------------------------------------------
template < typename Num_type > TOR_3vector_field<Num_type> &
TOR_embedding<Num_type>::
deriv_y( TOR_3vector_field<Num_type> & df_y ) const
{
  int ny(n_y());

  if ( TOR_deriv_order == 2 ) {
    for (int i = 0; i < ny; i++) {
      int i_pred(i-1), i_succ(i+1);
      if (i==0) i_pred = ny-1;
      else if (i==ny-1) i_succ = 0;

      df_y.m_(Range::all(), i) = substract( m_(Range::all(), i_succ),
					    m_(Range::all(), i_pred) );
    }
  }
  else if ( TOR_deriv_order == 4 ) {
    for (int i = 0; i < ny; i++) {
      int i_moins_2 = (ny + i - 2) % ny;
      int i_moins_1 = (i_moins_2 + 1) % ny;
      int i_plus_1 = (i + 1) % ny;
      int i_plus_2 = (i + 2) % ny;

      df_y.m_(Range::all(), i) = 
	(Num_type)8. * substract( m_(Range::all(), i_plus_1), 
				  m_(Range::all(), i_moins_1) ) + 
	substract( m_(Range::all(), i_moins_2), 
		   m_(Range::all(), i_plus_2) );
    }
  }
  return df_y;
}

//-------------------------------------------------------------
//  TOR_embedding<Num_type>::resample()
//
// 
// Resample the grid with nx' = floor(nx*ratio_x)
// ny' = floor(ny*ratio_y). The grid is not anymore squared !
//-------------------------------------------------------------
template <typename Num_type> void TOR_embedding<Num_type>::
resample( Num_type ratio_x, Num_type ratio_y )
{
  int nx = floor(n_x()*ratio_x);
  int ny = floor(n_y()*ratio_y);
  TOR_embedding resampled(nx, ny);

  int i;
  // #pragma omp parallel for shared(reranged) private(i)
  for (i = 0; i < nx; i++) {
    for (int j = 0; j <= ny; j++)
      resampled.mat(i, j) = 
	(*this)((Num_type)i/(Num_type)nx, (Num_type)j/(Num_type)ny);
  }
  *this = resampled;
}

#endif // _TOR_IMPLEMENTATION_H_
