//************************************************************************
//                                                    CYL_cylinder.h
//                                                  -------------------
//    written by         : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date                 : May 2009
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains class definitions for functions and embeddings
// on a cylinder. The parameter space is stored
// into a bidimentional matrix corresponding to a
// col*row uniform sampling of the canonical cylinder.
//
// =======================================================

#ifndef _CYL_CYLINDER_H_
#define _CYL_CYLINDER_H_

#include<iostream>
#include<fstream> 
#include <MAT_matrices.h>
#include <UTI_utils.h>
#include <TOR_torus.h>
const int CYL_interp_order = 3;  // 1 or 3
const int CYL_deriv_order = 3;  // 1 or 3
const int CYL_interp_order_S1 = 3;  // 1 or 3
const int CYL_deriv_order_S1 = 3;  // 1 or 3
#include "CYL_interp.h"
#include "CYL_interp_S1.h"

//-------------------------------------------------------------//
//              class CYL_function_S1<Num_type>
//
// stocks a function defined on a cylinder with values in 
// S1 = R/Z. Be cautious with the interpolation formula between
// samplings.
//-------------------------------------------------------------//
template <typename Num_type> class CYL_function_S1
{
  typedef MAT_matrix_2<Num_type> matrix_base;

  matrix_base m_;

public:
  typedef typename MAT_matrix_2<Num_type>::Array_base array_base;

  CYL_function_S1(): m_() {}
  CYL_function_S1( int nx, int ny ): m_(nx, ny + 1) {} // Cautious +1 for  
                                                    // boundaries of the cylinder !!
  CYL_function_S1( CYL_function_S1 const & other ) {m_ = other.m_;}
  CYL_function_S1( array_base const & other ) {m_ = other;}

  CYL_function_S1 & operator= ( CYL_function_S1 const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  CYL_function_S1 & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }
  CYL_function_S1 & operator= ( Num_type const c ) { // Affect a constant
    m_ = c;
    return *this;
  }

  int n_x() const { return m_.nb_row(); } 
  int n_y() const { return m_.nb_col() - 1; } // Caution : +1 for boundaries !!

  Num_type & mat( int const i, int const j ) {
    return m_(i,j);
  }     
  Num_type mat( int const i, int const j ) const {
    return m_(i,j);
  }

  Num_type operator() ( Num_type const x, Num_type const y ) const {
    return cyl_generic_interp_S1<Num_type, Num_type >( 
				     x - floor(x), y, m_, CYL_interp_order_S1);
  }     
  Num_type operator() (UTI_2vector<Num_type> const & param) const {
    return (*this)(param.x(), param.y());
  }

  Num_type deriv_x( Num_type const x, Num_type const y ) const {
    return cyl_generic_deriv_x_S1<Num_type, Num_type >(
                                     x - floor(x), y, m_, CYL_deriv_order_S1);
  }     
  Num_type deriv_x(UTI_2vector<Num_type> const & param) const {
    return (*this)(param.x(), param.y());
  }
  Num_type deriv_y( Num_type const x, Num_type const y ) const {
    return cyl_generic_deriv_y_S1<Num_type, Num_type >( 
                                      x - floor(x), y, m_, CYL_deriv_order_S1);
  }     
  Num_type deriv_y(UTI_2vector<Num_type> const & param) const {
    return (*this)(param.x(), param.y());
  }

  CYL_function_S1 & operator/= ( CYL_function_S1 const & quotient );
  CYL_function_S1 & operator/= (Num_type const & c) { // divide by 
    m_ /= c;                                          // a constant
      return *this;
  }
  CYL_function_S1 operator/ ( CYL_function_S1 const & quotient ) const;
  
  CYL_function_S1 & operator- ();
  CYL_function_S1 & operator+= (Num_type const & c) { // Add a constant
    m_ += c;
    return *this;
  }
  CYL_function_S1 & operator+= ( CYL_function_S1 const & other );

  void write_VTK(string const & nomFichier, 
		 UTI_2vector<Num_type> const & d1, 
		 UTI_2vector<Num_type> const & d2,
		 int step_x = 1) const ;

  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }
};

//-----------------------------------------------------------------------//
//              operator<< pour CYL_function_S1
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const CYL_function_S1 <Num_type> & f)
{
  return f.print(co);
}



//-------------------------------------------------------------//
//              class CYL_embedding<Num_type>
//-------------------------------------------------------------//
template <typename Num_type> class CYL_embedding 
{
  typedef MAT_matrix_2< UTI_3point<Num_type> > matrix_base;
  typedef typename matrix_base::Array_base array_base;

  matrix_base m_; // The matrix containing the points


public:

  CYL_embedding( int nx, int ny ): m_(nx, ny+1) {} // Caution : +1 for the 
                                                   // cylinder's boundaries !!
  CYL_embedding( CYL_embedding const & other ) {m_ = other.m_;}
  CYL_embedding( array_base const & other ) {m_ = other;}
  CYL_embedding & operator= ( CYL_embedding const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  CYL_embedding & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }

  int n_x() const { return m_.nb_row(); }
  int n_y() const { return m_.nb_col() - 1; } // Caution: +1 for the boundaries !!

  UTI_3point<Num_type> & mat ( int const i, int const j ) {
    return m_(i,j);
  }     
  UTI_3point<Num_type> mat ( int const i, int const j ) const {
    return m_(i,j);
  }

  UTI_3point<Num_type> operator() ( Num_type const x, 
				    Num_type const y ) const {
    return cyl_generic_interp<UTI_3point<Num_type>, UTI_3vector<Num_type>, 
      Num_type >( x - floor(x), y, m_, CYL_interp_order_S1 );
  }     
  UTI_3point<Num_type> operator() (UTI_2vector<Num_type> const & param) const {
    return (*this)(param.x(), param.y());
  }

  UTI_3vector<Num_type> deriv_x( Num_type const x, Num_type const y ) const {
    return cyl_generic_deriv_x<UTI_3point<Num_type>, UTI_3vector<Num_type>, 
      Num_type >( x - floor(x), y, m_, CYL_deriv_order_S1);
  }     
  UTI_3vector<Num_type> deriv_x(UTI_2vector<Num_type> const & param) const {
    return deriv_x(param.x(), param.y());
  }
  UTI_3vector<Num_type> deriv_y( Num_type const x, Num_type const y ) const {
    return cyl_generic_deriv_y<UTI_3point<Num_type>, UTI_3vector<Num_type>, 
      Num_type >( x - floor(x), y, m_, CYL_deriv_order_S1);
  }     
  UTI_3vector<Num_type> deriv_y(UTI_2vector<Num_type> const & param) const {
    return deriv_y(param.x(), param.y());
  }
  UTI_3vector<Num_type> normal( Num_type x, Num_type y ) const; // computes the 
                                                                // normal. 
  UTI_3vector<Num_type> normal( UTI_2vector<Num_type> const & param ) const
  { return normal(param.x(), param.y()); } 

  // pull_back returns the pull back at the point (x,y).
  UTI_sym_2form<Num_type> pull_back( Num_type x, Num_type y ) const; 
  UTI_sym_2form<Num_type> pull_back(UTI_2vector<Num_type> const & param ) const
  { return pull_back(param.x(), param.y()); } 

  void gluing( TOR_embedding<Num_type> const & f,
		CYL_function_S1<Num_type> const & phi_X,
		UTI_2vector<Num_type> const & d_2,
		UTI_2vector<Num_type> const & d_2_perp );

  TOR_embedding<Num_type> & cyl_to_torus( 
				       CYL_function_S1<Num_type> const & phi_X,
				       Num_type const lambda,
				       UTI_2vector<Num_type> const & d_2,
				       UTI_2vector<Num_type> const & d_3,
				       TOR_embedding<Num_type> & f) const;

  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }
  void write_VTK_local(string const & nomFichier, int nbx, int nby) const;
  void write_VTK(string const & nomFichier, 
		 int step_x = 1, 
		 int step_y = 1) const;
};

//-----------------------------------------------------------------------//
//                   operator<< pour CYL_embedding
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const CYL_embedding <Num_type> & f)
{
  return f.print(co);
}

#include "CYL_implementation.h"
#include "CYL_output.h"

#endif // _CYL_CYLINDER_H_
