//************************************************************************
//                                                     TOR_torus.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : May 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains class definitions for functions, vector fields,
// 2-forms, embedding... on a torus. The parameter space is stored
// into a bidimentional matrix corresponding to a
// col*row uniform sampling of the canonical torus R^2/Z^2.
//
// =======================================================

#ifndef _TOR_TORUS_H_
#define _TOR_TORUS_H_

#include<iostream>
#include<fstream> 
#include <MAT_matrices.h>
#include <UTI_utils.h>

const int TOR_interp_order = 3;  // 1 ou 3
const int TOR_deriv_order = 3;  // 1 ou 3

#include <TOR_interp.h>

template <typename Num_type> class TOR_3vector_field;
template <typename Num_type> class TOR_2vector_field;
template <typename Num_type> class TOR_function;
template <typename Num_type> class TOR_embedding;

template <typename Num_type> TOR_function<Num_type> 
J0_inv ( TOR_function<Num_type> const & f );

//--------------------------------------------------------
//  class TOR_function <Num_type>
//
//--------------------------------------------------------

template <typename Num_type> class TOR_function
{
  typedef MAT_matrix_2<Num_type> matrix_base;

  matrix_base m_; // The matrix storing the discrete values

  friend class TOR_3vector_field<Num_type>; // so that 3vector_field 
                                            // can access to m_ in norm()
  friend class TOR_2vector_field<Num_type>; // so that 2vector_field 
                                            // can access to m_ in operator*=()
  friend TOR_function<Num_type> 
  J0_inv<Num_type> ( TOR_function<Num_type> const & f ); // so that J0_inv
                                                         //access à m_

  //--------------------------------------------------------
  // interpolation formula for a point (x,y) in [0,1]^2 
  // replaced in a grid of size nx * ny.
  //--------------------------------------------------------

  Num_type interp_( Num_type const x, Num_type const y ) const {
    return tor_generic_interp< Num_type, Num_type, Num_type>( 
						  x, y, m_, TOR_interp_order );
  }

public:
  typedef typename MAT_matrix_2<Num_type>::Array_base array_base;

  TOR_function(): m_() {}
  TOR_function( int nx, int ny ): m_(nx, ny) {}
  TOR_function( TOR_function const & other ) {m_ = other.m_;}
  TOR_function( array_base const & other ) {m_ = other;}

  TOR_function & operator= ( TOR_function const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  TOR_function & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }
  TOR_function & operator= ( Num_type const c ) { // Affect a constant
    m_ = c;
    return *this;
  }

  int n_x() const { return m_.nb_row(); } 
  int n_y() const { return m_.nb_col(); } 
  Num_type & mat( int const i, int const j ) {
    return m_(i,j);
  }     
  Num_type mat( int const i, int const j ) const {
    return m_(i,j);
  }

  Num_type operator() ( Num_type const x, Num_type const y ) const {
    return interp_( x - floor(x) , y - floor(y) );
  }     
  Num_type operator() ( UTI_2vector<Num_type> const param ) const {
    return (*this)(param.x(), param.y());
  }

  TOR_function & operator*= ( TOR_function const & other ) {
    m_ *= other.m_;
    return *this;
  }
  TOR_function & operator/= ( TOR_function const & quotient );
  TOR_function & operator/= (Num_type const & c) { // divise par une constante
    m_ /= c;
      return *this;
  }
  TOR_function operator/ ( TOR_function const & quotient ) const;
  
  TOR_function & operator- ();
  TOR_function & operator-= ( TOR_function const & other ) {
    m_ -= other.m_;
    return *this;
  }
  TOR_function & operator+= (Num_type const & c) { // Add a constant
    m_ += c;
    return *this;
  }
  TOR_function & operator+= ( TOR_function const & other );

  //-----------------------------------------------------------------------//
  //                   operator J0_inv pour TOR_function
  //-----------------------------------------------------------------------//
  TOR_function J0_inv_func ( TOR_function const & f )
  {
    TOR_function j0_inv_f( f.n_x(), f.n_y() );
    
    j0_inv_f = (typename TOR_function<Num_type>::array_base) J0_inv( f.m_ );
    return j0_inv_f; //(typename TOR_function<Num_type>) J0_inv( f.m_ ); //
  }

  TOR_function & sqrt () // modifie this
  {
    m_ = array_base( blitz::sqrt( (array_base)m_ ) );
    return *this;  
  }  
  TOR_function cos ()
  {
    return array_base( blitz::cos( (array_base)m_ ) );
  }  
  TOR_function sin ()
  {
    return array_base( blitz::sin( (array_base)m_ ) );
  }  
  TOR_function J0 ()
  {
    return array_base( blitz::j0( (array_base)m_ ) );  
  }  
  Num_type max ()
  {
    return blitz::max( (array_base)m_ );  
  }  
  Num_type min ()
  {
    return blitz::min( (array_base)m_ );  
  }  
  Num_type mean ()
  {
    return blitz::mean( (array_base)m_ );  
  }  
  TOR_function abs ()
  {
    return array_base( blitz::abs( (array_base)m_ ) );  
  }  
  
  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }
};

//-----------------------------------------------------------------------//
//                          operator J0_inv pour TOR_function
//-----------------------------------------------------------------------//
template <typename Num_type> TOR_function<Num_type> 
J0_inv ( TOR_function<Num_type> const & f )
{
  TOR_function<Num_type> j0_inv_f( f.n_x(), f.n_y() );
  
//   j0_inv_f = (typename TOR_function<Num_type>::array_base) J0_inv( f.m_ );
//   return f; 
  return (typename TOR_function<Num_type>::array_base) ( J0_inv( f.m_ ) );
}

//-----------------------------------------------------------------------//
//                          operator<< pour TOR_function
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const TOR_function <Num_type> & f)
{
  return f.print(co);
}

//--------------------------------------------------------------------
//                  class TOR_3vector_field <Num_type>
//
//--------------------------------------------------------------------

template <typename Num_type> class TOR_3vector_field {

  typedef MAT_matrix_2< UTI_3vector<Num_type> > matrix_base;
  typedef typename matrix_base::Array_base array_base;
  typedef typename MAT_matrix_2<Num_type>::Array_base array_num_base;

  matrix_base m_; // The matrix of the vectors

  UTI_3vector<Num_type> interp_( Num_type const x, Num_type const y ) const {
    return tor_generic_interp< UTI_3vector<Num_type>, UTI_3vector<Num_type>, 
                           Num_type>( x, y, m_, TOR_interp_order );
  }

  friend class TOR_embedding<Num_type>; // so that TOR_embedding
                                        // can access to m_ in deriv_x,y
public:
  TOR_3vector_field( int nx, int ny ): m_(nx, ny) {}
  TOR_3vector_field( TOR_3vector_field const & other ) {m_ = other.m_; }
  TOR_3vector_field( array_base const & other ) {m_ = other;}

  TOR_3vector_field & operator= ( TOR_3vector_field const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  TOR_3vector_field & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }
  TOR_3vector_field & operator= ( UTI_3vector<Num_type> const & v ) { 
    m_ = v;                         // Affect a constant vector
    return *this;
  }

  UTI_3vector<Num_type> operator() ( Num_type const x, 
				     Num_type const y ) const {
    return interp_( x - floor(x) , y - floor(y) );
  }     
  UTI_3vector<Num_type> operator() (UTI_2vector<Num_type> const param ) const {
    return (*this)(param.x(), param.y());
  }
  UTI_3vector<Num_type> & mat ( int const i, int const j ) {
    return m_(i,j);
  }     
  UTI_3vector<Num_type> mat ( int const i, int const j ) const {
    return m_(i,j);
  }

  int n_x() const { return m_.nb_row(); } 
  int n_y() const { return m_.nb_col(); } 

  TOR_3vector_field & operator+= ( TOR_3vector_field const & other ) {
    m_ += other.m_;
    return *this;
  }
  TOR_3vector_field & operator*= ( TOR_function<Num_type> const & f ) {
    m_ *= f.m_;
    return *this;
  }
  TOR_3vector_field & operator*= ( Num_type const a ) {
    m_ *= a;
    return *this;
  }
  TOR_3vector_field & operator/= ( TOR_function<Num_type> const & f ) {
    m_ /= f.m_;
    return *this;
  }

  TOR_function<Num_type> operator* ( TOR_3vector_field const & other ) const
  {
    return (array_num_base) ( dot_prod( m_, other.m_ ) );
  }

  TOR_3vector_field cross_prod_mat ( TOR_3vector_field const & other ) const
  {
    return (array_base) ( cross_prod( m_, other.m_ ) );
  }

  TOR_function<Num_type> norm() const
  {
    TOR_function<Num_type> f = *this * *this;

    return (array_num_base) sqrt( f.m_ ); 
  }

  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }
};

//-----------------------------------------------------------------------//
//                          operator<< pour TOR_3vector_field
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const TOR_3vector_field <Num_type> & v)
{
  return v.print(co);
}

template <typename Num_type> class TOR_embedding;

//-------------------------------------------------------------------
//                 class TOR_2vector_field <Num_type>
//
//-------------------------------------------------------------------

template <typename Num_type> class TOR_2vector_field {

  typedef MAT_matrix_2< UTI_2vector<Num_type> > matrix_base;
  typedef typename matrix_base::Array_base array_base;

  matrix_base m_;       // The matrix containing the vectors...
  int n_x_, n_y_; // ... its dimensions. Only usefull for
                        // constant fields. 
  bool is_constant_;                  // So as to stock constant fields 
  UTI_2vector<Num_type> const_vec_;   // at low cost. However, a non-constant
                                      // field (! is_constant_) can have all 
                                      // its values equal.
  UTI_2vector<Num_type> interp_( Num_type const x, Num_type const y ) const {
    return tor_generic_interp< UTI_2vector<Num_type>, UTI_2vector<Num_type>, 
                           Num_type>( x, y, m_, TOR_interp_order );
  }


public:
  TOR_2vector_field( int nx, int ny ): m_(nx, ny), 
				       n_x_(nx),
				       n_y_(ny),
				       is_constant_(false) {}
  TOR_2vector_field( TOR_2vector_field const & other ): 
    is_constant_(other.is_constant_),  
    const_vec_(other.const_vec_), 
    n_x_(other.n_x_),
    n_y_(other.n_y_)
  {
    if (! is_constant_)  m_ = other.m_; 
  }
  TOR_2vector_field( UTI_2vector<Num_type> const & v, 
		     int nx = 1, 
		     int ny = 1): 
    is_constant_(true), 
    const_vec_(v), 
    n_x_(nx), 
    n_y_(ny) {}
  TOR_2vector_field(): is_constant_(true), // by default, a field is constant.
		       n_x_(0),
		       n_y_(0) {}
  TOR_2vector_field( array_base const & other ): is_constant_(false),
						 n_x_(other.rows()),
						 n_y_(other.columns()),
						 m_(other) {}
  TOR_2vector_field & operator= ( TOR_2vector_field const & other ) {
    if ( this != &other ) {
      if (! is_constant_)  {
	if ( other.is_constant_ ) m_ = other.const_vec_;
	else {
	  n_x_ = other.n_x_;
	  n_y_ = other.n_y_;
	  m_ = other.m_;
	}
      }
      else {
	is_constant_ = other.is_constant_;
	const_vec_ = other.const_vec_;
	n_x_ = other.n_x_;
	n_y_ = other.n_y_;
	m_ = other.m_;
      }
    }
    return *this;
  }
  TOR_2vector_field & operator= ( array_base const & other ) {
    is_constant_ = false;
    m_ = other;
    n_x_ = other.rows();
    n_y_ = other.columns();
    return *this;
  }
  TOR_2vector_field & operator= ( UTI_2vector<Num_type> const & v ) { 
    if (is_constant_)                         // Affect a constant.
      const_vec_ = v; 
    else
      m_ = v;
    return *this;
  }
  TOR_2vector_field & operator+= ( TOR_2vector_field const & other ) {
    assert( (n_x_ == other.n_x_) && (n_y_ == other.n_y_) );

    if (! is_constant_)  
      if ( other.is_constant_ ) 
	m_ += other.const_vec_;
      else
	m_ += other.m_;
    else 
      if ( other.is_constant_ ) 
	const_vec_ += other.const_vec_;
      else {
	is_constant_ = false;
	UTI_2vector<Num_type> temp_vec = const_vec_;
	m_ = other.m_;
	m_ += temp_vec;
      }
    return *this;
  }
//   TOR_2vector_field & operator+= (UTI_2vector<Num_type> const & v) { 
//     m_ += v;                       // Add a constant vector.
//     return *this;
//   }
  TOR_2vector_field & operator*= ( TOR_function<Num_type> const & f ) {
    assert( ! is_constant_ );
    m_ *= f.m_;
    return *this;
  }

  int n_x() const { assert( (is_constant_) || (n_x_ == m_.nb_row()) );
                    return n_x_; } 
  int n_y() const { assert( (is_constant_) || (n_y_ == m_.nb_col()) );
                    return n_y_; } 

  void set_constant_field( UTI_2vector<Num_type> const & v, 
			   int nx, 
			   int ny) 
  { 
    assert (is_constant_);
    const_vec_ = v;
    n_x_ = nx;
    n_y_ = ny;
  }

  TOR_3vector_field<Num_type> &           // apply a field to an embedding
  apply( TOR_embedding<Num_type> const & f, 
	 TOR_3vector_field<Num_type> & out_field ) const; 
  
  TOR_3vector_field<Num_type> &           // apply a field to an embedding
  operator() ( TOR_embedding<Num_type> const & f, 
	       TOR_3vector_field<Num_type> & out_field ) const { 
    return this->apply( f, out_field ); 
  }

  UTI_2vector<Num_type> operator() ( Num_type const x, 
				     Num_type const y ) const {
    if ( is_constant_ ) return const_vec_;
    else return interp_( x - floor(x) , y - floor(y) );
  }
  UTI_2vector<Num_type> operator() (UTI_2vector<Num_type> const param ) const {
    return *this(param.x(), param.y());
  }
  UTI_2vector<Num_type> & mat ( int const i, int const j ) {
    if ( is_constant_ ) return const_vec_; else return m_(i,j);
  }     
  UTI_2vector<Num_type> mat ( int const i, int const j ) const {
    if ( is_constant_ ) return const_vec_; else return m_(i,j);
  }

  std::ostream & print(std::ostream & co) const {
    if ( is_constant_ ) 
      return co << "champ constant " << const_vec_;
      else return co << m_;
  }
  void write_VTK(string const & nomFichier) const;
};

//-----------------------------------------------------------------------//
//                          operator<< pour TOR_2vector_field
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const TOR_2vector_field <Num_type> & v)
{
  return v.print(co);
}

//-----------------------------------------------------------------------//
//                          class TOR_primitive
//-----------------------------------------------------------------------//
template <typename Num_type> struct TOR_primitive:
public UTI_3vector<Num_type> {                // uniform primitive metric

  TOR_primitive(): UTI_3vector<Num_type>() {} // Constructor by default.
  TOR_primitive( Num_type x, Num_type y ):    // Constructor from a 
    UTI_3vector<Num_type>( x*x, x*y, y*y ) {} // 1-form.

  Num_type operator() ( UTI_2vector<Num_type> const & v1,
			UTI_2vector<Num_type> const & v2 ) const {
    UTI_2vector<Num_type> v_aux( (*this)[0]*v1.x() + (*this)[1]*v1.y(),
				 (*this)[1]*v1.x() + (*this)[2]*v1.y() );
    return v_aux * v2;
  }

  void set_from_1form( Num_type const x, Num_type const y ) {
    *(UTI_3vector<Num_type>*)this = UTI_3vector<Num_type>( x*x, x*y, y*y );
  }
  void set_from_1form( UTI_2vector<Num_type> const & l ) {
    set_from_1form( l.x(), l.y() );
  }
};


//-----------------------------------------------------------------------//
//                          class TOR_sym_2form_field
//-----------------------------------------------------------------------//
template <typename Num_type> class TOR_sym_2form_field {

  typedef MAT_matrix_2< UTI_3vector<Num_type> > matrix_base;
  typedef typename matrix_base::Array_base array_base;
  typedef typename MAT_matrix_2<Num_type>::Array_base array_num_base;

  matrix_base m_; // The E, F, G coeffs of the symetric 2-form

public:
  TOR_sym_2form_field( int nx, int ny ): m_(nx, ny) {}
  TOR_sym_2form_field( TOR_sym_2form_field const & other ) {
    m_ = other.m_; }
  TOR_sym_2form_field( array_base const & other ) {m_ = other;}
  TOR_sym_2form_field( int nx, 
		       int ny, 
		       UTI_3vector<Num_type> const v ): m_(nx, ny) {m_ = v;}

  TOR_sym_2form_field & operator= ( TOR_sym_2form_field const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  TOR_sym_2form_field & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }
  TOR_sym_2form_field & operator= ( UTI_3vector<Num_type> const v ) {
    m_ = v;
    return *this;
  }
  TOR_sym_2form_field & operator+= ( TOR_sym_2form_field const & other ) {
    m_ += other.m_;
    return *this;
  }
  TOR_sym_2form_field & operator-= ( TOR_sym_2form_field const & other ) {
    m_ -= other.m_;
    return *this;
  }
  TOR_sym_2form_field & operator*= ( Num_type c ) {
    m_ *= c;
    return *this;
  }

  int n_x() const { return m_.nb_row(); } 
  int n_y() const { return m_.nb_col(); }

  TOR_function<Num_type> & 
  decompose_alpha(            // returns the alpha coefficient of the decomposition in
	    UTI_3vector<Num_type> const & alpha_filter, // primitive metrics
	    TOR_function<Num_type> & alpha ) const; 

  Num_type & E ( int const i, int const j ) {
    return m_(i,j).x();
  }     
  Num_type E ( int const i, int const j ) const {
    return m_(i,j).x();
  }
  Num_type & F ( int const i, int const j ) {
    return m_(i,j).y();
  }     
  Num_type F ( int const i, int const j ) const {
    return m_(i,j).y();
  }
  Num_type & G ( int const i, int const j ) {
    return m_(i,j).z();
  }     
  Num_type G ( int const i, int const j ) const {
    return m_(i,j).z();
  }
  UTI_3vector<Num_type> & mat ( int const i, int const j ) {
    return m_(i,j);
  }     
  UTI_3vector<Num_type> mat ( int const i, int const j ) const {
    return m_(i,j);
  }
  
  TOR_function<Num_type> norm2() const
  {
    return (array_num_base) ( dot_prod( m_, m_ ) ); 
  }

  TOR_function<Num_type> frob_norm2() const
  {
    return (array_num_base) ( norm2_2form( m_ ) ); 
  }

  bool is_metrique();                         

  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }
};

//-----------------------------------------------------------------------//
//                          operator<< pour TOR_sym_2form_field
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const TOR_sym_2form_field <Num_type> & g)
{
  return g.print(co);
}

//-----------------------------------------------------------------------//
//                          operator* pour TOR_primitive
//-----------------------------------------------------------------------//
template <typename Num_type>  TOR_sym_2form_field<Num_type> 
operator* ( TOR_function<Num_type> const & f, 
	    TOR_primitive<Num_type> const & g )
{
  TOR_sym_2form_field<Num_type> f_g(f.n_x(), f.n_y()); 
  for( int i = 0; i < f.n_x(); i++ )
    for( int j = 0; j < f.n_y(); j++ )
      f_g.mat(i,j) = f.mat(i,j) * g;
  return f_g;
}


//-----------------------------------------------------------------------//
//                     class TOR_embedding <Num_type>
//-----------------------------------------------------------------------//
template <typename Num_type> class TOR_embedding {

  typedef MAT_matrix_2< UTI_3point<Num_type> > matrix_base;
  typedef typename matrix_base::Array_base array_base;

  matrix_base m_; // The matrix containing the points

  UTI_3point<Num_type> interp_( Num_type const x, Num_type const y ) const {
    return tor_generic_interp< UTI_3point<Num_type>, UTI_3vector<Num_type>, 
                               Num_type>( x, y, m_, TOR_interp_order );
  }
  UTI_3vector<Num_type> deriv_x_( Num_type const x, Num_type const y ) const {
    return tor_generic_deriv_x< UTI_3point<Num_type>, UTI_3vector<Num_type>, 
                                Num_type>( x, y, m_, TOR_deriv_order );
  }
  UTI_3vector<Num_type> deriv_y_( Num_type const x, Num_type const y ) const {
    return tor_generic_deriv_y< UTI_3point<Num_type>, UTI_3vector<Num_type>, 
                                Num_type>( x, y, m_, TOR_deriv_order );
  }
  void init_standard_(Num_type r1, Num_type r2);
  void init_standard_bis_(Num_type r1, Num_type r2);

public:

  TOR_embedding(int nx, 
		int ny, 
		TOR_embed_type e_type = TOR_EMPTY,
		Num_type r1 = 0,  
		Num_type r2 = 0): m_(nx, ny) {
    if (e_type == TOR_STANDARD) init_standard_(r1, r2); 
    else if (e_type == TOR_STANDARD_BIS) init_standard_bis_(r1, r2);
  }

  TOR_embedding( TOR_embedding const & other ) {
    m_ = other.m_; }
  TOR_embedding( array_base const & other ) {m_ = other;}
  TOR_embedding & operator= ( TOR_embedding const & other ) {
    if ( this != &other ) {
      m_ = other.m_;
    }
    return *this;
  }
  TOR_embedding & operator= ( array_base const & other ) {
    m_ = other;
    return *this;
  }

  int n_x() const { return m_.nb_row(); }
  int n_y() const { return m_.nb_col(); } 

  UTI_3point<Num_type> operator() ( Num_type const x, 
				    Num_type const y ) const {
    return interp_( x - floor(x) , y - floor(y) );
  }     
  UTI_3point<Num_type> operator() (UTI_2vector<Num_type> const & param) const {
    return (*this)(param.x(), param.y());
  }
  UTI_3vector<Num_type> deriv_x( Num_type const x, 
				 Num_type const y ) const {
    return deriv_x_( x - floor(x) , y - floor(y) );
  }     
  UTI_3vector<Num_type> deriv_x( UTI_2vector<Num_type> const & param ) const {
    return deriv_x(param.x(), param.y());
  }
  UTI_3vector<Num_type> deriv_y( Num_type const x, 
				 Num_type const y ) const {
    return deriv_y_( x - floor(x) , y - floor(y) );
  }     
  UTI_3vector<Num_type> deriv_y( UTI_2vector<Num_type> const & param ) const {
    return deriv_y(param.x(), param.y());
  }
  UTI_3vector<Num_type> deriv( UTI_2vector<Num_type> const & v,
			       Num_type const x, 
			       Num_type const y ) const {
    return v.x() * deriv_x(x,y) + v.y() * deriv_y(x,y);
  }     
  UTI_3vector<Num_type> deriv( UTI_2vector<Num_type> const & v,
			       UTI_2vector<Num_type> const & param ) const {
    return v.x() * deriv_x(param) + v.y() * deriv_y(param);
  }
  UTI_3point<Num_type> & mat ( int const i, int const j ) {
    return m_(i,j);
  }     
  UTI_3point<Num_type> mat ( int const i, int const j ) const {
    return m_(i,j);
  }

  TOR_3vector_field<Num_type> operator- ( TOR_embedding const & other ) const
  {
    return (typename TOR_3vector_field<Num_type>::array_base) 
      ( substract( m_, other.m_ ) );
  }

  TOR_3vector_field<Num_type> & normal_field( 
	  TOR_3vector_field<Num_type> & normals ) const; // computes the normal 
                                                         // fields.
  UTI_3vector<Num_type> normal( Num_type x, Num_type y ) const; // computes the
                                                                // normal. 
  UTI_3vector<Num_type> normal( UTI_2vector<Num_type> const & param ) const
  { return normal(param.x(), param.y()); } 
  
  TOR_sym_2form_field<Num_type> &
  pull_back( TOR_sym_2form_field<Num_type> & g ) const; // returns the pullback 
  UTI_sym_2form<Num_type> pull_back( Num_type x, Num_type y ) const; 
  UTI_sym_2form<Num_type> pull_back(UTI_2vector<Num_type> const & param ) const
  { return pull_back(param.x(), param.y()); } 

  std::ostream & print(std::ostream & co) const {
    return co << m_;
  }

  TOR_3vector_field<Num_type> &    // derivation in x
    deriv_x( TOR_3vector_field<Num_type> & df_x ) const;

  TOR_3vector_field<Num_type> &    // derivation in y
    deriv_y( TOR_3vector_field<Num_type> & df_y ) const;

  void resample( Num_type ratio_x, Num_type ratio_y );

  void write_VTK(string const & nomFichier, 
		 int step_x=1, 
		 int step_y=1) const;

  void write_VTK_local(string const & nomFichier, 
		       int nbre_x=0, 
		       int nbre_y=0) const;

  void write_VTK(string const & nomFichier, 
		 TOR_3vector_field<Num_type> const & champ, 
		 bool reduit=0, 
		 int nbre_x=0, 
		 int nbre_y=0) const;

  void write_Povray(string const & nomFichier, 
		    bool reduit=0, 
		    int nbre_x=0, 
		    int nbre_y=0) const;

  void write_OFF(string const & nomFichier, 
		    bool reduit=0, 
		    int nbre_x=0, 
		    int nbre_y=0) const;

  void write_VRML(string const & nomFichier, int n_max) const;

  void write_binary(string const & nomFichier) const;
  static TOR_embedding read_binary(string const & nomFichier);
};

//-----------------------------------------------------------------------//
//                          operator<< pour TOR_embedding
//-----------------------------------------------------------------------//
template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   const TOR_embedding <Num_type> & f)
{
  return f.print(co);
}

// //-----------------------------------------------------------------------//
// //                          operator>>
// //-----------------------------------------------------------------------//

// template <class Num_type>
// std::istream &operator>>(std::istream &is, 
// 			 TOR_embedding<Num_type>& f)
// {
// ...........
// }

#include "TOR_implementation.h"
#include "TOR_output.h"

#endif // _TOR_TORUS_H_
