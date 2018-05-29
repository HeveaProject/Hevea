//************************************************************************
//                                                      UTI_utils.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : May 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

/// =======================================================
//
// This file contains class definitions for points and vectors and various common functions.
//
// =======================================================
#ifndef _UTI_UTILS_H_
#define _UTI_UTILS_H_

#include <iostream>
#include <cassert>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_1_PI
#define M_1_PI	0.31830988618379067154	// 1/pi
#endif

#define UTI_PREC_J0	0.000000001	// pour test dans J0

enum TOR_embed_type {TOR_EMPTY = 0, TOR_STANDARD, TOR_STANDARD_BIS};

//----------------------------------------------------------------------//
//  Function interpolant()
//
// S-shaped interpolation function used for gluing cylinders.
// It satisfies : f(1) = 1 et f(0) = f'(0) = f'(1) = f''(0) = f''(1) = 0.
//----------------------------------------------------------------------//

template <typename Num_type> Num_type interpolant( Num_type x ) 
{
  return x*x*x*(x*(x*6 - 15) + 10);
  //  return x*x*x*x*(x*(x*(-20*x + 70) - 84) + 35);
}

//----------------------------------------------------------------------//
//  Function derive_interpolant()
//
// Dérivative of the previous function
//----------------------------------------------------------------------//

template <typename Num_type> Num_type derive_interpolant( Num_type x ) 
{
  return 30*x*x*(x-1)*(x-1);
  //  Num_type y = 1. - x;
  //  return 140*x*x*x*y*y*y;
}


//---------------------------------------------------------//
//       Function J0_inv defined on real numbers     
// Will be used by TOR_function<Num_type>  
//
// The template version does not seem to work with 
// Blitz matrices.                        
//---------------------------------------------------------//
//template <typename Num_type> Num_type J0_inv(Num_type  x)
double J0_inv(NB_TYPE  x)
{
  assert( x >= 0 && x <= 1 );
  NB_TYPE    a=2.404826;
  NB_TYPE    b;
  do
    {
      b=a;
      a=b+(j0(b)-x)/j1(b);
    }
  while   (abs(a-b) > UTI_PREC_J0);
  if  (a>0) return  a;
  else    return  0;
}

//---------------------------------------------------------//
//             class UTI_2vector
//
//---------------------------------------------------------//
template <typename Num_type> class UTI_2vector
{
  Num_type v_[2];

  typedef UTI_2vector Self;
  
public:
  UTI_2vector() {v_[0]=v_[1]=0.0;};
  UTI_2vector(Num_type const x, Num_type const y) {
    v_[0]=x; v_[1]=y;
  }
  UTI_2vector( UTI_2vector const & other) {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
  }   
  Self const & operator=( Self const & other) {
    if ( *this != other ) {
      v_[0] = other.v_[0];
      v_[1] = other.v_[1];
    }
    return *this;
  }      
  Num_type& operator[](int i)  {
    assert(i>=0 && i<2);
    return(v_[i]);
  }
  Num_type operator[](int i) const {
    assert(i>=0 && i<2);
    return(v_[i]);
  }
  operator Num_type *() {return(v_);}
  operator Num_type const *() const {return(v_);}

  Num_type& x() {return(v_[0]);}
  Num_type& y() {return(v_[1]);}
  Num_type x() const {return(v_[0]);}
  Num_type y() const {return(v_[1]);}

  bool operator==(Self const & other) const {
    return (v_[0] == other.v_[0]) && 
           (v_[1] == other.v_[1]);
  }
  bool operator!=(Self const & other) const {
    return !(*this == other);
  }

  Self operator-(Self const & other) const {
    return Self(v_[0] - other.v_[0], v_[1] - other.v_[1]);
  }
  Self operator+(Self const & other) const {
    return Self(v_[0] + other.v_[0], v_[1] + other.v_[1]);
  }
  Self & operator+= (Self const & other) {
    v_[0] += other.v_[0];
    v_[1] += other.v_[1];
    return *this;
  }
  Self operator/ (Num_type a) const {
    return Self(v_[0] / a, v_[1] / a);
  }
  Self & operator/= (Num_type a) {
    v_[0] /= a;
    v_[1] /= a;
    return *this;
  }
  Self & operator*= (Num_type a) {
    v_[0] *= a;
    v_[1] *= a;
    return *this;
  }
  Num_type operator* ( Self const & other ) const { // scalar product
    return v_[0]*other.v_[0] + v_[1]*other.v_[1];
  }  
  Num_type norm() const {
    return sqrt( (*this)*(*this) );
  }
};

template <typename Num_type> 
UTI_2vector<Num_type>
operator* (Num_type const & a, const UTI_2vector<Num_type> & v)
{
  return UTI_2vector<Num_type>(a*v[0], a*v[1]);
}



template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   UTI_2vector<Num_type> const & v)
{
  return co << "(" << v[0] << ", " << v[1] << ")";
}



//---------------------------------------------------------//
//             class UTI_3vector
//
//---------------------------------------------------------//
template <typename Num_type> class UTI_3vector
{
  Num_type v_[3];

  typedef UTI_3vector Self;
  
public:
  UTI_3vector() {v_[0]=v_[1]=v_[2]=0.0;};
  UTI_3vector(Num_type x, Num_type y, Num_type z) {
    v_[0]=x; v_[1]=y; v_[2]=z;
  }
  UTI_3vector( UTI_3vector const & other) {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
  }   
  Self const & operator=( Self const & other) {
    if ( *this != other ) {
      v_[0] = other.v_[0];
      v_[1] = other.v_[1];
      v_[2] = other.v_[2];
    }
    return *this;
  }      
  Num_type& operator[](int i)  {
    assert(i>=0 && i<3);
    return(v_[i]);
  }
  Num_type operator[](int i) const {
    assert(i>=0 && i<3);
    return(v_[i]);
  }
  operator Num_type *() {return(v_);}
  operator Num_type const *() const {return(v_);}

  Num_type& x() {return(v_[0]);}
  Num_type& y() {return(v_[1]);}
  Num_type& z() {return(v_[2]);}
  Num_type x() const {return(v_[0]);}
  Num_type y() const {return(v_[1]);}
  Num_type z() const {return(v_[2]);}

  bool operator==(Self const & other) const {
    return (v_[0] == other.v_[0]) && 
           (v_[1] == other.v_[1]) && 
           (v_[2] == other.v_[2]);
  }
  bool operator!=(Self const & other) const {
    return !(*this == other);
  }

  Self operator-(Self const & other) const {
    return Self(v_[0] - other.v_[0], v_[1] - other.v_[1], v_[2] - other.v_[2]);
  }
  Self & operator-= (Self const & other) {
    v_[0] -= other.v_[0];
    v_[1] -= other.v_[1];
    v_[2] -= other.v_[2];
    return *this;
  }
  Self operator+(Self const & other) const {
    return Self(v_[0] + other.v_[0], v_[1] + other.v_[1], v_[2] + other.v_[2]);
  }
  Self & operator+= (Self const & other) {
    v_[0] += other.v_[0];
    v_[1] += other.v_[1];
    v_[2] += other.v_[2];
    return *this;
  }

  Num_type operator* ( Self const & other ) const { // produit scalaire
    return v_[0]*other.v_[0] + v_[1]*other.v_[1] + v_[2]*other.v_[2];
  }
  Self & operator*= (Num_type a) {
    v_[0] *= a;
    v_[1] *= a;
    v_[2] *= a;

    return *this;
  }
  Self operator/ (Num_type a) const {
    return Self(v_[0] / a, v_[1] / a, v_[2] / a);
  }
  Self & operator/= (Num_type a) {
    v_[0] /= a;
    v_[1] /= a;
    v_[2] /= a;

    return *this;
  }

  Self & normalize() {
    Num_type norme = sqrt( (*this)*(*this) );
    return (*this /= norme);
  }

  Num_type norm() const {
    return sqrt( (*this)*(*this) );
  }
};

template <typename Num_type> 
UTI_3vector<Num_type>
operator* (Num_type const & a, UTI_3vector<Num_type> const & v)
{
  return UTI_3vector<Num_type>(a*v[0], a*v[1], a*v[2]);
}


template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   UTI_3vector<Num_type> const & v)
{
  return co << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
}


template <typename Num_type> Num_type 
dot_prod (  UTI_3vector< Num_type > const & left,
	    UTI_3vector< Num_type > const & right )
{
  return left * right;
}

template <typename Num_type> Num_type 
norm2_2form (  UTI_3vector< Num_type > const & f )
{
  return f[0]*f[0] + 2*f[1]*f[1] + f[2]*f[2];
}

template <typename Num_type> UTI_3vector< Num_type >
cross_prod( UTI_3vector< Num_type > const & left,
	    UTI_3vector< Num_type > const & right ) { // produit vectoriel
    
    return UTI_3vector<Num_type>( left[1]*right[2] - left[2]*right[1],
				  left[2]*right[0] - left[0]*right[2],
				  left[0]*right[1] - left[1]*right[0]);
}

template <typename Num_type> class UTI_3point
{
  Num_type p_[3];

  typedef UTI_3point Self;
  
public:
  UTI_3point() {p_[0]=p_[1]=p_[2]=0.0;};
  UTI_3point(Num_type x, Num_type y, Num_type z) {
    p_[0]=x; p_[1]=y; p_[2]=z;
  }
  UTI_3point( UTI_3point const & p) {
    p_[0] = p.p_[0];
    p_[1] = p.p_[1];
    p_[2] = p.p_[2];
  }   
  Self const & operator=( Self const & other) {
    if ( *this != other ) {
      p_[0] = other.p_[0];
      p_[1] = other.p_[1];
      p_[2] = other.p_[2];
    }
    return *this;
  }      
  Num_type& operator[](int i)  {
    assert(i>=0 && i<3);
    return(p_[i]);
  }
  Num_type operator[](int i) const {
    assert(i>=0 && i<3);
    return(p_[i]);
  }
  operator Num_type *() {return(p_);}
  operator Num_type const *() const {return(p_);}

  Num_type& x() {return(p_[0]);}
  Num_type& y() {return(p_[1]);}
  Num_type& z() {return(p_[2]);}
  Num_type x() const {return(p_[0]);}
  Num_type y() const {return(p_[1]);}
  Num_type z() const {return(p_[2]);}

  UTI_3vector<Num_type> operator-(Self const & other) const {
    return UTI_3vector<Num_type>(p_[0] - other.p_[0], 
				 p_[1] - other.p_[1], 
				 p_[2] - other.p_[2]);
  }
  Self & operator-=(UTI_3vector<Num_type> const & v) {
    p_[0] -= v[0];
    p_[1] -= v[1];
    p_[2] -= v[2];
    return *this;
  }
  Self operator+(UTI_3vector<Num_type> const & v) const {
    return UTI_3point(p_[0] + v[0], 
		      p_[1] + v[1], 
		      p_[2] + v[2]);
  }
  Self & operator+=(UTI_3vector<Num_type> const & v) {
    p_[0] += v[0];
    p_[1] += v[1];
    p_[2] += v[2];
    return *this;
  }

  bool operator==(UTI_3point const & q) const {
    return (p_[0] == q.p_[0]) && (p_[1] == q.p_[1]) && (p_[2] == q.p_[2]);
  }
  bool operator!=(UTI_3point const & q) const {
    return !(*this == q);
  }
  Self & scale(Num_type a) {
    p_[0] *= a;
    p_[1] *= a;
    p_[2] *= a;

    return *this;
  }
};


template <typename Num_type> 
std::ostream & operator<< (std::ostream & co, 
			   UTI_3point<Num_type> const & p)
{
  return co << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
}

template <typename NT> UTI_3vector<NT> substract(UTI_3point<NT> const & left, 
						 UTI_3point<NT> const & right ) 
{ 
  return  left - right; 
} 

template <typename Num_type> class UTI_sym_2form: 
public UTI_3vector<Num_type>
{
  typedef UTI_sym_2form Self;
  typedef UTI_3vector<Num_type> base;
public:
  UTI_sym_2form():
    UTI_3vector<Num_type>() {}

  UTI_sym_2form(Num_type x, Num_type y, Num_type z):
    UTI_3vector<Num_type>(x,y,z) {}

  UTI_sym_2form(UTI_3vector<Num_type> const & v):
    UTI_3vector<Num_type>(v) {}

  Num_type& E() {return this->x();}
  Num_type& F() {return this->y();}
  Num_type& G() {return this->z();}
  Num_type E() const {return this->x();}
  Num_type F() const {return this->y();}
  Num_type G() const {return this->z();}
  
  Self operator-(Self const & other) const {
    return (base)*this - (base)other;
  }
  Self operator+(Self const & other) const {
    return (base)*this + (base)other;
  }
  Self operator/ (Num_type a) const {
    return (base)*this / a;
  }
};

template <typename Num_type> 
UTI_sym_2form<Num_type>
operator* (Num_type const & a, UTI_sym_2form<Num_type> const & f)
{
  return UTI_sym_2form<Num_type>(a*f[0], a*f[1], a*f[2]);
}

template <typename Num_type> Num_type 
norm2_2form (  UTI_sym_2form< Num_type > const & f )
{
  return f[0]*f[0] + 2*f[1]*f[1] + f[2]*f[2];
}

#endif // _UTI_UTILS_H_
