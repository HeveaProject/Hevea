//************************************************************************
//                                                   MAT_matrices.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : May 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains class definitions for matrices. 
//
// =======================================================

#ifndef _MAT_MATRICES_H_
#define _MAT_MATRICES_H_

#include <blitz/array.h>

using namespace blitz;

template <typename NT> class MAT_matrix_2: 
public blitz::Array< NT, 2 >  {

  typedef MAT_matrix_2 Self; 

 public:
  typedef blitz::Array< NT, 2 > Array_base; 
  typedef NT Num_type;

  MAT_matrix_2(): Array_base() {}
  MAT_matrix_2(int i1, int i2): Array_base(i1,i2) {}
  //  MAT_matrix_2(MAT_matrix_2 const & other): Array_base(other) {}
  MAT_matrix_2(Array_base const & other): Array_base(other) {}

  int nb_col() const { return this->columns(); }
  int nb_row() const { return this->rows(); }

  Self & operator= ( Self const & other ) {
    return *this = (Array_base)other;
  }
  Self & operator= ( Array_base const & other ) {
    if ( (Array_base *) this != &other ) {
      this->resize(other.shape());
      (Array_base) *this = other;
    }
    return *this;
  }
  Self & operator= ( NT const & c ) { // Affectation to a constant
    (Array_base) *this = c;
    return *this;
  }

  Self & operator/= ( Self const & other ) {
    (Array_base) *this /= (Array_base)other;
    return *this;
  }
  template <typename other_NT> Self & operator/= ( 
				       MAT_matrix_2<other_NT> const & other ) {
    (Array_base) *this /= (Array< other_NT, 2 >)other;
    return *this;
  }
  template <typename coeff_type> Self & operator/= ( coeff_type c ) {
    (Array_base) *this /= c;
    return *this;
  }

  Self operator/ ( Self const & other ) const {
    return (Array_base)((Array_base) *this / (Array_base)other);
  }

  Self & operator*= ( Self const & other ) {
    (Array_base) *this *= (Array_base)other;
    return *this;
  }
  template <typename other_NT> Self & operator*= ( 
				       MAT_matrix_2<other_NT> const & other ) {
    (Array_base) *this *= (Array< other_NT, 2 >)other;
    return *this;
  }

  template <typename coeff_type> Self & operator*= ( coeff_type c ) {
    (Array_base) *this *= c;
    return *this;
  }

  Self operator+ ( Self const & other ) const {
    return (Array_base)((Array_base) *this + (Array_base)other);
  }

  Self & operator+= ( Self const & other ) {
    (Array_base) *this += (Array_base)other;
    return *this;
  }

  Self & operator+= ( NT const & c ) { // Add a constant
    (Array_base) *this += c;
    return *this;
  }
  Self & operator-= ( Self const & other ) {
    (Array_base) *this -= (Array_base)other;
    return *this;
  }

  Self operator- ( Self const & other ) const {
    return (Array_base)((Array_base) *this - (Array_base)other);
  }

  Self & operator- () {
    (Array_base) *this = -(Array_base)*this;
    return *this;
  }
};

template<typename NT1, typename NT2> MAT_matrix_2<NT2>
operator* (NT1 const a, MAT_matrix_2<NT2> const & mat )
{
  return (blitz::Array< NT2, 2 >)( a * (blitz::Array< NT2, 2 >)mat );
}


//========================
// Bi-periodic matrices
//========================

// template <class Num_type, size_t dim> MAT_matrix_bip {

//   typedef MAT_matrix_bip<Num_type, dim> Self 

//  public:
//   size_t nb_col();
//   size_t nb_row();

//   class circulator_matrix {
//     operator++();
//     operator--();
//     operator[]();
//   }

//   operator*(Num_type a)
//   operator+()
// }

#endif //_MAT_MATRICES_H_
