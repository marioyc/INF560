/*!
*  @file Vector.cpp
*  @brief source of class Vector
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

// project packages
#include "Vector.hpp"

// third-party packages


________________________________________________________________________________

//! @internal default constructor
template <class T, class U>
Vector<T, U>::Vector (
        void ) {

  // set size
  m_size = 0;
  // set pointer to NULL
  m_coef = NULL;

}

________________________________________________________________________________

//! @internal construct a vector from its size
template <class T, class U>
Vector<T, U>::Vector (
        const U size ) {

  // set size
  m_size = size;
  // allocate elements array
  m_coef = new T[size];

}

________________________________________________________________________________

//! @internal default destructor
template <class T, class U>
Vector<T, U>::~Vector (
        void ) {

  if ( m_coef != NULL ) {
    m_size = 0;
    delete [] m_coef;
    m_coef = NULL;
  }

}

________________________________________________________________________________

//! @internal get size of the vector
template <class T, class U>
U Vector<T, U>::GetSize (
        void ) const {

  return m_size;
}

________________________________________________________________________________

//! @internal get pointer to the elements of the vector
template <class T, class U>
T* Vector<T, U>::GetCoef (
        void ) const {

  return m_coef;
}

________________________________________________________________________________

//! @internal set data of the vector
template <class T, class U>
int Vector<T, U>::SetCoef (
        const T* coef ) {

  // copy elements
  for ( U i = 0; i < m_size; i++ ) {
    m_coef[i] = coef[i];
  }

  return 0;
}

________________________________________________________________________________

//! @internal get allocation status of the vector
template <class T, class U>
bool Vector<T, U>::Status (
        void ) const {

  return ( m_coef != NULL && m_size > 0 );
}

________________________________________________________________________________

//! @internal explicit allocation
template <class T, class U>
int Vector<T, U>::Allocate (
        const U size ) {

  // if already allocated, deallocate first
  this->Deallocate( );

  // set size
  m_size = size;
  // allocate elements array
  m_coef = new T[size];

  return 0;
}

________________________________________________________________________________

//! @internal explicit destructor
template <class T, class U>
int Vector<T, U>::Deallocate (
        void ) {

  if ( m_coef != NULL ) {
    m_size = 0;
    delete [] m_coef;
    m_coef = NULL;
  }

  return 0;
}


________________________________________________________________________________

//! @internal assign elements of the vector with a given value
template <class T, class U>
int Vector<T, U>::Assign (
        const U idx_begin,
        const U idx_end,
        const T value ) {

  // assign value to elements
  for ( U i = idx_begin; i <= idx_end; i++ ) {
    m_coef[i] = value;
  }

  return 0;
}

________________________________________________________________________________

//! @internal perform dot product (conj(*this) * y)
template <class T, class U>
T Vector<T, U>::Dot (
        const Vector<T,U>& x ) const {

  // perform dot product
  T res = static_cast<T>( 0 );
  for ( U i = 0; i < m_size; i++ ) {
    res = res + stdmrg::conj( m_coef[i] ) * x.m_coef[i];
  }

  return res;
}

________________________________________________________________________________

//! @internal vector p-norm (p = 2 | inf | ...)
template <class T, class U>
typename stdmrg::type_of<T>::value_type Vector<T, U>::NormLp (
        const int norm_type ) const {

  typedef typename stdmrg::type_of<T>::value_type R;
  R res_norm = R(0);

  if ( norm_type == mathmrg::norm::c_L1 || norm_type == 1 ) {
    // perform computation (norm 1)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += stdmrg::abs( m_coef[i] );
    }
  } else if ( norm_type == mathmrg::norm::c_L2 || norm_type == 2 ) {
    // perform computation (norm 2)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += stdmrg::norm( m_coef[i] );
    }
    res_norm = std::sqrt( res_norm );
  } else if ( norm_type == mathmrg::norm::c_LINF ) {
    // perform computation (norm inf)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm = std::max( res_norm, stdmrg::abs( m_coef[i] ) );
    }
//  } else if ( norm_type == mathmrg::norm::c_FRO ) {
//
  } else if ( norm_type > 0 ) {
    // perform computation (norm xx)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += std::pow( stdmrg::abs( m_coef[i] ), R(norm_type) );
    }
    res_norm = std::pow( res_norm, R(1) / R(norm_type) );
  }

  return res_norm;
}

________________________________________________________________________________

//! @internal vector p-norm (p = 2 | inf | ...), only sum
template <class T, class U>
typename stdmrg::type_of<T>::value_type Vector<T, U>::NormpLp (
        const int norm_type ) const {

  typedef typename stdmrg::type_of<T>::value_type R;
  R res_norm = R(0);

  if ( norm_type == mathmrg::norm::c_L1 || norm_type == 1 ) {
    // perform computation (norm 1)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += stdmrg::abs( m_coef[i] );
    }
  } else if ( norm_type == mathmrg::norm::c_L2 || norm_type == 2 ) {
    // perform computation (norm 2)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += stdmrg::norm( m_coef[i] );
    }
  } else if ( norm_type == mathmrg::norm::c_LINF ) {
    // perform computation (norm inf)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm = std::max( res_norm, stdmrg::abs( m_coef[i] ) );
    }
//  } else if ( norm_type == mathmrg::norm::c_FRO ) {
//
  } else if ( norm_type > 0 ) {
    // perform computation (norm xx)
    for ( U i = 0; i < m_size; i++ ) {
      res_norm += std::pow( stdmrg::abs( m_coef[i] ), R(norm_type) );
    }
  }

  return res_norm;
}


________________________________________________________________________________

//! @internal perform elementwise product of two vectors (x .* y)
template <class T, class U>
int Vector<T, U>::EWProduct (
        const Vector<T,U>& x,
        const Vector<T,U>& y ) {

  // perform computation
  for ( U i = 0; i < x.GetSize( ); i++ ) {
    m_coef[i] = x(i) * y(i);
  }

  return 0;
}

________________________________________________________________________________

//! @internal perform scaled vector addition (alpha * x + (*this))
template <class T, class U>
int Vector<T, U>::Saxpy (
        const T alpha,
        const Vector<T,U>& x ) {

  // perform computation
  for ( U i = 0; i < x.GetSize( ); i++ ) {
    m_coef[i] = alpha * x(i) + m_coef[i];
  }

  return 0;
}

________________________________________________________________________________

//! @internal perform linear combination of two vectors (alpha * x + beta * y)
template <class T, class U>
int Vector<T, U>::Saxpby (
        const T alpha,
        const Vector<T,U>& x,
        const T beta,
        const Vector<T,U>& y ) {

  // perform computation
  for ( U i = 0; i < x.GetSize( ); i++ ) {
    m_coef[i] = alpha * x(i) + beta * y(i);
  }

  return 0;
}

________________________________________________________________________________

//! @internal perform alpha * x
template <class T, class U>
int Vector<T, U>::Scal (
        const T alpha,
        const Vector<T,U>& x ) {

  // perform computation
  for ( U i = 0; i < x.GetSize( ); i++ ) {
    m_coef[i] = alpha * x(i);
  }

  return 0;
}

________________________________________________________________________________

//! @internal conjugate a vector
template <class T, class U>
Vector<T,U>& Vector<T, U>::Conjugate (
        void ) const {

  Vector<T,U>* conj_v = NULL;
  conj_v = new Vector<T,U> ( m_size );
  // perform conjugate of vector
  for ( U i = 0; i < m_size; i++ ) {
    // = conjugate of (*this) value
    (*conj_v)(i) = stdmrg::conj( m_coef[i] );
  }

  return *conj_v;
}

________________________________________________________________________________

//! @internal absolute a vector
template <class T, class U>
Vector<T,U>& Vector<T, U>::Abs (
        void ) const {

  Vector<T,U>* abs_v = NULL;
  abs_v = new Vector<T,U> ( m_size );
  // perform absolute of cpu vector
  for ( U i = 0; i < m_size; i++ ) {
    // = absolute of (*this) value
    (*abs_v)(i) = stdmrg::abs( m_coef[i] );
  }

  return *abs_v;
}

________________________________________________________________________________

//! @internal inverse a vector
template <class T, class U>
Vector<T,U>& Vector<T, U>::Inv (
        void ) const {

  Vector<T,U>* inv_v = NULL;
  inv_v = new Vector<T,U> ( m_size );
  // perform inverse of cpu vector
  for ( U i = 0; i < m_size; i++ ) {
    // = inverse of (*this) value
    (*inv_v)(i) = T(1) / m_coef[i];
  }

  return *inv_v;
}

________________________________________________________________________________

//! @internal overload operator "="
template <class T, class U>
Vector<T,U>& Vector<T, U>::operator= (
        const Vector<T,U>& copy_v ) {

  // if current object is different with the copy object
  if( this != &copy_v ) {
    this->Allocate( copy_v.GetSize( ) );

    // copy elements
    for ( U i = 0; i < copy_v.GetSize( ); i++ ) {
      m_coef[i] = copy_v.m_coef[i];
    }
  }

  return *this;
}

________________________________________________________________________________

//! @internal overload operator "( )"
template <class T, class U>
T& Vector<T, U>::operator() (
        const U idx ) {

  return m_coef[idx];
}

________________________________________________________________________________

//! @internal overload operator "( )" const
template <class T, class U>
T Vector<T, U>::operator() (
        const U idx ) const {

  return m_coef[idx];
}

________________________________________________________________________________

//! @internal overload operator "[ ]"
template <class T, class U>
T& Vector<T, U>::operator[] (
        const U idx ) {

  return m_coef[idx];
}

________________________________________________________________________________

//! @internal overload operator "[ ]" const
template <class T, class U>
T Vector<T, U>::operator[] (
        const U idx ) const {

  return m_coef[idx];
}

________________________________________________________________________________

//! @internal print vector on standard output
template <class T, class U>
int Vector<T, U>::WriteToStdout (
        const char separator,
        const U idx_begin,
        const U idx_end ) const {

  iomrg::printf( ".w. Writing on standard output (vector): \n" );

  // size of the vector
  U size = m_size;

  U _idx_begin = 0;
  U _idx_end = size - 1;
  if ( idx_end != 0 ) {
    size = idx_end - idx_begin + 1;
    _idx_begin = idx_begin;
    _idx_end = idx_end;
  }

  // -- write dimension
  std::cout << size << std::endl;
  // -- write elements
  for ( U i = _idx_begin; i < _idx_end; i++ ) {
    std::cout << m_coef[i] << separator;
  }
  std::cout << m_coef[_idx_end] << std::endl;

  return 0;
}

________________________________________________________________________________

//! @internal read vector from file
template <class T, class U>
int Vector<T, U>::ReadFromFileCsv (
        const char* file_name ) {

  iomrg::printf( ".r. Reading csv file (vector): %s \n", file_name );

  // -- open file
  std::ifstream file(file_name, std::ios::in);

  // read mrg information
  std::string current_line;

  // -- read new line
  std::getline( file, current_line );
  std::stringstream ss_current_line(current_line);
  // first line
  U size;
  ss_current_line >> size;
  // -- update size
  m_size = size;

  // -- read elements
  // allocate elements array
  m_coef = new T[m_size];
  // read elements
  for (U i = 0; i < m_size; i++) {
    file >> m_coef[i];
  }

  // -- close file
  file.close();

  return 0;
}

________________________________________________________________________________

//! @internal write vector into file
template <class T, class U>
int Vector<T, U>::WriteToFileCsv (
        const char* file_name,
        const char separator,
        const U idx_begin,
        const U idx_end ) const {

  iomrg::printf( ".w. Writing csv file (vector): %s \n", file_name );

  // -- open file
  std::ofstream file(file_name, std::ios::out | std::ios::trunc);

  // dimension of the array
  U size = m_size;

  U _idx_begin = 0;
  U _idx_end = size - 1;
  if ( idx_end != 0 ) {
    size = idx_end - idx_begin + 1;
    _idx_begin = idx_begin;
    _idx_end = idx_end;
  }

  // -- write dimension
  file << size << std::endl;

  // -- write elements
  if ( separator != '\0' ) {
    for ( U i = _idx_begin; i < _idx_end; i++ ) {
      file << m_coef[i] << separator;
    }
  } else {
    for ( U i = _idx_begin; i < _idx_end; i++ ) {
      file << m_coef[i];
    }
  }
  file << m_coef[_idx_end] << std::endl;

  // -- close file
  file.close();

  return 0;
}

________________________________________________________________________________

//! instantiate the class
INSTANTIATE_CLASS(Vector)

________________________________________________________________________________
