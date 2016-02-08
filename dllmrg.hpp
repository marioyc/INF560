/*!
*  @file dllmrg.hpp
*  @brief some tools and settings
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_DLLMRG_HPP_
#define GUARD_DLLMRG_HPP_

// basic packages
#include <stdio.h>
#include <stdarg.h>
#include <limits>
#include <complex>
#include <fstream>

// project packages

// third-party packages


#define \
________________________________________________________________________________ \
// _

#define TEMPLATE_OBJECT_T(name,T) template name<T>;

#define INSTANTIATE_CLASS(name) \
  template class name<double,int>;


//! @namespace stdmrg
//! @brief standard magoules research group zone
namespace stdmrg {

//! @struct type_of
//! @brief a component type used internally for the values of real and imag
//! @details programming example
//! T == float
//! typename type_of<T>::value_type x; ( == float x; )
//! typename type_of<T>::cplx_value_type y; ( == float y; )
template <class T>
struct type_of {
  typedef T value_type;
  typedef T cplx_value_type;
}; // struct type_of {

//! @struct type_of
//! @brief a component type used internally for the values of real and imag
//! @remarks kind of complex is std::complex
//! @details programming example
//! T == std::complex<float>
//! typename type_of<T>::value_type x; ( == float x; )
//! typename type_of<T>::cplx_value_type y; ( == std::complex<float> y; )
template <class T>
struct type_of<std::complex<T> > {
  typedef T value_type;
  typedef std::complex<T> cplx_value_type;
}; // struct type_of<std::complex<T> > {

// -----------------------------------------------------------------------------
// -- Real number
// -----------------------------------------------------------------------------

//! @brief absolute value of a real
//! @param [in] x = real number
//! @return absolute value of real number
template <class T>
T abs ( T x ) { return static_cast<T>( std::abs( static_cast<double>(x) ) ); }

//! @brief norm of a real
//! @param [in] x = real number
//! @return norm of real number
template <class T>
T norm ( T x ) { return x * x; }

//! @brief square root of a real
//! @param [in] x = real number
//! @return sqrt of real number
template <class T>
T sqrt ( T x ) { return std::sqrt(x); }

//! @brief conjugate of a real
//! @param [in] x = real number
//! @return real `conjugate' number
template <class T>
T conj ( T x ) { return x; }

// -----------------------------------------------------------------------------
// -- Complex number
// -----------------------------------------------------------------------------

//! @brief absolute value of a complex
//! @param [in] x = complex number
//! @return absolute value of complex number
template <class T>
T abs ( std::complex<T> x ) { return std::abs<T>(x); }

//! @brief norm of a complex
//! @param [in] x = complex number
//! @return norm of complex number
template <class T>
T norm ( std::complex<T> x ) { return std::norm(x); }

//! @brief square root of a complex
//! @param [in] x = complex number
//! @return sqrt of complex number
template <class T>
std::complex<T> sqrt ( std::complex<T> x ) { return std::sqrt(x); }

//! @brief conjugate of a complex
//! @param [in] x = complex number
//! @return complex conjugate number
template <class T>
std::complex<T> conj ( std::complex<T> x ) { return std::conj(x); }


} // namespace stdmrg {


//! @namespace mathmrg
//! @brief basic mathematics
namespace mathmrg {

//! @struct norm
//! @brief manage norm type
//! @remarks norm ranges:
struct norm {
  enum norm_enum {
    //! norm 1
    c_L1 = 1,
    //! norm 2
    c_L2 = 2,
    //! norm inf
    c_LINF = -1000,
    //! norm frobenius
    c_FRO = -2000,
    //! unknown norm
    c_UNKNOWN = -10000
  }  ; // enum norm_enum {

} ; // struct norm {


//! @brief find dim_1, dim_2; dim_1 * dim_2 = dim and abs(dim_1-dim_2) -> 0
//! @param [in,out] dim_1 = first dimension
//! @param [in,out] dim_2 = second dimension
//! @param [in] dim = dimension
//! @return error code
int CreateDims (
    int& dim_1,
    int& dim_2,
    int dim ) ;

} // namespace mathmrg {


//! @namespace iomrg
//! @brief basic io
namespace iomrg {


// -- printf
extern bool g_enabled_stdout;
extern bool g_enabled_logout;

// -- file to write log of codes
extern bool g_log_project;
extern int g_log_numb_procs;
extern int g_log_proc_numb;
extern int g_log_numb_entries;
extern bool g_log_file_reset;
extern bool g_log_close_file;
extern bool g_log_project_close;


//! @brief set log file name
//! @param [in] file_name  = file name
//! @return error code
int SetLogFileName (
        const char* file_name ) ;

//! @brief set log project name
//! @param [in] project_name  = project name
//! @return error code
int SetLogProjectName (
        const char* project_name ) ;

//! @brief finalize log
//! @return error code
int SetLogFinalize ( void ) ;

//! @brief custom mrg printf
//! @param [in] stream = file stream
//! @param [in] message = message to print
//! @param [in] args = list of arguments
//! @return error code
void mrgprintf (
        FILE* stream,
        const char* message,
        va_list args ) ;

//! @brief custom log mrg printf
//! @param [in] message = message to print
//! @param [in] args = list of arguments
//! @return error code
void mrgprintflog (
        const char* message,
        va_list args ) ;

//! @brief custom mrg printf
//! @param [in] message = message to print
//! @param [in] ... = more arguments
//! @return error code
inline void printf (
        const char* message,
        ... ) {

  va_list args;

  // -- print to stdout
  va_start( args, message );
  mrgprintf( stdout, message, args );
  va_end( args );

  // -- print to logfile
  va_start( args, message );
  mrgprintflog( message, args );
  va_end( args );

}

// -- example:
////  iomrg::g_log_project = true;
////  iomrg::SetLogProjectName("Test");

////  std::string log_file_name = "test.log";
////  iomrg::SetLogFileName( log_file_name.c_str( ) );

////  iomrg::g_log_numb_procs = numb_procs;
////  iomrg::g_log_proc_numb = proc_numb;
////  iomrg::g_log_close_file = false;
////  iomrg::g_log_file_reset = true;

////  // -- enable logout
////  iomrg::g_enabled_logout = false;
////  // -- print stdout only for the first processor (proc_numb = 0)
////  iomrg::g_enabled_stdout = false;
////  if ( proc_numb == 0 ) {
////    iomrg::g_enabled_logout = true;
////    iomrg::g_enabled_stdout = true;
////  }

////  iomrg::SetLogFinalize( );

//! @brief read mrg header from ifstream stream
//! @param [in] stream = stream (ifstream) file
//! @param [in] mrgheader_comments = mrgheader comments
//! @param [in] real_base = base of number inner the file
//! @param [in] standard_base = standard base of the type of file
//! @param [in] comment_symbol = comment symbol (first character)
//! @return error code
int ReadMrgHeader (
        std::ifstream& stream,
        std::string& mrgheader_comments,
        int& real_base,
        const int standard_base,
        const char comment_symbol = '%' ) ;


} // namespace iomrg {



#endif // #ifndef GUARD_DLLMRG_HPP_
