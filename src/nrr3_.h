#ifndef _NR3_H_
#define _NR3_H_

//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

/* JIMMY--
// all the system #include's we'll ever need
*/ 
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
//#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;

// Vector and Matrix Classes

//#ifdef _USESTDVECTOR_
//#define NRvector vector
//#else

template <class T>
class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(int n, const T &a);	//initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector & operator=(const NRvector &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~NRvector();
};

//#endif //ifdef _USESTDVECTOR_








// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type /// ERA DOUBLE (SIN LONG) ORIGINALMENTE (jimmy)
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;



// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

#endif /* _NR3_H_ */
//EOF
