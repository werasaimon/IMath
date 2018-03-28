#ifndef SRC_MATHS_RPCOMPLEX_H_
#define SRC_MATHS_RPCOMPLEX_H_

#include "func.hpp"
#include "rpMatrix2x2.hpp"

namespace CMath
{


// define Complex class
//=========================================================================
template<class T> class  rpComplex
{
 private:

	T real,imag;


 public:

	// default constructor
	SIMD_INLINE rpComplex(void)
    : real(0.0) ,
      imag(0.0)
    {}

	// constructor
	SIMD_INLINE rpComplex(T r)
	: real(r) , imag(0.0)
	{}

	// constructor
	SIMD_INLINE rpComplex(T r,T i)
	: real(r) , imag(i)
	{}



	// return the real part
	SIMD_INLINE T getReal(void) const { return(this->real); }

	// return the imagimary part
	SIMD_INLINE T getImag(void) const { return(this->imag); }


	// return the phase argument
	SIMD_INLINE T getArgument(void) const { return atan2(this->imag,this->real); }

	// return the magnitude square
	SIMD_INLINE T lengthSquare(void) const { return (real*real+imag*imag); }

	// return the magnitude
	SIMD_INLINE T length(void) const { return Sqrt(real*real+imag*imag); }


	// return the real matrix
	SIMD_INLINE rpMatrix2x2<T> getMatrix() const
	{
	  return rpMatrix2x2<T> ( real , -imag ,
			                  imag ,  real );
	}


	// return the complex conjugate
	SIMD_INLINE rpComplex<T> getConjugate(void) { return rpComplex<T>(this->real,-this->imag); }


	// overload the + operator to add 2 complex numbers
	SIMD_INLINE rpComplex<T> operator+(rpComplex<T> z) const
	{
		return rpComplex<T>(this->real+z.real,this->imag+z.imag);
	}

	// overload the + operator to add complex and a T
	SIMD_INLINE rpComplex<T> operator+(T a) const
	{
		return rpComplex<T>(this->real+a,this->imag);
	}

	// This is a friend function
	// overload the + operator to add a T and a complex
	friend SIMD_INLINE rpComplex<T> operator+(T a, rpComplex<T> z)
	{
		return rpComplex<T>(a+z.real,z.imag);
	}

	// overload the - operator to subtract 2 complex numbers
	SIMD_INLINE rpComplex<T> operator-(rpComplex<T> z) const
	{
		return rpComplex<T>(this->real-z.real,this->imag-z.imag);
	}

	// overload the - operator to subtract a T from a complex
	SIMD_INLINE rpComplex<T> operator-(T a) const
	{
		return rpComplex<T>(this->real-a,this->imag);
	}

	// overload the - operator to subtract a complex from a T
	friend SIMD_INLINE rpComplex<T> operator-(T a, rpComplex<T> z)
	{
		return rpComplex<T>(a-z.real,-z.imag);
	}

	// overload the - operator to take the negative of a complex
	friend SIMD_INLINE rpComplex<T> operator-(rpComplex<T> z)
	{
		return rpComplex<T>(-z.real,-z.imag);
	}

	// overload the * operator to multiply two complex numbers
	SIMD_INLINE rpComplex<T> operator*(rpComplex<T> z) const
	{
		return rpComplex<T>(this->real*z.real-this->imag*z.imag,
				            this->real*z.imag+this->imag*z.real);
	}

	// overload the * operator to multiply a complex by a T
	SIMD_INLINE rpComplex<T> operator*(T a) const
	{
		return rpComplex<T>(real*a,imag*a);
	}

	// overload the * operator to multiply a T by a complex
	friend SIMD_INLINE rpComplex<T> operator*(T a, rpComplex<T> z)
	{
		return rpComplex<T>(a*z.real,a*z.imag);
	}

	// overload the / operator to divide two complex numbers
	SIMD_INLINE rpComplex<T> operator/(rpComplex<T> z) const
	{
		rpComplex<T> top((*this)*z.getConjugate());
		T bottom(z.lengthSquare());
		rpComplex<T> quo(top/bottom);
		return quo;
	}

	// overload the / operator to divide a complex number by a T
	SIMD_INLINE rpComplex<T> operator/(T a) const
	{
		return rpComplex<T>(this->real/a,this->imag/a);
	}

	// overload the / operator to divide a T by a complex
	friend SIMD_INLINE rpComplex<T> operator/(T a, rpComplex<T> z)
	{
		rpComplex<T> top((a)*z.getConjugate());
		T bottom(z.lengthSquare());
		rpComplex<T> quo(top/bottom);
		return quo;
	}

	// overload the += operator
	SIMD_INLINE const rpComplex<T>& operator+=(const rpComplex<T>& z)
    {
		this->real+=z.real;
		this->imag+=z.imag;
		return *this;
	}

	// overload the -= operator
	SIMD_INLINE const rpComplex<T>& operator-=(const rpComplex<T>& z)
	{
		this->real-=z.real;
		this->imag-=z.imag;
		return *this;
	}

	// overload the == operator
	SIMD_INLINE bool operator==(rpComplex<T> z)
	{
		if (this->real == z.real &&
				this->imag == z.imag)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


    //--------------------------//

	static const rpComplex<T> j;
	static const rpComplex<T> i;

	//--------------------------//

	// take the square root of a complex number
	static SIMD_INLINE rpComplex<T> sqrt(rpComplex<T> z)
    {
		T zsqre,zsqim;

		zsqre = Sqrt(0.5*(z.length()+z.getReal()));
		zsqim = Sqrt(0.5*(z.length()-z.getImag()));

		if (z.getImag() >= 0.0)
		{
			return rpComplex<T>(zsqre,zsqim);
		}
		else
		{
			return rpComplex<T>(zsqre,-zsqim);
		}
	}

	// take the natural log of a complex number
	static SIMD_INLINE rpComplex<T> log(rpComplex<T> z)
	{
		if (z.real < 0 && z.imag == 0.0)
		{
			return Log(z.length())+j*M_PI;
		}
		else
		{
			return Log(z.length())+j*z.getArgument();
		}
	}

	//========= Triqonometriya Function ==========//

	// raise e to a complex number
	static SIMD_INLINE rpComplex<T> exp(rpComplex<T> z)
	{
		return Exp(z.real)*(Cos(z.imag)+j*Sin(z.imag));
	}

	// raise a complex number to a T
	static SIMD_INLINE rpComplex<T> pow(rpComplex<T> z, T c)
	{
		return exp(c*log(z));
	}


	// take the sin of a complex number
	static SIMD_INLINE rpComplex<T> sin(rpComplex<T> z)
	{
		return 0.5*(-j)*exp(j*z)-0.5*(-j)*exp(-j*z);
	}

	// take the cos of a complex number
	static SIMD_INLINE rpComplex<T> cos(rpComplex<T> z)
	{
		return 0.5*exp(j*z)+0.5*exp(-j*z);
	}

	// take the tan of a complex number
	static SIMD_INLINE rpComplex<T> tan(rpComplex<T> z) { return sin(z)/cos(z); }

	// take the sec of a complex number
	static SIMD_INLINE rpComplex<T> sec(rpComplex<T> z) { return 1/cos(z); }

	// take the csc of a complex number
	static SIMD_INLINE rpComplex<T> csc(rpComplex<T> z) { return 1/sin(z); }

	// take the cot of a complex number
	static SIMD_INLINE rpComplex<T> cot(rpComplex<T> z) { return cos(z)/sin(z); }

	// take the sinh of a complex number
	static SIMD_INLINE rpComplex<T> sinh(rpComplex<T> z) { return (exp(z)-exp(-z))/2.0; }

	// take the cosh of a complex number
	static SIMD_INLINE rpComplex<T> cosh(rpComplex<T> z) { return (exp(z)+exp(-z))/2.0; }

	// take the tanh of a complex number
	static SIMD_INLINE rpComplex<T> tanh(rpComplex<T> z) { return sinh(z)/cosh(z); }

	// take the sech of a complex number
	static SIMD_INLINE rpComplex<T> sech(rpComplex<T> z) { return 1/cosh(z); }

	// take the csch of a complex number
	static SIMD_INLINE rpComplex<T> csch(rpComplex<T> z) { return 1/sinh(z); }

	// take the coth of a complex number
	static SIMD_INLINE rpComplex<T> coth(rpComplex<T> z) { return cosh(z)/sinh(z); }

	// take the asin of a complex number
	static SIMD_INLINE rpComplex<T> asin(rpComplex<T> z) { return -j*log(j*z+sqrt(1.0-z*z)); }

	// take the acos of a complex number
	static SIMD_INLINE rpComplex<T> acos(rpComplex<T> z) { return -j*log(z+sqrt(z*z-1.0)); }

	// take the atan of a complex number
	static SIMD_INLINE rpComplex<T> atan(rpComplex<T> z) { return (0.5*j)*log((j+z)/(j-z)); }

	// take the asinh of a complex number
	static SIMD_INLINE rpComplex<T> asinh(rpComplex<T> z) { return log(z+sqrt(z*z+1.0)); }

	// take the acosh of a complex number
	static SIMD_INLINE rpComplex<T> acosh(rpComplex<T> z) { return log(z+sqrt(z*z-1.0)); }

	// take the atanh of a complex number
	static SIMD_INLINE rpComplex<T> atanh(rpComplex<T> z) { return 0.5*log((1.0+z)/(1.0-z)); }



	// round a complex number
	SIMD_INLINE rpComplex<T> rnd(int precision)
	{
		T rnum,inum;
		int tnum;

		rnum = this->real*Pow(10,precision);
		tnum = (int)(rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/Pow(10,precision);

		inum = this->imag*Pow(10,precision);
		tnum = (int)(inum < 0 ? inum-0.5 : inum + 0.5);
		inum = tnum/Pow(10,precision);


		return rpComplex<T>(rnum,inum);
	}

	static SIMD_INLINE rpComplex<T> rnd(rpComplex<T> z, int precision)
	{
		T rnum,inum;
		int tnum;

		rnum = z.real*Pow(10,precision);
		tnum = (int)(rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/Pow(10,precision);

		inum = z.imag*Pow(10,precision);
		tnum = (int)(inum < 0 ? inum-0.5 : inum + 0.5);
		inum = tnum/Pow(10,precision);


		return rpComplex<T>(rnum,inum);
	}

	// end of Complex class
	//=========================================================================

	// round a number
	static SIMD_INLINE T rnd(T num, int precision)
	{
		T rnum;
		int tnum;

		rnum = num*Pow(10,precision);
		tnum = (int)(rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/Pow(10,precision);

		return rnum;
	}


	// create an inserter function so
	// complex types work with cout
	friend std::ostream& operator<<(std::ostream& stream, rpComplex<T> z)
	{
        stream << "(" << "Re: " << z.real << " Im: " << z.imag << ")";

		return stream;
	}


	  /**
	  * Gets string representation.
	  */
	  std::string toString() const
	  {
	      std::ostringstream oss;
	      oss << *this;
	      return oss.str();
	  }

};


template<class T> const rpComplex<T> rpComplex<T>::i = rpComplex<T>(0.0f, 1.0f);
template<class T> const rpComplex<T> rpComplex<T>::j = rpComplex<T>(0.0f, 1.0f);



/// Three dimensional Complex of floats
typedef rpComplex<float> rpComplexf;
/// Three dimensional Complex of doubles
typedef rpComplex<double> rpComplexd;
/// Three dimensional Complex of ints
typedef rpComplex<int> rpComplexi;


} /* namespace */

#endif /* SRC_MATHS_RPCOMPLEX_H_ */
