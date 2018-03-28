#ifndef FUNC_H
#define FUNC_H



#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <limits>

#define LIGHT_MAX_VELOCITY_C 1.0

#ifndef M_PI
#define M_PI  acos(-1) // 3.14159265358979323846  /* pi */
#endif

#define DEG2RAD(x) ((x * M_PI) / 180.0)
//#define EPSILON (4.37114e-07)

const double epsilon = 4.37114e-05;
#define EPSILON epsilon

namespace CMath
{


//-------------------------------------------------------------------------------
//-- Typedefs, Structs ----------------------------------------------------------
//-------------------------------------------------------------------------------

#define kEpsilon    1.0e-5f

#define kPI         3.1415926535897932384626433832795f
#define kHalfPI     1.5707963267948966192313216916398f
#define kTwoPI      2.0f*kPI

#undef APPROXIMATION // commented out for the moment



#define SIMD_INLINE  inline

#define ATTRIBUTE_ALIGNED16(a)  __declspec(align(16)) a
#define ATTRIBUTE_ALIGNED64(a)  __declspec(align(64)) a
#define ATTRIBUTE_ALIGNED128(a) __declspec(align(128)) a


#define FLT_MIN 1.175494351e-38F /* min positive value */
#define FLT_MAX 3.402823466e+38F /* max value */


#define FLT_EPSILON 1.1920928955078125E-7f
#define MACHINE_EPSILON 10E-7f


#define btRecipSqrt(x) ((float)(float(1.0)/sqrt(float(x))))
#define SIMDSQRT12       float(0.7071067811865475244008443621048490)


const float TFSIMD_EPSILON = std::numeric_limits<float>::epsilon();

// FLOATING POINT VALIDITY
template<typename T> bool is_infinite(const T &value)
{
    T max_value = (std::numeric_limits<T>::max)();
    T min_value = -max_value;
    return !(min_value <= value && value <= max_value);
}

template<typename T> bool is_nan(const T &value)
{
    // True if NAN
    return value != value;
}

template<typename T> bool is_valid(const T &value)
{
    return !is_infinite(value) && !is_nan(value);
}




// ---------- Mathematics functions ---------- //

template<typename T> SIMD_INLINE T Min(T a, T b)
{
    return (a > b) ? b : a;
}

template<typename T> SIMD_INLINE T Max(T a, T b)
{
    return (a < b) ? b : a;
}

template<typename T> SIMD_INLINE T Min(T a, T b, T c)
{
    return Min<T>(Min<T>(a, b), c);
}

template<typename T> SIMD_INLINE T Max(T a, T b, T c)
{
    return Max<T>(Max<T>(a, b), c);
}

template<typename T> SIMD_INLINE T Clamp(T a, T min, T max)
{
    return Max<T>(Min<T>(a, max), min);
}

template<typename T> SIMD_INLINE T Wrap(T a, T min, T max)
{
    return (a < min) ? max - (min - a) : (a > max) ? min - (max - a) : a;
}

template<typename T> SIMD_INLINE void Swap(T& a, T& b)
{
    T c = a;
    a = b;
    b = c;
}

//-----------------------------------------------------//


template<typename T> SIMD_INLINE T Sign(T  x)
{
    return (x < 0.0f) ? -1.0f : 1.0f;
}


template<typename T> SIMD_INLINE T Pi(void)
{
    static const T  gPi = (T ) atan(1.0f) * 4.0f;
    return gPi;
}
template<typename T> SIMD_INLINE T TwoPi(void)
{
    static const T  gTwoPi = (T ) atan(1.0f) * 8.0f;
    return gTwoPi;
}
template<typename T> SIMD_INLINE T Modulo(T  x, T  div)
{
    return (T ) fmod((double) x, (double) div);
}
template<typename T> SIMD_INLINE T Abs(T  x)
{
    return (T ) fabs(x);
}

//-------------------------------------------------------------------------------
// Is this floating point value close to zero?
//-------------------------------------------------------------------------------
template<typename T>  bool IsZero( T a, T epsilon = kEpsilon )
{
    return (Abs(a) <= epsilon);
}


template<typename T> SIMD_INLINE T DegreesToRadians(T  Degrees)
{
    return Degrees * (M_PI / 180.0f);
}
template<typename T> SIMD_INLINE T RadiansToDegrees(T  Radians)
{
    return Radians * (180.0f / M_PI);
}


template<typename T> SIMD_INLINE T Sin(T  Radians)
{
    return (T ) sin(Radians);
}
template<typename T> SIMD_INLINE T Cos(T  Radians)
{
    return (T ) cos(Radians);
}

template<typename T> SIMD_INLINE T Sinh(T  Radians)
{
    return (T ) sinh(Radians);
}
template<typename T> SIMD_INLINE T Cosh(T  Radians)
{
    return (T ) cosh(Radians);
}

template<typename T> SIMD_INLINE T Sinc_pi(T x)
{
  return (T) sin(x) / x;
}

template<typename T> SIMD_INLINE T Cosc_pi(T x)
{
  return (T) cos(x) / x;
}

template<typename T> SIMD_INLINE T Sinhc_pi(T x)
{
	return (T) sinh(x) / x;
}

template<typename T> SIMD_INLINE T Coshc_pi(T x)
{
	return (T) cosh(x) / x;
}


template<typename T> SIMD_INLINE T Tangent(T  Radians)
{
    return (T ) tan(Radians);
}
template<typename T> SIMD_INLINE T ArcTangent(T  X)
{
    return (T ) atan(X);
}
template<typename T> SIMD_INLINE T ArcTang(T  X, T  Y)
{
    return (T ) atan(X);
}


template<typename T> SIMD_INLINE T ArcSin(T  X)
{
    return (T ) asin(X);
}
template<typename T> SIMD_INLINE T ArcCos(T  X)
{
    return (T ) acos(X);
}


template<typename T> SIMD_INLINE T Atan2(T  X , T  Y)
{
    return (T ) atan2(X , Y);
}

template<typename T> SIMD_INLINE T Log(T  X)
{
    return (T ) log(X);
}


template<typename T> SIMD_INLINE T Log2(T  X )
{
    return (T ) log2(X);
}

template<typename T> SIMD_INLINE T Exp(T  X )
{
    return (T ) exp(X);
}

template<typename T> SIMD_INLINE T Exp2(T  X )
{
    return (T ) exp2(X);
}

template<typename T> SIMD_INLINE T SquareRoot(T  In)
{
    return (T ) sqrt(In);
}
template<typename T> SIMD_INLINE T Sqrt(T  In)
{
    return (T ) sqrtf(In);
}

template<typename T> SIMD_INLINE T Rand(T  r = 1.0f)
{
    return rand() / ((T ) RAND_MAX) * r;
}
template<typename T> SIMD_INLINE T Rand(T  min, T  max)
{
    return min + Rand(max - min);
}

template<typename T> SIMD_INLINE T Pow(T  X , T  Y)
{
    return (T ) pow(X , Y);
}



template<typename T> SIMD_INLINE T min3(T a, T b, T c)
{
    return Min(Min(a, b), c);
}


template<typename T> SIMD_INLINE T max3(T a, T b, T c)
{
    return Max(Max(a, b), c);
}


/// Function to test if two real numbers are (almost) equal
/// We test if two numbers a and b are such that (a-b) are in [-EPSILON; EPSILON]
template<typename T> SIMD_INLINE bool approxEqual(T a, T b, T epsilon = MACHINE_EPSILON)
{
    return (Abs(a - b) < T(epsilon) );
}


template<typename T> SIMD_INLINE bool sameSign(T a, T b)
{
    return a * b >= T(0.0);
}


}  /* namespace */

#endif // FUNC_H
