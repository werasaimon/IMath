#ifndef RPVECTOR2_H
#define RPVECTOR2_H

#include "func.hpp"

namespace CMath
{



template<class T> class  rpMatrix2x2;

/**
 * Class for two dimensional vector.
 * There are three ways of accessing vector components.
 * Let's have <code>Vector2f v</code>, you can either:
 * <ul>
 * 	<li>access as position(x,y) &mdash; <code>v.x = v.y = 3;</code></li>
 * 	<li>access as texture coordinate (s,t) &mdash; <code>v.s = v.t = 3;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = 3;</code></li>
 * </ul>
 */

template<class T>
class rpVector2D
{
public:
    union
    {
        /**
         * First element of vector, alias for X-coordinate.
         */
        T x;

        /**
         * First element of vector, alias for U-coordinate.
         * For textures notation.
         */
        T u;
    };

    union
    {
        /**
         * Second element of vector, alias for Y-coordinate.
         */
        T y;

        /**
         * Second element of vector, alias for V-coordinate.
         * For textures notation.
         */
        T v;
    };

public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and sets to (0,0)
     */
    SIMD_INLINE rpVector2D()
            : x(0), y(0)
    {
    }

    /**
     * Creates and sets to (x,y)
     * @param nx initial x-coordinate value
     * @param ny initial y-coordinate value
     */
    SIMD_INLINE rpVector2D(T nx, T ny)
            : x(nx), y(ny)
    {
    }

//    /**
//     * Copy constructor.
//     * @param src Source of data for new created instance.
//     */
//    rpVector2D(const rpVector2D<T>& src)
//            : x(src.x), y(src.y)
//    {
//    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created instance.
     */
    template<class FromT>
    SIMD_INLINE rpVector2D(const rpVector2D<FromT>& src)
            : x(static_cast<T>(src.x)), y(static_cast<T>(src.y))
    {
    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void setToZero()
    {
      x = T(0);
      y = T(0);
    }


    SIMD_INLINE void setAllValues(T newX, T newY )
    {
      x = newX;
      y = newY;
    }

    SIMD_INLINE T getX() const { return x; }
    SIMD_INLINE T getY() const { return y; }

    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }


    //----------------[ access operators ]-------------------
    /**
     * Copy casting operator
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE rpVector2D<T>& operator=(const rpVector2D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        return *this;
    }

    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator=(const rpVector2D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, else reference to y
     * y coordinate.
     */
    SIMD_INLINE T& operator[](int n)
    {
    	static_assert(sizeof(*this) == sizeof(T[2]), "");
    	assert(n >= 0 && n < 2);
    	return (&x)[n];
    }

    /**
     * Constant array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, else reference to y
     * y coordinate.
     */
    SIMD_INLINE const T& operator[](int n) const
    {
    	static_assert(sizeof(*this) == sizeof(T[2]), "");
    	assert(n >= 0 && n < 2);
    	return (&x)[n];
    }

    //---------------[ vector aritmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator+(const rpVector2D<T>& rhs) const
    {
        return rpVector2D<T>(x + rhs.x, y + rhs.y);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator-(const rpVector2D<T>& rhs) const
    {
        return rpVector2D<T>(x - rhs.x, y - rhs.y);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator*(const rpVector2D<T>& rhs) const
    {
        return rpVector2D<T>(x * rhs.x, y * rhs.y);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator/(const rpVector2D<T>& rhs) const
    {
        return rpVector2D<T>(x / rhs.x, y / rhs.y);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator+=(const rpVector2D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    /**
     * Substraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator-=(const rpVector2D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator*=(const rpVector2D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator/=(const rpVector2D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        return *this;
    }

    //--------------[ scalar vector operator ]--------------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator+(T rhs) const
    {
        return rpVector2D<T>(x + rhs, y + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator-(T rhs) const
    {
        return rpVector2D<T>(x - rhs, y - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator*(T rhs) const
    {
        return rpVector2D<T>(x * rhs, y * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T> operator/(T rhs) const
    {
        return rpVector2D<T>(x / rhs, y / rhs);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector2D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }


    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  rpVector2D<T> operator*(T number, const rpVector2D<T>& vector)
    {
        return rpVector2D<T>(number * vector.x, number * vector.y);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  rpVector2D<T> operator/( T number , const rpVector2D<T>& vector )
    {
        return rpVector2D<T>(vector.x / number, vector.y / number);
    }


    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T dot(const rpVector2D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y;
    }

    /**
     * Cross product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T cross(const rpVector2D<T>& rhs) const
    {
        // just calculate the z-component
        return x * rhs.y - y * rhs.x;
    }

    //--------------[ equality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < EPSILON,
     * same for y-coordinate.
     */
    SIMD_INLINE bool operator==(const rpVector2D<T>& rhs) const
    {
        return (Abs(x - rhs.x) < EPSILON) && (Abs(y - rhs.y) < EPSILON);
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const rpVector2D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE rpVector2D<T> operator-() const
    {
        return rpVector2D<T>(-x, -y);
    }

    //-------------[ size operations ]---------------------------

    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length() method.
     */
    SIMD_INLINE T lengthSquare() const
    {
        return x * x + y * y;
    }


    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T length() const
    {
        return (T) Sqrt(x * x + y * y);
    }

    /**
     * Normalize vector
     */
    SIMD_INLINE void normalize()
    {
        T s = length();
        x /= s;
        y /= s;
    }


    /**
     * Normalize unit vector
     */
    SIMD_INLINE rpVector2D<T> getUnit() const
    {
        T lengthVector = length();
        if (lengthVector < EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return rpVector2D<T>( x * lengthInv,
                              y * lengthInv);
    }


    /**
     * Inverse vector
     */
    SIMD_INLINE rpVector2D<T> getInverse() const
    {
        return rpVector2D<T>( T(1.0/x) , T(1.0/y) );
    }





    //--------------[ misc. operations ]-----------------------
    /**
     * Linear interpolation of two vectors
     * @param fact Factor of interpolation. For translation from position
     * of this vector to vector r, values of factor goes from 0.0 to 1.0.
     * @param r Second Vector for interpolation
     * @note However values of fact parameter are reasonable only in interval
     * [0.0 , 1.0], you can pass also values outside of this interval and you
     * can get result (extrapolation?)
     */
    SIMD_INLINE rpVector2D<T> lerp(T fact, const rpVector2D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }


    SIMD_INLINE T getAngleBetween( const rpVector2D<T> &Vector2 ) const
    {
        rpVector2D<T> Vector1(*this);
        T dotProduct = Vector1.dot(Vector2);
        T vectorsMagnitude = (Vector1.length()) * (Vector2.length());
        T angle = acos(dotProduct / vectorsMagnitude);
        if( __isnan(angle)) return 0;
        return (angle);
    }

    //-------------[ conversion ]-----------------------------
    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class rpVector2D<T>)
     * used for passing rpVector2D<T> values to gl*2[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return (T*) this;
    }
    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class rpVector2D<T>)
     * used for passing rpVector2D<T> values to gl*2[fd] functions.
     */
    SIMD_INLINE operator const T*() const
    {
        return (const T*) this;
    }

    //-------------[ output operator ]------------------------
    /**
     * Output to stream operator
     * @param lhs Left hand side argument of operator (commonly ostream instance).
     * @param rhs Right hand side argument of operator.
     * @return Left hand side argument - the ostream object passed to operator.
     */
    friend std::ostream& operator<<(std::ostream& lhs, const rpVector2D<T>& rhs)
    {
        lhs << "[" << rhs[0] << "," << rhs[1] << "]";
        return lhs;
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


public:

    /**
     * The multiplicitive identity vector
     */
    static const rpVector2D<T> IDENTITY;
    /**
     * The additive identity vector.
     */
    static const rpVector2D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const rpVector2D<T> X;

    /**
     * The identity vector Y.
     */
    static const rpVector2D<T> Y;

};

template<class T> const rpVector2D<T> rpVector2D<T>::IDENTITY(1.0, 1.0);
template<class T> const rpVector2D<T> rpVector2D<T>::ZERO(0.0, 0.0);
template<class T> const rpVector2D<T> rpVector2D<T>::X(1.0, 0.0);
template<class T> const rpVector2D<T> rpVector2D<T>::Y(0.0, 1.0);


template<class T> const
static T cross(const rpVector2D<T>& a, const rpVector2D<T>& b)
{
    return a.cross(b);
}

template<class T> const
static T dot(const rpVector2D<T>& a, const rpVector2D<T>& b)
{
    return a.dot(b);
}

//--------------------------------------
// Typedef shortcuts for 2D vector
//-------------------------------------
/// Two dimensional Vector of floats
typedef class rpVector2D<float> Vector2f;
/// Two dimensional Vector of doubles
typedef class rpVector2D<double> Vector2d;
/// Two dimensional Vector of ints
typedef class rpVector2D<int> Vector2i;

} /* namespace */

#endif // RPVECTOR2_H
