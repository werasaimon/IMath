#ifndef RPVECTOR4D_H
#define RPVECTOR4D_H

#include "rpVector3D.hpp"


namespace CMath
{

  
/**
 * Class for four dimensional vector.
  * There are four ways of accessing vector components.
 * Let's have <code>Vector4f v</code>, you can either:
 * <ul>
 * 	<li>access as position in projective space (x,y,z,w) &mdash; <code>v.x = v.y = v.z = v.w = 1;</code></li>
 * 	<li>access as texture coordinate (s,t,u,v) &mdash; <code>v.s = v.t = v.u = v.v = 1;</code></li>
 * 	<li>access as color (r,g,b,a) &mdash; <code>v.r = v.g = v.b = v.a = 1;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = v[2] = v[3] = 1;</code></li>
 * </ul>
 */
template<class T>
class rpVector4D
{

public:

   //-------------------- Attributes --------------------//

    union
    {
        /**
         * First element of vector, alias for R-coordinate.
         * For color notation.
         */
        T r
        /**
         * First element of vector, alias for X-coordinate.
         */;
        T x;
    };

    union
    {
        /**
         * Second element of vector, alias for G-coordinate.
         * For color notation.
         */
        T g;
        /**
         * Second element of vector, alias for Y-coordinate.
         */
        T y;
    };

    union
    {
        /**
         * Third element of vector, alias for B-coordinate.
         * For color notation.
         */
        T b;
        /**
         * Third element of vector, alias for Z-coordinate.
         */
        T z;
    };

    union
    {
        /**
         * Fourth element of vector, alias for A-coordinate.
         * For color notation. This represnt aplha chanell
         */
        T a;
        /**
         * First element of vector, alias for W-coordinate.
         * @note For vectors (such as normals) should be set to 0.0
         * For vertices should be set to 1.0
         */
        T w;
    };

public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and sets to (0,0,0,0)
     */
    SIMD_INLINE rpVector4D()
       : x(0), y(0), z(0), w(0)
    {
    }

    /**
     * Creates and sets to (x,y,z,z)
     * @param nx initial x-coordinate value (R)
     * @param ny initial y-coordinate value (G)
     * @param nz initial z-coordinate value (B)
     * @param nw initial w-coordinate value (Alpha)
     */
    SIMD_INLINE rpVector4D(T nx, T ny, T nz, T nw)
       : x(nx), y(ny), z(nz), w(nw)
    {
    }

    /**
     * Copy constructor.
     * @param src Source of data for new created rpVector4D instance.
     */
    SIMD_INLINE rpVector4D(const rpVector4D<T>& src)
       : x(src.x), y(src.y), z(src.z), w(src.w)
    {
    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created rpVector4D instance.
     */
    template<class FromT>
    SIMD_INLINE rpVector4D(const rpVector4D<FromT>& src)
      : x(static_cast<T>(src.x)), y(static_cast<T>(src.y)), z(static_cast<T>(src.z)), w(static_cast<T>(src.w))
    {
    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void setToZero()
    {
      x = T(0);
      y = T(0);
      z = T(0);
      w = T(0);
    }


    SIMD_INLINE void setAllValues(T newX, T newY, T newZ, T newW)
    {
      x = newX;
      y = newY;
      z = newZ;
      w = newW;
    }


    SIMD_INLINE T getX() const { return x; }
    SIMD_INLINE T getY() const { return y; }
    SIMD_INLINE T getZ() const { return z; }
    SIMD_INLINE T getW() const { return w; }


    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }
    SIMD_INLINE void SetZ(T _z) { z = _z; }
    SIMD_INLINE void SetW(T _w) { w = _w; }

    //----------------[ access operators ]-------------------
    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator=(const rpVector4D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    /**
     * Copy casting operator
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE rpVector4D<T> operator=(const rpVector4D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        z = static_cast<T>(rhs.z);
        w = static_cast<T>(rhs.w);
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to w coordinate.
     */
    SIMD_INLINE T & operator[](int n)
    {
        static_assert(sizeof(*this) == sizeof(T[4]), "");
        assert(n >= 0 && n < 4);
        return (&x)[n];
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to w coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
    	static_assert(sizeof(*this) == sizeof(T[4]), "");
    	assert(n >= 0 && n < 4);
    	return (&x)[n];
    }

    //---------------[ vector aritmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator+(const rpVector4D<T>& rhs) const
    {
        return rpVector4D<T>(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator-(const rpVector4D<T>& rhs) const
    {
        return rpVector4D<T>(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator*(const rpVector4D<T> rhs) const
    {
        return rpVector4D<T>(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator/(const rpVector4D<T>& rhs) const
    {
        return rpVector4D<T>(x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator+=(const rpVector4D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator-=(const rpVector4D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator*=(const rpVector4D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator/=(const rpVector4D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        w /= rhs.w;
        return *this;
    }

    //--------------[ equiality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < EPSILON,
     * same for y-coordinate, z-coordinate, and w-coordinate.
     */
    SIMD_INLINE bool operator==(const rpVector4D<T>& rhs) const
    {
        return Abs(x - rhs.x) < EPSILON &&
               Abs(y - rhs.y) < EPSILON &&
               Abs(z - rhs.z) < EPSILON &&
               Abs(w - rhs.w) < EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const rpVector4D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE rpVector4D<T> operator-() const
    {
        return rpVector4D<T>(-x, -y, -z, -w);
    }

    //--------------[ scalar vector operator ]--------------------

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator+(T rhs) const
    {
        return rpVector4D<T>(x + rhs, y + rhs, z + rhs, w + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator-(T rhs) const
    {
        return rpVector4D<T>(x - rhs, y - rhs, z - rhs, w - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator*(T rhs) const
    {
        return rpVector4D<T>(x * rhs, y * rhs, z * rhs, w * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> operator/(T rhs) const
    {
        return rpVector4D<T>(x / rhs, y / rhs, z / rhs, w / rhs);
    }

    //------------------------------ Friends ----------------------------------------//

     /**
      * Multiplication operator
      * @param rhs Right hand side argument of binary operator.
      */
     friend SIMD_INLINE rpVector4D<T> operator*(T number, const rpVector4D<T>& vector)
     {
         return rpVector4D<T>(number * vector.x, number * vector.y, number * vector.z , number * vector.w);
     }


     /**
      * Division operator
      * @param rhs Right hand side argument of binary operator.
      */
     friend SIMD_INLINE  rpVector4D<T> operator/( T number , const rpVector4D<T>& vector )
     {
         return rpVector4D<T>(vector.x / number, vector.y / number, vector.z / number , vector.w / number);
     }


     //------------------------------ Dynamics ----------------------------------------//

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        w += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        w -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        w /= rhs;
        return *this;
    }

    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T dot(const rpVector4D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w;
    }

    /**
     * Cross tri product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpVector4D<T> cross(const rpVector4D<T>& b , const rpVector4D<T>& c) const
    {

        //Precompute some 2x2 matrix determinants for speed
         T Pxy = b.x*c.y - c.x*b.y;
         T Pxz = b.x*c.z - c.x*b.z;
         T Pxw = b.x*c.w - c.x*b.w;
         T Pyz = b.y*c.z - c.y*b.z;
         T Pyw = b.y*c.w - c.y*b.w;
         T Pzw = b.z*c.w - c.z*b.w;

          return rpVector4D<T>
          (
             y*Pzw - z*Pyw + w*Pyz,    //Note the lack of 'x' in this line
             z*Pxw - x*Pzw - w*Pxz,    //y, Etc.
             x*Pyw - y*Pxw + w*Pxy,
             y*Pxz - x*Pyz - z*Pxy
          );
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
        return x * x + y * y + z * z + w * w;
    }


    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T length() const
    {
        return (T) Sqrt(x * x + y * y + z * z + w * w);
    }

    /**
     * Normalize vector
     */
    SIMD_INLINE void normalize()
    {
        T s = length();
        x /= s;
        y /= s;
        z /= s;
        w /= s;
    }

    /**
     * Normalize Unit vector
     */
    SIMD_INLINE rpVector4D<T> getUnit() const
    {
        T lengthVector = length();
        if (lengthVector < EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return rpVector4D<T>( x * lengthInv,
                              y * lengthInv,
                              z * lengthInv,
                              w * lengthInv);
    }

    /**
    * Inverse vector
    */
    SIMD_INLINE rpVector4D<T> getInverse() const
    {
        return rpVector4D<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) , T(1.0/w));
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
    SIMD_INLINE rpVector4D<T> lerp(T fact, const rpVector4D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }

    //-------------[ conversion ]-----------------------------

    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class rpVector4D<T>)
     * used for passing rpVector4D<T> values to gl*4[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return (T*) this;
    }

    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class rpVector4D<T>)
     * used for passing rpVector4D<T> values to gl*4[fd] functions.
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
    friend std::ostream& operator<<(std::ostream& lhs, const rpVector4D<T>& rhs)
    {
        lhs << "[" << rhs[0] << "," << rhs[1] << "," << rhs[2] << "," << rhs[3] << "]";
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
    static const rpVector4D<T> IDENTITY;

    /**
     * The additive identity vector.
     */
    static const rpVector4D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const rpVector4D<T> X;
    /**
     * The identity vector Y.
     */
    static const rpVector4D<T> Y;

    /**
     * The identity vector Z.
     */
    static const rpVector4D<T> Z;

    /**
     * The identity vector W.
     */
    static const rpVector4D<T> W;

};


template<class T> const rpVector4D<T> rpVector4D<T>::IDENTITY(1.0, 1.0, 1.0, 1.0);
template<class T> const rpVector4D<T> rpVector4D<T>::ZERO(0.0, 0.0, 0.0, 0.0);
template<class T> const rpVector4D<T> rpVector4D<T>::X(1.0, 0.0, 0.0, 0.0);
template<class T> const rpVector4D<T> rpVector4D<T>::Y(0.0, 1.0, 0.0, 0.0);
template<class T> const rpVector4D<T> rpVector4D<T>::Z(0.0, 0.0, 1.0, 0.0);
template<class T> const rpVector4D<T> rpVector4D<T>::W(0.0, 0.0, 0.0, 1.0);


template<class T> const
static rpVector4D<T> cross(const rpVector4D<T>& a, const rpVector4D<T>& b , const rpVector4D<T>& c)
{
    return a.cross(b,c);
}

template<class T> const
static T dot(const rpVector4D<T>& a, const rpVector4D<T>& b)
{
    return a.dot(b);
}


/// Three dimensional Vector of floats
typedef rpVector4D<float> rpVector4Df;
/// Three dimensional Vector of doubles
typedef rpVector4D<double> rpVector4Dd;
/// Three dimensional Vector of ints
typedef rpVector4D<int> rpVector4Di;


} /* namespace */

#endif // RPVECTOR4D_H
