#ifndef IVECTOR3D_H_
#define IVECTOR3D_H_


#include "iFunc.h"

namespace IMath
{


template<class T> class  IMatrix3x3;

/**
 * Class for three dimensional vector.
 * There are four ways of accessing vector components.
 * Let's have <code>Vector3f v</code>, you can either:
 * <ul>
 * 	<li>access as position (x,y,z) &mdash; <code>v.x = v.y = v.z = 1;</code></li>
 * 	<li>access as texture coordinate (s,t,u) &mdash; <code>v.s = v.t = v.u = 1;</code></li>
 * 	<li>access as color (r,g,b) &mdash; <code>v.r = v.g = v.b = 1;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = v[2] = 1;</code></li>
 * </ul>
 **/
template<class T>
class IVector3D
{

public:

    //-------------------- Attributes --------------------//

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

        /**
         * First element of vector, alias for R-coordinate.
         * For color notation.
         */
        T r;
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
        /**
         * Second element of vector, alias for G-coordinate.
         * For color notation.
         */
        T g;
    };

    union
    {
        /**
         * Third element of vector, alias for Z-coordinate.
         */
        T z;

        /**
         * Third element of vector, alias for W-coordinate.
         * For textures notation.
         */
        T w;
        /**
         * Third element of vector, alias for B-coordinate.
         * For color notation.
         */
        T b;
    };


public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and sets to (0,0,0)
     */
    SIMD_INLINE IVector3D()
            : x(0), y(0), z(0)
    {
    }

    /**
     * Creates and sets to (x,y,z)
     * @param nx initial x-coordinate value
     * @param ny initial y-coordinate value
     * @param nz initial z-coordinate value
     */
    SIMD_INLINE IVector3D(T nx, T ny, T nz)
            : x(nx), y(ny), z(nz)
    {
    }

    /**
     * Copy constructor.
     * @param src Source of data for new created IVector3D instance.
     */
    SIMD_INLINE IVector3D(const IVector3D<T>& src)
            : x(src.x), y(src.y), z(src.z)
    {
    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created IVector3D instance.
     */
    template<class FromT>
    SIMD_INLINE IVector3D(const IVector3D<FromT>& src)
            : x(static_cast<T>(src.x)), y(static_cast<T>(src.y)), z(static_cast<T>(src.z))
    {
    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void setToZero()
    {
        x = T(0);
        y = T(0);
        z = T(0);
    }

    SIMD_INLINE bool isZero() const
    {
        return IApproxEqual<T>( lengthSquare(), T(0.0) );
    }


    SIMD_INLINE void setAllValues(T newX, T newY, T newZ)
    {
        x = newX;
        y = newY;
        z = newZ;
    }


    SIMD_INLINE T getX() const { return x; }
    SIMD_INLINE T getY() const { return y; }
    SIMD_INLINE T getZ() const { return z; }


    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }
    SIMD_INLINE void SetZ(T _z) { z = _z; }


    //----------------[ access operators ]-------------------
    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator=(const IVector3D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    /**
     * Copy casting operator.
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE IVector3D<T> operator=(const IVector3D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        z = static_cast<T>(rhs.z);
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y, else reference to z
     * y coordinate.
     */
    SIMD_INLINE T & operator[](int n)
    {
    	static_assert(sizeof(*this) == sizeof(T[3]), "");
    	assert(n >= 0 && n < 3);
    	return (&x)[n];
    }

    /**
     * Constant array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y, else reference to z
     * y coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
    	static_assert(sizeof(*this) == sizeof(T[3]), "");
    	assert(n >= 0 && n < 3);
    	return (&x)[n];
    }

    //---------------[ vector arithmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator+(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator-(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x * rhs.x, y * rhs.y, z * rhs.z);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator/(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x / rhs.x, y / rhs.y, z / rhs.z);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator+=(const IVector3D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator-=(const IVector3D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator*=(const IVector3D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator/=(const IVector3D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        return *this;
    }


    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T dot(const IVector3D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    /**
     * Cross product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> cross(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(y * rhs.z - rhs.y * z,
                             z * rhs.x - rhs.z * x,
                             x * rhs.y - rhs.x * y);
    }

    //--------------[ scalar vector operator ]--------------------
    /**
    * Addition operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE IVector3D<T> operator+(T rhs) const
    {
        return IVector3D<T>(x + rhs, y + rhs, z + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator-(T rhs) const
    {
        return IVector3D<T>(x - rhs, y - rhs, z - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator*(T rhs) const
    {
        return IVector3D<T>(x * rhs, y * rhs, z * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator/(T rhs) const
    {
        return IVector3D<T>(x / rhs, y / rhs, z / rhs);
    }


    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector3D<T> operator*(T number, const IVector3D<T>& vector)
    {
        return IVector3D<T>(number * vector.x, number * vector.y, number * vector.z);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector3D<T> operator/( T number , const IVector3D<T>& vector )
    {
        return IVector3D<T>(vector.x / number, vector.y / number, vector.z / number);
    }


    //------------------------------------------------------------------------------//

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        return *this;
    }

    //--------------[ Equality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < EPSILON,
     * same for y-coordinate, and z-coordinate.
     */
    SIMD_INLINE bool operator==(const IVector3D<T>& rhs) const
    {
        return IAbs(x - rhs.x) < MACHINE_EPSILON &&
               IAbs(y - rhs.y) < MACHINE_EPSILON &&
               IAbs(z - rhs.z) < MACHINE_EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const IVector3D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE IVector3D<T> operator-() const
    {
        return IVector3D<T>(-x, -y, -z);
    }

    //-------------[ size operations ]---------------------------
    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T length() const
    {
        return (T) ISqrt(x * x + y * y + z * z);
    }

    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length() method.
     */
    SIMD_INLINE T lengthSquare() const
    {
        return x * x + y * y + z * z;
    }


    SIMD_INLINE int getMinAxis() const
    {
        return (x < y ? (x < z ? 0 : 2) : (y < z ? 1 : 2));
    }

    SIMD_INLINE int getMaxAxis() const
    {
        return (x < y ? (y < z ? 2 : 1) : (x < z ? 2 : 0));
    }

    SIMD_INLINE T getMinValue() const
    {
      return IMin(IMin(x, y), z);
    }

    SIMD_INLINE T getMaxValue() const
    {
      return IMax(IMax(x, y), z);
    }

    //*********************************************//

    /**
     * Normalize vector
     */
    SIMD_INLINE void normalize()
    {
        T s = length();
        x /= s;
        y /= s;
        z /= s;
    }


    /**
     * Normalize unit vector
     */
    SIMD_INLINE IVector3D<T> getUnit() const
    {
        T lengthVector = length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector3D<T>( x * lengthInv,
                              y * lengthInv,
                              z * lengthInv);
    }



    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE IVector3D<T> normalized() const
    {
        T lengthVector = length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector3D<T>( x * lengthInv,
                              y * lengthInv,
                              z * lengthInv);
    }

    /**
    * Inverse vector
    */
    SIMD_INLINE IVector3D<T> getInverse() const
    {
        return IVector3D<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) );
    }

    /**
    * Orthogonal unit vector
    */
    SIMD_INLINE IVector3D<T> getOneUnitOrthogonalVector() const
    {
          assert(length() > MACHINE_EPSILON);

          // Get the minimum element of the vector
          IVector3D<T> vectorAbs(IAbs(x), IAbs(y), IAbs(z));
          int minElement = vectorAbs.getMinAxis();

          if (minElement == 0)
          {
             return IVector3D<T>(0.0, -z, y) / ISqrt(y*y + z*z);
          }
          else if (minElement == 1)
          {
             return IVector3D<T>(-z, 0.0, x) / ISqrt(x*x + z*z);
          }
          else
          {
             return IVector3D<T>(-y, x, 0.0) / ISqrt(x*x + y*y);
          }
    }


    //------------[ other operations ]---------------------------
    /**
     * Rotate vector around three axis.
     * @param ax Angle (in degrees) to be rotated around X-axis.
     * @param ay Angle (in degrees) to be rotated around Y-axis.
     * @param az Angle (in degrees) to be rotated around Z-axis.
     */
    SIMD_INLINE void rotate(T ax, T ay, T az)
    {
        T a = ICos(IDegreesToRadians(ax));
        T b = ISin(IDegreesToRadians(ax));
        T c = ICos(IDegreesToRadians(ay));
        T d = ISin(IDegreesToRadians(ay));
        T e = ICos(IDegreesToRadians(az));
        T f = ISin(IDegreesToRadians(az));
        T nx = c * e * x - c * f * y + d * z;
        T ny = (a * f + b * d * e) * x + (a * e - b * d * f) * y - b * c * z;
        T nz = (b * f - a * d * e) * x + (a * d * f + b * e) * y + a * c * z;
        x = nx;
        y = ny;
        z = nz;

    }

    /**
     * Linear interpolation of two vectors
     * @param fact Factor of interpolation. For translation from positon
     * of this vector to vector r, values of factor goes from 0.0 to 1.0.
     * @param r Second Vector for interpolation
     * @note However values of fact parameter are reasonable only in interval
     * [0.0 , 1.0], you can pass also values outside of this interval and you
     * can get result (extrapolation?)
     */
    SIMD_INLINE IVector3D<T> lerp(T fact, const IVector3D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }




    SIMD_INLINE T getAngleBetween( const IVector3D<T> &Vector2 ) const
    {
        IVector3D<T> Vector1(*this);
        T dotProduct = Vector1.dot(Vector2);
        T vectorsMagnitude = (Vector1.length()) * (Vector2.length());
        T angle = IArcCos(dotProduct / vectorsMagnitude);
        if( is_nan(angle)) return 0;
        return (angle);
    }

    //-------------[ conversion ]-----------------------------

    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class IVector3D<T>)
     * used for passing IVector3D<T> values to gl*3[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return (T*) this;
    }

    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class IVector3D<T>)
     * used for passing IVector3D<T> values to gl*3[fd] functions.
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
    friend std::ostream& operator<<(std::ostream& lhs, const IVector3D<T> rhs)
    {
        lhs << "[" << rhs[0] << "," << rhs[1] << "," << rhs[2] << "]";
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



    //========================== plugins =========================//


    /**
     * Compute the determinant of a matrix whose columns are three given vectors.
     * Useful property: det(a, b, c) = det(c, a, b) = det(b, c, a).
     */
    static SIMD_INLINE T Determinant(const IVector3D<T>& a, const IVector3D<T>& b, const IVector3D<T>& c)
    {
        return a.dot(b.cross(c));
    }


    static SIMD_INLINE IVector3D<T> clamp(const IVector3D<T>& vector, T maxLength)
    {
        if (vector.lengthSquare() > maxLength * maxLength)
        {
            return vector.getUnit() * maxLength;
        }

        return vector;
    }


    static T SIMD_INLINE AngleSigned(IVector3D<T> v1, IVector3D<T> v2, IVector3D<T> normal)
    {
        return IAtan2( normal.dot(v1.cross(v2)), v1.dot(v2));
    }


    /**
     * @brief triNormal
     * @param V0
     * @param V1
     * @param V2
     * @return Face_Normal
     */
    static  SIMD_INLINE  IVector3D<T> triNormal( const IVector3D<T>& V0,
                                                  const IVector3D<T>& V1,
                                                  const IVector3D<T>& V2)
    {
        IVector3D<T> Norm;
        IVector3D<T> E = V1;
        IVector3D<T> F = V2;
        E -= V0;
        F -= V1;
        Norm = E.cross(F);
        Norm.normalize();
        return Norm;
    }

    /**
     * @brief BiUnitGrammSchmidt
     * @param n -axis
     * @param p -ortogonal_left
     * @param q -ortogonal_up
     * Gramm schmidt process
     * https://www.math.hmc.edu/calculus/tutorials/gramschmidt/gramschmidt.pdf
     */
    static SIMD_INLINE void BiUnitGrammSchmidt(const IVector3D<T>& n_axis , IVector3D<T>& p, IVector3D<T>& q )
    {
        p = n_axis.getOneUnitOrthogonalVector();
        q = n_axis.cross(p);
    }


    #define btRecipSqrt(x)  (1.0/ISqrt(x))
    #define SIMDSQRT12      (0.7071067811865475244008443621048490)
    /// Bullet physics version (Gramm schmidt process)
    static SIMD_INLINE void BiUnitOrthogonalVector(const IVector3D<T>& n, IVector3D<T>& p, IVector3D<T>& q)
    {
        if (IAbs(n[2]) > SIMDSQRT12)
        {
            // choose p in y-z plane
            T a = n[1]*n[1] + n[2]*n[2];
            T k = btRecipSqrt(a);
            p[0] = T(0);
            p[1] = -n[2]*k;
            p[2] = n[1]*k;
            // set q = n x p
            q[0] = a*k;
            q[1] = -n[0]*p[2];
            q[2] = n[0]*p[1];
        }
        else
        {
            // choose p in x-y plane
            T a = n[0]*n[0] + n[1]*n[1];
            T k = btRecipSqrt(a);
            p[0] = -n[1]*k;
            p[1] = n[0]*k;
            p[2] = T(0);
            // set q = n x p
            q[0] = -n[2]*p[1];
            q[1] = n[2]*p[0];
            q[2] = a*k;
        }
    }



public:
    /**
     * The multiplicitive identity vector
     */
    static const IVector3D<T> IDENTITY;
    /**
     * The additive identity vector.
     */
    static const IVector3D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const IVector3D<T> X;
    /**
     * The identity vector Y.
     */
    static const IVector3D<T> Y;
    /**
     * The identity vector Z.
     */
    static const IVector3D<T> Z;


};

template<class T> const IVector3D<T> IVector3D<T>::IDENTITY(1.0, 1.0, 1.0);
template<class T> const IVector3D<T> IVector3D<T>::ZERO(0.0, 0.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::X(1.0, 0.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::Y(0.0, 1.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::Z(0.0, 0.0, 1.0);


template<class T> const
static IVector3D<T> cross(const IVector3D<T>& a, const IVector3D<T>& b)
{
    return a.cross(b);
}

template<class T> const
static T dot(const IVector3D<T>& a, const IVector3D<T>& b)
{
    return a.dot(b);
}

//--------------------------------------
// Typedef shortcuts for 3D vector
//-------------------------------------
/// Three dimensional Vector of floats
typedef IVector3D<float> IVector3f;
/// Three dimensional Vector of doubles
typedef IVector3D<double> IVector3d;
/// Three dimensional Vector of ints
typedef IVector3D<int> IVector3i;


} /* namespace */



#endif /* IVECTOR3D_H_ */
