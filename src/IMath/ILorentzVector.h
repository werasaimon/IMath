#ifndef ILORENTZVECTOR4D_H
#define ILORENTZVECTOR4D_H

#include "IVector3D.h"

namespace IMath
{



template<class T> class  ILorentzVector
{

  public:

    //---------------- attribute -------------------//

    // 3 vector component space
    T x;
    T y;
    T z;
    // 1 component time or energy of (x,y,z,t) or (px,py,pz,e)
    T t;


  public:

    // Constructor of the class Vector4D
    SIMD_INLINE ILorentzVector()
    : x(0.0), y(0.0), z(0.0), t(1.0)
    {

    }

    // Constructor with arguments
    SIMD_INLINE ILorentzVector( T newX, T newY, T newZ , T newT )
    : x(newX), y(newY), z(newZ) , t(newT)
    {

    }

    // Copy-constructor
    SIMD_INLINE ILorentzVector(const ILorentzVector<T>& vector)
    : x(vector.x), y(vector.y), z(vector.z) , t(vector.t)
    {

    }


    // Copy-constructor
    SIMD_INLINE ILorentzVector(const IVector3D<T>& vector , T time)
    : x(vector.x), y(vector.y), z(vector.z) , t(time)
    {

    }

    // Copy-constructor
    SIMD_INLINE ILorentzVector( T time , const IVector3D<T>& vector )
    : x(vector.x), y(vector.y), z(vector.z) , t(time)
    {

    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void setToZero()
    {
      x = T(0);
      y = T(0);
      z = T(0);
      t = T(0);
    }


    SIMD_INLINE void setAllValues(T newX, T newY, T newZ, T newT)
    {
      x = newX;
      y = newY;
      z = newZ;
      t = newT;
    }


    SIMD_INLINE T getX() const { return x; }
    SIMD_INLINE T getY() const { return y; }
    SIMD_INLINE T getZ() const { return z; }
    SIMD_INLINE T getT() const { return t; }


    SIMD_INLINE T getP()  const { return (x*x + y*y + z*z); }
    SIMD_INLINE T getE()  const { return t; }
    SIMD_INLINE T getEnergy()  const { return t; }


    SIMD_INLINE void setX(T _x) { x = _x; }
    SIMD_INLINE void setY(T _y) { y = _y; }
    SIMD_INLINE void setZ(T _z) { z = _z; }
    SIMD_INLINE void setW(T _t) { t = _t; }


    SIMD_INLINE void setXYZM(T x, T y, T z, T m)
    {
        if ( m  >= 0 )
        {
            setAllValues( x, y, z, ISqrt(x*x+y*y+z*z+m*m) );
        }
        else
        {
            setAllValues( x, y, z, ISqrt( IMax((x*x+y*y+z*z-m*m), 0. ) ) );
        }
    }


    SIMD_INLINE  void setPtEtaPhiM(T pt, T eta, T phi, T m)
    {
        pt = Abs(pt);
        SetXYZM(pt*Cos(phi), pt*ISin(phi), pt*ISinh(eta) ,m);
    }

    SIMD_INLINE void setPtEtaPhiE(T pt, T eta, T phi, T e)
    {
        pt = Abs(pt);
        SetXYZT(pt*Cos(phi), pt*ISin(phi), pt*ISinh(eta) ,e);
    }


    SIMD_INLINE void setVector3(const IVector3D<T> &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }


    /// System coordinate . 3 vector component
    SIMD_INLINE IVector3D<T> vect() const
    {
        return IVector3D<T>(x,y,z);
    }


    ///  Project on space-3D
    SIMD_INLINE IVector3D<T> getBoostVector() const
    {
       return IVector3D<T>(x/t, y/t, z/t);
    }

    /// Return the square of the length of the vector
    /// Metrices Minkowski Space
    SIMD_INLINE T lengthSquare() const
    {
    	const T c = LIGHT_MAX_VELOCITY_C;
         return (c*c) * (t*t) - (x*x + y*y + z*z);
    }

    /// Return the length of the vector
    SIMD_INLINE T length() const
    {
        T mm = lengthSquare();
        return mm < 0.0 ? -ISqrt(-mm) : ISqrt(mm);
    }


    SIMD_INLINE T getBeta() const
    {
       return ISqrt(x*x + y*y + z*z) / t;
    }

    SIMD_INLINE T getGamma() const
    {
       T b = getBeta();
       return 1.0/ISqrt(1- b*b);
    }


    SIMD_INLINE void Boost(T bx, T by, T bz)
    {
       /**
       //Boost this Lorentz vector
       T b2 = bx*bx + by*by + bz*bz;
       T gamma = 1.0 / ISqrt(1.0 - b2);
       T bp = bx*x + by*y + bz*z;
       T gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

       SetX(getX() + gamma2*bp*bx + gamma*bx*getT());
       SetY(getY() + gamma2*bp*by + gamma*by*getT());
       SetZ(getZ() + gamma2*bp*bz + gamma*bz*getT());
       SetT(gamma*(getT() + bp));
       **/

    	/// method Gerglocema
        *this = createBoost( *this , IVector3D<T>(bx,by,bz) );
    }

    SIMD_INLINE void Boost( const IVector3D<T> &b)
    {
    	/// method Gerglocema
    	*this = createBoost( *this , b );
    }


    SIMD_INLINE T Rapidity() const
    {
       //return rapidity
       return 0.5*log( (getE()+getZ()) / (getE()-getZ()) );
    }



    /**
     * Normalize Unit vector
     */
    SIMD_INLINE ILorentzVector<T> getUnit() const
    {
        T lengthVector = length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return ILorentzVector<T>( x * lengthInv ,
                                   y * lengthInv ,
                                   z * lengthInv ,
                                   t * lengthInv );
    }


    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE ILorentzVector<T> normalized() const
    {
        T lengthVector = length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return ILorentzVector<T>( x * lengthInv ,
                                   y * lengthInv ,
                                   z * lengthInv ,
                                   t * lengthInv );
    }


    /**
    * Inverse vector
    */
    SIMD_INLINE ILorentzVector<T> getInverse() const
    {
        return ILorentzVector<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) , T(1.0/t));
    }





    //---------------[ vector aritmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator+(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x + rhs.x, y + rhs.y, z + rhs.z, t + rhs.t);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator-(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x - rhs.x, y - rhs.y, z - rhs.z, t - rhs.t);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator*(const ILorentzVector<T> rhs) const
    {
        return ILorentzVector<T>(x * rhs.x, y * rhs.y, z * rhs.z, t * rhs.t);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator/(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x / rhs.x, y / rhs.y, z / rhs.z, t / rhs.t);
    }

    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE ILorentzVector<T> operator*(T number, const ILorentzVector<T>& vector)
    {
        return ILorentzVector<T>(number * vector.x, number * vector.y, number * vector.z , number * vector.t);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  ILorentzVector<T> operator/( T number , const ILorentzVector<T>& vector )
    {
        return ILorentzVector<T>(vector.x / number, vector.y / number, vector.z / number , vector.t / number);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator+=(const ILorentzVector<T>& rhs)
    {
    	x += rhs.x;
    	y += rhs.y;
    	z += rhs.z;
    	t += rhs.t;
    	return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator-=(const ILorentzVector<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        t -= rhs.t;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator*=(const ILorentzVector<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        t *= rhs.t;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator/=(const ILorentzVector<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        t /= rhs.t;
        return *this;
    }


    //--------------[ equiality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold MACHINE_EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < MACHINE_EPSILON,
     * same for y-coordinate, z-coordinate, and w-coordinate.
     */
    SIMD_INLINE bool operator==(const ILorentzVector<T>& rhs) const
    {
        return IAbs(x - rhs.x) < MACHINE_EPSILON &&
               IAbs(y - rhs.y) < MACHINE_EPSILON &&
               IAbs(z - rhs.z) < MACHINE_EPSILON &&
               IAbs(t - rhs.t) < MACHINE_EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const ILorentzVector<T>& rhs) const
    {
    	return !(*this == rhs);
    }


    //----------------[ access operators ]-------------------
     /**
      * Copy operator
      * @param rhs Right hand side argument of binary operator.
      */
    SIMD_INLINE ILorentzVector<T> operator=(const ILorentzVector<T>& rhs)
     {
         x = rhs.x;
         y = rhs.y;
         z = rhs.z;
         t = rhs.t;
         return *this;
     }

     /**
      * Copy casting operator
      * @param rhs Right hand side argument of binary operator.
      */
     template<class FromT>
     SIMD_INLINE ILorentzVector<T> operator=(const ILorentzVector<FromT>& rhs)
     {
         x = static_cast<T>(rhs.x);
         y = static_cast<T>(rhs.y);
         z = static_cast<T>(rhs.z);
         t = static_cast<T>(rhs.t);
         return *this;
     }


    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to t-time coordinate.
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
     * else reference to t time-coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
    	static_assert(sizeof(*this) == sizeof(T[4]), "");
    	assert(n >= 0 && n < 4);
    	return (&x)[n];
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE ILorentzVector<T> operator-() const
    {
        return ILorentzVector<T>(-x, -y, -z, -t);
    }

    //--------------[ scalar vector operator ]--------------------

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator+(T rhs) const
    {
        return ILorentzVector<T>(x + rhs, y + rhs, z + rhs, t + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator-(T rhs) const
    {
        return ILorentzVector<T>(x - rhs, y - rhs, z - rhs, t - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator*(T rhs) const
    {
        return ILorentzVector<T>(x * rhs, y * rhs, z * rhs, t * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator/(T rhs) const
    {
        return ILorentzVector<T>(x / rhs, y / rhs, z / rhs, t / rhs);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        t += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        t -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        t *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        t /= rhs;
        return *this;
    }


    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T dot(const ILorentzVector<T>& rhs) const
    {
         return t*rhs.t - z*rhs.z - y*rhs.y - x*rhs.x;
    }

    /**
     * Cross tri product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> cross(const ILorentzVector<T>& b , const ILorentzVector<T>& c) const
    {

        //Precompute some 2x2 matrix determinants for speed
         T Pxy = b.x*c.y - c.x*b.y;
         T Pxz = b.x*c.z - c.x*b.z;
         T Pxw = b.x*c.t - c.x*b.t;
         T Pyz = b.y*c.z - c.y*b.z;
         T Pyw = b.y*c.t - c.y*b.t;
         T Pzw = b.z*c.t - c.z*b.t;

          return ILorentzVector<T>
          (
             y*Pzw - z*Pyw + t*Pyz,    //Note the lack of 'x' in this line
             z*Pxw - x*Pzw - t*Pxz,    //y, Etc.
             x*Pyw - y*Pxw + t*Pxy,
             y*Pxz - x*Pyz - z*Pxy
          );
    }

    //================================ Plugin =====================================//

    /**
     *  Method Gerglocema
     */
    static SIMD_INLINE ILorentzVector<T> createBoost( const ILorentzVector<T> &pos , const IVector3D<T> &v)
    {
    	///Light speed
    	const T c = LIGHT_MAX_VELOCITY_C;

    	/// gamma factor
        float gamma = ISqrt( 1.0 - v.dot(v) / (c*c));

    	/// method Gerglocema
        ILorentzVector<T> res;
        IVector3D<T> ortogonalPredikat = (pos.vect().cross(v)).cross(v);
    	res.x = ((pos.x + v.x * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.x;
    	res.y = ((pos.y + v.y * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.y;
    	res.z = ((pos.z + v.z * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.z;
    	res.t =  (pos.t + (v.dot(pos.vect())/(c*c))) / gamma ;
    	return res;
    }


    //-------------[ output operator ]------------------------
     /**
     * Output to stream operator
     * @param lhs Left hand side argument of operator (commonly ostream instance).
     * @param rhs Right hand side argument of operator.
     * @return Left hand side argument - the ostream object passed to operator.
     */
     friend std::ostream& operator<<(std::ostream& lhs, const ILorentzVector<T>& rhs)
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

};


template<class T> const
static ILorentzVector<T> cross(const ILorentzVector<T>& a, const ILorentzVector<T>& b , const ILorentzVector<T>& c)
{
    return a.cross(b,c);
}

template<class T> const
static T dot(const ILorentzVector<T>& a, const ILorentzVector<T>& b)
{
    return a.dot(b);
}


//--------------------------------------
// Typedef shortcuts for 4D LorentzVector
//-------------------------------------
/// Three dimensional LorentzVector of floats
typedef ILorentzVector<float> ILorentzVectorf;
/// Three dimensional LorentzVector of doubles
typedef ILorentzVector<double> ILorentzVectord;
/// Three dimensional LorentzVector of ints
typedef ILorentzVector<int> ILorentzVectori;



} /* namespace */

#endif // ILORENTZVECTOR4D_H
