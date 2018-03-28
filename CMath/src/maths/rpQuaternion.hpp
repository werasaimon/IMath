
#ifndef SOURCE_REAL_PHYSICS_LINEARMATHS_RPQUATERNION_H_
#define SOURCE_REAL_PHYSICS_LINEARMATHS_RPQUATERNION_H_


// Libraries
#include <cmath>
#include "rpMatrix4x4.hpp"

namespace CMath
{
  
/**
 * Quaternion class implementing some quaternion algebra operations.
 * Quaternion is kind of complex number it consists of its real part (w)
 * and its complex part v. This complex part has three elements, so we
 * can express it as xi + yj + zk . Note that coordinates of (x,y,z) are
 * hold inside v field.
 */
template<class T> class rpQuaternion
{


private:

  // -------------------- Attributes -------------------- //


    /**
     * Real part of rpQuaternion.
     */
    T w;

    /**
     * Imaginary part of rpQuaternion.
     */
    rpVector3D<T> v;





public:

     //--------------------------[ constructors ]------------------------------- //

    /**
    * rpQuaternion constructor, sets rpQuaternion to (0 + 0i + 0j + 0k).
    */
    SIMD_INLINE rpQuaternion()
        : w(0), v(0, 0, 0)
    {
    }

    /**
    * Copy constructor.
    */
    SIMD_INLINE rpQuaternion(const rpQuaternion<T>& q)
        : w(q.w), v(q.v)
    {
    }

    /**
     * Copy casting constructor.
     */
    template<class FromT>
    SIMD_INLINE rpQuaternion(const rpQuaternion<FromT>& q)
        : w(static_cast<T>(q.w)), v(q.v)
    {
    }

    /**
    * Creates rpQuaternion object from real part w_ and complex part v_.
    * @param w_ Real part of rpQuaternion.
    * @param v_ Complex part of rpQuaternion (xi + yj + zk).
    */
    SIMD_INLINE rpQuaternion(T w_, const rpVector3D<T>& v_)
        : w(w_), v(v_)
    {
    }

    /**
    * Creates rpQuaternion object from real part w_ and complex part v_.
    * @param w_ Real part of rpQuaternion.
    * @param v_ Complex part of rpQuaternion (xi + yj + zk).
    */
    SIMD_INLINE rpQuaternion(const rpVector3D<T>& v_, T w_)
        : w(w_), v(v_)
    {
    }

    /**
    * Creates rpQuaternion object from value (w_ + xi + yj + zk).
    * @param x Complex coefficient for i complex constant.
    * @param y Complex coefficient for j complex constant.
    * @param z Complex coefficient for k complex constant.
    * @param w_ Real part of rpQuaternion.
    */
    SIMD_INLINE rpQuaternion(T _x, T _y, T _z, T _w)
     : w(_w), v(_x, _y, _z)
    {
    }


  /// Constructor which convert Euler angles (in radians) to a quaternion
  SIMD_INLINE rpQuaternion(T angleX, T angleY, T angleZ)
  {
     initWithEulerAngles(angleX, angleY, angleZ);
  }

  /// Constructor which convert Euler angles (in radians) to a quaternion
  SIMD_INLINE rpQuaternion(const rpVector3D<T>& eulerAngles)
  {
    initWithEulerAngles(eulerAngles.x, eulerAngles.y, eulerAngles.z);
  }

  /// Create a unit quaternion from a rotation matrix
  SIMD_INLINE rpQuaternion(const rpMatrix3x3<T>& matrix)
  {


	     T trace = matrix[0].x + matrix[1].y + matrix[2].z;
	     T temp[4];

	     if (trace > T(0.0))
	     {
	         T s = Sqrt(trace + T(1.0));
	         temp[3]=(s * T(0.5));
	         s = T(0.5) / s;

	         temp[0]=((matrix[2].y - matrix[1].z) * s);
	         temp[1]=((matrix[0].z - matrix[2].x) * s);
	         temp[2]=((matrix[1].x - matrix[0].y) * s);
	     }
	     else
	     {
	         int i = matrix[0].x < matrix[1].y ?
	                (matrix[1].y < matrix[2].z ? 2 : 1) :
	                (matrix[0].x < matrix[2].z ? 2 : 0);

	         int j = (i + 1) % 3;
	         int k = (i + 2) % 3;

	         T s = Sqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + T(1.0));
	         temp[i] = s * T(0.5);
	         s = T(0.5) / s;

	         temp[3] = (matrix[k][j] - matrix[j][k]) * s;
	         temp[j] = (matrix[j][i] + matrix[i][j]) * s;
	         temp[k] = (matrix[k][i] + matrix[i][k]) * s;
	     }

	     setAllValues(temp[0],temp[1],temp[2],temp[3]);

  }

  /// Create a unit quaternion from a transform matrix
  SIMD_INLINE rpQuaternion(const rpMatrix4x4<T>& matrix)
   : rpQuaternion( matrix.getRotMatrix() )
  {

  }

  //---------------------- Methods ---------------------//


  /// Set all the values
  SIMD_INLINE void setAllValues(T newX, T newY, T newZ, T newW)
  {
      v.x = newX;
      v.y = newY;
      v.z = newZ;
      w = newW;
  }

  /// Set the quaternion to zero
  SIMD_INLINE void setToZero()
  {
      v.x = 0;
      v.y = 0;
      v.z = 0;
      w = 0;
  }

  /// Set to the identity quaternion
  SIMD_INLINE void setToIdentity()
  {
      v.x = 0;
      v.y = 0;
      v.z = 0;
      w = 1;
  }

  SIMD_INLINE T R_component_1() const { return(w);   }
  SIMD_INLINE T R_component_2() const { return(v.x); }
  SIMD_INLINE T R_component_3() const { return(v.y); }
  SIMD_INLINE T R_component_4() const { return(v.z); }


  SIMD_INLINE rpVector3D<T> getV() const { return v; }
  SIMD_INLINE             T getW() const { return w; }


  /// Spinor operators : quantum mechanics
  SIMD_INLINE rpQuaternion<T> getConjugateSpinor()  const { return Quaternion(w,  v.x,  v.y,  v.z); }
  SIMD_INLINE rpQuaternion<T> getXConjugateSpinor() const { return Quaternion(w,  v.x, -v.y, -v.z); }
  SIMD_INLINE rpQuaternion<T> getYConjugateSpinor() const { return Quaternion(w, -v.x,  v.y, -v.z); }
  SIMD_INLINE rpQuaternion<T> getZConjugateSpinor() const { return Quaternion(w, -v.x, -v.y,  v.z); }


  SIMD_INLINE T getAngle() const { return T(2)*logRotor().length(); }


  /// Initialize the quaternion using Euler angles
  SIMD_INLINE void initWithEulerAngles(T angleX, T angleY, T angleZ)
  {
      T angle = angleX * T(0.5);
      const T sinX = Sin(angle);
      const T cosX = Cos(angle);

      angle = angleY * T(0.5);
      const T sinY = Sin(angle);
      const T cosY = Cos(angle);

      angle = angleZ * T(0.5);
      const T sinZ = Sin(angle);
      const T cosZ = Cos(angle);

      const T cosYcosZ = cosY * cosZ;
      const T sinYcosZ = sinY * cosZ;
      const T cosYsinZ = cosY * sinZ;
      const T sinYsinZ = sinY * sinZ;

      v.x = sinX * cosYcosZ - cosX * sinYsinZ;
      v.y = cosX * sinYcosZ + sinX * cosYsinZ;
      v.z = cosX * cosYsinZ - sinX * sinYcosZ;
      w   = cosX * cosYcosZ + sinX * sinYsinZ;

      // Normalize the quaternion
      normalize();
  }



  /// Quaternion using Euler angles
  SIMD_INLINE rpVector3D<T> getEulerAngles() const
  {
      // Store the Euler angles in radians
      rpVector3D<T> pitchYawRoll;

      rpQuaternion<T> q(*this);

      T sqw = q.w * q.w;
      T sqx = q.v.x * q.v.x;
      T sqy = q.v.y * q.v.y;
      T sqz = q.v.z * q.v.z;

      // If quaternion is normalised the unit is one, otherwise it is the correction factor
      T unit = sqx + sqy + sqz + sqw;
      T test = q.v.x * q.v.y + q.v.z * q.w;

      if (test > T(0.4999) * unit)                                // 0.4999f OR 0.5f - EPSILON
      {
          // Singularity at north pole
          pitchYawRoll.y = 2.f * Atan2(q.v.x, q.w);             // Yaw
          pitchYawRoll.x = M_PI * 0.5f;                         // Pitch
          pitchYawRoll.z = 0.f;                                 // Roll
          return pitchYawRoll;
      }
      else if (test < -T(0.4999) * unit)                          // -0.4999f OR -0.5f + EPSILON
      {
          // Singularity at south pole
          pitchYawRoll.y = -2.f * Atan2(q.v.x, q.w);            // Yaw
          pitchYawRoll.x = -M_PI * 0.5f;                        // Pitch
          pitchYawRoll.z = 0.f;                                 // Roll
          return pitchYawRoll;
      }
      else
      {
          pitchYawRoll.y = Atan2(2.f * q.v.y * q.w - 2.f * q.v.x * q.v.z,  sqx - sqy - sqz + sqw);      // Yaw
          pitchYawRoll.x = ArcSin(2.f * test / unit);                                                   // Pitch
          pitchYawRoll.z = Atan2(2.f * q.v.x * q.w - 2.f * q.v.y * q.v.z, -sqx + sqy - sqz + sqw);      // Roll
      }

      return pitchYawRoll;
  }


  /// Quaternion using camera look direction
  SIMD_INLINE rpQuaternion<T> lookRotation(rpVector3D<T>& lookAt, rpVector3D<T>& upDirection)
  {
      rpVector3D<T> forward = lookAt.getUnit();
      rpVector3D<T> right   = upDirection.getUnit().cross(forward);
      rpVector3D<T> up      = forward.cross(right);

      rpQuaternion<T> ret;
      ret.w = Sqrt(1.0f + right.x + up.y + forward.z) * 0.5f;
      T w4_recip = 1.0f / (4.0f * ret.w);
      ret.v.x = (forward.y - up.z) * w4_recip;
      ret.v.y = (right.z - forward.x) * w4_recip;
      ret.v.z = (up.x - right.y) * w4_recip;

      return ret;
  }


  //--------------------- operators -------------------------//

  /**
  * Copy operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator=(const rpQuaternion<T>& rhs)
  {
      v = rhs.v;
      w = rhs.w;
      return *this;
  }

  /**
  * Copy convert operator
  * @param rhs Right hand side argument of binary operator.
  */
  template<class FromT>
  SIMD_INLINE rpQuaternion<T>& operator=(const rpQuaternion<FromT>& rhs)
  {
      v = rhs.v;
      w = static_cast<T>(rhs.w);
      return *this;
  }


  /**
  * Addition operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T> operator+(const rpQuaternion<T>& rhs) const
  {
      const rpQuaternion<T>& lhs = *this;
      return rpQuaternion<T>(lhs.w + rhs.w, lhs.v + rhs.v);
  }


  /**
  * Subtraction operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T> operator-(const rpQuaternion<T>& rhs) const
  {
      const rpQuaternion<T>& lhs = *this;
      return rpQuaternion<T>(lhs.w - rhs.w, lhs.v - rhs.v);
  }


  /**
  * Multiplication operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T> operator*(T rhs) const
  {
      return rpQuaternion<T>(w * rhs, v * rhs);
  }


  /**
  * Multiplication invert operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T> operator/(T rhs) const
  {
      return rpQuaternion<T>(w / rhs, v / rhs);
  }



  /**
  * Multiplication operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T> operator*(const rpQuaternion<T>& rhs) const
  {
        //      const rpQuaternion<T>& lhs = *this;
        //      return rpQuaternion<T>(lhs.w * rhs.v.x + lhs.v.x * rhs.w   + lhs.v.y * rhs.v.z - lhs.v.z * rhs.v.y,
        //                             lhs.w * rhs.v.y - lhs.v.x * rhs.v.z + lhs.v.y * rhs.w   + lhs.v.z * rhs.v.x,
        //                             lhs.w * rhs.v.z + lhs.v.x * rhs.v.y - lhs.v.y * rhs.v.x + lhs.v.z * rhs.w  ,
        //                             lhs.w * rhs.w   - lhs.v.x * rhs.v.x - lhs.v.y * rhs.v.y - lhs.v.z * rhs.v.z);


       return rpQuaternion<T> (w * rhs.w - v.dot(rhs.v), w * rhs.v + rhs.w * v + v.cross(rhs.v));
  }


  SIMD_INLINE rpVector3D<T> nVidia_dot(const rpVector3D<T>& v) const
  {
          // nVidia SDK implementation
          rpVector3D<T> uv, uuv;
          rpVector3D<T> qvec(v.x, v.y, v.z);
          uv = qvec.cross(v);
          uuv = qvec.cross(uv);
          uv *= (2.0f * w);
          uuv *= 2.0f;
          return v + uv + uuv;
  }



  /**
  * Addition operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator+=(const rpQuaternion<T>& rhs)
  {
      w += rhs.w;
      v += rhs.v;
      return *this;
  }

  /**
  * Subtraction operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator-=(const rpQuaternion<T>& rhs)
  {
      w -= rhs.w;
      v -= rhs.v;
      return *this;
  }


  /**
  * Multiplication operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator*=(T rhs)
  {
      w *= rhs;
      v *= rhs;
      return *this;
  }


  /**
  * Multiplication invert operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator/=(T rhs)
  {
      w /= rhs;
      v /= rhs;
      return *this;
  }


  /**
  * Multiplication operator
  * @param rhs Right hand side argument of binary operator.
  */
  SIMD_INLINE rpQuaternion<T>& operator*=(const rpQuaternion<T>& rhs)
  {
      rpQuaternion q = (*this) * rhs;
      v = q.v;
      w = q.w;
      return *this;
  }



  /**
  * Equality test operator
  * @param rhs Right hand side argument of binary operator.
  * @note Test of equality is based of threshold EPSILON value. To be two
  * values equal, must satisfy this condition | lhs - rhs | < EPSILON,
  * for all rpQuaternion coordinates.
  */
  SIMD_INLINE bool operator==(const rpQuaternion<T>& rhs) const
  {
      const rpQuaternion<T>& lhs = *this;
      return (Abs(lhs.w - rhs.w) < MACHINE_EPSILON) && lhs.v == rhs.v;
  }

  /**
  * Inequality test operator
  * @param rhs Right hand side argument of binary operator.
  * @return not (lhs == rhs) :-P
  */
  SIMD_INLINE bool operator!=(const rpQuaternion<T>& rhs) const
  {
      return !(*this == rhs);
  }

  //---------------------- Methods ---------------------//

  SIMD_INLINE rpQuaternion<T> cross(const rpQuaternion<T>& Q) const
  {
	  return rpQuaternion<T>(0, -v.z*Q.v.y+v.y*Q.v.z,
						         v.z*Q.v.x-v.x*Q.v.z,
						        -v.y*Q.v.x+v.x*Q.v.y);
  }

  SIMD_INLINE rpQuaternion<T> dot(const rpQuaternion<T>& Q) const
  {
	  return v.x*Q.v.x+v.y*Q.v.y+v.z*Q.v.z;
  }


  SIMD_INLINE rpQuaternion<T> commutator(const rpQuaternion<T>& Q) const
  {
	  return rpQuaternion<T>(0, -2*v.z*Q.v.y+2*v.y*Q.v.z,
			                     2*v.z*Q.v.x-2*v.x*Q.v.z,
						        -2*v.y*Q.v.x+2*v.x*Q.v.y);
  }


  SIMD_INLINE rpQuaternion<T> pow(const T t) const
  {
	  return (this->log() * t).exp();
  }

  SIMD_INLINE rpQuaternion<T> pow(const rpQuaternion<T>& Q) const
  {
	  return (this->log() * Q).exp();
  }

  //-------------[ unary operations ]--------------------------

  /**
  * Unary negate operator
  * @return negated rpQuaternion
  */
  SIMD_INLINE rpQuaternion<T> operator-() const
  {
      return rpQuaternion<T>(-w, -v);
  }

  /**
  * Unary conjugate operator
  * @return conjugated rpQuaternion
  */
  SIMD_INLINE rpQuaternion<T> operator~() const
  {
      return rpQuaternion<T>(w, -v);
  }

  /**
  * Get lenght of rpQuaternion.
  * @return Length of rpQuaternion.
  */
  SIMD_INLINE T length() const
  {
      return (T) Sqrt(w * w + v.lengthSquare());
  }

  /**
  * Return square of length.
  * @return length ^ 2
  * @note This method is faster then length(). For comparison
  * of length of two rpQuaternion can be used just this value, instead
  * of more expensive length() method.
  */
  SIMD_INLINE T lengthSquare() const
  {
      return w * w + v.lengthSquare();
  }

  /**
  * Normalize rpQuaternion
  */
  SIMD_INLINE void normalize()
  {
      T len = length();
      w /= len;
      v /= len;
  }

  /**
  * Inverse rpQuaternion
  */
  SIMD_INLINE rpQuaternion<T> getInverse() const
  {
      T lengthSquareQuaternion = lengthSquare();

      assert (lengthSquareQuaternion > MACHINE_EPSILON);

      // Compute and return the inverse quaternion
      return rpQuaternion<T>(-v.x / lengthSquareQuaternion,
                             -v.y / lengthSquareQuaternion,
                             -v.z / lengthSquareQuaternion,
                              w   / lengthSquareQuaternion);
  }


  /**
  * Unit rpQuaternion
  */
  SIMD_INLINE rpQuaternion<T> getUnit() const
  {
    T lengthQuaternion = length();

    // Check if the length is not equal to zero
    assert (lengthQuaternion > MACHINE_EPSILON);

    // Compute and return the unit quaternion
    return rpQuaternion<T>(v.x / lengthQuaternion,
                           v.y / lengthQuaternion,
                           v.z / lengthQuaternion,
                           w / lengthQuaternion);
  }


  /**
  * Conjugate rpQuaternion
  */
  SIMD_INLINE rpQuaternion<T> getConjugate() const
  {
    return rpQuaternion<T>(-v, w);
  }


  /**
  * Creates rpQuaternion for eulers angles.
  * @param x Rotation around x axis (in degrees).
  * @param y Rotation around y axis (in degrees).
  * @param z Rotation around z axis (in degrees).
  * @return rpQuaternion object representing transformation.
  */
  static SIMD_INLINE rpQuaternion<T> fromEulerAngles(T x, T y, T z)
  {
      rpQuaternion<T> ret = fromAxisRot(rpVector3D<T>(1, 0, 0), x) *
                            fromAxisRot(rpVector3D<T>(0, 1, 0), y) *
                            fromAxisRot(rpVector3D<T>(0, 0, 1), z);
      return ret;
  }

  /**
  * Creates rpQuaternion as rotation around axis.
  * @param axis Unit vector expressing axis of rotation.
  * @param angleDeg Angle of rotation around axis (in degrees).
  */
  static SIMD_INLINE rpQuaternion<T> fromAxisRot(const rpVector3D<T> axis, T angleDeg)
  {
      T angleRad = DEG2RAD(angleDeg);
      T sa2 = Sin(angleRad / 2);
      T ca2 = Cos(angleRad / 2);
      return rpQuaternion<T>( ca2 , sa2 * axis );
  }

  /**
  * Converts rpQuaternion into rotation matrix.
  * @return Rotation matrix expressing this rpQuaternion.
  */
  SIMD_INLINE rpMatrix3x3<T> getMatrix() const
  {
      rpMatrix3x3<T> ret;


      T nQ = v.x*v.x + v.y*v.y + v.z*v.z + w*w;
      T s = 0.0;

      if(nQ > 0.0)
      {
          s = T(2.0) / nQ; /// normalize quaternion
      }

      // Computations used for optimization (less multiplications)
      T xs  = v.x*s;
      T ys  = v.y*s;
      T zs  = v.z*s;
      T wxs = w*xs;
      T wys = w*ys;
      T wzs = w*zs;
      T xxs = v.x*xs;
      T xys = v.x*ys;
      T xzs = v.x*zs;
      T yys = v.y*ys;
      T yzs = v.y*zs;
      T zzs = v.z*zs;

      // Create the matrix corresponding to the quaternion
      ret = rpMatrix3x3<T>( (1.0) - yys - zzs, xys-wzs, xzs + wys ,
                            xys + wzs, T(1.0) - xxs - zzs, yzs-wxs,
                            xzs-wys, yzs + wxs, T(1.0) - xxs - yys );


      return ret;
  }

  /**
  * Converts rpQuaternion into transformation matrix.
  * @note This method performs same operation as rotMatrix()
  * conversion method. But returns Matrix of 4x4 elements.
  * @return Transformation matrix expressing this rpQuaternion.
  */
  SIMD_INLINE rpMatrix4x4<T> transform() const
  {
      rpMatrix4x4<T> ret;

      T nQ = v.x*v.x + v.y*v.y + v.z*v.z + w*w;
      T s = 0.0;

      if(nQ > 0.0)
      {
          s = T(2.0) / nQ; /// normalize quaternion
      }

      // Computations used for optimization (less multiplications)
      T xs  = v.x*s;
      T ys  = v.y*s;
      T zs  = v.z*s;
      T wxs = w*xs;
      T wys = w*ys;
      T wzs = w*zs;
      T xxs = v.x*xs;
      T xys = v.x*ys;
      T xzs = v.x*zs;
      T yys = v.y*ys;
      T yzs = v.y*zs;
      T zzs = v.z*zs;

      // Create the matrix corresponding to the quaternion
      ret = rpMatrix4x4<T>( T(1.0) - yys - zzs, xys-wzs, xzs + wys, T(0.0),
                            xys + wzs, T(1.0) - xxs - zzs, yzs-wxs, T(0.0),
                            xzs-wys, yzs + wxs, T(1.0) - xxs - yys, T(0.0),
                            T(0.0),T(0.0),T(0.0),T(1.0));

      return ret;

  }


  /**
  * Linear interpolation of two rpQuaternions
  * @param fact Factor of interpolation. For translation from position
  * of this vector to rpQuaternion rhs, values of factor goes from 0.0 to 1.0.
  * @param rhs Second rpQuaternion for interpolation
  * @note However values of fact parameter are reasonable only in interval
  * [0.0 , 1.0], you can pass also values outside of this interval and you
  * can get result (extrapolation?)
  */
  SIMD_INLINE rpQuaternion<T> lerp(T fact, const rpQuaternion<T>& rhs) const
  {
      return rpQuaternion<T>((1 - fact) * w + fact * rhs.w, v.lerp(fact, rhs.v));
  }


  /**
  * Computes spherical interpolation between rpQuaternions (this, q2)
  * using coefficient of interpolation r (in [0, 1]).
  *
  * @param r The ratio of interpolation form this (r = 0) to q2 (r = 1).
  * @param q2 Second rpQuaternion for interpolation.
  * @return Result of interpolation.
  */
  SIMD_INLINE rpQuaternion<T> slerp(T r, const rpQuaternion<T>& q2) const
  {
      rpQuaternion<T> ret;
      T cosTheta = w * q2.w + v.x * q2.v.x + v.y * q2.v.y + v.z * q2.v.z;
      T theta = (T) ArcCos(cosTheta);
      if (Abs(theta) < epsilon)
      {
          ret = *this;
      }
      else
      {
          T sinTheta = (T) sqrt(1.0 - cosTheta * cosTheta);
          if (Abs(sinTheta) < epsilon)
          {
              ret.w = 0.5 * w + 0.5 * q2.w;
              ret.v = v.lerp(0.5, q2.v);
          }
          else
          {
              T rA = (T) Sin((1.0 - r) * theta) / sinTheta;
              T rB = (T) Sin(r * theta) / sinTheta;

              ret.w = w * rA + q2.w * rB;
              ret.v.x = v.x * rA + q2.v.x * rB;
              ret.v.y = v.y * rA + q2.v.y * rB;
              ret.v.z = v.z * rA + q2.v.z * rB;
          }
      }
      return ret;
  }



  /// Return logarithm of Quaternion.
  SIMD_INLINE rpQuaternion<T> log() const
  {
	rpQuaternion<T>  Result;
    const T b = Sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if(Abs(b) <= MACHINE_EPSILON*Abs(w))
    {
      if(w<0.0)
      {
        //INFOTOCERR << " Error: Infinitely many solutions for log of a negative scalar (w=" << w << ")." << endl;
        //throw(InfinitelyManySolutions);
      }
      Result.w = std::log(w);
    }
    else
    {
      const T _v = Atan2(b, w);
      const T _f = _v/b;
      Result.w = Log(w*w+b*b)/2.0;
      // Result.w = std::log(w/std::cos(v)); // Not nice for unit vectors [w=cos(v)=0]
      Result.v.x = _f*v.x;
      Result.v.y = _f*v.y;
      Result.v.z = _f*v.z;
    }
    return Result;
  }

  /// Return logarithm of a rotor.
  SIMD_INLINE rpQuaternion<T> logRotor() const
  {
	  /// This function is just like the standard log function, except
	  /// that a negative scalar does not raise an exception; instead, the
	  /// scalar pi is returned.
	  rpQuaternion<T> Result;
	  const T b = Sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

	  if(Abs(b) <= MACHINE_EPSILON*Abs(w))
	  {
		  if(w<0.0)
		  {
			  //        INFOTOCERR
			  //             << "\nWarning: Infinitely many solutions for log of a negative scalar (w=" << w << "); "
			  //             << "arbitrarily returning the one in the x direction." << endl;
			  Result.v.x = M_PI;
			  if(Abs(w+1)>MACHINE_EPSILON)
			  {
				  Result.w = Log(-w);
			  }
		  }
		  else
		  {
			  Result.w = Log(w);
		  }
	  }
	  else
	  {
		  const T _v = Atan2(b, w);
		  const T _f = _v/b;
		  Result.w = Log(w*w+b*b)/2.0;
		  // Result.w = std::log(w/std::cos(v)); // Not nice for unit vectors [w=cos(v)=0]
		  Result.v.x = _f*v.x;
		  Result.v.y = _f*v.y;
		  Result.v.z = _f*v.z;
	  }

	  return Result;
  }

  /// Return exponent of Quaternion.
  SIMD_INLINE rpQuaternion<T>  exp() const
  {
	rpQuaternion<T>  Result;
    const T b = Sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if(Abs(b)<=MACHINE_EPSILON*Abs(w))
    {
      Result.w = Exp(w);
    }
    else
    {
      const T e = Exp(w);
      const T f = Sin(b)/b; // Note: b is never 0.0 at this point
      Result.w = e*Cos(b);
      Result.v.x = e*f*v.x;
      Result.v.y = e*f*v.y;
      Result.v.z = e*f*v.z;
    }
    return Result;
  }

  //----------[ output operator ]----------------------------
  /**
  * Provides output to standard output stream.
  */
  friend std::ostream& operator <<(std::ostream& oss, const rpQuaternion<T>& q)
  {
	  oss << "(" << "Re: " << q.w << " Im: " << q.v << ")";
      return oss;
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


  //---------------------------------- Plugins ----------------------------------------//

  /// The angle between the vectors is simple to find: the dot product gives its cosine.
  /// The needed axis is also simple to find: it’s the cross product of the two vectors.
  static SIMD_INLINE rpQuaternion<T> rotationBetweenSlerp(  rpVector3D<T> start ,  rpVector3D<T> dest )
  {

      //start.normalize();
      //dest.normalize();

      T cosTheta = start.dot(dest);
      rpVector3D<T> rotationAxis;

      if (cosTheta < -1 + 0.0001f)
      {
          // special case when vectors in opposite directions:
          // there is no "ideal" rotation axis
          // So guess one; any will do as long as it's perpendicular to start
          rotationAxis = rpVector3D<T>::Z.cross(start);
          if (rotationAxis.length2() < T(0.0001) )
          {
              // bad luck, they were parallel, try again!
              rotationAxis = rpVector3D<T>::X.cross(start);
          }

          rotationAxis.normalize();
          static rpQuaternion<T> q;
          return q.fromAxisRot(rotationAxis,T(180.0f));
      }

      rotationAxis = start.cross(dest);
      //rotationAxis = rotationAxis.normalize();



      T s = Sqrt( (1+cosTheta)*2 );
      T invs = 1.0 / s;

      return rpQuaternion<T>
      (
        rotationAxis.x * invs,
        rotationAxis.y * invs,
        rotationAxis.z * invs,
        s * 0.5f
      );

  }



public:
    /**
     * The multiplicitive identity quaternion
     */
    static const rpQuaternion<T> IDENTITY;
    /**
     * The additive identity quaternion.
     */
    static const rpQuaternion<T> ZERO;

};

template<class T> const rpQuaternion<T> rpQuaternion<T>::IDENTITY(0.0, 0.0, 0.0, 1.0);
template<class T> const rpQuaternion<T> rpQuaternion<T>::ZERO(0.0, 0.0, 0.0, 0.0);


/// Octonion of floats
typedef rpQuaternion<float> rpQuaternionf;
/// Octonionof doubles
typedef rpQuaternion<double> rpQuaterniond;
/// Octonion of int
typedef rpQuaternion<double> rpQuaternioni;

} /* namespace  */



#endif /* SOURCE_REAL_PHYSICS_LINEARMATHS_RPQUATERNION_H_ */
