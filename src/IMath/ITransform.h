 /********************************************************************************
 *
 * ITransform.h
 *
 * IMath : 3d_math library,
 * Copyright (c)  *
 * Created on: 3 July. 2018 г.
 * Author: werasaimon                                     *
 *********************************************************************************
 *                                                                               *
 * This software is provided 'as-is', without any express or implied warranty.   *
 * In no event will the authors be held liable for any damages arising from the  *
 * use of this software.                                                         *
 *                                                                               *
 * Permission is granted to anyone to use this software for any purpose,         *
 * including commercial applications, and to alter it and redistribute it        *
 * freely, subject to the following restrictions:                                *
 *                                                                               *
 * 1. The origin of this software must not be misrepresented; you must not claim *
 *    that you wrote the original software. If you use this software in a        *
 *    product, an acknowledgment in the product documentation would be           *
 *    appreciated but is not required.                                           *
 *                                                                               *
 * 2. Altered source versions must be plainly marked as such, and must not be    *
 *    misrepresented as being the original software.                             *
 *                                                                               *
 * 3. This notice may not be removed or altered from any source distribution.    *
 *                                                                               *
 ********************************************************************************/

#ifndef ITRANSFORM_H_
#define ITRANSFORM_H_


// Libraries
#include "IMatrix3x3.h"
#include "IVector3D.h"
#include "IQuaternion.h"


namespace IMath
{

namespace
{
#define SIMD_2_PI (6.283185307179586232)
#define SIMD_PI (SIMD_2_PI * (0.5))
#define SIMD_HALF_PI (SIMD_2_PI * (0.25))
#define ANGULAR_MOTION_THRESHOLD (0.5)*SIMD_HALF_PI
}

template<class T> class  ITransform
{


public:


  //---------------------------------- Attributes ---------------------------------//

  T                    mTime;
  IVector3D<T>         mPosition;
  IMatrix3x3<T>        mBasis;



  //--------------------------------  Constructor ---------------------------------//

  // Constructor
  SIMD_INLINE ITransform()
   : mPosition(IVector3D<T>::ZERO) , mBasis(IMatrix3x3<T>::IDENTITY) , mTime(1.0)
  {

  }
  // Constructor
  SIMD_INLINE ITransform(const IVector3D<T>& position , const IMatrix3x3<T>& _basis = IMatrix3x3<T>::IDENTITY ,T  _time = 1.0)
  : mPosition(position) , mBasis(_basis) , mTime(_time)
  {

  }


  // Constructor
  SIMD_INLINE ITransform(const IVector3D<T>& position , const IQuaternion<T>& _quaternion = IQuaternion<T>::IDENTITY ,T  _time = 1.0)
  : mPosition(position) , mTime(_time)
  {
      mBasis = IMatrix3x3<T>::createRotation(_quaternion);
  }


  // Copy-constructor
  SIMD_INLINE ITransform(const ITransform<T>& transform)
  : mPosition(transform.mPosition) , mBasis(transform.mBasis) ,  mTime(transform.mTime)
  {

  }

  //-------------------------------   Methods ---------------------------------//

  /// Return the identity transform
  static SIMD_INLINE ITransform<T> identity()
  {
      return ITransform<T>();
  }


  ///**@brief Return a quaternion representing the rotation */
   IQuaternion<T> getRotation() const
   {
       IQuaternion<T> q(mBasis);
       return q;
   }


  SIMD_INLINE T getTime() const
  {
     return mTime;
  }

  SIMD_INLINE void setTime(const T &time)
  {
     mTime = time;
  }


  SIMD_INLINE IVector3D<T> getPosition() const
  {
      return mPosition;
  }

  SIMD_INLINE void setPosition(const IVector3D<T> &position)
  {
      mPosition = position;
  }


  SIMD_INLINE IMatrix3x3<T> getBasis() const
  {
      return mBasis;
  }

  SIMD_INLINE void setBasis(const IMatrix3x3<T> &basis)
  {
      mBasis = basis;
  }



  /**@brief Return the inverse of this transform */
  SIMD_INLINE ITransform<T> getInverse() const
  {
      IMatrix3x3<T> inv = mBasis.getInverse();
      return ITransform<T>(inv * -mPosition , inv , T(1.0) / mTime);
  }


  SIMD_INLINE IVector3D<T> invXform(const IVector3D<T>& inVec) const
  {
      IVector3D<T> v = inVec - mPosition;
      return (mBasis.getTranspose() * v);
  }


  /**@brief Return the transform of the vector */
  SIMD_INLINE IVector3D<T> operator()(const IVector3D<T>& x) const
  {
      return mBasis * x + mPosition;
  }

  /**@brief Return the transform of the vector */
  SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& x) const
  {
      return (*this)(x);
  }

  //**@brief Return the transform of the btQuaternion */
  SIMD_INLINE IQuaternion<T> operator*(const IQuaternion<T>& q) const
  {
      return getRotation() * q;
  }


  SIMD_INLINE ITransform<T> operator*(const ITransform<T>& t) const
  {
      return ITransform<T>( (*this)(t.mPosition) , mBasis * t.mBasis , mTime * t.mTime);
  }


  SIMD_INLINE bool operator == (const ITransform<T>& transform2) const
  {
    return  (mTime     == transform2.mTime)     &&
            (mPosition == transform2.mPosition) &&
            (mBasis    == transform2.mBasis);
  }


  SIMD_INLINE bool operator !=(const ITransform<T>& transform2) const
  {
    return !(*this == transform2);
  }

  //--------------------------------- function OpenGL -----------------------------------//


  /// Set the transform from an OpenGL transform matrix
  SIMD_INLINE void setFromOpenGL(T* openglMatrix)
  {
      IMatrix3x3<T> matrix(openglMatrix[0], openglMatrix[4], openglMatrix[8],
                            openglMatrix[1], openglMatrix[5], openglMatrix[9],
                            openglMatrix[2], openglMatrix[6], openglMatrix[10]);

      mBasis = (matrix.getTranspose());

      IVector3D<T> pos( openglMatrix[12],
                         openglMatrix[13],
                         openglMatrix[14]);

      mPosition = pos;


  }


  /// Get the OpenGL matrix of the transform
  SIMD_INLINE void getOpenGLMatrix(T* openglMatrix) const
  {

      const IMatrix3x3<T>& matrix = getBasis();

      openglMatrix[0]  = matrix[0][0];
      openglMatrix[1]  = matrix[1][0];
      openglMatrix[2]  = matrix[2][0];
      openglMatrix[3]  = 0.0;

      openglMatrix[4]  = matrix[0][1];
      openglMatrix[5]  = matrix[1][1];
      openglMatrix[6]  = matrix[2][1];
      openglMatrix[7]  = 0.0;

      openglMatrix[8]  = matrix[0][2];
      openglMatrix[9]  = matrix[1][2];
      openglMatrix[10] = matrix[2][2];
      openglMatrix[11] = 0.0;

      openglMatrix[12] = mPosition.x;
      openglMatrix[13] = mPosition.y;
      openglMatrix[14] = mPosition.z;
      openglMatrix[15] = 1.0;

  }



  ///========================== Plugin ==============================///

  static ITransform<T> integrateTransform( const ITransform<T>& curTrans , const IVector3D<T>& linvel, const IVector3D<T>& angvel, T timeStep)
  {
	  //Exponential map
	  //google for "Practical Parameterization of Rotations Using the Exponential Map", F. Sebastian Grassia
      IVector3D<T> axis;
	  T fAngle = angvel.length();

	  if (fAngle < T(0.001))
	  {
		  // use Taylor's expansions of sync function
		  axis = angvel * (T(0.5) * timeStep -
				  (timeStep * timeStep * timeStep) *
				  (T(0.020833333333)) * fAngle * fAngle);
	  }
	  else
	  {
		  axis = angvel * (Sin(T(0.5) * fAngle * timeStep) / fAngle);
	  }



      IQuaternion<T> dorn(axis, ICos(fAngle * timeStep * T(0.5)));
      IQuaternion<T> predictedOrn = dorn * IQuaternion<T>(curTrans.getBasis());
	  predictedOrn.normalize();


          return ITransform<T>( curTrans.getPosition() + linvel * timeStep ,  predictedOrn.getMatrix() , curTrans.getTime() );

  }


  static IQuaternion<T> integrateAngular( const IQuaternion<T>& curQuaternion , const IVector3D<T>& angvel, T timeStep )
  {
      //Exponential map
      //google for "Practical Parameterization of Rotations Using the Exponential Map", F. Sebastian Grassia
      IVector3D<T> axis;
      T fAngle = angvel.length();


      if (fAngle < (0.001))
      {
          // use Taylor's expansions of sync function
          axis = angvel * ((0.5) * timeStep -
                          (timeStep * timeStep * timeStep) *
                          ((0.020833333333)) * fAngle * fAngle);
      }
      else
      {
          axis = angvel * (ISin( (0.5) * fAngle * timeStep) / fAngle);
      }



      IQuaternion<T> dorn(axis, ICos(fAngle * timeStep * (0.5)));
      IQuaternion<T> predictedOrn = dorn * curQuaternion;
      predictedOrn.normalize();

      return predictedOrn;

  }


  static IVector3D<T> integrateLinear( const IVector3D<T>& curPosition , const IVector3D<T>& linvel, T timeStep )
  {
      return curPosition + linvel * timeStep;
  }

};





}



#endif /* ITRANSFORM_H_ */
