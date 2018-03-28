/*
 * rpTransformComposite.h
 *
 *  Created on: 3 февр. 2018 г.
 *      Author: werqa
 */

#ifndef SRC_MATHS_RPTRANSFORM_H_
#define SRC_MATHS_RPTRANSFORM_H_

// Libraries
#include "rpMatrix3x3.hpp"
#include "rpVector3D.hpp"
#include "rpQuaternion.hpp"


namespace CMath
{


namespace
{
#define SIMD_2_PI (6.283185307179586232)
#define SIMD_PI (SIMD_2_PI * (0.5))
#define SIMD_HALF_PI (SIMD_2_PI * (0.25))
#define ANGULAR_MOTION_THRESHOLD (0.5)*SIMD_HALF_PI
}

template<class T> class  rpTransform
{


public:


  //---------------------------------- Attributes ---------------------------------//

  T                     mTime;
  rpVector3D<T>         mPosition;
  rpMatrix3x3<T>        mBasis;



  //--------------------------------  Constructor ---------------------------------//

  // Constructor
  SIMD_INLINE rpTransform()
   : mPosition(rpVector3D<T>::ZERO) , mBasis(rpMatrix3x3<T>::IDENTITY) , mTime(1.0)
  {

  }
  // Constructor
  SIMD_INLINE rpTransform(const rpVector3D<T>& position , const rpMatrix3x3<T>& _basis = rpMatrix3x3<T>::IDENTITY ,T  _time = 1.0)
  : mPosition(position) , mBasis(_basis) , mTime(_time)
  {

  }


  // Copy-constructor
  SIMD_INLINE rpTransform(const rpTransform<T>& transform)
  : mPosition(transform.mPosition) , mBasis(transform.mBasis) ,  mTime(transform.mTime)
  {

  }

  //-------------------------------   Methods ---------------------------------//

  /// Return the identity transform
  static SIMD_INLINE rpTransform<T> identity()
  {
      return rpTransform<T>();
  }



  SIMD_INLINE T getTime() const
  {
     return mTime;
  }

  SIMD_INLINE void setTime(const T &time)
  {
     mTime = time;
  }


  SIMD_INLINE rpVector3D<T> getPosition() const
  {
      return mPosition;
  }

  SIMD_INLINE void setPosition(const rpVector3D<T> &position)
  {
      mPosition = position;
  }


  SIMD_INLINE rpMatrix3x3<T> getBasis() const
  {
      return mBasis;
  }

  SIMD_INLINE void setBasis(const rpMatrix3x3<T> &basis)
  {
      mBasis = basis;
  }


  /// Return the inverse of the transform
  SIMD_INLINE rpTransform<T> getInverse() const
  {
      const rpMatrix3x3<T>  &invMatrix     = getBasis().getInverse();
      return rpTransform<T>( invMatrix * (-mPosition)  , invMatrix , T(1.0/ mTime));
  }


  SIMD_INLINE rpVector3D<T> operator * (const rpVector3D<T>& vector) const
  {
    return  (getBasis() * vector) + mPosition;
  }

  SIMD_INLINE rpTransform<T> operator * (const rpTransform<T>& transform2) const
  {
    return rpTransform<T>(mPosition + (getBasis() * transform2.mPosition) ,  /*transform2.mBasis * mBasis*/ mBasis * transform2.mBasis , mTime * transform2.mTime);
  }


  SIMD_INLINE bool operator == (const rpTransform<T>& transform2) const
  {
    return  (mTime     == transform2.mTime)     &&
            (mPosition == transform2.mPosition) &&
            (mBasis    == transform2.mBasis);
  }


  SIMD_INLINE bool operator !=(const rpTransform<T>& transform2) const
  {
    return !(*this == transform2);
  }

  //--------------------------------- function OpenGL -----------------------------------//


  /// Set the transform from an OpenGL transform matrix
  SIMD_INLINE void setFromOpenGL(T* openglMatrix)
  {
      rpMatrix3x3<T> matrix(openglMatrix[0], openglMatrix[4], openglMatrix[8],
                            openglMatrix[1], openglMatrix[5], openglMatrix[9],
                            openglMatrix[2], openglMatrix[6], openglMatrix[10]);

      mBasis = (matrix.getTranspose());

      rpVector3D<T> pos( openglMatrix[12],
                         openglMatrix[13],
                         openglMatrix[14]);

      mPosition = pos;


  }


  /// Get the OpenGL matrix of the transform
  SIMD_INLINE void getOpenGLMatrix(T* openglMatrix) const
  {

      const rpMatrix3x3<T>& matrix = getBasis();

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

  static rpTransform<T> integrateTransform( const rpTransform<T>& curTrans , const rpVector3D<T>& linvel, const rpVector3D<T>& angvel, T timeStep)
  {
	  //Exponential map
	  //google for "Practical Parameterization of Rotations Using the Exponential Map", F. Sebastian Grassia
	  rpVector3D<T> axis;
	  T fAngle = angvel.length();

	  //limit the angular motion
	  if (fAngle*timeStep > ANGULAR_MOTION_THRESHOLD)
	  {
		  fAngle = ANGULAR_MOTION_THRESHOLD / timeStep;
	  }

	  if (fAngle < T(0.001))
	  {
		  // use Taylor's expansions of sync function
		  axis = angvel * (T(0.5) * timeStep -
				  (timeStep * timeStep * timeStep) *
				  (T(0.020833333333)) * fAngle * fAngle);
	  }
	  else
	  {
		  // sync(fAngle) = sin(c*fAngle)/t
		  axis = angvel * (Sin(T(0.5) * fAngle * timeStep) / fAngle);
	  }



	  rpQuaternion<T> dorn(axis, Cos(fAngle * timeStep * T(0.5)));
	  rpQuaternion<T> predictedOrn = dorn * rpQuaternion<T>(curTrans.getBasis());
	  predictedOrn.normalize();


	  return rpTransform<T>(curTrans.getPosition() + linvel * timeStep , predictedOrn , curTrans.getTime() );

  }




};





}



#endif /* SRC_MATHS_RPTRANSFORM_H_ */
