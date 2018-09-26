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
    : mPosition(IVector3D<T>::ZERO) ,
      mBasis(IMatrix3x3<T>::IDENTITY) ,
      mTime(1.0)
    {

    }
    // Constructor
    SIMD_INLINE ITransform(const IVector3D<T>& position , const IMatrix3x3<T>& _basis = IMatrix3x3<T>::IDENTITY ,T  _time = 1.0)
    : mPosition(position) ,
      mBasis(_basis) ,
      mTime(_time)
    {

    }


    // Constructor
    SIMD_INLINE ITransform(const IVector3D<T>& position , const IQuaternion<T>& _quaternion = IQuaternion<T>::IDENTITY ,T  _time = 1.0)
    : mPosition(position) ,
      mTime(_time),
      mBasis(_quaternion.GetRotMatrix())
    {

    }


    // Copy-constructor
    SIMD_INLINE ITransform(const ITransform<T>& transform)
    : mPosition(transform.mPosition) ,
      mBasis(transform.mBasis) ,
      mTime(transform.mTime)
    {

    }

    //-------------------------------   Methods ---------------------------------//

    /// Return the identity transform
    static SIMD_INLINE ITransform<T> Identity()
    {
        return ITransform<T>();
    }


    ///**@brief Return a quaternion representing the rotation */
    IQuaternion<T> GetRotation() const
    {
        IQuaternion<T> q(mBasis);
        return q;
    }


    SIMD_INLINE T GetTime() const
    {
        return mTime;
    }

    SIMD_INLINE void SetTime(const T &time)
    {
        mTime = time;
    }


    SIMD_INLINE IVector3D<T> GetPosition() const
    {
        return mPosition;
    }

    SIMD_INLINE void SetPosition(const IVector3D<T> &position)
    {
        mPosition = position;
    }


    SIMD_INLINE IMatrix3x3<T> GetBasis() const
    {
        return mBasis;
    }

    SIMD_INLINE void SetBasis(const IMatrix3x3<T> &basis)
    {
        mBasis = basis;
    }



    /**@brief Return the inverse of this transform */
    SIMD_INLINE ITransform<T> GetInverse() const
    {
        IMatrix3x3<T> inv = mBasis.GetInverse();
        return ITransform<T>(inv * -mPosition , inv , T(1.0) / mTime);
    }


    SIMD_INLINE IVector3D<T> InvXform(const IVector3D<T>& inVec) const
    {
        IVector3D<T> v = inVec - mPosition;
        return (mBasis.GetTranspose() * v);
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
        return GetRotation() * q;
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

        mBasis = (matrix.GetTranspose());

        IVector3D<T> pos( openglMatrix[12],
                          openglMatrix[13],
                          openglMatrix[14]);

        mPosition = pos;


    }


    /// Get the OpenGL matrix of the transform
    SIMD_INLINE void GetOpenGLMatrix(T* openglMatrix) const
    {
        const IMatrix3x3<T>& matrix = GetBasis();

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

};





}



#endif /* ITRANSFORM_H_ */
