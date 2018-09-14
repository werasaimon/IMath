#ifndef ISCALARTYPE_H
#define ISCALARTYPE_H



#include "IVector.h"
#include "IVector2D.h"
#include "IVector3D.h"
#include "IVector4D.h"
#include "IMatrix3x3.h"
#include "IMatrix4x4.h"
#include "IQuaternion.h"


namespace IMath
{


//! Provides the scalar type of the scalar, vector, or matrix type specified by 'T'.
template <typename T>
struct ScalarType
{
    using Type = T;
};

template <typename T, std::size_t N>
struct ScalarType<IVector<T, N>>
{
    using Type = T;
};

//  template <typename T, std::size_t Rows, std::size_t Cols>
//  struct ScalarType<Matrix<T, Rows, Cols>>
//  {
//      using Type = T;
//  };

template <typename T>
struct ScalarType<IMatrix3x3<T>>
{
    using Type = T;
};

template <typename T>
struct ScalarType<IMatrix4x4<T>>
{
    using Type = T;
};

template <typename T>
struct ScalarType<IQuaternion<T>>
{
    using Type = T;
};

}

#endif // ISCALARTYPE_H
