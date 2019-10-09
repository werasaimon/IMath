#ifndef IVECTORTYPE_H
#define IVECTORTYPE_H

#include "IVector2D.h"
#include "IVector3D.h"
#include "IVector4D.h"

#include <vector>


namespace IMath
{

template <typename T, std::size_t N>
struct VectorType
{
    using Vector = std::vector<T>;
};

template <typename T>
struct VectorType<T, 2u>
{
    using Vector = IVector2D<T>;
};

template <typename T>
struct VectorType<T, 3u>
{
    using Vector = IVector3D<T>;
};

template <typename T>
struct VectorType<T, 4u>
{
    using Vector = IVector4D<T>;
};

}

#endif // IVECTORTYPE_H
