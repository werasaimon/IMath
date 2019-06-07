#ifndef IVECTOR_H
#define IVECTOR_H

#include "IFunc.h"

#include <algorithm>
#include <iterator>
#include <cstdint>

namespace IMath
{


/**
  \brief Base vector class with N components.
  \tparam T Specifies the data type of the vector components.
  This should be a primitive data type such as float, double, int etc.
  \tparam N Specifies the number of components. There are specialized templates for N = 2, 3, and 4.
  */
template <typename T, std::size_t N>
class IVector
{

private:

    T v_[N];


public:

    //! Specifies the typename of the scalar components.
    using ScalarType = T;

    //! Specifies the number of vector components.
    static const std::size_t components = N;

#ifndef I_DISABLE_AUTO_INIT
    IVector()
    {
        std::fill(std::begin(v_), std::end(v_), T(0));
    }
#else
    IVector() = default;
#endif

    IVector(const IVector<T, N>& rhs)
    {
        std::copy(std::begin(rhs.v_), std::end(rhs.v_), v_);
    }

    IVector(const int& deciminal)
    {
        std::fill(std::begin(v_), std::end(v_), deciminal);
    }

    explicit IVector(const T& scalar)
    {
        std::fill(std::begin(v_), std::end(v_), scalar);
    }


    IVector<T, N>& operator += (const IVector<T, N>& rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] += rhs[i];
        return *this;
    }

    IVector<T, N>& operator -= (const IVector<T, N>& rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] -= rhs[i];
        return *this;
    }

    IVector<T, N>& operator *= (const IVector<T, N>& rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] *= rhs[i];
        return *this;
    }

    IVector<T, N>& operator /= (const IVector<T, N>& rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] /= rhs[i];
        return *this;
    }

    IVector<T, N>& operator *= (const T rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] *= rhs;
        return *this;
    }

    IVector<T, N>& operator /= (const T rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            v_[i] /= rhs;
        return *this;
    }


    IVector<T, N>& Dot(const IVector<T, N>& rhs) const
    {
        T result = T(0);
        for (std::size_t i = 0; i < N; ++i)
        {
            result += (*this)[i]*rhs[i];
        }

        return result;
    }


    /**
     \brief Returns the specified vector component.
     \param[in] component Specifies the vector component index. This must be in the range [0, N).
    */
    T& operator [] (std::size_t component)
    {
        assert(component < N);
        return v_[component];
    }

    /**
     \brief Returns the specified vector component.
     \param[in] component Specifies the vector component index. This must be in the range [0, N).
     */
    const T& operator [] (std::size_t component) const
    {
        assert(component < N);
        return v_[component];
    }



    /**
     Returns a type casted instance of this vector.
     \tparam C Specifies the static cast type.
    */
    template <typename C>
    IVector<C, N> Cast() const
    {
        IVector<C, N> result;

        for (std::size_t i = 0; i < N; ++i)
            result[i] = static_cast<C>(v_[i]);

        return result;
    }

    //! Returns a pointer to the first element of this vector.
    T* Ptr()
    {
        return v_;
    }

    //! Returns a constant pointer to the first element of this vector.
    const T* Ptr() const
    {
        return v_;
    }


    //-------------[ conversion ]-----------------------------
    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class IVector2D<T>)
     * used for passing IVector<T,N> values to gl*n[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return &v_[0];
    }
    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class IVector2D<T>)
     * used for passing IVector<T,N> values to gl*n[fd] functions.
     */
    SIMD_INLINE operator const T*() const
    {
        return &v_[0];
    }


    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE  IVector<T, N> operator - () const
    {
        T result = T(0);
        for (std::size_t i = 0; i < N; ++i)
        {
            result += (*this)[i]*(*this)[i];
        }
        return result;
    }

    //-------------[ size operations ]---------------------------

    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length()^ 2 method.
     */
    SIMD_INLINE T LengthSquare() const
    {
        auto result = T(0);
        for (std::size_t i = 0; i < N; ++i)
        {
            result = result + (*this)[i] * (*this)[i];
        }
        return result;
    }

    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T Length() const
    {
        return ISqrt(LengthSquare());
    }


    /**
     * Normalize vector
     */
    SIMD_INLINE void Normalize()
    {
        const T s = Length();
        auto result = *this;
        for (std::size_t i = 0; i < N; ++i)
        {
            (*this)[i] = result[i] / s;
        }
    }



    /**
     * Inverse vector
     */
    SIMD_INLINE IVector<T, N> GetInverse() const
    {
        T s = Length();
        auto result = *this;
        for (std::size_t i = 0; i < N; ++i)
        {
            result[i] = T(1.0)/result[i];
        }
        return result;
    }


    /**
     * Normalize unit vector
     */
    SIMD_INLINE IVector<T, N> GetUnit() const
    {
        auto  v = (*this);
        v.Normalize();
        return v;
    }

    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE IVector<T, N> Normalized() const
    {
        auto  v = (*this);
        v.Normalize();
        return v;
    }


    //-------------[ output operator ]------------------------
    /**
     * Output to stream operator
     * @param lhs Left hand side argument of operator (commonly ostream instance).
     * @param rhs Right hand side argument of operator.
     * @return Left hand side argument - the ostream object passed to operator.
     */
    friend std::ostream& operator<<(std::ostream& lhs, const IVector<T, N>& rhs)
    {
        lhs << "[";
        for (std::size_t i = 0; i < N; ++i)
        {
            lhs << rhs[i];
            if(i < N-1)
                lhs << ",";
        }

        lhs << "]";
        return lhs;
    }

    /**
           * Gets string representation.
           */
    std::string ToString() const
    {
        std::ostringstream oss;
        oss << *this;
        return oss.str();
    }


};


/* --- Global Operators --- */

template <typename T, std::size_t N>
IVector<T, N> operator + (const IVector<T, N>& lhs, const IVector<T, N>& rhs)
{
    auto result = lhs;
    result += rhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator - (const IVector<T, N>& lhs, const IVector<T, N>& rhs)
{
    auto result = lhs;
    result -= rhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator * (const IVector<T, N>& lhs, const IVector<T, N>& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator / (const IVector<T, N>& lhs, const IVector<T, N>& rhs)
{
    auto result = lhs;
    result /= rhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator * (const IVector<T, N>& lhs, const T& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

//! \note This implementation is equivavlent to (rhs * lhs) for optimization purposes.
template <typename T, std::size_t N>
IVector<T, N> operator * (const T& lhs, const IVector<T, N>& rhs)
{
    auto result = rhs;
    result *= lhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator / (const IVector<T, N>& lhs, const T& rhs)
{
    auto result = lhs;
    result /= rhs;
    return result;
}

template <typename T, std::size_t N>
IVector<T, N> operator / (const T& lhs, const IVector<T, N>& rhs)
{
    auto result = IVector<T, N> { lhs };
    result /= rhs;
    return result;
}


}

#endif // IVECTOR_H
