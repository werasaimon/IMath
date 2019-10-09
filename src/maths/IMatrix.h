/********************************************************************************
*
* IMatrix.h
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

#ifndef IMATRIX_H
#define IMATRIX_H


#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <assert.h>

#include "IReal.h"
#include "IFunc.h"

namespace IMath
{

template <typename T, std::size_t Rows, std::size_t Cols>
class IMatrix;

template <typename T, std::size_t N>
T Determinantt(const IMatrix<T, N, N>& m);

template <typename T, std::size_t N>
bool Inverse(IMatrix<T, N, N>& inv, const IMatrix<T, N, N>& m);


namespace
{

template <typename T, std::size_t Rows, std::size_t Cols>
struct MatrixDefaultInitializer
{
    static void Initialize(IMatrix<T, Rows, Cols>& matrix)
    {
        matrix.Reset();
    }
};

template <typename T, std::size_t N>
struct MatrixDefaultInitializer<T, N, N>
{
    static void Initialize(IMatrix<T, N, N>& matrix)
    {
        matrix.LoadIdentity();
    }
};



#ifdef GS_ROW_MAJOR_STORAGE
#   define GS_FOREACH_ROW_COL(r, c)             \
        for (std::size_t r = 0; r < Rows; ++r)  \
        for (std::size_t c = 0; c < Cols; ++c)
#else
#   define GS_FOREACH_ROW_COL(r, c)             \
        for (std::size_t c = 0; c < Cols; ++c)  \
        for (std::size_t r = 0; r < Rows; ++r)
#endif



//! Internal class for implementation details.
template <template <typename, std::size_t, std::size_t>
          class M,  typename T,
          std::size_t Rows,
          std::size_t Cols>
class MatrixHelper
{

        MatrixHelper() = delete;

    public:


        //friend bool Gs::Inverse<T, Rows>(Matrix<T, Rows, Cols>&, const Matrix<T, Rows, Cols>&);
        static std::vector<T> MatrixToArray(const M<T, Rows, Cols>& mat)
        {
            std::vector<T> vec(Rows*Cols);

            for (std::size_t r = 0, i = 0; r < Rows; ++r)
            {
                for (std::size_t c = 0; c < Cols; ++c)
                    vec[i++] = mat(r, c);
            }

            return vec;
        }

        static T OrderedDeterminant(const std::vector<T>& mat, std::size_t order)
        {
            if (order == 1) return mat[0];
            std::vector<T> minorMat((order - 1)*(order - 1), T());

            T det = T(0.0);
            for (std::size_t i = 0; i < order; ++i)
            {
                GetMinorMatrix(mat, minorMat, 0, i, order);
                if (i % 2 == 1)
                {
                    det -= mat[i] * OrderedDeterminant(minorMat, order - 1);
                }
                else
                {
                    det += mat[i] * OrderedDeterminant(minorMat, order - 1);
                }
            }

            return det;
        }

        static void OrderedInverse(const std::vector<T>& mat, std::size_t order)
        {
            T det = OrderedDeterminant(mat, order);

            //todo...
        }

    private:

        static void GetMinorMatrix( const std::vector<T>& mat, std::vector<T>& minorMat,
                                    std::size_t row, std::size_t column,
                                    std::size_t order)
        {
            for (std::size_t r = 1, i = 0; r < order; ++r)
            {
                if (r != row)
                {
                    for (std::size_t c = 0, j = 0; c < order; ++c)
                    {
                        if (c != column)
                        {
                            minorMat[i*(order - 1) + j] = mat[r*order + c];
                            ++j;
                        }
                    }
                  ++i;
                }
            }
        }

};



}






/**
\brief Base matrix class.
\tparam T Specifies the data type of the matrix components.
This should be a primitive data type such as float, double, int etc.
\remarks The macro GS_ROW_MAJOR_STORAGE can be defined, to use row-major storage layout.
By default column-major storage layout is used.
The macro GS_ROW_VECTORS can be defined, to use row vectors. By default column vectors are used.
Here is an example, how a 4x4 matrix is laid-out with column- and row vectors:
\code
// 4x4 matrix with column vectors:
// / x1 y1 z1 w1 \
// | x2 y2 z2 w2 |
// | x3 y3 z3 w3 |
// \ x4 y4 z4 w4 /
// 4x4 matrix with row vectors:
// / x1 x2 x3 x4 \
// | y1 y2 y3 y4 |
// | z1 z2 z3 z4 |
// \ w1 w2 w3 w4 /
// In both cases, (w1, w2, w3, w4) stores the position in an affine transformation.
\endcode
Matrix elements can be accessed by the bracket operator:
\code
Gs::Matrix4 A;
A(0, 0) = row0column0;
A(2, 1) = row2column1;
\endcode
This is independent of the matrix storage layout and independent of the usage of row- or column vectors.
But the following function is dependent of the usage of row- or column vectors:
\code
// For column vectors:
A.At(2, 1) = row2column1;
// For row vectors:
A.At(2, 1) = row1column2;
\endcode
This function is used for easier support between row- and column vectors.
*/
template <typename T, std::size_t Rows, std::size_t Cols>
class IMatrix
{

    public:

        static_assert(Rows*Cols > 0, "matrices must consist of at least 1x1 elements");

        /* ----- Static members ----- */

        //! Number of rows of this matrix type.
        static const std::size_t rows       = Rows;

        //! Number of columns of this matrix type.
        static const std::size_t columns    = Cols;

        //! Number of scalar elements of this matrix type.
        static const std::size_t elements   = Rows*Cols;

        /* ----- Typenames ----- */

        //! Specifies the typename of the scalar components.
        using ScalarType        = T;

        //! Typename of this matrix type.
        using ThisType          = IMatrix<T, Rows, Cols>;

        //! Typename of the transposed of this matrix type.
        using TransposedType    = IMatrix<T, Cols, Rows>;

        /* ----- Functions ----- */

        /**
        \brief Default constructor.
        \remarks If the 'GS_DISABLE_AUTO_INIT' is NOT defined, the matrix elements will be initialized. Otherwise, the matrix is in an uninitialized state.
        */
        IMatrix()
        {
            #ifndef GS_DISABLE_AUTO_INIT
            MatrixDefaultInitializer<T, Rows, Cols>::Initialize(*this);
            #endif
        }

        //! Copy constructor.
        IMatrix(const ThisType& rhs)
        {
            *this = rhs;
        }

        //! Initializes this matrix with the specified values (row by row, and column by column).
        IMatrix(const std::initializer_list<T>& values)
        {
            std::size_t i = 0, n = values.size();
            for (auto it = values.begin(); i < n; ++i, ++it)
            {
                (*this)(i / columns, i % columns) = *it;
            }

            for (; i < elements; ++i)
            {
                (*this)(i / columns, i % columns) = T(0);
            }
        }

//        /**
//        \brief Explicitly uninitialization constructor.
//        \remarks With this constructor, the matrix is always in an uninitialized state.
//        */
//        explicit Matrix(UninitializeTag)
//        {
//            // do nothing
//        }

        /**
        \brief Returns a reference to a single matrix element at the specified location.
        \param[in] row Specifies the zero-based row index. Must be in the range [0, Rows).
        \param[in] col Specifies the zero-based column index. Must be in the range [0, Cols).
        \throws std::runtime_error If the macro 'GS_ENABLE_ASSERT' and the macro 'GS_ASSERT_EXCEPTION' are defined,
        and either the row or the column is out of range.
        */
        T& operator () (std::size_t row, std::size_t col)
        {
            assert(row < Rows);
            assert(col < Cols);
            #ifdef GS_ROW_MAJOR_STORAGE
            return m_[row*Cols + col];
            #else
            return m_[col*Rows + row];
            #endif
        }

        /**
        \brief Returns a constant reference to a single matrix element at the specified location.
        \param[in] row Specifies the zero-based row index. Must be in the range [0, Rows).
        \param[in] col Specifies the zero-based column index. Must be in the range [0, Cols).
        \throws std::runtime_error If the macro 'GS_ENABLE_ASSERT' and the macro 'GS_ASSERT_EXCEPTION' are defined,
        and either the row or the column is out of range.
        */
        const T& operator () (std::size_t row, std::size_t col) const
        {
            assert(row < Rows);
            assert(col < Cols);
            #ifdef GS_ROW_MAJOR_STORAGE
            return m_[row*Cols + col];
            #else
            return m_[col*Rows + row];
            #endif
        }

        T& operator [] (std::size_t element)
        {
            assert(element < ThisType::elements);
            return m_[element];
        }

        const T& operator [] (std::size_t element) const
        {
            assert(element < ThisType::elements);
            return m_[element];
        }

        ThisType& operator += (const ThisType& rhs)
        {
            for (std::size_t i = 0; i < ThisType::elements; ++i)
                m_[i] += rhs.m_[i];
            return *this;
        }

        ThisType& operator -= (const ThisType& rhs)
        {
            for (std::size_t i = 0; i < ThisType::elements; ++i)
                m_[i] -= rhs.m_[i];
            return *this;
        }

        ThisType& operator *= (const ThisType& rhs)
        {
            *this = (*this * rhs);
            return *this;
        }

        ThisType& operator *= (const T& rhs)
        {
            for (std::size_t i = 0; i < ThisType::elements; ++i)
            {
                m_[i] *= rhs;
            }
            return *this;
        }

        ThisType& operator = (const ThisType& rhs)
        {
            for (std::size_t i = 0; i < ThisType::elements; ++i)
            {
                m_[i] = rhs.m_[i];
            }
            return *this;
        }

        #ifdef GS_ROW_VECTORS

        T& At(std::size_t col, std::size_t row)
        {
            return (*this)(row, col);
        }

        const T& At(std::size_t col, std::size_t row) const
        {
            return (*this)(row, col);
        }

        #else

        T& At(std::size_t row, std::size_t col)
        {
            return (*this)(row, col);
        }

        const T& At(std::size_t row, std::size_t col) const
        {
            return (*this)(row, col);
        }

        #endif

        //! Restes all matrix elements to zero.
        void Reset()
        {
            for (std::size_t i = 0; i < ThisType::elements; ++i)
            {
                m_[i] = T(0);
            }
        }

        //! Loads the identity for this matrix.
        void LoadIdentity()
        {
            GS_FOREACH_ROW_COL(r, c)
            {
                (*this)(r, c) = (r == c ? T(1.0) : T(0.0));
            }
        }

        //! Returns an identity matrix.
        static ThisType Identity()
        {
            ThisType result;
            result.LoadIdentity();
            return result;
        }

        //! Returns a transposed copy of this matrix.
        TransposedType GetTransposed() const
        {
            TransposedType result;
            GS_FOREACH_ROW_COL(r, c)
            {
                result(c, r) = (*this)(r, c);
            }
            return result;
        }

        //! Transposes this matrix.
        void GetTranspose()
        {
            for (std::size_t i = 0; i + 1 < Cols; ++i)
            {
                for (std::size_t j = 1; j + i < Cols; ++j)
                {
                    std::swap(
                        m_[i*(Cols + 1) + j],
                        m_[(j + i)*Cols + i]
                    );
                }
            }
        }

        /**
        \brief Returns the determinant of this matrix.
        \see Gs::Determinant
        */
        T GetDeterminant() const
        {
            return Determinantt(*this);
        }

        /**
        Returns the trace of this matrix: M(0, 0) + M(1, 1) + ... + M(N - 1, N - 1).
        \note This can only be used for squared matrices!
        */
        T GetTrace() const
        {
            static_assert(Rows == Cols, "traces can only be computed for squared matrices");

            T trace = T(0);

            for (std::size_t i = 0; i < Rows; ++i)
                trace += (*this)(i, i);

            return trace;
        }

        IMatrix<T, Rows, Cols> GetInverse() const
        {
            IMatrix<T, Rows, Cols> inv { *this };
            inv.MakeInverse();
            return inv;
        }

        bool MakeInverse()
        {
            IMatrix<T, Rows, Cols> in { *this };
            return Inverse(*this, in);
        }

        //! Returns a pointer to the first element of this matrix.
        T* Ptr()
        {
            return &(m_[0]);
        }

        //! Returns a constant pointer to the first element of this matrix.
        const T* Ptr() const
        {
            return &(m_[0]);
        }

        /**
        Returns a type casted instance of this matrix.
        \tparam C Specifies the static cast type.
        */
        template <typename C> IMatrix<C, Rows, Cols> Cast() const
        {
            IMatrix<C, Rows, Cols> result;

            for (std::size_t i = 0; i < ThisType::elements; ++i)
            {
                result[i] = static_cast<C>(m_[i]);
            }

            return result;
        }


#ifdef ENABLE_STL_SUPPORT

        //----------[ output operator ]----------------------------
        /**
        * Output to stream operator
        * @param lhs Left hand side argument of operator (commonly ostream instance).
        * @param rhs Right hand side argument of operator.
        * @return Left hand side argument - the ostream object passed to operator.
        */
        friend std::ostream& operator <<(std::ostream& lhs, const IMatrix<T,Rows,Cols>& rhs)
        {
            for (int i = 0; i < Rows; i++)
            {
                lhs << "|\t";
                for (int j = 0; j < Cols; j++)
                {
                    lhs << rhs(i,j)  << "\t";
                }
                lhs << "|" << std::endl;
            }
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

#endif


    private:

        T m_[ThisType::elements];

};


/* --- Global Operators --- */

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator + (const IMatrix<T, Rows, Cols>& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = lhs;
    result += rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator - (const IMatrix<T, Rows, Cols>& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = lhs;
    result -= rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const IMatrix<T, Rows, Cols>& lhs, const T& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const T& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = rhs;
    result *= lhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t ColsRows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const IMatrix<T, Rows, ColsRows>& lhs,
                                   const IMatrix<T, ColsRows, Cols>& rhs)
{
    IMatrix<T, Rows, Cols> result;

    GS_FOREACH_ROW_COL(r, c)
    {
        result(c, r) = T(0.0);
        for (std::size_t i = 0; i < ColsRows; ++i)
        {
            result(c, r) += lhs(i, r)*rhs(c, i);
        }
    }

    return result;
}



/**
\brief Computes the determinant of an arbitrary NxN matrix.
\tparam M Specifies the matrix type. This should be "Matrix".
\tparam T Specifies the data type. This should be float or double.
\tparam Rows Specifies the rows of the matrix.
\tparam Cols Specifies the columns of the matrix.
\remarks The template arguments 'Rows' and 'Cols' must be equal, otherwise a compile time error will occur,
since a determinant is only defined for squared matrices.
\param[in] m Specifies the squared matrix for which the determinant is to be computed.
*/
template <typename T, std::size_t N>
T Determinantt(const IMatrix<T, N, N>& m)
{
   using Helper = MatrixHelper<IMatrix, T, N, N>;
   return T(Helper::OrderedDeterminant(Helper::MatrixToArray(m), N));
}


template <typename T, std::size_t N>
bool Inverse(IMatrix<T, N, N>& inv, const IMatrix<T, N, N>& m)
{
    using Helper = MatrixHelper<IMatrix, T, N, N>;
    return Helper::OrderedInverse(Helper::MatrixToArray(m) , N);
    return false;//!!!
}

/* --- Type Alias --- */

#undef GS_FOREACH_ROW_COL

} // /namespace Gs


#endif // IMATRIX_H
