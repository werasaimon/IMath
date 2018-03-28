#ifndef RPMATRIX2X2_H
#define RPMATRIX2X2_H



#include "rpVector2D.hpp"

namespace CMath
{


template<class T> class rpMatrix2x2
{

    private:
	 //-------------------- Attributes --------------------//

	    union
	    {
	        T                mData[4];
	        rpVector2D<T>    mRows[2];
	    };



    public:
        /**
         * Default constructor.
         *
         * Constructs the identity matrix.
         */
	    SIMD_INLINE rpMatrix2x2()
        : mData{1.0f, 0.0f,
        	    0.0f, 1.0f}
        {
                // Nothing to do.
        }


        /// Constructor
        SIMD_INLINE rpMatrix2x2(T value)
        {
            setAllValues(value, value,
                         value, value);
        }


        /**
         * Constructor.
         *
         * Constructs the matrix from the passed array.
         *
         * @param arr Array of floating point values in row-major order.
         */
        SIMD_INLINE rpMatrix2x2(const T arr[4])
        {
            std::memcpy(mData, arr, 4 * sizeof(T));
        }


        /**
         * Constructor.
         *
         * Constructs the matrix from the specified entries.
         *
         * @param entry00 Entry at row 0 column 0.
         * @param entry01 Entry at row 0 column 1.
         * @param entry10 Entry at row 1 column 0.
         * @param entry11 Entry at row 1 column 1.
         * @param entry20 Entry at row 2 column 0.
         */
        SIMD_INLINE rpMatrix2x2(T entry00, T entry01, T entry10, T entry11)
          : mData{entry00, entry01,
                  entry10, entry11}
        {
                // Nothing to do.
        }

        /**
         * Copy constructor.
         *
         * @param other The other matrix to copy.
         */
        SIMD_INLINE rpMatrix2x2(const rpMatrix2x2<T>& other) = default;



        /**
        * Copy casting constructor.
        * @param src Data source for new created instance of rpMatrix3x3
        */
        template<class FromT>
        SIMD_INLINE rpMatrix2x2(const rpMatrix2x2<FromT>& src)
        {
            for (int i = 0; i < 4; i++)
            {
                mData[i] = static_cast<T>(src.mData[i]);
            }
        }



        //---------------------- Methods ---------------------//

        /// Set the matrix to the identity matrix
        SIMD_INLINE void setToIdentity()
        {
        	// Initialize all values in the matrix to identity
        	setAllValues(1.0, 0.0,
        			     0.0, 1.0);
        }


        /// Set the matrix to zero
        SIMD_INLINE void setToZero()
        {
        	mRows[0].setToZero();
        	mRows[1].setToZero();
        }


        /// Set all the values in the matrix
        SIMD_INLINE void setAllValues(T a1, T a2,
                                      T b1, T b2)
        {
        	mRows[0][0] = a1; mRows[0][1] = a2;
     	    mRows[1][0] = b1; mRows[1][1] = b2;

        }



        /**
         * Matrix entry accessor operator.
         *
         * @note Entry indicies are in the range 0 <= `index` <= 8.
         *
         * @param index Index for the entry to return.
         * @return Entry at position `index`.
         */
        SIMD_INLINE T operator[](std::size_t index) const
        {
            assert(index < 4);
            return mData[index];
        }

        /**
         * Equality operator.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return True if the two supplied matrices are equal. False otherwise.
         */
        friend SIMD_INLINE bool operator==(const rpMatrix2x2<T>& A, const rpMatrix2x2<T>& B)
        {
            const T epsilon = MACHINE_EPSILON;
            for (int i = 0; i < 4; ++i)
            {
                if (Abs(A[i] - B[i]) > epsilon) return false;
            }

            return true;
        }

        /**
         * Non-equality operator.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return True if the two supplied matrices are not equal. False
         * otherwise.
         */
        friend SIMD_INLINE  bool operator!=(const rpMatrix2x2<T>& A, const rpMatrix2x2<T>& B)
        {
            return !(A == B);
        }



        //---------------------[ assignment operations ]---------------------------------

        /**
         * Copy operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE rpMatrix2x2<T>& operator=(const rpMatrix2x2<T>& rhs)
        {
        	std::memcpy(mData, rhs.mData, sizeof(T) * 4);
        	return *this;
        }

        /**
         * Copy casting operator
         * @param rhs Right hand side argument of binary operator.
         */
        template<class FromT>
        SIMD_INLINE rpMatrix2x2<T>& operator=(const rpMatrix2x2<FromT>& rhs)
        {
        	for (int i = 0; i < 4; i++)
        	{
        		mData[i] = static_cast<T>(rhs.mData[i]);
        	}
        	return *this;
        }

        /**
         * Copy operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE rpMatrix2x2<T>& operator=(const T* rhs)
        {
        	std::memcpy(mData, rhs, sizeof(T) * 4);
        	return *this;
        }


        //-------------[ conversion data ]-----------------------------

        /**
         * Conversion to pointer operator
         * @return Pointer to internally stored (in management of class rpMatrix2x2<T>)
         * used for passing rpMatrix2x2<T> values to gl*[fd]v functions.
         */
        SIMD_INLINE operator T*()
        {
        	return (T*) &mRows[0][0];
        }

        /**
         * Conversion to pointer operator
         * @return Constant Pointer to internally stored (in management of class rpMatrix2x2<T>)
         * used for passing rpMatrix2x2<T> values to gl*[fd]v functions.
         */
        SIMD_INLINE operator const T*() const
        {
        	return (const T*) &mRows[0][0];
        }


        /**
         * Conversion to pointer operator
         */
        SIMD_INLINE const T* getData() const
        {
        	return &mRows[0][0];
        }

        /**
         * Conversion to pointer operator
         */
        SIMD_INLINE T* getData()
        {
        	return &mRows[0][0];
        }




        /// Overloaded operator to read element of the matrix.
        SIMD_INLINE const rpVector2D<T>& operator[](int row) const
        {
        	return mRows[row];
        }

        /// Overloaded operator to read/write element of the matrix.
        SIMD_INLINE rpVector2D<T>& operator[](int row)
        {
        	return mRows[row];
        }




        /// Return a column
        SIMD_INLINE rpVector2D<T> getColumn(int i) const
        {
        	assert(i>= 0 && i<2);
        	return rpVector2D<T> (mRows[0][i], mRows[1][i] );
        }

        /// Return a row
        SIMD_INLINE rpVector2D<T> getRow(int i) const
        {
        	assert(i>= 0 && i<2);
        	return mRows[i];
        }


        /**
         * Matrix addition operator.
         *
         * @note Matrix addition is commutative.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return The matrix equal to the sum of `A` and `B`.
         */
        SIMD_INLINE rpMatrix2x2<T> operator+(const rpMatrix2x2<T>& rhs) const
        {
        	rpMatrix2x2<T> ret;
        	for (int i = 0; i < 2; i++)  ret.mRows[i] = mRows[i] + rhs.mRows[i];
        	return ret;

        }



        /**
         * Matrix subtraction operator.
         *
         * @note Matrix subtraction is not commutative.
         *
         * @param lhs The left hand side matrix.
         * @param rhs The right hand side matrix.
         * @return The matrix equal to the `rhs` matrix subtracted from the `lhs`
         * matrix.
         */
        SIMD_INLINE rpMatrix2x2<T> operator-(const rpMatrix2x2<T>& rhs) const
        {
        	rpMatrix2x2<T> ret;
        	for (int i = 0; i < 2; i++) ret.mRows[i] = mRows[i] - rhs.mRows[i];
        	return ret;
        }



        /**
         * Matrix negation operator.
         *
         * @param A The matrix to negate.
         * @return The additive inverse of the matrix `A`.
         */
        friend SIMD_INLINE rpMatrix2x2<T> operator-(const rpMatrix2x2<T> &A)
        {
        	return rpMatrix2x2<T>( -A[0], -A[1],
        			               -A[2], -A[3]
        	);
        }

        /**
         * Scalar multiplication operator.
         *
         * Multiplies each entry of a matrix by a given scalar value.
         *
         * @param A The matrix to be multiplied by the given scalar.
         * @param s The scalar value.
         * @return The matrix `A` multiplied by the scalar `s`.
         */
        friend SIMD_INLINE rpMatrix2x2<T> operator*(const rpMatrix2x2<T>& matrix, const T nb)
        {
        	return rpMatrix2x2<T>(matrix.mRows[0][0] * nb, matrix.mRows[0][1] * nb,
        	        			  matrix.mRows[1][0] * nb, matrix.mRows[1][1] * nb);
        }


        /**
         * Scalar multiplication operator.
         *
         * Multiplies each entry of a matrix by a given scalar value.
         *
         * @param s The scalar value.
         * @param A The matrix to be multiplied by the given scalar.
         * @return The matrix `A` multiplied by the scalar `s`.
         */
        friend SIMD_INLINE rpMatrix2x2<T> operator*(const T s, const rpMatrix2x2<T>& A)
        {
        	return A * s;
        }

        /// Overloaded operator for inveret multiplication with a number
        friend SIMD_INLINE rpMatrix2x2<T>  operator/(T nb, const rpMatrix2x2<T>& matrix)
        {
        	return rpMatrix2x2<T>(matrix.mRows[0][0] / nb, matrix.mRows[0][1] / nb,
        			              matrix.mRows[1][0] / nb, matrix.mRows[1][1] / nb);
        }

        /// Overloaded operator for inveret multiplication with a matrix
        friend SIMD_INLINE rpMatrix2x2<T>  operator/(const rpMatrix2x2<T>& matrix, T nb)
        {
        	return nb / matrix;
        }

        /**
         * Vector multiplication operator.
         *
         * Multiplies the matrix `rhs` on the left by the row vector `lhs`,
         * returning the resulting vector.
         *
         * @param lhs The matrix.
         * @param rhs The row vector.
         * @return The vector `lhs` multiplied by the matrix `rhs` on the right.
         */
        friend SIMD_INLINE rpVector2D<T> operator*(const rpVector2D<T>& rhs, const rpMatrix2x2<T>& lhs)
        {
        	return lhs * rhs;
        }

        /**
         * Vector multiplication operator.
         *
         * Multiplies the column vector `rhs` on the left by the matrix `lhs`,
         * returning the resulting vector.
         *
         * @param lhs The matrix.
         * @param rhs The column vector.
         * @return The vector `rhs` multiplied by the matrix `lhs` on the left.
         */
        SIMD_INLINE rpVector2D<T> operator*(const rpVector2D<T>& rhs) const
        {
        	T fX = rhs.x;
        	T fY = rhs.y;

        	rpVector2D<T> Point;
        	Point.x = ( fX * mRows[0][0] + fY * mRows[0][1] );
        	Point.y = ( fX * mRows[1][0] + fY * mRows[1][1] );

        	return Point;

        }



        /**
         * Matrix multiplication operator.
         *
         * @note Matrix multiplication is not commutative.
         *
         * @param lhs The left hand side matrix.
         * @param rhs The right hand side matrix.
         * @return The matrix equal to the product `lhs` x `rhs`.
         */
        SIMD_INLINE rpMatrix2x2<T> operator*(rpMatrix2x2<T> rhs) const
        {
        	rpMatrix2x2<T> w;
        	for(int i = 0; i < 2; i++)
        	{
        		for (int j = 0; j < 2; j++)
        		{
        			T n = 0;
        			for (int k = 0; k < 2; k++)
        			{
        				n += rhs.mRows[i][k] * mRows[k][j];
        			}
        			w.mRows[i][j] = n;
        		}
        	}
        	return w;
        }


        //---------------------------[ misc operations ]----------------------------


        /**
         * Returns the determinant of the matrix.
         *
         * @note A square matrix is invertable if and only if its determinant is
         * nonzero.
         *
         * @return The determinant.
         */
        SIMD_INLINE T getDeterminant() const
        {
        	return mData[0] * mData[3] - mData[1] * mData[2];
        }

        /**
         * Returns a copy of the multiplicitive inverse of this matrix.
         *
         * @return The multiplicitive inverse of this matrix.
         */
        SIMD_INLINE rpMatrix2x2<T> getInverse() const
        {
        	// Ensure that the matrix is not singular.
        	const T det = getDeterminant();
        	assert(det != 0.0f);

        	// Return a copy of the inverse of this matrix.
        	const T invDet = 1.0f / det;
        	return rpMatrix2x2<T>( mData[3] * invDet, -mData[1] * invDet,
        			              -mData[2] * invDet,  mData[0] * invDet );
        }

        /**
         * Returns a copy of this matrix transposed so that the rows now form
         * columns.
         *
         * @return Transposed copy of this matrix.
         */
        SIMD_INLINE rpMatrix2x2<T> getTranspose() const
        {
        	rpMatrix2x2<T> ret;
        	for (int i = 0; i < 2; i++)
        	{
        		for (int j = 0; j < 2; j++)
        		{
        			ret.mRows[i][j] = mRows[j][i];
        		}
        	}
        	return ret;
        }


        /**
         * Return the matrix with absolute values
         */
        SIMD_INLINE rpMatrix2x2<T> getAbsoluteMatrix() const
        {
        	return rpMatrix2x2<T>(Abs(mRows[0][0]), Abs(mRows[0][1]),
        			              Abs(mRows[1][0]), Abs(mRows[1][1]));
        }


        /**
         * Return the trace of the matrix
         */
        SIMD_INLINE T getTrace() const
        {
        	// Compute and return the trace
        	return (mRows[0][0] + mRows[1][1]);
        }


        //--------------------------------------------------------------------//


        /// Overloaded operator for addition with assignment
        SIMD_INLINE rpMatrix2x2<T>& operator+=(const rpMatrix2x2<T>& matrix)
        {
        	mRows[0] += matrix.mRows[0];
        	mRows[1] += matrix.mRows[1];
        	return *this;
        }

        /// Overloaded operator for substraction with assignment
        SIMD_INLINE rpMatrix2x2<T>& operator-=(const rpMatrix2x2<T>& matrix)
        {
        	mRows[0] -= matrix.mRows[0];
        	mRows[1] -= matrix.mRows[1];
        	return *this;
        }

        /// Overloaded operator for multiplication with a number with assignment
        SIMD_INLINE rpMatrix2x2<T>& operator*=(T nb)
        {
        	mRows[0] *= nb;
        	mRows[1] *= nb;
        	return *this;
        }

        /// Overloaded operator for invert multiplication with a number with assignment
        SIMD_INLINE rpMatrix2x2<T> &operator/=(T nb)
        {
        	mRows[0] /= nb;
        	mRows[1] /= nb;
        	return *this;
        }



        //-------------------------- Plugins ----------------------------//


        /// Return a skew-symmetric matrix using a given vector that can be used
        /// to compute cross product with another vector using matrix multiplication
        static SIMD_INLINE rpMatrix2x2<T> computeSkewSymmetricMatrixForCrossProduct(const rpVector2D<T>& vector)
        {
            return rpMatrix2x2<T>(0 , -vector.x,
                                  vector.y , 0 );
        }


        /// Return a symmetric matrix using a given vector that can be used
        /// to compute dot product with another vector using matrix multiplication
         static SIMD_INLINE rpMatrix2x2<T> computeSymmetricMatrix(const rpVector2D<T>& vector)
         {
             return rpMatrix2x2<T>(0 , vector.x,
                                   vector.y , 0 );
         }

        /**
         * Returns a scaling matrix that scales by `scaleFactors.x` and
         * 'scaleFactors.y' in the x and y axes respectively.
         *
         * @param scaleFactors Scale factors.
         * @return Scaling matrix.
         */
        static SIMD_INLINE rpMatrix2x2<T> createScaling(const rpVector2D<T>& scaleFactors)
        {
            return rpMatrix2x2<T>(  scaleFactors.x, 0.0f,
                                    0.0f, scaleFactors.y);
        }

        /**
         * Returns a scaling matrix that scales by `factor` uniformly.
         *
         * @param scale Uniform scale factor.
         * @return Scaling matrix.
         */
        static SIMD_INLINE rpMatrix2x2<T> createScaling(const T factor)
        {
            return rpMatrix2x2<T>( factor, 0.0f,
                                   0.0f, factor);
        }

        /**
         * Returns a rotation matrix that rotates by `angle` radians.
         *
         * @param angle Angle (in radians) for the rotation.
         * @return Rotation matrix that rotates `angle` radians
         * counter-clockwise.
         */
        static SIMD_INLINE rpMatrix2x2<T> createRotation(const T angle)
        {

            const T cosTheta = Cos(angle);
            const T sinTheta = Sin(angle);
            return rpMatrix2x2<T>( cosTheta, -sinTheta,
                                   sinTheta,  cosTheta);
        }

        /**
         * Returns a rotation matrix that represents the sortest rotation from
         * the `fromDirection` to the `toDirection`.
         *
         * @param fromDirection Direction for the matrix to rotate from.
         * @param toDirection Direction for the matrix to rotate to.
         * @return Rotation matrix corresponding to the rotation from
         * `fromDirection` to `toDirection`.
         */
        static SIMD_INLINE rpMatrix2x2<T> createFromToRotation(const rpVector2D<T>& fromDirection, const rpVector2D<T>& toDirection)
        {
            assert(fromDirection.sqrMagnitude() > 0.0f && toDirection.sqrMagnitude() > 0.0f);

            // Compute the angle between the two vectors.
            const T theta = angle(fromDirection, toDirection);

            // Return the rotation matrix.
            return angleRotation(theta);
        }



        //----------[ output operator ]----------------------------
        /**
        * Output to stream operator
        * @param lhs Left hand side argument of operator (commonly ostream instance).
        * @param rhs Right hand side argument of operator.
        * @return Left hand side argument - the ostream object passed to operator.
        */
        friend std::ostream& operator <<(std::ostream& lhs, const rpMatrix2x2<T>& rhs)
        {
            for (int i = 0; i < 2; i++)
            {
                lhs << "|\t";
                for (int j = 0; j < 2; j++)
                {
                    lhs << rhs[i][j]  << "\t";
                }
                lhs << "|" << std::endl;
            }
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


    public:
        /**
         * The multiplicitive identity matrix.
         */
        static const rpMatrix2x2<T> IDENTITY;

        /**
         * The additive identity matrix.
         */
        static const rpMatrix2x2<T> ZERO;

};

template<class T> const rpMatrix2x2<T> rpMatrix2x2<T>::IDENTITY = rpMatrix2x2<T>(1.0, 0.0,
																				 0.0, 1.0);

template<class T> const rpMatrix2x2<T> rpMatrix2x2<T>::ZERO = rpMatrix2x2<T>(0.0, 0.0,
																			 0.0, 0.0);


template<class T>  SIMD_INLINE rpMatrix2x2<T> operator ^ (const rpVector2D<T> lhs , const rpVector2D<T> rhs )
{
	return rpMatrix2x2<T>( lhs.x * rhs.x , lhs.x * rhs.y,
			               lhs.y * rhs.x , lhs.y * rhs.y);
}

/// Matrix 3x3 of floats
typedef rpMatrix2x2<float> rpMatrix2x2f;
/// Matrix 3x3 of doubles
typedef rpMatrix2x2<double> rpMatrix2x2d;
/// Matrix 3x3 of int
typedef rpMatrix2x2<int> rpMatrix2x2i;



} /* namespace */


#endif // RPMATRIX2X2_H
