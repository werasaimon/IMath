
#ifndef SRC_PHYSICS_ENGINE_LINEARMATHS_RPMATRIX4X4_H_
#define SRC_PHYSICS_ENGINE_LINEARMATHS_RPMATRIX4X4_H_

#include "rpLorentzVector.hpp"
#include "rpMatrix3x3.hpp"


namespace CMath
{

/**
 * Class for matrix 4x4
 * @note Data stored in this matrix are in column major order. This arrangement suits OpenGL.
 * If you're using row major matrix, consider using fromRowMajorArray as way for construction
 * Matrix4<T> instance.
 */
template<class T>
class rpMatrix4x4
{
 private:

    //-------------------- Attributes --------------------//

    union
    {
        T             mData[16];
        rpVector4D<T> mRows[4];
    };

 public:

    //--------------------------[ constructors ]-------------------------------
    /**
     *Creates identity matrix
     */
    SIMD_INLINE rpMatrix4x4()
    {
        // Initialize all values in the matrix to zero
        setAllValues(1.0, 0.0, 0.0, 0.0 ,
                     0.0, 1.0, 0.0, 0.0 ,
                     0.0, 0.0, 1.0, 0.0 ,
                     0.0, 0.0, 0.0, 1.0);
    }

    /**
     * Copy matrix values from array (these data must be in column
     * major order!)
     */
    SIMD_INLINE rpMatrix4x4(const T * dt)
    {
        std::memcpy(mData, dt, sizeof(T) * 16);
    }

    /**
     * Copy constructor.
     * @param src Data source for new created instance of rpMatrix4x4.
     */
    SIMD_INLINE rpMatrix4x4(const rpMatrix4x4<T>& src)
    {
        std::memcpy(mData, src.mData, sizeof(T) * 16);
    }

    // Constructor
    SIMD_INLINE rpMatrix4x4<T>(T value)
    {
        setAllValues(value, value, value, value,
                     value, value, value, value,
                     value, value, value, value,
                     value, value, value, value);
    }

    // Constructor
    SIMD_INLINE rpMatrix4x4<T>(T a1, T a2, T a3, T a4,
                               T b1, T b2, T b3, T b4,
                               T c1, T c2, T c3, T c4,
                               T d1, T d2, T d3, T d4)
    {
        // Initialize the matrix with the values
        setAllValues(a1, a2, a3, a4 ,
                     b1, b2, b3, b4 ,
                     c1, c2, c3, c4 ,
                     d1, d2, d3, d4);
    }


    // Constructor with arguments
    SIMD_INLINE rpMatrix4x4( T data[4][4] )
    {
        setAllValues(data[0][0], data[0][1], data[0][2], data[0][3],
                     data[1][0], data[1][1], data[1][2], data[1][3],
                     data[2][0], data[2][1], data[2][2], data[2][3],
                     data[3][0], data[3][1], data[3][2], data[3][3]);
    }


    /**
     * Copy casting constructor.
     * @param src Data source for new created instance of rpMatrix4x4.
     */
    template<class FromT>
    SIMD_INLINE rpMatrix4x4(const rpMatrix4x4<FromT>& src)
    {
        for (int i = 0; i < 16; i++)
        {
            mData[i] = static_cast<T>(src.mData[i]);
        }
    }



   //---------------------- Methods ---------------------//

    /**
    * Resets matrix to be identity matrix
    */
    SIMD_INLINE void setToIdentity()
    {
    	// Initialize all values in the matrix to identity
    	setAllValues(1.0, 0.0, 0.0, 0.0,
    			     0.0, 1.0, 0.0, 0.0,
				     0.0, 0.0, 1.0, 0.0,
				     0.0, 0.0, 0.0, 1.0);
    }



    /**
    * Resets matrix to zero
    */
    SIMD_INLINE void setToZero()
    {
      // Initialize all values in the matrix to zero
      mRows[0].setToZero();
      mRows[1].setToZero();
      mRows[2].setToZero();
      mRows[3].setToZero();
    }



   /**
    * Resets matrix to be value matrix
    */
    SIMD_INLINE void setAllValues(T a1, T a2, T a3, T a4,
                                  T b1, T b2, T b3, T b4,
                                  T c1, T c2, T c3, T c4,
                                  T d1, T d2, T d3, T d4)
    {
        mRows[0][0] = a1; mRows[0][1] = a2; mRows[0][2] = a3; mRows[0][3] = a4;
        mRows[1][0] = b1; mRows[1][1] = b2; mRows[1][2] = b3; mRows[1][3] = b4;
        mRows[2][0] = c1; mRows[2][1] = c2; mRows[2][2] = c3; mRows[2][3] = c4;
        mRows[3][0] = d1; mRows[3][1] = d2; mRows[3][2] = d3; mRows[3][3] = d4;
    }



    //---------------------[ Equality operators ]------------------------------

    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold EPSILON value. To be two
     * values equal, must satisfy this condition all elements of matrix
     * | lhs[i] - rhs[i] | < EPSILON,
     * same for y-coordinate, z-coordinate, and w-coordinate.
     */
    SIMD_INLINE bool operator==(const rpMatrix4x4<T>& rhs) const
    {
        for (int i = 0; i < 16; i++)
        {
            if (Abs(mData[i] - rhs.mData[i]) >= EPSILON )
                return false;
        }
        return true;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const rpMatrix4x4<T>& rhs) const
    {
        return !(*this == rhs);
    }





    //-------------[ conversion data ]-----------------------------

    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class rpMatrix4x4<T>)
     * used for passing rpMatrix4x4<T> values to gl*[fd]v functions.
     */
    SIMD_INLINE operator T*()
    {
        return (T*) &mRows[0][0];
    }

    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class rpMatrix4x4<T>)
     * used for passing rpMatrix4x4<T> values to gl*[fd]v functions.
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
    SIMD_INLINE const rpVector4D<T>& operator[](int row) const
    {
       return mRows[row];
    }

    /// Overloaded operator to read/write element of the matrix.
    SIMD_INLINE rpVector4D<T>& operator[](int row)
    {
       return mRows[row];
    }



    /// Return a column
    SIMD_INLINE rpVector4D<T> getColumn(int i) const
    {
        assert(i>= 0 && i<4);
        return rpVector4D<T> (mRows[0][i], mRows[1][i], mRows[2][i] , mRows[3][i]);
    }

    /// Return a row
    SIMD_INLINE rpVector4D<T> getRow(int i) const
    {
        assert(i>= 0 && i<4);
        return mRows[i];
    }



    //---------------------[ assignment operations ]---------------------------------

    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpMatrix4x4<T>& operator=(const rpMatrix4x4<T>& rhs)
    {
        std::memcpy(mData, rhs.mData, sizeof(T) * 16);
        return *this;
    }

    /**
     * Copy casting operator
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE rpMatrix4x4<T>& operator=(const rpMatrix4x4<FromT>& rhs)
    {
        for (int i = 0; i < 16; i++)
        {
            mData[i] = static_cast<T>(rhs.mData[i]);
        }
        return *this;
    }

    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE rpMatrix4x4<T>& operator=(const T* rhs)
    {
        std::memcpy(mData, rhs, sizeof(T) * 16);
        return *this;
    }


    /**
     * Sets rotation part (matrix 3x3) of matrix.
     *
     * @param m Rotation part of matrix
     */
    SIMD_INLINE void setRotation(const rpMatrix3x3<T>& m)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                mRows[i][j]  = m.mRows[i][j];
            }
        }
    }



    //--------------------[ matrix with matrix operations ]---------------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
     SIMD_INLINE rpMatrix4x4<T> operator+(const rpMatrix4x4<T>& rhs) const
     {
        rpMatrix4x4<T> ret;
        for (int i = 0; i < 4; i++)   ret.mRows[i] = mRows[i] + rhs.mRows[i];
        return ret;
     }

     /**
      * Subtraction operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE rpMatrix4x4<T> operator-(const rpMatrix4x4<T>& rhs) const
      {
         rpMatrix4x4<T> ret;
         for (int i = 0; i < 4; i++)   ret.mRows[i] = mRows[i] - rhs.mRows[i];
         return ret;
      }






    //--------------------[ matrix with scalar operations ]---------------------

    /**
    * Addition operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE rpMatrix4x4<T> operator+(T rhs) const
    {
       rpMatrix4x4<T> ret;
       for (int i = 0; i < 4; i++) ret.mRows[i] = mRows[i] + rhs;
       return ret;
    }

    /**
    * Subtraction operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE rpMatrix4x4<T> operator-(T rhs) const
    {
      rpMatrix4x4<T> ret;
      for (int i = 0; i < 4; i++) ret.mRows[i] = mRows[i] - rhs;
      return ret;
    }


    /// Overloaded operator for multiplication with a number
    friend SIMD_INLINE rpMatrix4x4<T>  operator*(T nb, const rpMatrix4x4<T>& matrix)
    {
        return rpMatrix4x4<T>(matrix.mRows[0][0] * nb, matrix.mRows[0][1] * nb, matrix.mRows[0][2] * nb, matrix.mRows[0][3] * nb,
                              matrix.mRows[1][0] * nb, matrix.mRows[1][1] * nb, matrix.mRows[1][2] * nb, matrix.mRows[1][3] * nb,
                              matrix.mRows[2][0] * nb, matrix.mRows[2][1] * nb, matrix.mRows[2][2] * nb, matrix.mRows[2][3] * nb,
                              matrix.mRows[3][0] * nb, matrix.mRows[3][1] * nb, matrix.mRows[3][2] * nb, matrix.mRows[3][3] * nb);
    }

    /// Overloaded operator for multiplication with a matrix
    friend SIMD_INLINE rpMatrix4x4<T>  operator*(const rpMatrix4x4<T>& matrix, T nb)
    {
        return nb * matrix;
    }


    /// Overloaded operator for invert multiplication with a number
    friend SIMD_INLINE rpMatrix4x4<T>  operator/(T nb, const rpMatrix4x4<T>& matrix)
    {
        return rpMatrix4x4<T>(matrix.mRows[0][0] / nb, matrix.mRows[0][1] / nb, matrix.mRows[0][2] / nb, matrix.mRows[0][3] / nb,
                              matrix.mRows[1][0] / nb, matrix.mRows[1][1] / nb, matrix.mRows[1][2] / nb, matrix.mRows[1][3] / nb,
                              matrix.mRows[2][0] / nb, matrix.mRows[2][1] / nb, matrix.mRows[2][2] / nb, matrix.mRows[2][3] / nb,
                              matrix.mRows[3][0] / nb, matrix.mRows[3][1] / nb, matrix.mRows[3][2] / nb, matrix.mRows[3][3] / nb);
    }

    /// Overloaded operator for invert multiplication with a matrix
    friend SIMD_INLINE rpMatrix4x4<T>  operator/(const rpMatrix4x4<T>& matrix, T nb)
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
    friend SIMD_INLINE rpVector3D<T> operator*(const rpVector3D<T>& rhs , const rpMatrix4x4<T>& lhs)
    {
    	return lhs * rhs;
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
    friend SIMD_INLINE rpVector4D<T> operator*(const rpVector4D<T>& rhs , const rpMatrix4x4<T>& lhs)
    {
    	return lhs * rhs;
    }

    //--------------------[ multiply operators ]--------------------------------
    /**
    * Multiplication operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE rpVector4D<T> operator*(const rpVector4D<T>& rhs) const
    {
        rpVector4D<T> u =rpVector4D<T>(mRows[0][0]*rhs.x + mRows[0][1]*rhs.y + mRows[0][2]*rhs.z + rhs.w*mRows[0][3],
                                       mRows[1][0]*rhs.x + mRows[1][1]*rhs.y + mRows[1][2]*rhs.z + rhs.w*mRows[1][3],
                                       mRows[2][0]*rhs.x + mRows[2][1]*rhs.y + mRows[2][2]*rhs.z + rhs.w*mRows[2][3],
                                       mRows[3][0]*rhs.x + mRows[3][1]*rhs.y + mRows[3][2]*rhs.z + rhs.w*mRows[3][3]);

        return u;
    }


    /**
    * Multiplication operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE rpLorentzVector<T> operator*(const rpLorentzVector<T>& rhs) const
    {
       rpLorentzVector<T> u = rpLorentzVector<T>(mRows[0][0]*rhs.x + mRows[0][1]*rhs.y + mRows[0][2]*rhs.z + rhs.t*mRows[0][3],
                                                 mRows[1][0]*rhs.x + mRows[1][1]*rhs.y + mRows[1][2]*rhs.z + rhs.t*mRows[1][3],
                                                 mRows[2][0]*rhs.x + mRows[2][1]*rhs.y + mRows[2][2]*rhs.z + rhs.t*mRows[2][3],
                                                 mRows[3][0]*rhs.x + mRows[3][1]*rhs.y + mRows[3][2]*rhs.z + rhs.t*mRows[3][3]);

        return u;
    }


    /**
    * Multiplication operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE rpVector3D<T> operator*(const rpVector3D<T>& rhs) const
    {

    	rpVector3D<T> Point;
    	Point.x = ( rhs.x * mRows[0][0] + rhs.y * mRows[0][1] + rhs.z * mRows[0][2] + mRows[3][0]);
    	Point.y = ( rhs.x * mRows[1][0] + rhs.y * mRows[1][1] + rhs.z * mRows[1][2] + mRows[3][1]);
    	Point.z = ( rhs.x * mRows[2][0] + rhs.y * mRows[2][1] + rhs.z * mRows[2][2] + mRows[3][2]);

    	T w = mRows[0][3]*rhs.x +
    		  mRows[1][3]*rhs.y +
			  mRows[2][3]*rhs.z +
			  mRows[3][3];

    	return Point / w;

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
    SIMD_INLINE rpMatrix4x4<T> operator*(rpMatrix4x4<T> rhs) const
    {
    	rpMatrix4x4<T> w;
    	for(int i = 0; i < 4; i++)
    	{
    		for (int j = 0; j < 4; j++)
    		{
    			T n = 0;
    			for (int k = 0; k < 4; k++)
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
    * Computes of matrix_3x3
    * @return recombinate of matrix_3x3
    */
    SIMD_INLINE rpMatrix3x3<T> getRotMatrix() const
    {
        rpMatrix3x3<T> M;

        M[0][0] = mRows[0][0];
        M[1][0] = mRows[1][0];
        M[2][0] = mRows[2][0];

        M[0][1] = mRows[0][1];
        M[1][1] = mRows[1][1];
        M[2][1] = mRows[2][1];

        M[0][2] = mRows[0][2];
        M[1][2] = mRows[1][2];
        M[2][2] = mRows[2][2];

        return M;
    }


    /**
     * Computes determinant of matrix
     * @return Determinant of matrix
     * @note This function does 3 * 4 * 6 mul, 3 * 6 add.
     */
    SIMD_INLINE T getDeterminant() const
    {

        return  + mRows[3][0] * mRows[2][1] * mRows[1][2] * mRows[0][3] - mRows[2][0] * mRows[3][1] * mRows[1][2] * mRows[0][3]
                - mRows[3][0] * mRows[1][1] * mRows[2][2] * mRows[0][3] + mRows[1][0] * mRows[3][1] * mRows[2][2] * mRows[0][3]

                + mRows[2][0] * mRows[1][1] * mRows[3][2] * mRows[0][3] - mRows[1][0] * mRows[2][1] * mRows[3][2] * mRows[0][3]
                - mRows[3][0] * mRows[2][1] * mRows[0][2] * mRows[1][3] + mRows[2][0] * mRows[3][1] * mRows[0][2] * mRows[1][3]

                + mRows[3][0] * mRows[0][1] * mRows[2][2] * mRows[1][3] - mRows[0][0] * mRows[3][1] * mRows[2][2] * mRows[1][3]
                - mRows[2][0] * mRows[0][1] * mRows[3][2] * mRows[1][3] + mRows[0][0] * mRows[2][1] * mRows[3][2] * mRows[1][3]

                + mRows[3][0] * mRows[1][1] * mRows[0][2] * mRows[2][3] - mRows[1][0] * mRows[3][1] * mRows[0][2] * mRows[2][3]
                - mRows[3][0] * mRows[0][1] * mRows[1][2] * mRows[2][3] + mRows[0][0] * mRows[3][1] * mRows[1][2] * mRows[2][3]

                + mRows[1][0] * mRows[0][1] * mRows[3][2] * mRows[2][3] - mRows[0][0] * mRows[1][1] * mRows[3][2] * mRows[2][3]
                - mRows[2][0] * mRows[1][1] * mRows[0][2] * mRows[3][3] + mRows[1][0] * mRows[2][1] * mRows[0][2] * mRows[3][3]

                + mRows[2][0] * mRows[0][1] * mRows[1][2] * mRows[3][3] - mRows[0][0] * mRows[2][1] * mRows[1][2] * mRows[3][3]
                - mRows[1][0] * mRows[0][1] * mRows[2][2] * mRows[3][3] + mRows[0][0] * mRows[1][1] * mRows[2][2] * mRows[3][3];

    }

    /**
     * Computes inverse matrix
     * @return Inverse matrix of this matrix.
     * @note This is a little bit time consuming operation
     * (16 * 6 * 3 mul, 16 * 5 add + det() + mul() functions)
     */
    SIMD_INLINE rpMatrix4x4<T> getInverse() const
    {
        // Compute the determinant of the matrix
        T determinant = getDeterminant();

        // Check if the determinant is equal to zero
        assert(Abs(determinant) > MACHINE_EPSILON);

        rpMatrix4x4<T> ret;

                  ret.mRows[0][0] = + mRows[2][1] * mRows[3][2] * mRows[1][3] - mRows[3][1] * mRows[2][2] * mRows[1][3] + mRows[3][1] * mRows[1][2] * mRows[2][3]
                                    - mRows[1][1] * mRows[3][2] * mRows[2][3] - mRows[2][1] * mRows[1][2] * mRows[3][3] + mRows[1][1] * mRows[2][2] * mRows[3][3];

                  ret.mRows[1][0] = + mRows[3][0] * mRows[2][2] * mRows[1][3] - mRows[2][0] * mRows[3][2] * mRows[1][3] - mRows[3][0] * mRows[1][2] * mRows[2][3]
                                    + mRows[1][0] * mRows[3][2] * mRows[2][3] + mRows[2][0] * mRows[1][2] * mRows[3][3] - mRows[1][0] * mRows[2][2] * mRows[3][3];

                  ret.mRows[2][0] = + mRows[2][0] * mRows[3][1] * mRows[1][3] - mRows[3][0] * mRows[2][1] * mRows[1][3] + mRows[3][0] * mRows[1][1] * mRows[2][3]
                                    - mRows[1][0] * mRows[3][1] * mRows[2][3] - mRows[2][0] * mRows[1][1] * mRows[3][3] + mRows[1][0] * mRows[2][1] * mRows[3][3];

                  ret.mRows[3][0] = + mRows[3][0] * mRows[2][1] * mRows[1][2] - mRows[2][0] * mRows[3][1] * mRows[1][2] - mRows[3][0] * mRows[1][1] * mRows[2][2]
                                    + mRows[1][0] * mRows[3][1] * mRows[2][2] + mRows[2][0] * mRows[1][1] * mRows[3][2] - mRows[1][0] * mRows[2][1] * mRows[3][2];

                  ret.mRows[0][1] = + mRows[3][1] * mRows[2][2] * mRows[0][3] - mRows[2][1] * mRows[3][2] * mRows[0][3] - mRows[3][1] * mRows[0][2] * mRows[2][3]
                                    + mRows[0][1] * mRows[3][2] * mRows[2][3] + mRows[2][1] * mRows[0][2] * mRows[3][3] - mRows[0][1] * mRows[2][2] * mRows[3][3];

                  ret.mRows[1][1] = + mRows[2][0] * mRows[3][2] * mRows[0][3] - mRows[3][0] * mRows[2][2] * mRows[0][3] + mRows[3][0] * mRows[0][2] * mRows[2][3]
                                    - mRows[0][0] * mRows[3][2] * mRows[2][3] - mRows[2][0] * mRows[0][2] * mRows[3][3] + mRows[0][0] * mRows[2][2] * mRows[3][3];

                  ret.mRows[2][1] = + mRows[3][0] * mRows[2][1] * mRows[0][3] - mRows[2][0] * mRows[3][1] * mRows[0][3] - mRows[3][0] * mRows[0][1] * mRows[2][3]
                                    + mRows[0][0] * mRows[3][1] * mRows[2][3] + mRows[2][0] * mRows[0][1] * mRows[3][3] - mRows[0][0] * mRows[2][1] * mRows[3][3];

                  ret.mRows[3][1] = + mRows[2][0] * mRows[3][1] * mRows[0][2] - mRows[3][0] * mRows[2][1] * mRows[0][2] + mRows[3][0] * mRows[0][1] * mRows[2][2]
                                    - mRows[0][0] * mRows[3][1] * mRows[2][2] - mRows[2][0] * mRows[0][1] * mRows[3][2] + mRows[0][0] * mRows[2][1] * mRows[3][2];

                  ret.mRows[0][2] = + mRows[1][1] * mRows[3][2] * mRows[0][3] - mRows[3][1] * mRows[1][2] * mRows[0][3] + mRows[3][1] * mRows[0][2] * mRows[1][3]
                                    - mRows[0][1] * mRows[3][2] * mRows[1][3] - mRows[1][1] * mRows[0][2] * mRows[3][3] + mRows[0][1] * mRows[1][2] * mRows[3][3];

                  ret.mRows[1][2] = + mRows[3][0] * mRows[1][2] * mRows[0][3] - mRows[1][0] * mRows[3][2] * mRows[0][3] - mRows[3][0] * mRows[0][2] * mRows[1][3]
                                    + mRows[0][0] * mRows[3][2] * mRows[1][3] + mRows[1][0] * mRows[0][2] * mRows[3][3] - mRows[0][0] * mRows[1][2] * mRows[3][3];

                  ret.mRows[2][2] = + mRows[1][0] * mRows[3][1] * mRows[0][3] - mRows[3][0] * mRows[1][1] * mRows[0][3] + mRows[3][0] * mRows[0][1] * mRows[1][3]
                                    - mRows[0][0] * mRows[3][1] * mRows[1][3] - mRows[1][0] * mRows[0][1] * mRows[3][3] + mRows[0][0] * mRows[1][1] * mRows[3][3];

                  ret.mRows[3][2] = + mRows[3][0] * mRows[1][1] * mRows[0][2] - mRows[1][0] * mRows[3][1] * mRows[0][2] - mRows[3][0] * mRows[0][1] * mRows[1][2]
                                    + mRows[0][0] * mRows[3][1] * mRows[1][2] + mRows[1][0] * mRows[0][1] * mRows[3][2] - mRows[0][0] * mRows[1][1] * mRows[3][2];

                  ret.mRows[0][3] = + mRows[2][1] * mRows[1][2] * mRows[0][3] - mRows[1][1] * mRows[2][2] * mRows[0][3] - mRows[2][1] * mRows[0][2] * mRows[1][3]
                                    + mRows[0][1] * mRows[2][2] * mRows[1][3] + mRows[1][1] * mRows[0][2] * mRows[2][3] - mRows[0][1] * mRows[1][2] * mRows[2][3];

                  ret.mRows[1][3] = + mRows[1][0] * mRows[2][2] * mRows[0][3] - mRows[2][0] * mRows[1][2] * mRows[0][3] + mRows[2][0] * mRows[0][2] * mRows[1][3]
                                    - mRows[0][0] * mRows[2][2] * mRows[1][3] - mRows[1][0] * mRows[0][2] * mRows[2][3] + mRows[0][0] * mRows[1][2] * mRows[2][3];

                  ret.mRows[2][3] = + mRows[2][0] * mRows[1][1] * mRows[0][3] - mRows[1][0] * mRows[2][1] * mRows[0][3] - mRows[2][0] * mRows[0][1] * mRows[1][3]
                                    + mRows[0][0] * mRows[2][1] * mRows[1][3] + mRows[1][0] * mRows[0][1] * mRows[2][3] - mRows[0][0] * mRows[1][1] * mRows[2][3];

                  ret.mRows[3][3] = + mRows[1][0] * mRows[2][1] * mRows[0][2] - mRows[2][0] * mRows[1][1] * mRows[0][2] + mRows[2][0] * mRows[0][1] * mRows[1][2]
                                    - mRows[0][0] * mRows[2][1] * mRows[1][2] - mRows[1][0] * mRows[0][1] * mRows[2][2] + mRows[0][0] * mRows[1][1] * mRows[2][2];

        return ret / determinant;
    }

    /**
    * Transpose matrix.
    */
    SIMD_INLINE rpMatrix4x4<T> getTranspose() const
    {
        rpMatrix4x4<T> ret;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                ret.mRows[i][j] = mRows[j][i];
            }
        }
        return ret;
    }


    /// Return the matrix with absolute values
    SIMD_INLINE rpMatrix4x4<T> getAbsoluteMatrix() const
    {
        return rpMatrix4x4<T>(Abs(mRows[0][0]), Abs(mRows[0][1]), Abs(mRows[0][2]), Abs(mRows[0][3]),
                              Abs(mRows[1][0]), Abs(mRows[1][1]), Abs(mRows[1][2]), Abs(mRows[1][3]),
                              Abs(mRows[2][0]), Abs(mRows[2][1]), Abs(mRows[2][2]), Abs(mRows[2][3]),
                              Abs(mRows[3][0]), Abs(mRows[3][1]), Abs(mRows[3][2]), Abs(mRows[3][3]));
    }


    /// Return the trace of the matrix
    SIMD_INLINE T getTrace() const
    {
        // Compute and return the trace
        return (mRows[0][0] + mRows[1][1] + mRows[2][2] + mRows[3][3]);
    }



    /**
     * Linear interpolation of two matrices
     * @param fact Factor of interpolation. For translation from positon
     * of this matrix (lhs) to matrix rhs, values of factor goes from 0.0 to 1.0.
     * @param rhs Second Matrix for interpolation
     * @note However values of fact parameter are reasonable only in interval
     * [0.0 , 1.0], you can pass also values outside of this interval and you
     * can get result (extrapolation?)
     */
    SIMD_INLINE rpMatrix4x4<T> lerp(T fact, const rpMatrix4x4<T>& rhs) const
    {
        rpMatrix4x4<T> ret = (*this) + (rhs - (*this)) * fact;
        return ret;
    }




    //--------------------------------------------------------------------//


    /// Overloaded operator for addition with assignment
    SIMD_INLINE rpMatrix4x4<T>& operator+=(const rpMatrix4x4<T>& matrix)
    {
        mRows[0] += matrix.mRows[0];
        mRows[1] += matrix.mRows[1];
        mRows[2] += matrix.mRows[2];
        mRows[3] += matrix.mRows[3];

        return *this;
    }

    /// Overloaded operator for substraction with assignment
    SIMD_INLINE rpMatrix4x4<T>& operator-=(const rpMatrix4x4<T>& matrix)
    {
        mRows[0] -= matrix.mRows[0];
        mRows[1] -= matrix.mRows[1];
        mRows[2] -= matrix.mRows[2];
        mRows[3] -= matrix.mRows[3];

        return *this;
    }

    /// Overloaded operator for multiplication with a number with assignment
    SIMD_INLINE rpMatrix4x4<T>& operator*=(T nb)
    {

        mRows[0] *= nb;
        mRows[1] *= nb;
        mRows[2] *= nb;
        mRows[3] *= nb;

        return *this;
    }

    /// Overloaded operator for invert multiplication with a number with assignment
    SIMD_INLINE rpMatrix4x4<T> &operator/=(T nb)
    {
        mRows[0] /= nb;
        mRows[1] /= nb;
        mRows[2] /= nb;
        mRows[3] /= nb;

        return *this;
    }





    ///---------------------------[ Pulgins ] -----------------------------------///



    /// Return a skew-symmetric matrix using a given vector that can be used
    /// to compute cross product with another vector using matrix multiplication
    static SIMD_INLINE rpMatrix4x4<T> computeSkewSymmetricMatrixForCrossProduct(const rpVector4D<T>& vector)
    {
        return rpMatrix4x4<T>(0       , -vector.z,  vector.y,  vector.t,
                             vector.z , 0        , -vector.x,  vector.y,
                            -vector.y , vector.x , 0        , -vector.x,
                            -vector.t ,-vector.y , vector.x , 0);
    }




    /// Return a symmetric matrix using a given vector that can be used
    /// to compute dot product with another vector using matrix multiplication
     static SIMD_INLINE rpMatrix4x4<T> computeSymmetricMatrix(const rpVector4D<T>& vector)
     {
         return rpMatrix4x4<T>(0       ,  vector.z,  vector.y,  vector.t,
                              vector.z , 0        ,  vector.x,  vector.y,
                              vector.y , vector.x , 0        ,  vector.x,
                              vector.t , vector.y , vector.x , 0);
     }



     /*****************************************************
      *  Help info to web site:  https://arxiv.org/pdf/1103.0156.pdf
      *****************************************************/

     /**
      * class TLorentzRotation
             \ingroup Physics
         The TLorentzRotation class describes Lorentz transformations including
         Lorentz boosts and rotations (see TRotation)
         ~~~
                     | xx  xy  xz  xt |
                     |                |
                     | yx  yy  yz  yt |
            lambda = |                |
                     | zx  zy  zz  zt |
                     |                |
                     | tx  ty  tz  tt |
         ~~~
         ### Declaration
         By default it is initialized to the identity matrix, but it may also be
         intialized by an other TLorentzRotation,
         by a pure TRotation or by a boost:
          TLorentzRotation l; // l is
         initialized as identity
          TLorentzRotation m(l); // m = l
          TRotation r;
          TLorentzRotation lr(r);
          TLorentzRotation lb1(bx,by,bz);
          TVector3 b;
          TLorentzRotation lb2(b);
         The Matrix for a Lorentz boosts is:
         ~~~
          | 1+gamma'*bx*bx  gamma'*bx*by   gamma'*bx*bz  gamma*bx |
          |  gamma'*by*bx  1+gamma'*by*by  gamma'*by*bz  gamma*by |
          |  gamma'*bz*bx   gamma'*bz*by  1+gamma'*bz*bz gamma*bz |
          |    gamma*bx       gamma*by       gamma*bz     gamma   |
         ~~~
         with the boost vector b=(bx,by,bz) and gamma=1/Sqrt(1-beta*beta)
         and gamma'=(gamma-1)/beta*beta.
         ### Access to the matrix components/Comparisons
         Access to the matrix components is possible through the member functions
         XX(), XY() .. TT(),
         through the operator (int,int):
      */

     static SIMD_INLINE rpMatrix4x4<T> createLorentzBoost( const rpVector3D<T> &vel )
     {


        static rpMatrix4x4<T> M;


        const rpVector3D<T> n = vel.getUnit();
        const T             v = vel.length();

        const T c = LIGHT_MAX_VELOCITY_C;

         //boost this Lorentz vector

          T gamma = 1.0 * Sqrt( 1.0 - (v*v) / (c*c) );

        // T bgamma = gamma * gamma / (1.0 + gamma);
          T bgamma = (gamma - 1.0);


        M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
        M[1][0] =     ((bgamma)*((n.y * n.x)));
        M[2][0] =     ((bgamma)*((n.z * n.x)));
        M[3][0] = (v*n.x*gamma);

        M[0][1] =     ((bgamma)*((n.x * n.y)));
        M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
        M[2][1] =     ((bgamma)*((n.z * n.y)));
        M[3][1] = (v*n.y*gamma);

        M[0][2] =      ((bgamma)*((n.x * n.z)));
        M[1][2] =      ((bgamma)*((n.y * n.z)));
        M[2][2] =  1.0+((bgamma)*((n.z * n.z)));
        M[3][2] = (v*n.z*gamma);

        M[0][3] =  v*n.x*gamma/(c*c);
        M[1][3] =  v*n.y*gamma/(c*c);
        M[2][3] =  v*n.z*gamma/(c*c);
        M[3][3] =  gamma;


         return M;
    }



     static SIMD_INLINE rpMatrix4x4<T> createLorentzBoost( const rpVector3D<T> vel , T gamma )
     {

    	 static rpMatrix4x4<T> M;

    	 const rpVector3D<T> n = vel.getUnit();
    	 const T             v = vel.length();

    	 // T bgamma = gamma * gamma / (1.0 + gamma);
    	 T bgamma = (gamma - 1.0);

    	 const T c = LIGHT_MAX_VELOCITY_C;

    	 M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
    	 M[1][0] =     ((bgamma)*((n.y * n.x)));
    	 M[2][0] =     ((bgamma)*((n.z * n.x)));
    	 M[3][0] = (v*n.x*gamma);

    	 M[0][1] =     ((bgamma)*((n.x * n.y)));
    	 M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
    	 M[2][1] =     ((bgamma)*((n.z * n.y)));
    	 M[3][1] = (v*n.y*gamma);

    	 M[0][2] =      ((bgamma)*((n.x * n.z)));
    	 M[1][2] =      ((bgamma)*((n.y * n.z)));
    	 M[2][2] =  1.0+((bgamma)*((n.z * n.z)));
    	 M[3][2] = (v*n.z*gamma);

    	 M[0][3] =  v*n.x*gamma/(c*c);
    	 M[1][3] =  v*n.y*gamma/(c*c);
    	 M[2][3] =  v*n.z*gamma/(c*c);
    	 M[3][3] =  gamma;


    	 return M;
     }



       /// Creates translation matrix
       /**
        * Creates translation matrix.
        * @param x X-direction translation
        * @param y Y-direction translation
        * @param z Z-direction translation
        * @param w for W-coordinate translation (implicitly set to 1)
        */
       static SIMD_INLINE rpMatrix4x4<T> createTranslation(T x, T y, T z, T w = 1)
       {
           rpMatrix4x4 ret;
           ret.mRows[3][0] = x;
           ret.mRows[3][1] = y;
           ret.mRows[3][2] = z;
           ret.mRows[3][3] = w;

           return ret;
       }


       /// Creates translation matrix
       /**
        * Creates translation matrix.
        * @param x X-direction translation
        * @param y Y-direction translation
        * @param z Z-direction translation
        * @param w for W-coordinate translation (implicitly set to 1)
        */
       static SIMD_INLINE rpMatrix4x4<T> createTranslation(const rpVector3D<T> & v, T w = 1)
       {
           rpMatrix4x4 ret;
           ret.mRows[3][0] = v.x;
           ret.mRows[3][1] = v.y;
           ret.mRows[3][2] = v.z;
           ret.mRows[3][3] = w;

           return ret;
       }


       /**
        * Creates rotation matrix by rotation around axis.
        * @param xDeg Angle (in degrees) of rotation around axis X.
        * @param yDeg Angle (in degrees) of rotation around axis Y.
        * @param zDeg Angle (in degrees) of rotation around axis Z.
        */
       static SIMD_INLINE rpMatrix4x4<T> createRotationAroundAxis(T xDeg, T yDeg, T zDeg)
       {
           T xRads(DEG2RAD(xDeg));
           T yRads(DEG2RAD(yDeg));
           T zRads(DEG2RAD(zDeg));

           rpMatrix4x4<T> ma, mb, mc;
           float ac = Cos(xRads);
           float as = Sin(xRads);
           float bc = Cos(yRads);
           float bs = Sin(yRads);
           float cc = Cos(zRads);
           float cs = Sin(zRads);

           ma.mRows[1][1] = ac;
           ma.mRows[2][1] = as;
           ma.mRows[1][2] = -as;
           ma.mRows[2][2] = ac;

           mb.mRows[0][0] = bc;
           mb.mRows[2][0] = -bs;
           mb.mRows[0][2] = bs;
           mb.mRows[2][2] = bc;

           mc.mRows[0][0] = cc;
           mc.mRows[1][0] = cs;
           mc.mRows[0][1] = -cs;
           mc.mRows[1][1] = cc;

           rpMatrix4x4<T> ret = ma * mb * mc;

           return ret;
       }


       // Return a 4x4 rotation of axis to matrix
       static SIMD_INLINE rpMatrix4x4<T> createRotationAxis(const rpVector3D<T>& axis, T angle)
       {

           //angle = angle / 180.0f * (float)M_PI;

           T cosA = Cos(angle);
           T sinA = Sin(angle);
           rpMatrix4x4<T> rotationMatrix;

           rotationMatrix.mRows[0][0] = cosA + (1-cosA) * axis.x * axis.x;
           rotationMatrix.mRows[1][0] = (1-cosA) * axis.x * axis.y - axis.z * sinA;
           rotationMatrix.mRows[2][0] = (1-cosA) * axis.x * axis.z + axis.y * sinA;
           rotationMatrix.mRows[3][0] = 0.f;

           rotationMatrix.mRows[0][1] = (1-cosA) * axis.x * axis.y + axis.z * sinA;
           rotationMatrix.mRows[1][1] = cosA + (1-cosA) * axis.y * axis.y;
           rotationMatrix.mRows[2][1] = (1-cosA) * axis.y * axis.z - axis.x * sinA;
           rotationMatrix.mRows[3][1] = 0.f;

           rotationMatrix.mRows[0][2] = (1-cosA) * axis.x * axis.z - axis.y * sinA;
           rotationMatrix.mRows[1][2] = (1-cosA) * axis.y * axis.z + axis.x * sinA;
           rotationMatrix.mRows[2][2] = cosA + (1-cosA) * axis.z * axis.z;
           rotationMatrix.mRows[3][2] = 0.f;

           rotationMatrix.mRows[0][3] = 0.f;
           rotationMatrix.mRows[1][3] = 0.f;
           rotationMatrix.mRows[2][3] = 0.f;
           rotationMatrix.mRows[3][3] = 1.f;

           return rotationMatrix;
       }


        // Return a 4x4 rotation of quaternion to matrix
       static  SIMD_INLINE rpMatrix4x4<T>& createRotation(const rpQuaternion<T>& Quat)
       {
           static rpMatrix4x4<T> M;

           T D1, D2, D3, D4, D5, D6, D7, D8, D9; //Dummy variables to hold precalcs

           D1 = (Quat.x * Quat.x) * 2.0f;
           D2 = (Quat.y * Quat.y) * 2.0f;
           D3 = (Quat.z * Quat.z) * 2.0f;

           T RTimesTwo = Quat.w * 2.0f;
           D4 = Quat.x * RTimesTwo;
           D5 = Quat.y * RTimesTwo;
           D6 = Quat.z * RTimesTwo;

           D7 = (Quat.x * Quat.y) * 2.0f;
           D8 = (Quat.x * Quat.z) * 2.0f;
           D9 = (Quat.y * Quat.z) * 2.0f;

           M.mRows[0][0] = 1.0f - D2 - D3;
           M.mRows[1][0] = D7 - D6;
           M.mRows[2][0] = D8 + D5;
           M.mRows[3][0] = 0;

           M.mRows[0][1] = D7 + D6;
           M.mRows[1][1] = 1.0f - D1 - D3;
           M.mRows[2][1] = D9 - D4;
           M.mRows[3][1] = 0;

           M.mRows[0][2] = D8 - D5;
           M.mRows[1][2] = D9 + D4;
           M.mRows[2][2] = 1.0f - D1 - D2;
           M.mRows[3][2] = 0;

           M.mRows[0][3] = 0;
           M.mRows[1][3] = 0;
           M.mRows[2][3] = 0;
           M.mRows[3][3] = 1;

           return M;
       }


       /**
        * Creates new view matrix to look from specified position @a eyePos to specified position @a centerPos
        * @param eyePos A position of camera
        * @param centerPos A position where camera looks-at
        * @param upDir Direction of up vector
        * @return Resulting view matrix that looks from and at specific position.
        */
       static SIMD_INLINE rpMatrix4x4<T> createLookAt(const rpVector3D<T>& eyePos, const rpVector3D<T>& centerPos, const rpVector3D<T>& upDir)
       {
           rpVector3D<T> forward, side, up;
           rpMatrix4x4<T> m;

           forward = centerPos - eyePos;
           up = upDir;

           forward.normalize();

           // Side = forward x up
           side = forward.cross(up);
           side.normalize();

           // Recompute up as: up = side x forward
           up = side.cross(forward);

           m.mRows[0][0] = side.x;
           m.mRows[1][0] = side.y;
           m.mRows[2][0] = side.z;

           m.mRows[0][1] = up.x;
           m.mRows[1][1] = up.y;
           m.mRows[2][1] = up.z;

           m.mRows[0][2] = -forward.x;
           m.mRows[1][2] = -forward.y;
           m.mRows[2][2] = -forward.z;

           m = m * rpMatrix4x4<T>::createTranslation(-eyePos.x, -eyePos.y, -eyePos.z);
           return m;
       }



       /*!
           Multiplies this matrix by another that applies a perspective
           projection. The vertical field of view will be \a verticalAngle degrees
           within a window with a given \a aspectRatio that determines the horizontal
           field of view.
           The projection will have the specified \a nearPlane and \a farPlane clipping
           planes which are the distances from the viewer to the corresponding planes.
           \sa ortho(), frustum()
       */
       static SIMD_INLINE rpMatrix4x4<T> createPerspective(T verticalAngle, T aspectRatio, T nearPlane, T farPlane)
       {
           // Bail out if the projection volume is zero-sized.
           //if (nearPlane == farPlane || aspectRatio == 0.0f) *this->identity();

           // Construct the projection.
           rpMatrix4x4 m;
           T radians = DegreesToRadians(verticalAngle / 2.0f);
           T sine = Sin(radians);

           if (sine == 0.0f) m.setToIdentity();

           T cotan = Cos(radians) / sine;
           T clip = farPlane - nearPlane;

           m.mRows[0][0] = cotan / aspectRatio;
           m.mRows[1][0] = 0.0f;
           m.mRows[2][0] = 0.0f;
           m.mRows[3][0] = 0.0f;

           m.mRows[0][1] = 0.0f;
           m.mRows[1][1] = cotan;
           m.mRows[2][1] = 0.0f;
           m.mRows[3][1] = 0.0f;

           m.mRows[0][2] = 0.0f;
           m.mRows[1][2] = 0.0f;
           m.mRows[2][2] = -(nearPlane + farPlane) / clip;
           m.mRows[3][2] = -(2.0f * nearPlane * farPlane) / clip;

           m.mRows[0][3] = 0.0f;
           m.mRows[1][3] = 0.0f;
           m.mRows[2][3] = -1.0f;
           m.mRows[3][3] = 0.0f;

           return m;
       }


       /*!
           Multiplies this matrix by another that performs the scale and bias
           transformation used by OpenGL to transform from normalized device
           coordinates (NDC) to viewport (window) coordinates. That is it maps
           points from the cube ranging over [-1, 1] in each dimension to the
           viewport with it's near-lower-left corner at (\a left, \a bottom, \a nearPlane)
           and with size (\a width, \a height, \a farPlane - \a nearPlane).
           This matches the transform used by the fixed function OpenGL viewport
           transform controlled by the functions glViewport() and glDepthRange().
        */
       static SIMD_INLINE rpMatrix4x4<T> createViewport(T left, T bottom, T width, T height, T nearPlane, T farPlane)
       {
           const T w2 = width / 2.0f;
           const T h2 = height / 2.0f;

           rpMatrix4x4 m;
           m.mRows[0][0] = w2;
           m.mRows[1][0] = 0.0f;
           m.mRows[2][0] = 0.0f;
           m.mRows[3][0] = left + w2;
           m.mRows[0][1] = 0.0f;
           m.mRows[1][1] = h2;
           m.mRows[2][1] = 0.0f;
           m.mRows[3][1] = bottom + h2;
           m.mRows[0][2] = 0.0f;
           m.mRows[1][2] = 0.0f;
           m.mRows[2][2] = (farPlane - nearPlane) / 2.0f;
           m.mRows[3][2] = (nearPlane + farPlane) / 2.0f;
           m.mRows[0][3] = 0.0f;
           m.mRows[1][3] = 0.0f;
           m.mRows[2][3] = 0.0f;
           m.mRows[3][3] = 1.0f;

           return m;
       }

       /**
        * Creates OpenGL compatible perspective projection according specified frustum parameters.
        *
        * @param left Specify the coordinate for the left vertical clipping plane,
        * @param right Specify the coordinate for the right vertical clipping plane.
        * @param bottom Specify the coordinate for the bottom horizontal clipping plane,
        * @param top Specify the coordinate for the top horizontal clipping plane.
        * @param zNear Specify the distance to the near clipping plane.  Distance must be positive.
        * @param zFar Specify the distance to the far depth clipping plane.  Distance must be positive.
        *
        * @return Projection matrix for specified frustum.
        */
       static SIMD_INLINE rpMatrix4x4<T> createFrustum(T left, T right, T bottom, T top, T zNear, T zFar)
       {
           rpMatrix4x4<T> ret;

           const T invWidth = 1.0 / (right - left);
           const T invHeight = 1.0 / (top - bottom);
           const T invDepth = 1.0 / (zFar - zNear);

           const T twoZNear = 2 * zNear;

           ret.mRows[0,0] = twoZNear * invWidth;
           ret.mRows[1,1] = twoZNear * invHeight;

           ret.mRows[2,0] = (right + left) * invWidth;
           ret.mRows[2,1] = (top + bottom) * invHeight;
           ret.mRows[2,2] = - (zFar + zNear) * invDepth;
           ret.mRows[2,3] = -1;

           ret.mRows[3,2] = - twoZNear * zFar * invDepth;

           return ret;
       }

       /**
        * Creates OpenGL compatible orthographic projection matrix.
        * @param left Specify the coordinate for the left vertical clipping plane,
        * @param right Specify the coordinate for the right vertical clipping plane.
        * @param bottom Specify the coordinate for the bottom horizontal clipping plane,
        * @param top Specify the coordinate for the top horizontal clipping plane.
        * @param zNear Specify the distance to the nearer depth clipping plane.
        *       This value is negative if the plane is to be behind the viewer,
        * @param zFar Specify the distance to the farther depth clipping plane.
        *       This value is negative if the plane is to be behind the viewer.
        * @return Othrographic projection matrix.
        */
       static SIMD_INLINE rpMatrix4x4<T> createOrtho(T left, T right, T bottom, T top, T zNear, T zFar)
       {
           const T invWidth = 1.0 / (right  - left);
           const T invHeight = 1.0 / (top - bottom);
           const T invDepth = 1.0 / (zFar - zNear);

           rpMatrix4x4<T> ret;

           ret.mRows[0,0] =  2 * invWidth;
           ret.mRows[1,1] =  2 * invHeight;
           ret.mRows[2,2] = -2 * invDepth;

           ret.mRows[3,0] = -(right + left) * invWidth;
           ret.mRows[3,1] = -(top + bottom) * invHeight;
           ret.mRows[3,2] = -(zFar + zNear) * invDepth;

           return ret;
       }





    //----------[ output operator ]----------------------------
    /**
    * Output to stream operator
    * @param lhs Left hand side argument of operator (commonly ostream instance).
    * @param rhs Right hand side argument of operator.
    * @return Left hand side argument - the ostream object passed to operator.
    */
    friend std::ostream& operator <<(std::ostream& lhs, const rpMatrix4x4<T>& rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            lhs << "|\t";
            for (int j = 0; j < 4; j++)
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
    static const rpMatrix4x4<T> IDENTITY;

    /**
     * The additive identity matrix.
     */
    static const rpMatrix4x4<T> ZERO;

};

template<class T> const rpMatrix4x4<T> rpMatrix4x4<T>::IDENTITY = rpMatrix4x4<T>(1.0, 0.0, 0.0, 0.0,
																				  0.0, 1.0, 0.0, 0.0,
																				  0.0, 0.0, 1.0, 0.0,
																				  0.0, 0.0, 0.0, 1.0);

template<class T> const rpMatrix4x4<T> rpMatrix4x4<T>::ZERO = rpMatrix4x4<T>(0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0);





template<class T>  SIMD_INLINE rpMatrix4x4<T> operator ^ (const rpVector4D<T> lhs , const rpVector4D<T> rhs )
{
	return rpMatrix4x4<T>( lhs.x * rhs.x , lhs.x * rhs.y, lhs.x * rhs.z, lhs.x * rhs.w ,
			               lhs.y * rhs.x , lhs.y * rhs.y, lhs.y * rhs.z, lhs.y * rhs.w ,
						   lhs.z * rhs.x , lhs.z * rhs.y, lhs.z * rhs.z, lhs.z * rhs.w ,
						   lhs.w * rhs.x , lhs.w * rhs.y, lhs.w * rhs.z, lhs.w * rhs.w);
}

template<class T>  SIMD_INLINE rpMatrix4x4<T> operator ^ (const rpLorentzVector<T> lhs , const rpLorentzVector<T> rhs )
{
	return rpMatrix4x4<T>( lhs.x * rhs.x , lhs.x * rhs.y, lhs.x * rhs.z, lhs.x * rhs.t ,
			               lhs.y * rhs.x , lhs.y * rhs.y, lhs.y * rhs.z, lhs.y * rhs.t ,
						   lhs.z * rhs.x , lhs.z * rhs.y, lhs.z * rhs.z, lhs.z * rhs.t ,
						   lhs.t * rhs.x , lhs.t * rhs.y, lhs.t * rhs.z, lhs.t * rhs.t);
}

/// Matrix 4x4 of floats
typedef rpMatrix4x4<float> rpMatrix4x4f;
/// Matrix 4x4 of doubles
typedef rpMatrix4x4<double> rpMatrix4x4d;
/// Matrix 4x4 of int
typedef rpMatrix4x4<int> rpMatrix4x4i;



} /* namespace  */

#endif /* SRC_PHYSICS_ENGINE_LINEARMATHS_RPMATRIX4X4_H_ */
