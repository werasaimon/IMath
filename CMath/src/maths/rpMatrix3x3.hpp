
#ifndef SOURCE_REAL_PHYSICS_LINEARMATHS_RPMATRIX3X3_H_
#define SOURCE_REAL_PHYSICS_LINEARMATHS_RPMATRIX3X3_H_


// Libraries
#include <cassert>
#include "rpVector4D.hpp"



namespace CMath
{



template<class T> class rpQuaternion;

//**
// * Class for matrix 3x3.
// * @note Data stored in this matrix are in column major order. This arrangement suits OpenGL.
// * If you're using row major matrix, consider using fromRowMajorArray as way for construction
// * Matrix3<T> instance.
// */
template<class T> class rpMatrix3x3
{

   private:

   //-------------------- Attributes --------------------//

    union
    {
        T             mData[9];
        rpVector3D<T> mRows[3];
    };


    public:


     //--------------------------[ constructors ]-------------------------------

      // Constructor of the class Matrix3x3
     SIMD_INLINE rpMatrix3x3()
     {
          // Initialize all values in the matrix to identity
          setAllValues(1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0);
     }

      // Constructor
     SIMD_INLINE rpMatrix3x3(T value)
      {
          setAllValues(value, value, value,
                       value, value, value,
                       value, value, value);
      }

      /**
      * Copy casting constructor.
      * @param src Data source for new created instance of rpMatrix3x3
      */
      SIMD_INLINE rpMatrix3x3(T a1, T a2, T a3,
                              T b1, T b2, T b3,
                              T c1, T c2, T c3)
      {
                  setAllValues(a1, a2, a3,
                               b1, b2, b3,
                               c1, c2, c3);
      }


      // Constructor with arguments
      SIMD_INLINE rpMatrix3x3( T data[3][3] )
      {
          setAllValues(data[0][0], data[0][1], data[0][2],
                       data[1][0], data[1][1], data[1][2],
                       data[2][0], data[2][1], data[2][2]);
      }


      // Copy-constructor
      SIMD_INLINE rpMatrix3x3(const rpMatrix3x3<T>& matrix)
      {
          setAllValues(matrix.mRows[0][0], matrix.mRows[0][1], matrix.mRows[0][2],
                       matrix.mRows[1][0], matrix.mRows[1][1], matrix.mRows[1][2],
                       matrix.mRows[2][0], matrix.mRows[2][1], matrix.mRows[2][2]);
      }


      /**
      * Copy matrix values from array (these data must be in column
      * major order!)
      */
      SIMD_INLINE rpMatrix3x3(const T * dt)
      {
          std::memcpy(mData, dt, sizeof(T) * 9);
      }


      /**
      * Copy casting constructor.
      * @param src Data source for new created instance of rpMatrix3x3
      */
      template<class FromT>
      SIMD_INLINE rpMatrix3x3(const rpMatrix3x3<FromT>& src)
      {
          for (int i = 0; i < 9; i++)
          {
              mData[i] = static_cast<T>(src.mData[i]);
          }
      }


      //---------------------- Methods ---------------------//

      /// Set the matrix to the identity matrix
      SIMD_INLINE void setToIdentity()
      {
    	  // Initialize all values in the matrix to identity
    	  setAllValues(1.0, 0.0, 0.0,
    			       0.0, 1.0, 0.0,
				       0.0, 0.0, 1.0);
      }


      /// Set the matrix to zero
      SIMD_INLINE void setToZero()
      {
          mRows[0].setToZero();
          mRows[1].setToZero();
          mRows[2].setToZero();
      }


      /// Set all the values in the matrix
      SIMD_INLINE void setAllValues(T a1, T a2, T a3,
                                    T b1, T b2, T b3,
                                    T c1, T c2, T c3)
      {
          mRows[0][0] = a1; mRows[0][1] = a2; mRows[0][2] = a3;
          mRows[1][0] = b1; mRows[1][1] = b2; mRows[1][2] = b3;
          mRows[2][0] = c1; mRows[2][1] = c2; mRows[2][2] = c3;
      }



      //---------------------[ assignment operations ]---------------------------------

      /**
      * Copy operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE rpMatrix3x3<T>& operator=(const rpMatrix3x3<T>& rhs)
      {
          std::memcpy(mData, rhs.mData, sizeof(T) * 9);
          return *this;
      }

      /**
      * Copy casting operator
      * @param rhs Right hand side argument of binary operator.
      */
      template<class FromT>
      SIMD_INLINE rpMatrix3x3<T>& operator=(const rpMatrix3x3<FromT>& rhs)
      {
          for (int i = 0; i < 9; i++)
          {
              mData[i] = static_cast<T>(rhs.mData[i]);
          }
          return *this;
      }

      /**
      * Copy operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE rpMatrix3x3<T>& operator=(const T* rhs)
      {
          std::memcpy(mData, rhs, sizeof(T) * 9);
          return *this;
      }


      /// Overloaded operator for equality condition
      SIMD_INLINE bool operator == (const rpMatrix3x3<T>& matrix) const
      {
          return (mRows[0] == matrix.mRows[0] &&
                  mRows[1] == matrix.mRows[1] &&
                  mRows[2] == matrix.mRows[2] );
      }

      /// Overloaded operator for the is different condition
      SIMD_INLINE bool operator != (const rpMatrix3x3<T>& matrix) const
      {
            return !(*this == matrix);
      }




      //-------------[ conversion data ]-----------------------------

      /**
       * Conversion to pointer operator
       * @return Pointer to internally stored (in management of class rpMatrix3x3<T>)
       * used for passing rpMatrix3x3<T> values to gl*[fd]v functions.
       */
      SIMD_INLINE operator T*()
      {
          return (T*) &mRows[0][0];
      }

      /**
       * Conversion to pointer operator
       * @return Constant Pointer to internally stored (in management of class rpMatrix3x3<T>)
       * used for passing rpMatrix3x3<T> values to gl*[fd]v functions.
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
      SIMD_INLINE const rpVector3D<T>& operator[](int row) const
      {
           return mRows[row];
      }

      /// Overloaded operator to read/write element of the matrix.
      SIMD_INLINE rpVector3D<T>& operator[](int row)
      {
           return mRows[row];
      }




      /// Return a column
      SIMD_INLINE rpVector3D<T> getColumn(int i) const
      {
          assert(i>= 0 && i<3);
          return rpVector3D<T> (mRows[0][i], mRows[1][i], mRows[2][i]);
      }

      /// Return a row
      SIMD_INLINE rpVector3D<T> getRow(int i) const
      {
          assert(i>= 0 && i<3);
          return mRows[i];
      }




      //--------------------[ matrix with scalar operations ]---------------------
      /**
      * Addition operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE rpMatrix3x3<T> operator+(T rhs) const
      {
          rpMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++) ret.mRows[i] = mRows[i] + rhs;
          return ret;
      }

      /**
      * Subtraction operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE rpMatrix3x3<T> operator-(T rhs) const
      {
          rpMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++) ret.mRows[i] = mRows[i] - rhs;
          return ret;
      }


      //--------------------[ matrix with matrix operations ]---------------------



      /**
       * Addition operator
       * @param rhs Right hand side argument of binary operator.
       */
      SIMD_INLINE rpMatrix3x3<T> operator+(const rpMatrix3x3<T>& rhs) const
      {
    	  rpMatrix3x3<T> ret;
    	  for (int i = 0; i < 3; i++)   ret.mRows[i] = mRows[i] + rhs.mRows[i];
    	  return ret;
      }

      /**
       * Subtraction operator
       * @param rhs Right hand side argument of binary operator.
       */
      SIMD_INLINE rpMatrix3x3<T> operator-(const rpMatrix3x3<T>& rhs) const
      {
    	  rpMatrix3x3<T> ret;
    	  for (int i = 0; i < 3; i++)   ret.mRows[i] = mRows[i] - rhs.mRows[i];
    	  return ret;
      }

      //---------------------------- Friend method -------------------------//

      /// Overloaded operator for the negative of the matrix
      friend SIMD_INLINE rpMatrix3x3<T>  operator-(const rpMatrix3x3<T> & matrix)
      {
          return rpMatrix3x3<T>(-matrix.mRows[0][0], -matrix.mRows[0][1], -matrix.mRows[0][2],
                                -matrix.mRows[1][0], -matrix.mRows[1][1], -matrix.mRows[1][2],
                                -matrix.mRows[2][0], -matrix.mRows[2][1], -matrix.mRows[2][2]);
      }

      /// Overloaded operator for multiplication with a number
      friend SIMD_INLINE rpMatrix3x3<T>  operator*(T nb, const rpMatrix3x3<T>& matrix)
      {
          return rpMatrix3x3<T>(matrix.mRows[0][0] * nb, matrix.mRows[0][1] * nb, matrix.mRows[0][2] * nb,
                                matrix.mRows[1][0] * nb, matrix.mRows[1][1] * nb, matrix.mRows[1][2] * nb,
                                matrix.mRows[2][0] * nb, matrix.mRows[2][1] * nb, matrix.mRows[2][2] * nb);
      }

      /// Overloaded operator for multiplication with a matrix
      friend SIMD_INLINE rpMatrix3x3<T>  operator*(const rpMatrix3x3<T>& matrix, T nb)
      {
          return nb * matrix;
      }

      /// Overloaded operator for inveret multiplication with a number
      friend SIMD_INLINE rpMatrix3x3<T>  operator/(T nb, const rpMatrix3x3<T>& matrix)
      {
          return rpMatrix3x3<T>(matrix.mRows[0][0] / nb, matrix.mRows[0][1] / nb, matrix.mRows[0][2] / nb,
                                matrix.mRows[1][0] / nb, matrix.mRows[1][1] / nb, matrix.mRows[1][2] / nb,
                                matrix.mRows[2][0] / nb, matrix.mRows[2][1] / nb, matrix.mRows[2][2] / nb);
      }

      /// Overloaded operator for inveret multiplication with a matrix
      friend SIMD_INLINE rpMatrix3x3<T>  operator/(const rpMatrix3x3<T>& matrix, T nb)
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
      friend SIMD_INLINE rpVector3D<T> operator*(const rpVector3D<T>& rhs , const rpMatrix3x3<T>& lhs)
      {
    	  return lhs * rhs;
      }


       //--------------------[ multiply operators ]--------------------------------

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
      SIMD_INLINE rpVector3D<T> operator*(const rpVector3D<T>& rhs) const
      {
    	  T fX = rhs.x;
    	  T fY = rhs.y;
    	  T fZ = rhs.z;

    	  rpVector3D<T> Point;
    	  Point.x = ( fX * mRows[0][0] + fY * mRows[0][1] + fZ * mRows[0][2]);
    	  Point.y = ( fX * mRows[1][0] + fY * mRows[1][1] + fZ * mRows[1][2]);
    	  Point.z = ( fX * mRows[2][0] + fY * mRows[2][1] + fZ * mRows[2][2]);

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
      SIMD_INLINE rpMatrix3x3<T> operator*(rpMatrix3x3<T> rhs) const
      {
    	  rpMatrix3x3<T> w;
    	  for(int i = 0; i < 3; i++)
    	  {
    		  for (int j = 0; j < 3; j++)
    		  {
    			  T n = 0;
    			  for (int k = 0; k < 3; k++)
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
      * Return the determinant minor of the matrix
      */
      SIMD_INLINE T getDeterminantOfMinor( int  theRowHeightY , int  theColumnWidthX ) const
      {
          int x1 = theColumnWidthX == 0 ? 1 : 0;  /* always either 0 or 1 */
          int x2 = theColumnWidthX == 2 ? 1 : 2;  /* always either 1 or 2 */
          int y1 = theRowHeightY   == 0 ? 1 : 0;  /* always either 0 or 1 */
          int y2 = theRowHeightY   == 2 ? 1 : 2;  /* always either 1 or 2 */

          return ( mRows[y1][x1]  *  mRows[y2][x2] )
              -  ( mRows[y1][x2]  *  mRows[y2][x1] );
      }

      /**
       * Computes determinant of matrix
       * @return Determinant of matrix
       * @note This function does.
       */
      SIMD_INLINE T getDeterminant() const
      {
          return ( mRows[0][0] * getDeterminantOfMinor(0,0) )
              -  ( mRows[1][0] * getDeterminantOfMinor(1,0) )
              +  ( mRows[2][0] * getDeterminantOfMinor(2,0) );
      }

      /**
      * Computes inverse matrix
      * @return Inverse matrix of this matrix.
      * @note This is a little bit time consuming operation
      */
      SIMD_INLINE rpMatrix3x3<T> getInverse() const
      {
          // Compute the determinant of the matrix
            T determinant = getDeterminant();

            // Check if the determinant is equal to zero
             assert(Abs(determinant) > MACHINE_EPSILON);

            T invDeterminant = T(1.0) / determinant;

            rpMatrix3x3<T> tempMatrix((mRows[1][1]*mRows[2][2]-mRows[2][1]*mRows[1][2]),
                                     -(mRows[0][1]*mRows[2][2]-mRows[2][1]*mRows[0][2]),
                                      (mRows[0][1]*mRows[1][2]-mRows[0][2]*mRows[1][1]),
                                     -(mRows[1][0]*mRows[2][2]-mRows[2][0]*mRows[1][2]),
                                      (mRows[0][0]*mRows[2][2]-mRows[2][0]*mRows[0][2]),
                                     -(mRows[0][0]*mRows[1][2]-mRows[1][0]*mRows[0][2]),
                                      (mRows[1][0]*mRows[2][1]-mRows[2][0]*mRows[1][1]),
                                     -(mRows[0][0]*mRows[2][1]-mRows[2][0]*mRows[0][1]),
                                      (mRows[0][0]*mRows[1][1]-mRows[0][1]*mRows[1][0]));

            // Return the inverse matrix
        return (invDeterminant * tempMatrix);
      }


      /**
      * Transpose matrix.
      */
      SIMD_INLINE rpMatrix3x3<T> getTranspose() const
      {
          rpMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++)
          {
              for (int j = 0; j < 3; j++)
              {
                  ret.mRows[i][j] = mRows[j][i];
              }
          }
          return ret;
      }


      /**
      * Return the matrix with absolute values
      */
      SIMD_INLINE rpMatrix3x3<T> getAbsoluteMatrix() const
      {
          return rpMatrix3x3<T>(Abs(mRows[0][0]), Abs(mRows[0][1]), Abs(mRows[0][2]),
                                Abs(mRows[1][0]), Abs(mRows[1][1]), Abs(mRows[1][2]),
                                Abs(mRows[2][0]), Abs(mRows[2][1]), Abs(mRows[2][2]));
      }


      /**
      * Return the trace of the matrix
      */
      SIMD_INLINE T getTrace() const
      {
          // Compute and return the trace
          return (mRows[0][0] + mRows[1][1] + mRows[2][2]);
      }

      /**
      * Return the diagonalize of the matrix
      */
      SIMD_INLINE rpMatrix3x3<T> getDiagonalize( T threshold, int maxSteps )
      {
    	  rpMatrix3x3<T> rot;
    	  rot.setToIdentity();
    	  for (int step = maxSteps; step > 0; step--)
    	  {
    		  // find off-diagonal element [p][q] with largest magnitude
    		  int p = 0;
    		  int q = 1;
    		  int r = 2;
    		  T max = Abs(mRows[0][1]);
    		  T v = Abs(mRows[0][2]);
    		  if (v > max)
    		  {
    			  q = 2;
    			  r = 1;
    			  max = v;
    		  }
    		  v = Abs(mRows[1][2]);
    		  if (v > max)
    		  {
    			  p = 1;
    			  q = 2;
    			  r = 0;
    			  max = v;
    		  }

    		  T t = threshold * (Abs(mRows[0][0]) + Abs(mRows[1][1]) + Abs(mRows[2][2]));
    		  if (max <= t)
    		  {
    			  if (max <= EPSILON * t)
    			  {
    				  break;
    			  }
    			  step = 1;
    		  }

    		  // compute Jacobi rotation J which leads to a zero for element [p][q]
		      T mpq = mRows[p][q];
    		  T theta = (mRows[q][q] - mRows[p][p]) / (2 * mpq);
    		  T theta2 = theta * theta;
    		  T cos;
    		  T sin;
    		  if (theta2 * theta2 < T(10 / EPSILON))
    		  {
    			  t = (theta >= 0) ? 1 / (theta + Sqrt(1 + theta2))
    					  : 1 / (theta - Sqrt(1 + theta2));
    			  cos = 1 / Sqrt(1 + t * t);
    			  sin = cos * t;
    		  }
    		  else
    		  {
    			  // approximation for large theta-value, i.e., a nearly diagonal matrix
				  t = 1 / (theta * (2 + T(0.5) / theta2));
				  cos = 1 - T(0.5) * t * t;
				  sin = cos * t;
    		  }

    		  // apply rotation to matrix (this = J^T * this * J)
    		  mRows[p][q] = mRows[q][p] = 0;
    		  mRows[p][p] -= t * mpq;
    		  mRows[q][q] += t * mpq;
    		  T mrp = mRows[r][p];
    		  T mrq = mRows[r][q];
    		  mRows[r][p] = mRows[p][r] = cos * mrp - sin * mrq;
    		  mRows[r][q] = mRows[q][r] = cos * mrq + sin * mrp;
    		  // apply rotation to rot (rot = rot * J)
    		  for (int i = 0; i < 3; i++)
    		  {
    			  rpVector3D<T>& row = rot[i];
    			  mrp = row[p];
    			  mrq = row[q];
    			  row[p] = cos * mrp - sin * mrq;
    			  row[q] = cos * mrq + sin * mrp;
    		  }

    	  }

    	  return rot;
      }


      //--------------------------------------------------------------------//


      /// Overloaded operator for addition with assignment
      SIMD_INLINE rpMatrix3x3<T>& operator+=(const rpMatrix3x3<T>& matrix)
      {
          mRows[0] += matrix.mRows[0];
          mRows[1] += matrix.mRows[1];
          mRows[2] += matrix.mRows[2];
          return *this;
      }

      /// Overloaded operator for substraction with assignment
      SIMD_INLINE rpMatrix3x3<T>& operator-=(const rpMatrix3x3<T>& matrix)
      {
          mRows[0] -= matrix.mRows[0];
          mRows[1] -= matrix.mRows[1];
          mRows[2] -= matrix.mRows[2];
          return *this;
      }

      /// Overloaded operator for multiplication with a number with assignment
      SIMD_INLINE rpMatrix3x3<T>& operator*=(T nb)
      {
          mRows[0] *= nb;
          mRows[1] *= nb;
          mRows[2] *= nb;
          return *this;
      }

      /// Overloaded operator for invert multiplication with a number with assignment
      SIMD_INLINE rpMatrix3x3<T> &operator/=(T nb)
      {
          mRows[0] /= nb;
          mRows[1] /= nb;
          mRows[2] /= nb;
          return *this;
      }





      ///---------------------------[ Pulgins ] -----------------------------------///



      /// Return a skew-symmetric matrix using a given vector that can be used
      /// to compute cross product with another vector using matrix multiplication
      static SIMD_INLINE rpMatrix3x3<T> computeSkewSymmetricMatrixForCrossProduct(const rpVector3D<T>& vector)
      {
          return rpMatrix3x3<T>(0, -vector.z, vector.y,
                                vector.z, 0, -vector.x,
                               -vector.y, vector.x, 0);
      }


      /// Return a symmetric matrix using a given vector that can be used
      /// to compute dot product with another vector using matrix multiplication
       static SIMD_INLINE rpMatrix3x3<T> computeSymmetricMatrix(const rpVector3D<T>& vector)
       {
           return rpMatrix3x3<T>(0, vector.z, vector.y,
                                 vector.z, 0, vector.x,
                                 vector.y, vector.x, 0);
       }



       /**
        * Returns a scaling matrix that scales by `scaleFactors.x` and
        * 'scaleFactors.y' in the x and y axes respectively.
        *
        * @param scaleFactors Scale factors.
        * @return Scaling matrix.
        */
       static SIMD_INLINE rpMatrix3x3<T> createScaling(const rpVector3D<T>& scaleFactors)
       {
           return rpMatrix3x3<T>(  scaleFactors.x, 0.0f , 0.0f ,
                                   0.0f, scaleFactors.y , 0.0f ,
                                   0.0f , 0.0f , scaleFactors.z);
       }

       /**
        * Returns a scaling matrix that scales by `factor` uniformly.
        *
        * @param scale Uniform scale factor.
        * @return Scaling matrix.
        */
       static SIMD_INLINE rpMatrix3x3<T> createScaling(const T factor)
       {
           return rpMatrix3x3<T>( factor, 0.0f  , 0.0f ,
                                  0.0f  , factor, 0.0f ,
                                  0.0f  , 0.0f  , factor);
       }



      /*****************************************************
       *  Help info to web site:  https://arxiv.org/pdf/1103.0156.pdf
       *****************************************************/
       /// Return lorentz demission distance world
      static SIMD_INLINE rpMatrix3x3<T> createLorentzRotationBoost(  rpVector3D<T> vel )
      {
          const rpVector3D<T> n = vel.getUnit();
          const T             v = vel.length();

          static rpMatrix3x3<T> M;

          const T c = LIGHT_MAX_VELOCITY_C;

          //boost this Lorentz vector
          T gamma = 1.0 * Sqrt( 1.0 - (v*v) / (c*c) );

         // T bgamma = gamma * gamma / (1.0 + gamma);
          T bgamma = (gamma - 1.0);


          M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
          M[1][0] =     ((bgamma)*((n.y * n.x)));
          M[2][0] =     ((bgamma)*((n.z * n.x)));


          M[0][1] =     ((bgamma)*((n.x * n.y)));
          M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
          M[2][1] =     ((bgamma)*((n.z * n.y)));


          M[0][2] =     ((bgamma)*((n.x * n.z)));
          M[1][2] =     ((bgamma)*((n.y * n.z)));
          M[2][2] = 1.0+((bgamma)*((n.z * n.z)));

          return M;
      }


      static SIMD_INLINE rpMatrix3x3<T> createLorentzRotationBoost(  const rpVector3D<T> n ,  T gamma )
      {

              static rpMatrix3x3<T> M;

             // T bgamma = gamma * gamma / (1.0 + gamma);
              T bgamma = (gamma - 1.0);


              M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
              M[1][0] =     ((bgamma)*((n.y * n.x)));
              M[2][0] =     ((bgamma)*((n.z * n.x)));


              M[0][1] =     ((bgamma)*((n.x * n.y)));
              M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
              M[2][1] =     ((bgamma)*((n.z * n.y)));


              M[0][2] =     ((bgamma)*((n.x * n.z)));
              M[1][2] =     ((bgamma)*((n.y * n.z)));
              M[2][2] = 1.0+((bgamma)*((n.z * n.z)));

              return M;
       }


      // Return a 4x4 rotation of axis to matrix
      static SIMD_INLINE rpMatrix3x3<T> createRotationAxis(const rpVector3D<T>& axis, T angle)
      {

          //angle = angle / 180.0f * (float)M_PI;

          T cosA = Cos(angle);
          T sinA = Sin(angle);
          rpMatrix3x3<T> rotationMatrix;

          rotationMatrix.mRows[0][0] = cosA + (1-cosA) * axis.x * axis.x;
          rotationMatrix.mRows[1][0] = (1-cosA) * axis.x * axis.y - axis.z * sinA;
          rotationMatrix.mRows[2][0] = (1-cosA) * axis.x * axis.z + axis.y * sinA;

          rotationMatrix.mRows[0][1] = (1-cosA) * axis.x * axis.y + axis.z * sinA;
          rotationMatrix.mRows[1][1] = cosA + (1-cosA) * axis.y * axis.y;
          rotationMatrix.mRows[2][1] = (1-cosA) * axis.y * axis.z - axis.x * sinA;


          rotationMatrix.mRows[0][2] = (1-cosA) * axis.x * axis.z - axis.y * sinA;
          rotationMatrix.mRows[1][2] = (1-cosA) * axis.y * axis.z + axis.x * sinA;
          rotationMatrix.mRows[2][2] = cosA + (1-cosA) * axis.z * axis.z;


          return rotationMatrix;
      }



      static  SIMD_INLINE rpMatrix3x3<T>& createRotation(const rpQuaternion<T>& Quat)
      {
          static rpMatrix3x3<T> M;

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

          M.mRows[0][1] = D7 + D6;
          M.mRows[1][1] = 1.0f - D1 - D3;
          M.mRows[2][1] = D9 - D4;

          M.mRows[0][2] = D8 - D5;
          M.mRows[1][2] = D9 + D4;
          M.mRows[2][2] = 1.0f - D1 - D2;

          return M;
      }




      //----------[ output operator ]----------------------------
      /**
      * Output to stream operator
      * @param lhs Left hand side argument of operator (commonly ostream instance).
      * @param rhs Right hand side argument of operator.
      * @return Left hand side argument - the ostream object passed to operator.
      */
      friend std::ostream& operator <<(std::ostream& lhs, const rpMatrix3x3<T>& rhs)
      {
          for (int i = 0; i < 3; i++)
          {
              lhs << "|\t";
              for (int j = 0; j < 3; j++)
              {
                  lhs << rhs[i][j] << "\t";
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
    static const rpMatrix3x3<T> IDENTITY;

    /**
     * The additive identity matrix.
     */
    static const rpMatrix3x3<T> ZERO;
};



template<class T> const rpMatrix3x3<T> rpMatrix3x3<T>::IDENTITY = rpMatrix3x3<T>(1.0, 0.0, 0.0,
																				 0.0, 1.0, 0.0,
																				 0.0, 0.0, 1.0);

template<class T> const rpMatrix3x3<T> rpMatrix3x3<T>::ZERO = rpMatrix3x3<T>(0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0);


template<class T>  SIMD_INLINE rpMatrix3x3<T> operator ^ (const rpVector3D<T> lhs , const rpVector3D<T> rhs )
{
	return rpMatrix3x3<T>( lhs.x * rhs.x , lhs.x * rhs.y, lhs.x * rhs.z,
			               lhs.y * rhs.x , lhs.y * rhs.y, lhs.y * rhs.z,
						   lhs.z * rhs.x , lhs.z * rhs.y, lhs.z * rhs.z);
}

/// Matrix 3x3 of floats
typedef rpMatrix3x3<float> rpMatrix3x3f;
/// Matrix 3x3 of doubles
typedef rpMatrix3x3<double> rpMatrix3x3d;
/// Matrix 3x3 of int
typedef rpMatrix3x3<int> rpMatrix3x3i;



} /* namespace */


#endif /* SOURCE_REAL_PHYSICS_LINEARMATHS_RPMATRIX3X3_H_ */
