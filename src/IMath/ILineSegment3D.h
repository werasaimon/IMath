 /********************************************************************************
 *
 * ILineSegment3D.h
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


#ifndef ILINESEGMENT3D_H
#define ILINESEGMENT3D_H

#include "ILine3D.h"

namespace IMath
{

template<class T> class ILineSegment3D;

template<class T>
T DistanceSquared( const ILineSegment3D<T>& seg0, const ILineSegment3D<T>& seg1, T& s_c, T& t_c );

template<class T>
T DistanceSquared( const ILineSegment3D<T>& segment, const ILine3D<T>& line,  T& s_c, T& t_c );

//-------------------------------------------------------------------------------
//-- Classes --------------------------------------------------------------------
//-------------------------------------------------------------------------------
template<class T>
class ILineSegment3D
{
public:
    // constructor/destructor
    ILineSegment3D()
     : mOrigin( 0.0f, 0.0f, 0.0f ),
       mDirection( 1.0f, 0.0f, 0.0f )
    {

    }

    ILineSegment3D( const IVector3D<T>& endpoint0, const IVector3D<T>& endpoint1 )
     : mOrigin( endpoint0 ),
       mDirection( endpoint1-endpoint0 )
    {

    }

    SIMD_INLINE ~ILineSegment3D() {}

    // copy operations
    ILineSegment3D(const ILineSegment3D<T>& other)
     : mOrigin( other.mOrigin ),
       mDirection( other.mDirection )
    {

    }

    // ---------------------------------------------------------------------------
    // Assigment operator
    //-----------------------------------------------------------------------------
    ILineSegment3D& operator=(const ILineSegment3D<T>& other)
    {
        // if same object
         if ( this == &other )
             return *this;

         mOrigin = other.mOrigin;
         mDirection = other.mDirection;

         return *this;
    }


    // accessors
    T& operator()(unsigned int i, unsigned int j);
    T  operator()(unsigned int i, unsigned int j) const;

    SIMD_INLINE IVector3D<T> GetOrigin() const { return mOrigin; }
    SIMD_INLINE IVector3D<T> GetDirection() const { return mDirection; }
    SIMD_INLINE IVector3D<T> GetEndpoint0() const { return mOrigin; }
    SIMD_INLINE IVector3D<T> GetEndpoint1() const { return mOrigin + mDirection; }
    SIMD_INLINE IVector3D<T> GetCenter() const { return mOrigin + T(0.5f)*mDirection; }

    // ---------------------------------------------------------------------------
    // Returns the two endpoints
    //-----------------------------------------------------------------------------
    void Get( IVector3D<T>& end0, IVector3D<T>& end1 ) const
    {
        end0 = mOrigin;
        end1 = mOrigin + mDirection;
    }

    // manipulators
    void Set( const IVector3D<T>& end0, const IVector3D<T>& end1 )
    {
        mOrigin = end0;
        mDirection = end1-end0;
    }


    // ---------------------------------------------------------------------------
    // Returns the distance between two endpoints
    //-----------------------------------------------------------------------------
    T Length() const
    {
         return mDirection.length();
    }

    // ---------------------------------------------------------------------------
    // Returns the squared distance between two endpoints
    //-----------------------------------------------------------------------------
    T LengthSquared() const
    {
        return mDirection.lengthSquare();
    }

    // ---------------------------------------------------------------------------
    // Are two IvLineSegment3's equal?
    //----------------------------------------------------------------------------
    bool operator==( const ILineSegment3D<T>& segment ) const
    {
        return ((segment.mOrigin == mOrigin && segment.mDirection == mDirection) ||
                (segment.mOrigin == mOrigin+mDirection && segment.mDirection == -mDirection));
    }

    // ---------------------------------------------------------------------------
    // Are two IvLineSegment3's not equal?
    //----------------------------------------------------------------------------
    bool operator!=( const ILineSegment3D<T>& segment ) const
    {
        return !((segment.mOrigin == mOrigin && segment.mDirection == mDirection) ||
                 (segment.mOrigin == mOrigin+mDirection && segment.mDirection == -mDirection));
    }



    // ---------------------------------------------------------------------------
    // Transforms segment into new space
    //-----------------------------------------------------------------------------
    ILineSegment3D Transform( T scale, const IQuaternion<T>& quat, const IVector3D<T>& translate ) const
    {
        ILineSegment3D<T> segment;
        IMatrix3x3<T>  rotate = quat.GetRotMatrix();

        segment.mDirection = rotate * mDirection;
        segment.mDirection *= scale;

        segment.mOrigin = rotate * mOrigin;
        segment.mOrigin *= scale;
        segment.mOrigin += translate;

        return segment;
    }

    // ---------------------------------------------------------------------------
    // Transforms segment into new space
    //-----------------------------------------------------------------------------
    ILineSegment3D Transform( T scale, const  IMatrix3x3<T>& rotate, const IVector3D<T>& translate ) const
    {
        ILineSegment3D<T> segment;

        segment.mDirection = rotate * mDirection;
        segment.mDirection *= scale;

        segment.mOrigin = rotate * mOrigin;
        segment.mOrigin *= scale;
        segment.mOrigin += translate;

        return segment;
    }



    // ---------------------------------------------------------------------------
    // Returns the distance squared between two line segments.
    // Based on article and code by Dan Sunday at www.geometryalgorithms.com
    //-----------------------------------------------------------------------------
    static T DistanceSquared( const ILineSegment3D<T>& segment0, const ILineSegment3D<T>& segment1,  T& s_c, T& t_c )
    {
            // compute intermediate parameters
            IVector3D<T> w0 = segment0.mOrigin - segment1.mOrigin;
            T a = segment0.mDirection.Dot( segment0.mDirection );
            T b = segment0.mDirection.Dot( segment1.mDirection );
            T c = segment1.mDirection.Dot( segment1.mDirection );
            T d = segment0.mDirection.Dot( w0 );
            T e = segment1.mDirection.Dot( w0 );

            T denom = a*c - b*b;
            // parameters to compute s_c, t_c
            T sn, sd, tn, td;

            // if denom is zero, try finding closest point on segment1 to origin0
            if ( IsZero(denom) )
            {
                // clamp s_c to 0
                sd = td = c;
                sn = 0.0f;
                tn = e;
            }
            else
            {
                // clamp s_c within [0,1]
                sd = td = denom;
                sn = b*e - c*d;
                tn = a*e - b*d;

                // clamp s_c to 0
                if (sn < 0.0f)
                {
                    sn = 0.0f;
                    tn = e;
                    td = c;
                }
                // clamp s_c to 1
                else if (sn > sd)
                {
                    sn = sd;
                    tn = e + b;
                    td = c;
                }
            }

            // clamp t_c within [0,1]
            // clamp t_c to 0
            if (tn < 0.0f)
            {
                t_c = 0.0f;
                // clamp s_c to 0
                if ( -d < 0.0f )
                {
                    s_c = 0.0f;
                }
                // clamp s_c to 1
                else if ( -d > a )
                {
                    s_c = 1.0f;
                }
                else
                {
                    s_c = -d/a;
                }
            }
            // clamp t_c to 1
            else if (tn > td)
            {
                t_c = 1.0f;
                // clamp s_c to 0
                if ( (-d+b) < 0.0f )
                {
                    s_c = 0.0f;
                }
                // clamp s_c to 1
                else if ( (-d+b) > a )
                {
                    s_c = 1.0f;
                }
                else
                {
                    s_c = (-d+b)/a;
                }
            }
            else
            {
                t_c = tn/td;
                s_c = sn/sd;
            }

            // compute difference vector and distance squared
            IVector3D<T> wc = w0 + s_c*segment0.mDirection - t_c*segment1.mDirection;
            return wc.Dot(wc);
    }

    // ---------------------------------------------------------------------------
    // Returns the distance squared between line segment and line.
    // Based on article and code by Dan Sunday at www.geometryalgorithms.com
    //-----------------------------------------------------------------------------
    static T DistanceSquared( const ILineSegment3D<T>& segment,  const ILine3D<T>& line,  T& s_c, T& t_c )
    {
        // compute intermediate parameters
        IVector3D<T> w0 = segment.mOrigin - line.GetOrigin();
        T a = segment.mDirection.Dot( segment.mDirection );
        T b = segment.mDirection.Dot( line.GetDirection() );
        T c = line.GetDirection().Dot( line.GetDirection() );
        T d = segment.mDirection.Dot( w0 );
        T e = line.GetDirection().Dot( w0 );

        T denom = a*c - b*b;

        // if denom is zero, try finding closest point on segment1 to origin0
        if ( IsZero(denom) )
        {
            s_c = 0.0f;
            t_c = e/c;
            // compute difference vector and distance squared
            IVector3D<T> wc = w0 - t_c*line.GetDirection();
            return wc.Dot(wc);
        }
        else
        {
            // parameters to compute s_c, t_c
            T sn;

            // clamp s_c within [0,1]
            sn = b*e - c*d;

            // clamp s_c to 0
            if (sn < 0.0f)
            {
                s_c = 0.0f;
                t_c = e/c;
            }
            // clamp s_c to 1
            else if (sn > denom)
            {
                s_c = 1.0f;
                t_c = (e+b)/c;
            }
            else
            {
                s_c = sn/denom;
                t_c = (a*e - b*d)/denom;
            }

            // compute difference vector and distance squared
            IVector3D<T> wc = w0 + s_c*segment.mDirection - t_c*line.GetDirection();
            return wc.Dot(wc);
        }
    }

    static T DistanceSquared( const ILine3D<T>& line,  const ILineSegment3D<T>& segment, T& s_c, T& t_c )
    {
        return DistanceSquared( segment, line, t_c, s_c );
    }

    // ---------------------------------------------------------------------------
    // Returns the distance squared between line segment and point.
    //-----------------------------------------------------------------------------
    static T DistanceSquared( const ILineSegment3D<T>& segment,  const IVector3D<T>& point,  T& t_c )
    {
        IVector3D<T> w = point - segment.mOrigin;
        T proj = w.Dot(segment.mDirection);
        // endpoint 0 is closest point
        if ( proj <= 0 )
        {
            t_c = 0.0f;
            return w.Dot(w);
        }
        else
        {
            T vsq = segment.mDirection.Dot(segment.mDirection);
            // endpoint 1 is closest point
            if ( proj >= vsq )
            {
                t_c = 1.0f;
                return w.Dot(w) - 2.0f*proj + vsq;
            }
            // otherwise somewhere else in segment
            else
            {
                t_c = proj/vsq;
                return w.Dot(w) - t_c*proj;
            }
        }
    }



    // ---------------------------------------------------------------------------
    // Returns the closest points between two line segments.
    //-----------------------------------------------------------------------------
    static void ClosestPoints( const ILineSegment3D<T>& segment0,  const ILineSegment3D<T>& segment1 , IVector3D<T>& point0, IVector3D<T>& point1 )
    {
            // compute intermediate parameters
            IVector3D<T> w0 = segment0.mOrigin - segment1.mOrigin;
            T a = segment0.mDirection.Dot( segment0.mDirection );
            T b = segment0.mDirection.Dot( segment1.mDirection );
            T c = segment1.mDirection.Dot( segment1.mDirection );
            T d = segment0.mDirection.Dot( w0 );
            T e = segment1.mDirection.Dot( w0 );

            T denom = a*c - b*b;
            // parameters to compute s_c, t_c
            T s_c, t_c;
            T sn, sd, tn, td;

            // if denom is zero, try finding closest point on segment1 to origin0
            if ( IsZero(denom) )
            {
                // clamp s_c to 0
                sd = td = c;
                sn = 0.0f;
                tn = e;
            }
            else
            {
                // clamp s_c within [0,1]
                sd = td = denom;
                sn = b*e - c*d;
                tn = a*e - b*d;

                // clamp s_c to 0
                if (sn < 0.0f)
                {
                    sn = 0.0f;
                    tn = e;
                    td = c;
                }
                // clamp s_c to 1
                else if (sn > sd)
                {
                    sn = sd;
                    tn = e + b;
                    td = c;
                }
            }

            // clamp t_c within [0,1]
            // clamp t_c to 0
            if (tn < 0.0f)
            {
                t_c = 0.0f;
                // clamp s_c to 0
                if ( -d < 0.0f )
                {
                    s_c = 0.0f;
                }
                // clamp s_c to 1
                else if ( -d > a )
                {
                    s_c = 1.0f;
                }
                else
                {
                    s_c = -d/a;
                }
            }
            // clamp t_c to 1
            else if (tn > td)
            {
                t_c = 1.0f;
                // clamp s_c to 0
                if ( (-d+b) < 0.0f )
                {
                    s_c = 0.0f;
                }
                // clamp s_c to 1
                else if ( (-d+b) > a )
                {
                    s_c = 1.0f;
                }
                else
                {
                    s_c = (-d+b)/a;
                }
            }
            else
            {
                t_c = tn/td;
                s_c = sn/sd;
            }

            // compute closest points
            point0 = segment0.mOrigin + s_c*segment0.mDirection;
            point1 = segment1.mOrigin + t_c*segment1.mDirection;

    }


    // ---------------------------------------------------------------------------
    // Returns the closest points between line segment and line.
    //-----------------------------------------------------------------------------
    static void ClosestPoints( const ILineSegment3D& segment,  const ILine3D<T>& line , IVector3D<T>& point0, IVector3D<T>& point1 )
    {

        // compute intermediate parameters
        IVector3D<T> w0 = segment.mOrigin - line.GetOrigin();
        T a = segment.mDirection.Dot( segment.mDirection );
        T b = segment.mDirection.Dot( line.GetDirection() );
        T c = line.GetDirection().Dot( line.GetDirection() );
        T d = segment.mDirection.Dot( w0 );
        T e = line.GetDirection().Dot( w0 );

        T denom = a*c - b*b;

        // if denom is zero, try finding closest point on line to segment origin
        if ( IsZero(denom) )
        {
            // compute closest points
            point0 = segment.mOrigin;
            point1 = line.GetOrigin() + (e/c)*line.GetDirection();
        }
        else
        {
            // parameters to compute s_c, t_c
            T s_c, t_c;
            T sn;

            // clamp s_c within [0,1]
            sn = b*e - c*d;

            // clamp s_c to 0
            if (sn < 0.0f)
            {
                s_c = 0.0f;
                t_c = e/c;
            }
            // clamp s_c to 1
            else if (sn > denom)
            {
                s_c = 1.0f;
                t_c = (e+b)/c;
            }
            else
            {
                s_c = sn/denom;
                t_c = (a*e - b*d)/denom;
            }

            // compute closest points
            point0 = segment.mOrigin + s_c*segment.mDirection;
            point1 = line.GetOrigin() + t_c*line.GetDirection();
        }

    }

    // ---------------------------------------------------------------------------
    // Returns the closest point on line segment to point
    //-----------------------------------------------------------------------------
    IVector3D<T> ClosestPoint( const IVector3D<T>& point ) const
    {
        IVector3D<T> w = point - mOrigin;
        T proj = w.Dot(mDirection);
        // endpoint 0 is closest point
        if ( proj <= 0.0f )
            return mOrigin;
        else
        {
            T vsq = mDirection.Dot(mDirection);
            // endpoint 1 is closest point
            if ( proj >= vsq )
                return mOrigin + mDirection;
            // else somewhere else in segment
            else
                return mOrigin + (proj/vsq)*mDirection;
        }
    }




    bool IsPointOnLine(const IVector3D<T> p) // returns true if (p) is on lineSegment [GetEndpoint0(), GetEndpoint1()]
    {
        return (p - GetEndpoint0()).Dot(p - GetEndpoint1()) <= 0;
    }

    //----------[ output operator ]----------------------------
     /**
     * Provides output to standard output stream.
     */
     friend std::ostream& operator <<(std::ostream& oss, const ILineSegment3D<T>& rhs)
     {
   	     oss << "(" << "origin: " << rhs.mOrigin << " direction: " << rhs.mDirection << ")";
         return oss;
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


protected:

    IVector3D<T> mOrigin;
    IVector3D<T> mDirection;

};



//--------------------------------------
// Typedef shortcuts for LineSegment3D
//-------------------------------------

using ILineSegment3r    = ILineSegment3D<Real>;
using ILineSegment3f    = ILineSegment3D<float>;
using ILineSegment3d    = ILineSegment3D<double>;
using ILineSegment3i    = ILineSegment3D<std::int32_t>;
using ILineSegment3ui   = ILineSegment3D<std::uint32_t>;
using ILineSegment3b    = ILineSegment3D<std::int8_t>;
using ILineSegment3ub   = ILineSegment3D<std::uint8_t>;


}

#endif // ILINESEGMENT3D_H
