#ifndef RPPLANE_H
#define RPPLANE_H

#include "rpVector3D.hpp"
#include "rpQuaternion.hpp"
#include "rpLine3D.hpp"
#include "rpLineSegment3D.hpp"

namespace CMath
{
  
//-------------------------------------------------------------------------------
//-- Classes --------------------------------------------------------------------
//-------------------------------------------------------------------------------
template<class T>
class rpPlane
{
public:
    // constructor/destructor
    rpPlane( const rpVector3D<T>& normal = rpVector3D<T>(1,0,0) , T offset = T(0) )
      : mNormal(normal) ,
        mOffset(offset)
    {

    }

    rpPlane( T a, T b, T c, T d )
    {
       Set( a, b, c, d );
    }

    rpPlane( const rpVector3D<T>& p0,
             const rpVector3D<T>& p1,
             const rpVector3D<T>& p2 )
    {
       Set( p0, p1, p2 );
    }



    // copy operations
    rpPlane(const rpPlane& other)
     : mNormal( other.mNormal ),
       mOffset( other.mOffset )
    {

    }

    SIMD_INLINE ~rpPlane() {}


    // ---------------------------------------------------------------------------
    // Assigment operator
    //-----------------------------------------------------------------------------
    rpPlane& operator=(const rpPlane& other)
    {
        // if same object
        if ( this == &other )
            return *this;

        mNormal = other.mNormal;
        mOffset = other.mOffset;

        return *this;
    }


    // accessors
    SIMD_INLINE const rpVector3D<T>& GetNormal() const { return mNormal; }
    SIMD_INLINE T GetOffset() const { return mOffset; }

    // ---------------------------------------------------------------------------
    // Returns the two endpoints
    //-----------------------------------------------------------------------------
    void Get( rpVector3D<T>& normal, T& offset ) const
    {
        normal = mNormal;
        offset = mOffset;
    }

    // ---------------------------------------------------------------------------
    // Are two rpPlane's equal?
    //----------------------------------------------------------------------------
    bool operator==( const rpPlane&  plane ) const
    {
        return (plane.mNormal == mNormal &&
                plane.mOffset == mOffset);
    }

    // ---------------------------------------------------------------------------
    // Are two rpPlane's not equal?
    //----------------------------------------------------------------------------
    bool operator!=( const rpPlane&  plane ) const
    {
        return !(plane.mNormal == mNormal &&
                 plane.mOffset == mOffset);
    }

    // manipulators
    SIMD_INLINE void Set( const rpVector3D<T>& n, T d )
    {
        Set( n.x, n.y, n.z, d );
    }

    // ---------------------------------------------------------------------------
    // Sets the parameters
    //-----------------------------------------------------------------------------
    void Set( T a, T b, T c, T d )
    {
        // normalize for cheap distance checks
        T lensq = a*a + b*b + c*c;
        // length of normal had better not be zero
        assert( !IsZero( lensq ) );

        // recover gracefully
        if ( IsZero( lensq ) )
        {
            mNormal = rpVector3D<T>::X;
            mOffset = 0.0f;
        }
        else
        {
            T recip = 1.0/Sqrt(lensq);
            mNormal.setAllValues( a*recip, b*recip, c*recip );
            mOffset = d*recip;
        }
    }



    // ---------------------------------------------------------------------------
    // Sets the parameters
    //-----------------------------------------------------------------------------
    void Set( const rpVector3D<T>& p0, const rpVector3D<T>& p1, const rpVector3D<T>& p2 )
    {
        // get plane vectors
        rpVector3D<T> u = p1 - p0;
        rpVector3D<T> v = p2 - p0;
        rpVector3D<T> w = u.cross(v);

        // normalize for cheap distance checks
        T lensq = w.x*w.x + w.y*w.y + w.z*w.z;
        // length of normal had better not be zero
        assert( !IsZero( lensq ) );

        // recover gracefully
        if ( IsZero( lensq ) )
        {
            mNormal = rpVector3D<T>::X;
            mOffset = 0.0f;
        }
        else
        {
            T recip = 1.0f/lensq;
            mNormal.setAllValues( w.x*recip, w.y*recip, w.z*recip );
            mOffset = -mNormal.dot(p0);
        }
    }

    // ---------------------------------------------------------------------------
    // Transforms plane into new space
    //-----------------------------------------------------------------------------
    rpPlane Transform( T scale, const rpQuaternion<T>& rotate, const rpVector3D<T>& translate ) const
    {
        rpPlane<T> plane;

        // get rotation matrix
        rpMatrix3x3<T>  rotmatrix = rotate.getMatrix();

        // transform to get normal
        plane.mNormal = rotmatrix*mNormal/scale;

        // transform to get offset
        rpVector3D<T> newTrans = translate*rotmatrix;
        plane.mOffset = -newTrans.dot( mNormal )/scale + mOffset;

        return plane;
    }

    // ---------------------------------------------------------------------------
    // Transforms plane into new space
    //-----------------------------------------------------------------------------
    rpPlane Transform( T scale, const rpMatrix3x3<T>& rotmatrix, const rpVector3D<T>& translate ) const
    {
        rpPlane<T> plane;

        // transform to get normal
        plane.mNormal = rotmatrix*mNormal/scale;

        // transform to get offset
        rpVector3D<T> newTrans = translate*rotmatrix;
        plane.mOffset = -newTrans.dot( mNormal )/scale + mOffset;

        return plane;
    }

    // distance
    static T Distance( const rpPlane& plane, const rpVector3D<T>& point )
    {
        return Abs( plane.Test( point ) );
    }

    // ---------------------------------------------------------------------------
    // Returns the closest point on plane to point
    //-----------------------------------------------------------------------------
    rpVector3D<T> ClosestPoint( const rpVector3D<T>& point ) const
    {
        return point - Test(point)*mNormal;
    }

    // result of plane test
    SIMD_INLINE T Test( const rpVector3D<T>& point ) const
    {
        return mNormal.dot(point) - mOffset;
    }



    // intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
    //    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
    //    Output: *I0 = the intersect point (when it exists)
    //    Return: 0 = disjoint (no intersection)
    //            1 =  intersection in the unique point *I0
    //            2 = the  segment lies in the plane
    int IntersectSegmentToPlane( const rpLineSegment3D<T>& line , rpVector3D<T>* I )
    {
        #define SMALL_NUM  0.00001// anything that avoids division overflow

        rpVector3D<T>    u = line.GetEndpoint1() - line.GetEndpoint0();

        T D =  mNormal.dot(u);
        T N = -Test(line.GetEndpoint0());

        if (Abs(D) < SMALL_NUM) {            // segment is parallel to plane
            if (N == 0)                      // segment lies in plane
                return 2;
            else
                return 0;                    // no intersection
        }
        // they are not parallel
        // compute intersect param
        T sI = T(N / D);
        if (sI < 0 || sI > 1)
            return 0;                        // no intersection

        *I = line.GetEndpoint0() + sI * u;               // compute segment intersect point
        return 1;
    }


    //----------[ output operator ]----------------------------
     /**
     * Provides output to standard output stream.
     */
     friend std::ostream& operator <<(std::ostream& oss, const rpPlane<T>& rhs)
     {
   	     oss << "(" << "normal: " << rhs.mNormal << " offset: " << rhs.mOffset << ")";
         return oss;
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

protected:
    rpVector3D<T> mNormal;
    T             mOffset;

private:
};


}


#endif // RPPLANE_H
