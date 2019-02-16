#ifndef IAXISALIGNEDBOX3D_H
#define IAXISALIGNEDBOX3D_H

#include "IVector3D.h"
#include "IMatrix4x4.h"
#include "IRay.h"

namespace IMath
{


   template<class T> class IAxisAlignedBox3D
   {
       IVector3D<T> mMin;
       IVector3D<T> mMax;


     public:

       static const IAxisAlignedBox3D<T> Empty;
       static const IAxisAlignedBox3D<T> Zero;
       static const IAxisAlignedBox3D<T> UnitPositive;


       SIMD_INLINE IAxisAlignedBox3D()
       {
           mMin = IVector3D<T>( 1,  1,  1);
           mMax = IVector3D<T>(-1, -1, -1);
       }

       SIMD_INLINE IAxisAlignedBox3D(T xmin, T ymin, T zmin,
                                     T xmax, T ymax, T zmax)
       {
           mMin = IVector3D<T>(xmin, ymin, zmin);
           mMax = IVector3D<T>(xmax, ymax, zmax);
       }

       /// <summary>
       /// init box [0,size] x [0,size] x [0,size]
       /// </summary>
       SIMD_INLINE IAxisAlignedBox3D(T fCubeSize)
       {
           mMin = IVector3D<T>(0, 0, 0);
           mMax = IVector3D<T>(fCubeSize, fCubeSize, fCubeSize);
       }

       /// <summary>
       /// Init box [0,width] x [0,height] x [0,depth]
       /// </summary>
       SIMD_INLINE IAxisAlignedBox3D(T fWidth, T fHeight, T fDepth)
       {
           mMin = IVector3D<T>(0, 0, 0);
           mMax = IVector3D<T>(fWidth, fHeight, fDepth);
       }

       SIMD_INLINE IAxisAlignedBox3D(IVector3D<T>& vMin, IVector3D<T>& vMax)
       {
           mMin = IVector3D<T>(IMin(vMin.x, vMax.x), IMin(vMin.y, vMax.y), IMin(vMin.z, vMax.z));
           mMax = IVector3D<T>(IMax(vMin.x, vMax.x), IMax(vMin.y, vMax.y), IMax(vMin.z, vMax.z));
       }

       SIMD_INLINE IAxisAlignedBox3D(const IVector3D<T>& vMin, const IVector3D<T>& vMax)
       {
           mMin = IVector3D<T>(IMin(vMin.x, vMax.x), IMin(vMin.y, vMax.y), IMin(vMin.z, vMax.z));
           mMax = IVector3D<T>(IMax(vMin.x, vMax.x), IMax(vMin.y, vMax.y), IMax(vMin.z, vMax.z));
       }

       SIMD_INLINE IAxisAlignedBox3D(const IVector3D<T>& vCenter, T fHalfWidth, T fHalfHeight, T fHalfDepth)
       {
           mMin = IVector3D<T>(vCenter.x - fHalfWidth, vCenter.y - fHalfHeight, vCenter.z - fHalfDepth);
           mMax = IVector3D<T>(vCenter.x + fHalfWidth, vCenter.y + fHalfHeight, vCenter.z + fHalfDepth);
       }

       SIMD_INLINE IAxisAlignedBox3D(IVector3D<T>& vCenter, T fHalfWidth, T fHalfHeight, T fHalfDepth)
       {
           mMin = IVector3D<T>(vCenter.x - fHalfWidth, vCenter.y - fHalfHeight, vCenter.z - fHalfDepth);
           mMax = IVector3D<T>(vCenter.x + fHalfWidth, vCenter.y + fHalfHeight, vCenter.z + fHalfDepth);
       }

       SIMD_INLINE IAxisAlignedBox3D(const IVector3D<T>& vCenter, T fHalfSize)
       {
           mMin = IVector3D<T>(vCenter.x - fHalfSize, vCenter.y - fHalfSize, vCenter.z - fHalfSize);
           mMax = IVector3D<T>(vCenter.x + fHalfSize, vCenter.y + fHalfSize, vCenter.z + fHalfSize);
       }

       SIMD_INLINE IVector3D<T> GetMax() const
       {
           return mMax;
       }

       SIMD_INLINE IVector3D<T> GetMin() const
       {
           return mMin;
       }

       SIMD_INLINE IVector3D<T> HalfSize() const
       {
           IVector3D<T> size;
           for (std::size_t i = 0; i < IVector3D<T>::components; ++i)
           {
               size[i] = IAbs(mMax[i] - mMin[i]);
           }
           return size;
       }

       SIMD_INLINE IAxisAlignedBox3D(const IVector3D<T>& vCenter)
       {
           mMin = mMax = vCenter;
       }

       SIMD_INLINE T Width()  const
       {
            return IMax(mMax.x - mMin.x, T(0));
       }

       SIMD_INLINE T Height()  const
       {
            return IMax(mMax.y - mMin.y, 0);
       }

       SIMD_INLINE T Depth()  const
       {
            return IMax(mMax.z - mMin.z, 0);
       }

       SIMD_INLINE T Volume()  const
       {
            return Width() * Height() * Depth();
       }

       SIMD_INLINE T DiagonalLength()  const
       {
            return ISqrt((mMax.x - mMin.x) * (mMax.x - mMin.x) +
                         (mMax.y - mMin.y) * (mMax.y - mMin.y) +
                         (mMax.z - mMin.z) * (mMax.z - mMin.z));
       }

       SIMD_INLINE T MaxDim()  const
       {
            return IMax(Width(), IMax(Height(), Depth()));
       }

       SIMD_INLINE IVector3D<T> Diagonal()  const
       {
            return IVector3D<T>(mMax.x - mMin.x, mMax.y - mMin.y, mMax.z-mMin.z);
       }

       SIMD_INLINE IVector3D<T> Extents()  const
       {
            return IVector3D<T>((mMax.x - mMin.x)*0.5,
                                (mMax.y - mMin.y)*0.5,
                                (mMax.z - mMin.z)*0.5);
       }

       SIMD_INLINE IVector3D<T> Center()  const
       {
            return IVector3D<T>(0.5 * (mMin.x + mMax.x),
                                0.5 * (mMin.y + mMax.y),
                                0.5 * (mMin.z + mMax.z));
       }

       SIMD_INLINE bool operator ==(const IAxisAlignedBox3D& b) const
       {
           return mMin == b.mMin && mMax == b.mMax;
       }

       SIMD_INLINE bool operator != (const IAxisAlignedBox3D& b)
       {
           return mMin != b.mMin || mMax != b.mMax;
       }

       SIMD_INLINE bool Equals(IAxisAlignedBox3D other)
       {
           return this == other;
       }

//       int CompareTo(IAxisAlignedBox3D other)
//       {
//           int c = this.Min.CompareTo(other.Min);
//           if (c == 0)
//               return this.Max.CompareTo(other.Max);
//           return c;
//       }

       SIMD_INLINE int GetHashCode()  const
       {
           // Overflow is fine, just wrap
           int hash = (int) 2166136261;
           hash = (hash * 16777619) ^ mMin.GetHashCode();
           hash = (hash * 16777619) ^ mMax.GetHashCode();
           return hash;
       }


       // See Box3.Corner for details on which corner is which
       SIMD_INLINE IVector3D<T> Corner(int i)  const
       {
           T x = (  ((i&1) != 0) ^ ((i&2) != 0) ) ? (mMax.x) : (mMin.x);
           T y = ( (i / 2) % 2 == 0 ) ? (mMin.y) : (mMax.y);
           T z = (i < 4) ? (mMin.z) : (mMax.z);
           return IVector3D<T>(x, y, z);
       }

       /// <summary>
       /// Returns point on face/edge/corner. For each coord value neg==min, 0==center, pos==max
       /// </summary>
       SIMD_INLINE IVector3D<T> Point(int xi, int yi, int zi)  const
       {
           T x = (xi < 0) ? mMin.x : ((xi == 0) ? (0.5*(mMin.x + mMax.x)) : mMax.x);
           T y = (yi < 0) ? mMin.y : ((yi == 0) ? (0.5*(mMin.y + mMax.y)) : mMax.y);
           T z = (zi < 0) ? mMin.z : ((zi == 0) ? (0.5*(mMin.z + mMax.z)) : mMax.z);
           return IVector3D<T>(x, y, z);
       }


       // TODO
       ////! 0 == bottom-left, 1 = bottom-right, 2 == top-right, 3 == top-left
       //IVector3D<T> GetCorner(int i) {
       //    return IVector3D<T>((i % 3 == 0) ? Min.x : Max.x, (i < 2) ? Min.y : Max.y);
       //}

       //! value is subtracted from min and added to max
       SIMD_INLINE void Expand(T fRadius)   const
       {
           mMin.x -= fRadius; mMin.y -= fRadius; mMin.z -= fRadius;
           mMax.x += fRadius; mMax.y += fRadius; mMax.z += fRadius;
       }

       //! return this box expanded by radius
       SIMD_INLINE IAxisAlignedBox3D Expanded(T fRadius) const
       {
           return IAxisAlignedBox3D(
               mMin.x - fRadius, mMin.y - fRadius, mMin.z - fRadius,
               mMax.x + fRadius, mMax.y + fRadius, mMax.z + fRadius);
       }

       //! value is added to min and subtracted from max
       SIMD_INLINE void Contract(T fRadius)  const
       {
           T w = 2 * fRadius;
           if ( w > mMax.x-mMin.x )
           {
               mMin.x = mMax.x = 0.5 * (mMin.x + mMax.x);
           }
           else
           {
               mMin.x += fRadius; mMax.x -= fRadius;
           }
           if ( w > mMax.y-mMin.y )
           {
               mMin.y = mMax.y = 0.5 * (mMin.y + mMax.y);
           }
           else
           {
               mMin.y += fRadius; mMax.y -= fRadius;
           }
           if ( w > mMax.z-mMin.z )
           {
               mMin.z = mMax.z = 0.5 * (mMin.z + mMax.z);
           }
           else
           {
               mMin.z += fRadius; mMax.z -= fRadius;
           }
       }

       //! return this box expanded by radius
       SIMD_INLINE IAxisAlignedBox3D Contracted(T fRadius)  const
       {
           IAxisAlignedBox3D result = IAxisAlignedBox3D(
               mMin.x + fRadius, mMin.y + fRadius, mMin.z + fRadius,
               mMax.x - fRadius, mMax.y - fRadius, mMax.z - fRadius);
           if (result.Min.x > result.Max.x) { result.Min.x = result.Max.x = 0.5 * (mMin.x + mMax.x); }
           if (result.Min.y > result.Max.y) { result.Min.y = result.Max.y = 0.5 * (mMin.y + mMax.y); }
           if (result.Min.z > result.Max.z) { result.Min.z = result.Max.z = 0.5 * (mMin.z + mMax.z); }
           return result;
       }


       SIMD_INLINE void Scale(T sx, T sy, T sz)
       {
           IVector3D<T> c = Center();
           IVector3D<T> e = Extents(); e.x *= sx; e.y *= sy; e.z *= sz;
           mMin = IVector3D<T>(c.x - e.x, c.y - e.y, c.z - e.z);
           mMax = IVector3D<T>(c.x + e.x, c.y + e.y, c.z + e.z);
       }

       SIMD_INLINE void MergeWithAABB(IVector3D<T>& v)
       {
           mMin.x = IMin(mMin.x, v.x);
           mMin.y = IMin(mMin.y, v.y);
           mMin.z = IMin(mMin.z, v.z);
           mMax.x = IMax(mMax.x, v.x);
           mMax.y = IMax(mMax.y, v.y);
           mMax.z = IMax(mMax.z, v.z);
       }

       SIMD_INLINE void MergeWithAABB(const IAxisAlignedBox3D& box)
       {
           mMin.x = IMin(mMin.x, box.Min.x);
           mMin.y = IMin(mMin.y, box.Min.y);
           mMin.z = IMin(mMin.z, box.Min.z);
           mMax.x = IMax(mMax.x, box.Max.x);
           mMax.y = IMax(mMax.y, box.Max.y);
           mMax.z = IMax(mMax.z, box.Max.z);
       }


       // Replace the current AABB with a new AABB that is the union of two AABBs in parameters
       SIMD_INLINE void ContainsTwoAABBs(const IAxisAlignedBox3D& aabb1, const IAxisAlignedBox3D& aabb2)
       {
           mMin.x = IMin(aabb1.mMin.x, aabb2.mMin.x);
           mMin.y = IMin(aabb1.mMin.y, aabb2.mMin.y);
           mMin.z = IMin(aabb1.mMin.z, aabb2.mMin.z);

           mMax.x = IMax(aabb1.mMax.x, aabb2.mMax.x);
           mMax.y = IMax(aabb1.mMax.y, aabb2.mMax.y);
           mMax.z = IMax(aabb1.mMax.z, aabb2.mMax.z);
       }

       SIMD_INLINE IAxisAlignedBox3D Intersect(const IAxisAlignedBox3D& box)  const
       {
               IAxisAlignedBox3D intersect = IAxisAlignedBox3D(
               IMax(mMin.x, box.Min.x), IMax(mMin.y, box.Min.y), IMax(mMin.z, box.Min.z),
               IMin(mMax.x, box.Max.x), IMin(mMax.y, box.Max.y), IMin(mMax.z, box.Max.z));
               if (intersect.Height <= 0 || intersect.Width <= 0 || intersect.Depth <= 0)
               {
                   return Empty;
               }
               else
               {
                   return intersect;
               }
       }

       // Return true if the IAABB of a triangle intersects the IAABB
       SIMD_INLINE bool TestCollisionTriangleAABB(const IVector3D<T>* trianglePoints) const
       {

           if (IMin3(trianglePoints[0].x, trianglePoints[1].x, trianglePoints[2].x) > mMax.x) return false;
           if (IMin3(trianglePoints[0].y, trianglePoints[1].y, trianglePoints[2].y) > mMax.y) return false;
           if (IMin3(trianglePoints[0].z, trianglePoints[1].z, trianglePoints[2].z) > mMax.z) return false;

           if (IMax3(trianglePoints[0].x, trianglePoints[1].x, trianglePoints[2].x) < mMin.x) return false;
           if (IMax3(trianglePoints[0].y, trianglePoints[1].y, trianglePoints[2].y) < mMin.y) return false;
           if (IMax3(trianglePoints[0].z, trianglePoints[1].z, trianglePoints[2].z) < mMin.z) return false;

           return true;
       }


       SIMD_INLINE bool Contains(IVector3D<T>& v)  const
       {
           return (mMin.x <= v.x) && (mMin.y <= v.y) && (mMin.z <= v.z)
               && (mMax.x >= v.x) && (mMax.y >= v.y) && (mMax.z >= v.z);
       }

       SIMD_INLINE bool Contains(const IVector3D<T>& v)  const
       {
           return (mMin.x <= v.x) && (mMin.y <= v.y) && (mMin.z <= v.z) &&
                  (mMax.x >= v.x) && (mMax.y >= v.y) && (mMax.z >= v.z);
       }

       SIMD_INLINE bool Contains(IAxisAlignedBox3D box2)  const
       {
           return Contains(box2.Min) && Contains(box2.Max);
       }

       SIMD_INLINE bool Contains(const IAxisAlignedBox3D& box2)  const
       {
           return Contains(box2.Min) && Contains(box2.Max);
       }


       SIMD_INLINE bool Intersects(const IAxisAlignedBox3D& box)  const
       {
           return !((box.Max.x <= mMin.x) || (box.Min.x >= mMax.x)
                 || (box.Max.y <= mMin.y) || (box.Min.y >= mMax.y)
                 || (box.Max.z <= mMin.z) || (box.Min.z >= mMax.z) );
       }


       SIMD_INLINE bool IntersectsRay( IRay<T> &r )
       {

           T enter = T(0.0f);
           T exit  = T(1.0f);

           T tx1 = (mMin[0] - r.Origin[0]) / r.Direction[0];
           T tx2 = (mMax[0] - r.Origin[0]) / r.Direction[0];

           T tmin = IMath::IMin(tx1, tx2);
           T tmax = IMath::IMax(tx1, tx2);

           for (std::size_t i = 0; i < IVector3D<T>::components; ++i)
           {
               //if(IMath::IAbs(r.direction[i]) < T(0.00001) ) continue;


               if (IMath::IAbs(r.Direction[i]) < T(0.0001))
               {
                   // The segment is parallel to slab. No hit if origin not within slab.
                   if (r.Origin[i] < mMin[i] || r.Origin[i] > mMax[i])
                   {
                       return false;
                   }
               }
               else
               {
                   T tx1 = (mMin[i] - r.Origin[i]) / r.Direction[i];
                   T tx2 = (mMax[i] - r.Origin[i]) / r.Direction[i];

                   // Make t1 be intersection with near plane, t2 with far plane.
                   if (tx1 > tx2)
                   {
                       IMath::ISwap(tx1, tx2);
                   }

                   tmin = IMath::IMax(tmin, IMath::IMin(tx1, tx2));
                   tmax = IMath::IMin(tmax, IMath::IMax(tx1, tx2));

                   if(tmin > tmax)
                   {
                     //IMath::ISwap(tmin,tmax);
                     return false;
                   }

                   //Reduce interval based on intersection
                   if(tmin > enter) enter = tmin;
                   if(tmax < exit)  exit  = tmax;
               }
           }

           r.maxFraction = enter;

           return true;
       }



       SIMD_INLINE T DistanceSquared(const IVector3D<T>& v)  const
       {
           T dx = (v.x < mMin.x) ? mMin.x - v.x : (v.x > mMax.x ? v.x - mMax.x : 0);
           T dy = (v.y < mMin.y) ? mMin.y - v.y : (v.y > mMax.y ? v.y - mMax.y : 0);
           T dz = (v.z < mMin.z) ? mMin.z - v.z : (v.z > mMax.z ? v.z - mMax.z : 0);
           return dx * dx + dy * dy + dz * dz;
       }

       SIMD_INLINE T Distance(const IVector3D<T>& v)  const
       {
           return ISqrt(DistanceSquared(v));
       }

       SIMD_INLINE T SignedDistance(const IVector3D<T>& v)  const
       {
           if ( Contains(v) == false ) {
               return Distance(v);
           } else {
               T dx = IMin(IAbs(v.x - mMin.x), IAbs(v.x - mMax.x));
               T dy = IMin(IAbs(v.y - mMin.y), IAbs(v.y - mMax.y));
               T dz = IMin(IAbs(v.z - mMin.z), IAbs(v.z - mMax.z));
               return -IMin(dx, dy, dz);
           }
       }


       SIMD_INLINE T DistanceSquared(const IAxisAlignedBox3D& box2)  const
       {
           // compute lensqr( max(0, abs(center1-center2) - (extent1+extent2)) )
           T delta_x = IAbs((box2.Min.x + box2.Max.x) - (mMin.x + mMax.x)) - ((mMax.x - mMin.x) +
                            (box2.Max.x - box2.Min.x));
           if ( delta_x < 0 )
                delta_x = 0;
           T delta_y = IAbs((box2.Min.y + box2.Max.y) - (mMin.y + mMax.y))
                         - ((mMax.y - mMin.y) + (box2.Max.y - box2.Min.y));
           if (delta_y < 0)
               delta_y = 0;
           T delta_z = IAbs((box2.Min.z + box2.Max.z) - (mMin.z + mMax.z))
                         - ((mMax.z - mMin.z) + (box2.Max.z - box2.Min.z));
           if (delta_z < 0)
               delta_z = 0;
           return 0.25 * (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
       }


       // [TODO] we have handled corner cases, but not edge cases!
       //   those are 2D, so it would be like (dx > width && dy > height)
       //T Distance(const IVector3D<T>& v)
       //{
       //    T dx = (double)IAbs(v.x - Center.x);
       //    T dy = (double)IAbs(v.y - Center.y);
       //    T dz = (double)IAbs(v.z - Center.z);
       //    T fWidth = Width * 0.5;
       //    T fHeight = Height * 0.5;
       //    T fDepth = Depth * 0.5;
       //    if (dx < fWidth && dy < fHeight && dz < Depth)
       //        return 0.0f;
       //    else if (dx > fWidth && dy > fHeight && dz > fDepth)
       //        return (double)Math.Sqrt((dx - fWidth) * (dx - fWidth) + (dy - fHeight) * (dy - fHeight) + (dz - fDepth) * (dz - fDepth));
       //    else if (dx > fWidth)
       //        return dx - fWidth;
       //    else if (dy > fHeight)
       //        return dy - fHeight;
       //    else if (dz > fDepth)
       //        return dz - fDepth;
       //    return 0.0f;
       //}


       SIMD_INLINE IAxisAlignedBox3D GetAxisAlignedBoxTransform(const IMatrix4x4<T> &_TransformMatrix) const
       {
           IAxisAlignedBox3D box;

           // Get the local bounds in x,y and z direction
           IVector3D<T> minBounds = mMin;
           IVector3D<T> maxBounds = mMax;

           // Scale size AABB local
           minBounds *= 1.04;
           maxBounds *= 1.04;


           // Rotate the local bounds according to the orientation of the body
           IMatrix3x3<T> worldAxis =  (_TransformMatrix.GetRotMatrix()).GetAbsoluteMatrix();

           IVector3D<T>   worldMinBounds(worldAxis.GetColumn(0).Dot(minBounds),
                                         worldAxis.GetColumn(1).Dot(minBounds),
                                         worldAxis.GetColumn(2).Dot(minBounds));

           IVector3D<T>   worldMaxBounds(worldAxis.GetColumn(0).Dot(maxBounds),
                                         worldAxis.GetColumn(1).Dot(maxBounds),
                                         worldAxis.GetColumn(2).Dot(maxBounds));


           // Compute the minimum and maximum coordinates of the rotated extents
           IVector3D<T> minCoordinates = _TransformMatrix.GetTranslation() + worldMinBounds;
           IVector3D<T> maxCoordinates = _TransformMatrix.GetTranslation() + worldMaxBounds;


           // Update the AABB with the new minimum and maximum coordinates
           return IAxisAlignedBox3D(minCoordinates , maxCoordinates);

       }



       //! relative translation
       SIMD_INLINE void Translate(const IVector3D<T>& vTranslate)
       {
           mMin += vTranslate;
           mMax += vTranslate;
       }

       SIMD_INLINE void MoveMin(const IVector3D<T>& vNewMin)
       {
           mMax.x = vNewMin.x + (mMax.x - mMin.x);
           mMax.y = vNewMin.y + (mMax.y - mMin.y);
           mMax.z = vNewMin.z + (mMax.z - mMin.z);
           mMin.SetAllValues(vNewMin);
       }

       SIMD_INLINE void MoveMin(T fNewX, T fNewY, T fNewZ)
       {
           mMax.x = fNewX + (mMax.x - mMin.x);
           mMax.y = fNewY + (mMax.y - mMin.y);
           mMax.z = fNewZ + (mMax.z - mMin.z);
           mMin.SetAllValues(fNewX, fNewY, fNewZ);
       }


       SIMD_INLINE void Insert(const IVector3D<T>& point)
       {
           for (std::size_t i = 0; i < IVector3D<T>::components; ++i)
           {
               if (mMin[i] > point[i]) mMin[i] = point[i];
               if (mMax[i] < point[i]) mMax[i] = point[i];
           }
       }

       SIMD_INLINE void Insert(const IAxisAlignedBox3D& aabb)
       {
           for (std::size_t i = 0; i < IVector3D<T>::components; ++i)
           {
               if (mMin[i] > aabb.mMin[i]) mMin[i] = aabb.mMin[i];
               if (mMax[i] < aabb.mMax[i]) mMax[i] = aabb.mMax[i];
           }
       }

       SIMD_INLINE void Repair()
       {
           for (std::size_t i = 0; i < IVector3D<T>::components; ++i)
           {
               if (mMin[i] > mMax[i])
               {
                   ISwap(mMin[i], mMax[i]);
               }
           }
       }

//        explicit operator IAxisAlignedBox3D(IAxisAlignedBox3D v)
//       {
//           return IAxisAlignedBox3D(v.Min, v.Max);
//       }
//       static explicit operator AxisAlignedBox3f(IAxisAlignedBox3D v)
//       {
//           return AxisAlignedBox3f((Vector3f)v.Min, (Vector3f)v.Max);
//       }




       //-------------[ output operator ]------------------------
       /**
            * Output to stream operator
            * @param lhs Left hand side argument of operator (commonly ostream instance).
            * @param rhs Right hand side argument of operator.
            * @return Left hand side argument - the ostream object passed to operator.
            */
       friend std::ostream& operator<<(std::ostream& lhs, const IAxisAlignedBox3D& rhs)
       {
           lhs << "Min[" << rhs.Min[0] << "," << rhs.Min[1] << "," << rhs.Min[2] << "] ";
           lhs << "Max[" << rhs.Max[0] << "," << rhs.Max[1] << "," << rhs.Max[2] << "] ";
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



template<class T> const IAxisAlignedBox3D<T> IAxisAlignedBox3D<T>::Empty = IAxisAlignedBox3D(false);
template<class T> const IAxisAlignedBox3D<T> IAxisAlignedBox3D<T>::Zero  = IAxisAlignedBox3D(0);
template<class T> const IAxisAlignedBox3D<T> IAxisAlignedBox3D<T>::UnitPositive = IAxisAlignedBox3D(1);



//--------------------------------------
// Typedef shortcuts for IAxisAlignedBox3D
//-------------------------------------
using IAxisAlignedBox3Dr    = IAxisAlignedBox3D<Real>;
using IAxisAlignedBox3Df    = IAxisAlignedBox3D<float>;
using IAxisAlignedBox3Dd    = IAxisAlignedBox3D<double>;
using IAxisAlignedBox3Di    = IAxisAlignedBox3D<std::int32_t>;
using IAxisAlignedBox3Dui   = IAxisAlignedBox3D<std::uint32_t>;
using IAxisAlignedBox3Db    = IAxisAlignedBox3D<std::int8_t>;
using IAxisAlignedBox3Dub   = IAxisAlignedBox3D<std::uint8_t>;




}




#endif // IAXISALIGNEDBOX3D_H
