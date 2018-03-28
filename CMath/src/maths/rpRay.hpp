#ifndef RPRAY_H
#define RPRAY_H


#include "rpVector3D.hpp"

namespace CMath
{


template<class T> class rpRay
{

  
    public:

        // -------------------- Attributes -------------------- //

        /// First point of the ray (origin)
        rpVector3D<T> point1;

        /// Second point of the ray
        rpVector3D<T> point2;

        /// Maximum fraction value
        T maxFraction;


        // -------------------- Methods -------------------- //


        /// Constructor default
        SIMD_INLINE rpRay(){}

        /// Constructor with arguments
        SIMD_INLINE rpRay(const rpVector3D<T>& p1, const rpVector3D<T>& p2, T maxFrac = T(1.0))
            : point1(p1), point2(p2), maxFraction(maxFrac)
        {
        }

        /// Copy-constructor
        SIMD_INLINE rpRay(const rpRay<T>& ray)
            : point1(ray.point1), point2(ray.point2), maxFraction(ray.maxFraction)
        {

        }


        /// Overloaded assignment operator
        SIMD_INLINE rpRay& operator=(const rpRay& ray)
        {
            if (&ray != this)
            {
                point1 = ray.point1;
                point2 = ray.point2;
                maxFraction = ray.maxFraction;
            }
            return *this;
        }

};

} /* namespace */



#endif // RPRAY_H
