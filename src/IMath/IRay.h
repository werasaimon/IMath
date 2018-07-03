#ifndef IRAY_H
#define IRAY_H

#include "IVector3D.h"

namespace IMath
{

template<class T> class IRay
{
    public:

        // -------------------- Attributes -------------------- //

        /// First point of the ray (origin)
        IVector3D<T> point1;

        /// Second point of the ray
        IVector3D<T> point2;

        /// Maximum fraction value
        T maxFraction;


        // -------------------- Methods -------------------- //


        /// Constructor default
        SIMD_INLINE IRay(){}

        /// Constructor with arguments
        SIMD_INLINE IRay(const IVector3D<T>& p1, const IVector3D<T>& p2, T maxFrac = T(1.0))
            : point1(p1), point2(p2), maxFraction(maxFrac)
        {
        }

        /// Copy-constructor
        SIMD_INLINE IRay(const IRay<T>& ray)
            : point1(ray.point1), point2(ray.point2), maxFraction(ray.maxFraction)
        {

        }


        /// Overloaded assignment operator
        SIMD_INLINE IRay& operator=(const IRay& ray)
        {
            if (&ray != this)
            {
                point1 = ray.point1;
                point2 = ray.point2;
                maxFraction = ray.maxFraction;
            }
            return *this;
        }



        //-------------[ output operator ]------------------------
        /**
        * Output to stream operator
        * @param lhs Left hand side argument of operator (commonly ostream instance).
        * @param rhs Right hand side argument of operator.
        * @return Left hand side argument - the ostream object passed to operator.
        */
        friend std::ostream& operator<<(std::ostream& lhs, const IRay<T>& rhs)
        {
            lhs << "point1:" << rhs.point1 << " , " <<  "point2:" << rhs.point2 << " , " << "maxFriction:" << rhs.maxFraction;
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

};



//--------------------------------------
// Typedef shortcuts for Ray
//-------------------------------------
/// Three dimensional Ray of floats
typedef IRay<float> IRayf;
/// Three dimensional Ray of doubles
typedef IRay<double> IRayd;
/// Three dimensional Ray of ints
typedef IRay<int> IRayi;

} /* namespace */



#endif // IRAY_H
