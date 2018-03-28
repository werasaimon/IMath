#ifndef MATHS_H
#define MATHS_H

#include "func.hpp"
#include "rpVector2D.hpp"
#include "rpVector3D.hpp"
#include "rpVector4D.hpp"
#include "rpLorentzVector.hpp"
#include "rpMatrix2x2.hpp"
#include "rpMatrix3x3.hpp"
#include "rpMatrix4x4.hpp"
#include "rpComplex.hpp"
#include "rpQuaternion.hpp"
#include "rpOctonion.hpp"
#include "rpRay.hpp"
#include "rpTransform.hpp"
#include "rpLine3D.hpp"
#include "rpLineSegment3D.hpp"
#include "rpPlane.hpp"

#include <limits>

namespace CMath
{

typedef float scalar;


typedef rpVector2D<scalar>       Vector2;
typedef rpVector3D<scalar>       Vector3;
typedef rpVector4D<scalar>       Vector4;
typedef rpLorentzVector<scalar>  LorentzVector;
typedef rpMatrix2x2<scalar>      Matrix2;
typedef rpMatrix3x3<scalar>      Matrix3;
typedef rpMatrix4x4<scalar>      Matrix4;
typedef rpComplex<scalar>        Complex;
typedef rpQuaternion<scalar>     Quaternion;
typedef rpOctonion<scalar>       Octonion;
typedef rpRay<scalar>            Ray;
typedef rpTransform<scalar>      Transform;

typedef rpLine3D<scalar>         Line3;
typedef rpLineSegment3D<scalar>  LineSegment3;
typedef rpPlane<scalar>          Plane;

const scalar DECIMAL_SMALLEST = -std::numeric_limits<scalar>::max();
const scalar DECIMAL_LARGEST  =  std::numeric_limits<scalar>::max();

} /* namespace */

#endif // MATHS_H
