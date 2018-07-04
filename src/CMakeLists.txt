cmake_minimum_required(VERSION 2.8)

project(Math)
add_executable(Math
               main.cpp
               IMath/iFunc.h
               IMath/ILine3D.h
               IMath/ILineSegment3D.h
               IMath/ILorentzVector.h
               IMath/IMaths.h
               IMath/IMatrix2x2.h
               IMath/IMatrix3x3.h
               IMath/IMatrix4x4.h
               IMath/IPlane.h
               IMath/IQuaternion.h
               IMath/IRay.h
               IMath/ITransform.h
               IMath/IVector2D.h
               IMath/IVector3D.h
               IMath/IVector4D.h
               IMath/IComplex.h
               IMath/IOctonion)


