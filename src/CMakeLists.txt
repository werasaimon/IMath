cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "-std=c++11")
SET(CMAKE_CXX_STANDARD 11)

project(IMath)
add_executable(${PROJECT_NAME}
    main.cpp
    IMath/IReal.h
    IMath/IFunc.h
    IMath/IAlgebra.h
    IMath/ILine3D.h
    IMath/ILineSegment3D.h
    IMath/ILorentzVector.h
    IMath/IMaths.h
    IMath/IMatrix.h
    IMath/IMatrix2x2.h
    IMath/IMatrix3x3.h
    IMath/IMatrix4x4.h
    IMath/IPlane.h
    IMath/IQuaternion.h
    IMath/IRay.h
    IMath/ITransform.h
    IMath/IVector.h
    IMath/IVector2D.h
    IMath/IVector3D.h
    IMath/IVector4D.h
    IMath/IComplex.h
    IMath/IOctonion.h
    IMath/ISpherical.h
    IMath/IScalarType.h
    IMath/IVectorType.h)


