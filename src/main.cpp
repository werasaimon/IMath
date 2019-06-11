#include <iostream>

#include "IMath/IMaths.h"

using namespace std;

typedef IMath::Real scalar;

typedef IMath::IVector2D<scalar>       Vector2;
typedef IMath::IVector3D<scalar>       Vector3;
typedef IMath::IVector4D<scalar>       Vector4;
typedef IMath::ILorentzVector<scalar>  LorentzVector;
typedef IMath::IMatrix2x2<scalar>      Matrix2;
typedef IMath::IMatrix3x3<scalar>      Matrix3;
typedef IMath::IMatrix4x4<scalar>      Matrix4;
typedef IMath::IQuaternion<scalar>     Quaternion;
typedef IMath::IRay<scalar>            Ray;
typedef IMath::ITransform<scalar>      Transform;
typedef IMath::IComplex<scalar>        Complex;
typedef IMath::IOctonion<scalar>       Octonion;
typedef IMath::ILine3D<scalar>         Line3;
typedef IMath::ILineSegment3D<scalar>  LineSegment3;
typedef IMath::IPlane<scalar>          Plane;

int main()
{
    cout << "Hello World!" << endl;



    /**/

        {
            cout<< "*********************************" << endl;
            cout<< "---- Unit test Vector2 -----" <<endl<<endl;

            Vector2 v0(7 , 3);
            Vector2 v1(1 , 6);

            cout<< "----constant operator----" <<endl;
            cout << " Vector2:  " <<  v0 << " + " << v1 << " = " << (v0 + v1) <<endl;
            cout << " Vector2:  " <<  v0 << " - " << v1 << " = " << (v0 - v1) <<endl;
            cout << " Vector2:  " <<  v0 << " * " << v1 << " = " << (v0 * v1) <<endl;
            cout << " Vector2:  " <<  v0 << " / " << v1 << " = " << (v0 / v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector2:  " <<  v0 << " += " << v1 << " = " << (v0 += v1) <<endl;
            cout << " Vector2:  " <<  v0 << " -= " << v1 << " = " << (v0 -= v1) <<endl;
            cout << " Vector2:  " <<  v0 << " *= " << v1 << " = " << (v0 *= v1) <<endl;
            cout << " Vector2:  " <<  v0 << " /= " << v1 << " = " << (v0 /= v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector2:  " <<  v0 << " .dot " << v1 << " = " << (v0.Dot(v1)) <<endl;
            cout << " Vector2:  " <<  v0 << " .cross " << v1 << " = " << (v0.Cross(v1)) <<endl;
            cout << " Vector2:  " <<  v0 << " .angle " << v1 << " = " << (v0.GetAngleBetween(v1)) <<endl;
            cout << "----function operator----" <<endl;
            cout << " Vector2:  " <<  v0 << ".lengthSquare" << " = " << v0.LengthSquare() <<endl;
            cout << " Vector2:  " <<  v0 << ".length" << " = " << v0.Length() <<endl;
            cout << " Vector2:  " <<  v0 << ".inverse" << " = " << v0.GetInverse() <<endl;
            cout << " Vector2:  " <<  v0 << ".unit" << " = " << v0.GetUnit() <<endl;
            Vector2 v = v0; v.Normalize();
            cout << " Vector2:  " <<  v0 << ".normalize" << " = " << v <<endl;
            cout<<endl<<endl;
        }


        {

            cout<< "*********************************" << endl;
            cout<< "---- Unit test Vector3 -----" <<endl<<endl;

            Vector3 v0(7 , 3 , 1);
            Vector3 v1(1 , 6 , 8);

            cout<< "----constant operator----" <<endl;
            cout << " Vector3:  " <<  v0 << " + " << v1 << " = " << (v0 + v1) <<endl;
            cout << " Vector3:  " <<  v0 << " - " << v1 << " = " << (v0 - v1) <<endl;
            cout << " Vector3:  " <<  v0 << " * " << v1 << " = " << (v0 * v1) <<endl;
            cout << " Vector3:  " <<  v0 << " / " << v1 << " = " << (v0 / v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector3:  " <<  v0 << " += " << v1 << " = " << (v0 += v1) <<endl;
            cout << " Vector3:  " <<  v0 << " -= " << v1 << " = " << (v0 -= v1) <<endl;
            cout << " Vector3:  " <<  v0 << " *= " << v1 << " = " << (v0 *= v1) <<endl;
            cout << " Vector3:  " <<  v0 << " /= " << v1 << " = " << (v0 /= v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector3:  " <<  v0 << " .dot " << v1 << " = " << (v0.Dot(v1)) <<endl;
            cout << " Vector3:  " <<  v0 << " .cross " << v1 << " = " << (v0.Cross(v1)) <<endl;
            cout << " Vector3:  " <<  v0 << " .angle " << v1 << " = " << (v0.GetAngleBetween(v1)) <<endl;
            cout << "----function operator----" <<endl;
            cout << " Vector3:  " <<  v0 << ".lengthSquare" << " = " << v0.LengthSquare() <<endl;
            cout << " Vector3:  " <<  v0 << ".length" << " = " << v0.Length() <<endl;
            cout << " Vector3:  " <<  v0 << ".inverse" << " = " << v0.GetInverse() <<endl;
            cout << " Vector3:  " <<  v0 << ".unit" << " = " << v0.GetUnit() <<endl;
            Vector3 v = v0; v.Normalize();
            cout << " Vector3:  " <<  v0 << ".normalize" << " = " << v <<endl;
            cout<<endl<<endl;

        }


        {

            cout<< "*********************************" << endl;
            cout<< "---- Unit test Vector4 -----" <<endl<<endl;

            Vector4 v0(7 , 3 , 1, 2);
            Vector4 v1(1 , 6 , 8, 9);

            cout<< "----constant operator----" <<endl;
            cout << " Vector4:  " <<  v0 << " + " << v1 << " = " << (v0 + v1) <<endl;
            cout << " Vector4:  " <<  v0 << " - " << v1 << " = " << (v0 - v1) <<endl;
            cout << " Vector4:  " <<  v0 << " * " << v1 << " = " << (v0 * v1) <<endl;
            cout << " Vector4:  " <<  v0 << " / " << v1 << " = " << (v0 / v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector4:  " <<  v0 << " += " << v1 << " = " << (v0 += v1) <<endl;
            cout << " Vector4:  " <<  v0 << " -= " << v1 << " = " << (v0 -= v1) <<endl;
            cout << " Vector4:  " <<  v0 << " *= " << v1 << " = " << (v0 *= v1) <<endl;
            cout << " Vector4:  " <<  v0 << " /= " << v1 << " = " << (v0 /= v1) <<endl;
            cout << "----dynamics operator----" <<endl;
            cout << " Vector4:  " <<  v0 << " .dot " << v1 << " = " << (v0.Dot(v1)) <<endl;
            //cout << " Vector4:  " <<  v0 << " .cross " << v1 << " = " << (v0.cross(v1)) <<endl;
            //cout << " Vector4:  " <<  v0 << " .angle " << v1 << " = " << (v0.getAngleBetween(v1)) <<endl;
            cout << "----function operator----" <<endl;
            cout << " Vector4:  " <<  v0 << ".lengthSquare" << " = " << v0.LengthSquare() <<endl;
            cout << " Vector4:  " <<  v0 << ".length" << " = " << v0.Length() <<endl;
            cout << " Vector4:  " <<  v0 << ".inverse" << " = " << v0.GetInverse() <<endl;
            cout << " Vector4:  " <<  v0 << ".unit" << " = " << v0.GetUnit() <<endl;
            Vector4 v = v0; v.Normalize();
            cout << " Vector4:  " <<  v0 << ".normalize" << " = " << v <<endl;
            cout<<endl<<endl;

        }


        {

                cout<< "*********************************" << endl;
                cout<< "---- Unit test Complex -----" <<endl<<endl;

                Complex c0(0.5,0.1);
                Complex c1(0.1,2.0);

                cout<< c0 << " + " << c1 << " = " << c0 + c1 <<endl;
                cout<< c0 << " - " << c1 << " = " << c0 - c1 <<endl;

                cout<< c0 << " * " << c1 << " = " << c0 * c1 <<endl;
                cout<< c0 << " / " << c1 << " = " << c0 / c1 <<endl;

                cout<< c0 << " += " << c1 << " = "; (c0 += c1);
                cout<< c0 <<endl;
                cout<< c0 << " -= " << c1 << " = "; (c0 -= c1);
                cout<< c0 <<endl;

                scalar mT = 2.8;

                cout<< c0 << " + " << mT << " = " <<  c0 + mT <<endl;
                cout<< c0 << " - " << mT << " = " <<  c0 - mT <<endl;

                cout<< c0 << " * " << mT << " = " <<  c0 * mT <<endl;
                cout<< c0 << " / " << mT << " = " <<  c0 / mT <<endl;

                cout<< "-----------Function-----------" <<endl;
                cout<< c0 << ".lengthSquare = " << c0.LengthSquare() <<endl;
                cout<< c0 << ".length = " << c0.Length() <<endl;
                cout<< c0 << ".argument =" << c0.GetArgument() <<endl;
                cout<< c0 << ".conjugate =" << c0.GetConjugate() <<endl;

                cout<< "-----------Complex Function-----------" <<endl;

                cout<< c0 << ".exp = " << Complex::Exp(c0) <<endl;
                cout<< c0 << ".log = " << Complex::Log(c0) <<endl;
                cout<< c0 << ".sqrt = " << Complex::Sqrt(c0) <<endl;


                cout<< "-----------Triqonometriya Function-----------" <<endl;

                cout<< c0 << ".sin = " << Complex::Sin(c0) <<endl;
                cout<< c0 << ".cos = " << Complex::Cos(c0) <<endl;
                cout<< c0 << ".tan = " << Complex::Tan(c0) <<endl;
                cout<< c0 << ".sec = " << Complex::Sec(c0) <<endl;
                cout<< c0 << ".csc = " << Complex::Csc(c0) <<endl;
                cout<< c0 << ".cot = " << Complex::Cot(c0) <<endl;
                cout<< c0 << ".sinh = " << Complex::Sinh(c0) <<endl;
                cout<< c0 << ".cosh = " << Complex::Cosh(c0) <<endl;
                cout<< c0 << ".tanh = " << Complex::Tanh(c0) <<endl;
                cout<< c0 << ".sech = " << Complex::Sech(c0) <<endl;
                cout<< c0 << ".csch = " << Complex::Csch(c0) <<endl;
                cout<< c0 << ".coth = " << Complex::Coth(c0) <<endl;
                cout<< c0 << ".asin = " << Complex::Asin(c0) <<endl;
                cout<< c0 << ".acos = " << Complex::Acos(c0) <<endl;
                cout<< c0 << ".atan = " << Complex::Atan(c0) <<endl;
                cout<< c0 << ".asinh = " << Complex::Asinh(c0) <<endl;
                cout<< c0 << ".acosh = " << Complex::Acosh(c0) <<endl;
                cout<< c0 << ".atanh = " << Complex::Atanh(c0) <<endl;
                cout<<endl;

            }



            {

                cout<< "*********************************" << endl;
                cout<< "---- Unit test Quaternion -----" <<endl<<endl;

                Quaternion q0( 0 , 0 , 1.0 , 0.5 );
                Quaternion q1( 0.1 , 0 , 1.5 , 0.025 );

                cout<< "-----------operators-----------" <<endl;

                cout<< q0 << " +  " << q1 << " = " << (q0 + q1) <<endl;
                cout<< q0 << " -  " << q1 << " = " << (q0 - q1) <<endl;


                cout<< q0 << " += " << q1 << " = "; (q0 += q1);
                cout << q0 <<endl;
                cout<< q0 << " -= " << q1 << " = "; (q0 -= q1);
                cout << q0 <<endl;
                cout<< q0 << " *  " << q1 << " = " << (q0 * q1) <<endl;

                cout<< q0 << " *=  " << q1 << " = "; (q0 *= q1);
                cout<< q0 <<endl;

                scalar mT = 9.8;

                cout<< q0 << " *  " << mT << " = " << (q0 * mT) <<endl;
                cout<< q0 << " /  " << mT << " = " << (q0 / mT) <<endl;

                cout<< q0 << " *=  " << mT << " = "; (q0 *= mT);
                cout<< q0 <<endl;
                cout<< q0 << " /=  " << mT << " = "; (q0 /= mT);
                cout<< q0 <<endl;
                cout<<endl <<endl;


                cout<< "-----------Function-----------" <<endl;
                cout<< q0 << ".lengthSquare = " << q0.LengthSquare() <<endl;
                cout<< q0 << ".length = " << q0.Length() <<endl;
                cout<< q0 << ".conjugate =" << q0.GetConjugate() <<endl;
                cout<< q0 << ".inverse = " << q0.GetInverse() <<endl;
                cout<< q0 << ".unit = " << q0.GetUnit() <<endl;
                cout<< q0 << ".Angle =" << q0.GetAngle() <<endl;
                cout<< q0 << ".EulerAngle =" << q0.GetEulerAngles() <<endl;
                cout<< q0 << ".exponent = " << q0.Exp() <<endl;
                cout<< q0 << ".logarithm = " << q0.Log() <<endl;
                cout<< q0 << ".Matrix3x3 = " << endl;
                cout<< q0.GetRotMatrix() <<endl;
                cout<<endl;
            }



            {

                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Octonion -----" <<endl<<endl;

                    Octonion OC0(1.0 , 2.0 , 1.1 , 2.2 ,
                                 1.4 , 4.1 , 4.1 , 4.1);

                    Octonion OC1(1.4 , 1.0 , 2.1 , 2.2 ,
                                 2.4 , 1.1 , 3.1 , 1.1);


                    scalar t = 2.3;

                    Complex c( 1.2 , 2.1);

                    Quaternion q( 1.0 , 1.1 , 1.4 , 4.4 );

                    cout<< "-----------------operators-------------------" <<endl;

                    cout<< OC0 << " + " << t << " = " << (OC0 + t) << endl;
                    cout<< OC0 << " + " << c << " = " << (OC0 + c) << endl;
                    cout<< OC0 << " + " << q << " = " << (OC0 + q) << endl;
                    cout<< OC0 << " + " << OC1 << " = " << (OC0 + OC1) << endl;

                    cout<< " ---------------------------------------- " <<endl;

                    cout<< OC0 << " - " << t << " = " << (OC0 - t) << endl;
                    cout<< OC0 << " - " << c << " = " << (OC0 - c) << endl;
                    cout<< OC0 << " - " << q << " = " << (OC0 - q) << endl;
                    cout<< OC0 << " - " << OC1 << " = " << (OC0 - OC1) << endl;

                    cout<< " ---------------------------------------- " <<endl;

                    cout<< OC0 << " * " << t << " = " << (OC0 * t) << endl;
                    cout<< OC0 << " * " << c << " = " << (OC0 * c) << endl;
                    cout<< OC0 << " * " << q << " = " << (OC0 * q) << endl;
                    cout<< OC0 << " * " << OC1 << " = " << (OC0 * OC1) << endl;


                    cout<< " ---------------------------------------- " <<endl;

                    cout<< OC0 << " / " << t << " = " << (OC0 / t) << endl;
                    cout<< OC0 << " / " << c << " = " << (OC0 / c) << endl;
                    cout<< OC0 << " / " << q << " = " << (OC0 / q) << endl;
                    cout<< OC0 << " / " << OC1 << " = " << (OC0 / OC1) << endl;

                    cout<< "-----------------Function-------------------" <<endl;

                    cout << OC0 << " .lenght = " << OC0.Length() <<endl;
                    cout << OC0 << " .lenghtSquare = " << OC0.LengthSquare() <<endl;
                    cout << OC0 << " .norm = " << OC0.Norm() <<endl;
                    cout << OC0 << " .abs = " << OC0.Abs() <<endl;
                    cout << OC0 << " .conjugate = " << OC0.GetConjugate() <<endl;
                    cout << OC0 << " .exponenta = " << Octonion::Exp(OC0) <<endl;

                    cout << OC0 << " .cosinus = " << Octonion::Cos(OC0) <<endl;
                    cout << OC0 << " .sininus = " << Octonion::Sin(OC0) <<endl;
                    cout << OC0 << " .tanget  = " << Octonion::Tan(OC0) <<endl;

                    cout << OC0 << " .hyperbolic cosinus = " << Octonion::Cosh(OC0) <<endl;
                    cout << OC0 << " .hyperbolic sininus = " << Octonion::Sinh(OC0) <<endl;
                    cout << OC0 << " .hyperbolic tanget  = " << Octonion::Tanh(OC0) <<endl;

                }




                {

                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Matrix_2x2 -----" <<endl<<endl;

                    Matrix2 m0( -2.0 , 7.0 ,
                                -4.0 , 2.0);
                    Matrix2 m1( 3.0 , 1.0 ,
                                2.0 , 4.0);



                    cout << "------ operator + --------" <<endl<<endl;
                    cout<< m0 << " + " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 + m1) <<endl<<endl;

                    cout << "------ operator - --------" <<endl<<endl;
                    cout<< m0 << " - " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 - m1) <<endl<<endl;

                    cout << "------ operator += -------" <<endl<<endl;
                    cout<< m0 << " += " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 += m1) <<endl<<endl;

                    cout << "------ operator -= -------" <<endl<<endl;
                    cout<< m0 << " -= " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 -= m1) <<endl<<endl;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1.GetInverse()) <<endl<<endl;


                    scalar mT = 8.9;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mT) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 / mT) <<endl<<endl;


                    cout << "------ operator *= -------" <<endl<<endl;
                    cout<< m0 << " *= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 *= mT) <<endl<<endl;


                    cout << "------ operator /= -------" <<endl<<endl;
                    cout<< m0 << " /= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 /= mT) <<endl<<endl;


                    Vector2 mV(2,4);

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mV) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0.GetInverse() * mV) <<endl<<endl;


                    cout << "------ Matrix Function -------" <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .determinat = " << (m0.GetDeterminant()) <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .trace = " << (m0.GetTrace()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .inverse = " <<endl;
                    cout<< (m0.GetInverse()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .transpose = " <<endl;
                    cout<< (m0.GetTranspose()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .absolute = " <<endl;
                    cout<< (m0.GetAbsoluteMatrix()) <<endl<<endl;

                }


                {

                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Matrix_3x3 -----" <<endl<<endl;

                    Matrix3 m0( -2.0 , 7.0 , 1.1 ,
                                -4.0 , 2.0 , 2.5 ,
                                2.0 , 1.0 , 1.0);

                    Matrix3 m1( 3.0 , 1.0 , 1.1 ,
                                2.0 , 4.0 , 2.1 ,
                                1.0 , 5.0 , 1.0);



                    cout << "------ operator + --------" <<endl<<endl;
                    cout<< m0 << " + " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 + m1) <<endl<<endl;

                    cout << "------ operator - --------" <<endl<<endl;
                    cout<< m0 << " - " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 - m1) <<endl<<endl;

                    cout << "------ operator += -------" <<endl<<endl;
                    cout<< m0 << " += " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 += m1) <<endl<<endl;

                    cout << "------ operator -= -------" <<endl<<endl;
                    cout<< m0 << " -= " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 -= m1) <<endl<<endl;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1.GetInverse()) <<endl<<endl;


                    scalar mT = 8.9;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mT) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 / mT) <<endl<<endl;


                    cout << "------ operator *= -------" <<endl<<endl;
                    cout<< m0 << " *= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 *= mT) <<endl<<endl;


                    cout << "------ operator /= -------" <<endl<<endl;
                    cout<< m0 << " /= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 /= mT) <<endl<<endl;


                    Vector3 mV(2,4,4);

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mV) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0.GetInverse() * mV) <<endl<<endl;


                    cout << "------ Matrix Function -------" <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .determinat = " << (m0.GetDeterminant()) <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .trace = " << (m0.GetTrace()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .inverse = " <<endl;
                    cout<< (m0.GetInverse()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .transpose = " <<endl;
                    cout<< (m0.GetTranspose()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .absolute = " <<endl;
                    cout<< (m0.GetAbsoluteMatrix()) <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .diagonalize = " <<endl;
                    cout<< (m0.GetDiagonalize(0.01 , 10)) <<endl<<endl;

                }



                {

                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Matrix_4x4 -----" <<endl<<endl;

                    Matrix4 m0( -2.0 , 7.0 , 1.1 , 2.2 ,
                                -4.0 , 2.0 , 2.5 , 3.3 ,
                                2.0 , 1.0 , 1.0 , 2.3 ,
                                1.0 , 2.0 , 4.0 , 8.8);

                    Matrix4 m1( 3.0 , 1.0 , 1.1 , 2.0 ,
                                2.0 , 4.0 , 2.1 , 1.0 ,
                                1.0 , 5.0 , 1.0 , 2.2 ,
                                1.0 , 2.0 , 2.0 , 3.3);



                    cout << "------ operator + --------" <<endl<<endl;
                    cout<< m0 << " + "  <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 + m1) <<endl<<endl;

                    cout << "------ operator - --------" <<endl<<endl;
                    cout<< m0 << " - " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 - m1) <<endl<<endl;

                    cout << "------ operator += -------" <<endl<<endl;
                    cout<< m0 << " += " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 += m1) <<endl<<endl;

                    cout << "------ operator -= -------" <<endl<<endl;
                    cout<< m0 << " -= " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 -= m1) <<endl<<endl;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< m1 << " = "  <<endl;
                    cout<< (m0 * m1.GetInverse()) <<endl<<endl;


                    scalar mT = 8.9;

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mT) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 / mT) <<endl<<endl;


                    cout << "------ operator *= -------" <<endl<<endl;
                    cout<< m0 << " *= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 *= mT) <<endl<<endl;


                    cout << "------ operator /= -------" <<endl<<endl;
                    cout<< m0 << " /= " <<endl ;
                    cout<< mT <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 /= mT) <<endl<<endl;


                    Vector3 mV(2,4,4);

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * mV) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0.GetInverse() * mV) <<endl<<endl;


                    Vector4 _mV(2,4,4 ,2.3);

                    cout << "------ operator * -------" <<endl<<endl;
                    cout<< m0 << " * " <<endl ;
                    cout<< _mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0 * _mV) <<endl<<endl;


                    cout << "------ operator / -------" <<endl<<endl;
                    cout<< m0 << " / " <<endl ;
                    cout<< _mV <<endl;
                    cout<<  " = "  <<endl;
                    cout<< (m0.GetInverse() * _mV) <<endl<<endl;


                    cout << "------ Matrix Function -------" <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .determinat = " << (m0.GetDeterminant()) <<endl<<endl;

                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .trace = " << (m0.GetTrace()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .inverse = " <<endl;
                    cout<< (m0.GetInverse()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .transpose = " <<endl;
                    cout<< (m0.GetTranspose()) <<endl<<endl;


                    cout << "------------------------------" <<endl<<endl;
                    cout<< m0  <<endl ;
                    cout<<  " .absolute = " <<endl;
                    cout<< (m0.GetAbsoluteMatrix()) <<endl<<endl;

                }



                {
                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Line3 -----" <<endl<<endl;

                    Vector3 A(1,2,3);
                    Vector3 B(9,6,3);
                    Vector3 C(0,9,4);
                    Vector3 D(0,0,1);


                    Line3 line_a(A,B-A);
                    Line3 line_b(C,D-C);

                    cout<<  " line_line_a " << line_a <<endl;
                    cout<<  " line_line_b " << line_b <<endl;

                    Vector3 cross_a;
                    Vector3 cross_b;
                    Line3::ClosestPoints(line_a , line_b , cross_a , cross_b);

                    cout<< endl;
                    cout<< " cross pointA in lineA to lineB: = " << cross_a <<endl;
                    cout<< " cross pointB in lineA to lineB: = " << cross_b <<endl;

                    cout<<endl;
                    cout<< " closset point B to point in lineA " << line_a.ClosestPoint(B) << endl;
                }



                {
                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test LineSegment3 -----" <<endl<<endl;

                    Vector3 A(1,2,3);
                    Vector3 B(9,6,3);
                    Vector3 C(0,9,4);
                    Vector3 D(0,0,1);

                    LineSegment3 line_a(A,B);
                    LineSegment3 line_b(C,D);

                    cout<<  " line_segment_a " << line_a <<endl;
                    cout<<  " line_segment_b " << line_b <<endl;

                    Vector3 cross_a;
                    Vector3 cross_b;
                    LineSegment3::ClosestPoints(line_a , line_b , cross_a , cross_b);

                    cout<< endl;
                    cout<< " cross pointA in segmetnA to segmentB: = " << cross_a <<endl;
                    cout<< " cross pointB in segmetnA to segmentB: = " << cross_b <<endl;

                    cout<<endl;
                    cout<< " closset point B to point in segmentA " << line_a.ClosestPoint(B) << endl;
                }


                {
                    cout<< "*********************************" << endl;
                    cout<< "---- Unit test Euler Angle in Matrix3x3 -----" <<endl<<endl;

                    Vector3 eulerAngle( 0.5 , 0.7 , 0.1);

                    cout<< " ---------- push angle -------- " << endl;
                    cout<< "Set Angle X: " << eulerAngle.x << endl;
                    cout<< "Set Angle Y: " << eulerAngle.y << endl;
                    cout<< "Set Angle Z: " << eulerAngle.z << endl;
                    cout<< endl;
                    cout<< " ---------- extract angle -------- " << endl;


                    Matrix3 M;
                    M = M.CreateRotationEulerAngle( eulerAngle.x ,
                                                    eulerAngle.y ,
                                                    eulerAngle.z );


                    // eulerAngle =  Q.GetEulerAngleGimbalLock(Quaternion::RotSeq::xyz);
                    eulerAngle = Quaternion(M).GetRotMatrix().GetEulerAngles();

                    cout<< "Get Angle X: " << eulerAngle.x << endl;
                    cout<< "Get Angle Y: " << eulerAngle.y << endl;
                    cout<< "Get Angle Z: " << eulerAngle.z << endl;
                }


    return 0;
}
