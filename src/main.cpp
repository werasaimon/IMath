#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include "maths/maths.h"

using namespace CMath;


int main(void)
{
	/**

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
		cout << " Vector2:  " <<  v0 << " .dot " << v1 << " = " << (v0.dot(v1)) <<endl;
		cout << " Vector2:  " <<  v0 << " .cross " << v1 << " = " << (v0.cross(v1)) <<endl;
		cout << " Vector2:  " <<  v0 << " .angle " << v1 << " = " << (v0.getAngleBetween(v1)) <<endl;
		cout << "----function operator----" <<endl;
		cout << " Vector2:  " <<  v0 << ".lengthSquare" << " = " << v0.lengthSquare() <<endl;
		cout << " Vector2:  " <<  v0 << ".length" << " = " << v0.length() <<endl;
		cout << " Vector2:  " <<  v0 << ".inverse" << " = " << v0.getInverse() <<endl;
		cout << " Vector2:  " <<  v0 << ".unit" << " = " << v0.getUnit() <<endl;
		Vector2 v = v0; v.normalize();
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
		cout << " Vector3:  " <<  v0 << " .dot " << v1 << " = " << (v0.dot(v1)) <<endl;
		cout << " Vector3:  " <<  v0 << " .cross " << v1 << " = " << (v0.cross(v1)) <<endl;
		cout << " Vector3:  " <<  v0 << " .angle " << v1 << " = " << (v0.getAngleBetween(v1)) <<endl;
		cout << "----function operator----" <<endl;
		cout << " Vector3:  " <<  v0 << ".lengthSquare" << " = " << v0.lengthSquare() <<endl;
		cout << " Vector3:  " <<  v0 << ".length" << " = " << v0.length() <<endl;
		cout << " Vector3:  " <<  v0 << ".inverse" << " = " << v0.getInverse() <<endl;
		cout << " Vector3:  " <<  v0 << ".unit" << " = " << v0.getUnit() <<endl;
		Vector3 v = v0; v.normalize();
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
		cout << " Vector4:  " <<  v0 << " .dot " << v1 << " = " << (v0.dot(v1)) <<endl;
		//cout << " Vector4:  " <<  v0 << " .cross " << v1 << " = " << (v0.cross(v1)) <<endl;
		//cout << " Vector4:  " <<  v0 << " .angle " << v1 << " = " << (v0.getAngleBetween(v1)) <<endl;
		cout << "----function operator----" <<endl;
		cout << " Vector4:  " <<  v0 << ".lengthSquare" << " = " << v0.lengthSquare() <<endl;
		cout << " Vector4:  " <<  v0 << ".length" << " = " << v0.length() <<endl;
		cout << " Vector4:  " <<  v0 << ".inverse" << " = " << v0.getInverse() <<endl;
		cout << " Vector4:  " <<  v0 << ".unit" << " = " << v0.getUnit() <<endl;
		Vector4 v = v0; v.normalize();
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
		cout<< c0 << ".lengthSquare = " << c0.lengthSquare() <<endl;
		cout<< c0 << ".length = " << c0.length() <<endl;
		cout<< c0 << ".argument =" << c0.getArgument() <<endl;
		cout<< c0 << ".conjugate =" << c0.getConjugate() <<endl;

		cout<< "-----------Complex Function-----------" <<endl;

		cout<< c0 << ".exp = " << Complex::exp(c0) <<endl;
		cout<< c0 << ".log = " << Complex::log(c0) <<endl;
		cout<< c0 << ".sqrt = " << Complex::sqrt(c0) <<endl;


		cout<< "-----------Triqonometriya Function-----------" <<endl;

		cout<< c0 << ".sin = " << Complex::sin(c0) <<endl;
		cout<< c0 << ".cos = " << Complex::cos(c0) <<endl;
		cout<< c0 << ".tan = " << Complex::tan(c0) <<endl;
		cout<< c0 << ".sec = " << Complex::sec(c0) <<endl;
		cout<< c0 << ".csc = " << Complex::csc(c0) <<endl;
		cout<< c0 << ".cot = " << Complex::cot(c0) <<endl;
		cout<< c0 << ".sinh = " << Complex::sinh(c0) <<endl;
		cout<< c0 << ".cosh = " << Complex::cosh(c0) <<endl;
		cout<< c0 << ".tanh = " << Complex::tanh(c0) <<endl;
		cout<< c0 << ".sech = " << Complex::sech(c0) <<endl;
		cout<< c0 << ".csch = " << Complex::csch(c0) <<endl;
		cout<< c0 << ".coth = " << Complex::coth(c0) <<endl;
		cout<< c0 << ".asin = " << Complex::asin(c0) <<endl;
		cout<< c0 << ".acos = " << Complex::acos(c0) <<endl;
		cout<< c0 << ".atan = " << Complex::atan(c0) <<endl;
		cout<< c0 << ".asinh = " << Complex::asinh(c0) <<endl;
		cout<< c0 << ".acosh = " << Complex::acosh(c0) <<endl;
		cout<< c0 << ".atanh = " << Complex::atanh(c0) <<endl;
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
			cout<< q0 << ".lengthSquare = " << q0.lengthSquare() <<endl;
			cout<< q0 << ".length = " << q0.length() <<endl;
			cout<< q0 << ".conjugate =" << q0.getConjugate() <<endl;
			cout<< q0 << ".inverse = " << q0.getInverse() <<endl;
			cout<< q0 << ".unit = " << q0.getUnit() <<endl;
			cout<< q0 << ".Angle =" << q0.getAngle() <<endl;
			cout<< q0 << ".EulerAngle =" << q0.getEulerAngles() <<endl;
			cout<< q0 << ".exponent = " << q0.exp() <<endl;
			cout<< q0 << ".logarithm = " << q0.log() <<endl;
			cout<< q0 << ".Matrix3x3 = " << endl;
			cout<< q0.getMatrix() <<endl;
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

			cout << OC0 << " .lenght = " << OC0.length() <<endl;
			cout << OC0 << " .lenghtSquare = " << OC0.lengthSquare() <<endl;
			cout << OC0 << " .norm = " << OC0.norm() <<endl;
			cout << OC0 << " .abs = " << OC0.abs() <<endl;
			cout << OC0 << " .conjugate = " << OC0.getConjugate() <<endl;
			cout << OC0 << " .exponenta = " << Octonion::exp(OC0) <<endl;

			cout << OC0 << " .cosinus = " << Octonion::cos(OC0) <<endl;
			cout << OC0 << " .sininus = " << Octonion::sin(OC0) <<endl;
			cout << OC0 << " .tanget  = " << Octonion::tan(OC0) <<endl;

			cout << OC0 << " .hyperbolic cosinus = " << Octonion::cosh(OC0) <<endl;
			cout << OC0 << " .hyperbolic sininus = " << Octonion::sinh(OC0) <<endl;
			cout << OC0 << " .hyperbolic tanget  = " << Octonion::tanh(OC0) <<endl;

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
		cout<< (m0 * m1.getInverse()) <<endl<<endl;


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
		cout<< (m0.getInverse() * mV) <<endl<<endl;


		cout << "------ Matrix Function -------" <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .determinat = " << (m0.getDeterminant()) <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .trace = " << (m0.getTrace()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .inverse = " <<endl;
		cout<< (m0.getInverse()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .transpose = " <<endl;
		cout<< (m0.getTranspose()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .absolute = " <<endl;
		cout<< (m0.getAbsoluteMatrix()) <<endl<<endl;


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
		cout<< (m0 * m1.getInverse()) <<endl<<endl;


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
		cout<< (m0.getInverse() * mV) <<endl<<endl;


		cout << "------ Matrix Function -------" <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .determinat = " << (m0.getDeterminant()) <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .trace = " << (m0.getTrace()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .inverse = " <<endl;
		cout<< (m0.getInverse()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .transpose = " <<endl;
		cout<< (m0.getTranspose()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .absolute = " <<endl;
		cout<< (m0.getAbsoluteMatrix()) <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .diagonalize = " <<endl;
		cout<< (m0.getDiagonalize(0.01 , 10)) <<endl<<endl;


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
		cout<< (m0 * m1.getInverse()) <<endl<<endl;


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
		cout<< (m0.getInverse() * mV) <<endl<<endl;


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
		cout<< (m0.getInverse() * _mV) <<endl<<endl;


		cout << "------ Matrix Function -------" <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .determinat = " << (m0.getDeterminant()) <<endl<<endl;

		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .trace = " << (m0.getTrace()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .inverse = " <<endl;
		cout<< (m0.getInverse()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .transpose = " <<endl;
		cout<< (m0.getTranspose()) <<endl<<endl;


		cout << "------------------------------" <<endl<<endl;
		cout<< m0  <<endl ;
		cout<<  " .absolute = " <<endl;
		cout<< (m0.getAbsoluteMatrix()) <<endl<<endl;


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


    /**/


	return EXIT_SUCCESS;
}
