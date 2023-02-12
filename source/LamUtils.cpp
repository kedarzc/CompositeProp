#include "LamUtils.h"


// Computes the direction cosine of the lamina system
void LamUtils::CalcLamTransfMat(double angle, Matrix3d& a)
{
	// Declare the value of pi
	const double pi = 3.141592653589793238;
		
	// Set all elements to zero
	a.setZero();
	
	// Fill in all the non-zero components
	a(0, 0) = cos(angle  * (pi / 180.0));
	a(0, 1) =  sin(angle  * (pi / 180.0));
	a(1, 0) = -sin(angle  * (pi / 180.0));
	a(1, 1) =  cos(angle  * (pi / 180.0));
	a(2, 2) = 1.0;
}

void LamUtils::CalcRotMat(const Matrix3d a, Matrix <double, 6, 6>& T)
{	
	// First row
	T(0, 0) = a(0, 0) * a(0, 0);
	T(0, 1) = a(0, 1) * a(0, 1);
	T(0, 2) = a(0, 2) * a(0, 2);
	T(0, 3) = 2.0 * a(0, 1) * a(0, 2);
	T(0, 4) = 2.0 * a(0, 0) * a(0, 2);
	T(0, 5) = 2.0 * a(0, 0) * a(0, 1);

	// Second row
	T(1, 0) = a(1, 0) * a(1, 0);
	T(1, 1) = a(1, 1) * a(1, 1);
	T(1, 2) = a(1, 2) * a(1, 2);
	T(1, 3) = 2.0 * a(1, 1) * a(1, 2);
	T(1, 4) = 2.0 * a(1, 0) * a(1, 2);
	T(1, 5) = 2.0 * a(1, 0) * a(1, 1);
	
	// Third row
	T(2, 0) = a(2, 0) * a(2, 0);
	T(2, 1) = a(2, 1) * a(2, 1);
	T(2, 2) = a(2, 2) * a(2, 2);
	T(2, 3) = 2.0 * a(2, 1) * a(2, 2);
	T(2, 4) = 2.0 * a(2, 0) * a(2, 2);
	T(2, 5) = 2.0 * a(2, 0) * a(2, 1);
	
	// Fourth row
	T(3, 0) = a(1, 0) * a(2, 0);
	T(3, 1) = a(1, 1) * a(2, 1);
	T(3, 2) = a(1, 2) * a(2, 2);
	T(3, 3) = a(1, 1) * a(2, 2) + a(1, 2) * a(2, 1);
	T(3, 4) = a(1, 0) * a(2, 2) + a(1, 2) * a(2, 0);
	T(3, 5) = a(1, 0) * a(2, 1) + a(1, 1) * a(2, 0);
	
	// Fifth row
	T(4, 0) = a(0, 0) * a(2, 0);
	T(4, 1) = a(0, 1) * a(2, 1);
	T(4, 2) = a(0, 2) * a(2, 2);
	T(4, 3) = a(0, 1) * a(2, 2) + a(0, 2) * a(2, 1);
	T(4, 4) = a(0, 0) * a(2, 2) + a(0, 2) * a(2, 0);
	T(4, 5) = a(0, 0) * a(2, 1) + a(0, 1) * a(2, 0);
	//
	// Sixth row
	T(5, 0) = a(0, 0) * a(1, 0);
	T(5, 1) = a(0, 1) * a(1, 1);
	T(5, 2) = a(0, 2) * a(1, 2);
	T(5, 3) = a(0, 1) * a(1, 2) + a(0, 2) * a(1, 1);
	T(5, 4) = a(0, 0) * a(1, 2) + a(0, 2) * a(1, 0);
	T(5, 5) = a(0, 0) * a(1, 1) + a(0, 1) * a(1, 0);
}

void LamUtils::CalcComplianceMat(const double E1, const double E2, const double E3, const double G12, const double G13, const double G23, 
	const double nu12, const double nu13, const double nu23, Matrix <double, 6, 6>& S)
{
	// First initialize the compliance matrix
	S.setZero();

	// Dialgonal Terms
	S(0, 0) = 1.0 / E1;
	S(1, 1) = 1.0 / E2;
	S(2, 2) = 1.0 / E3;
	S(3, 3) = 1.0 / G23;
	S(4, 4) = 1.0 / G13;
	S(5, 5) = 1.0 / G12;

	// Non-diagonal terms
	S(0, 1) = -nu12 / E2;
	S(0, 2) = -nu13 / E3;
	S(1, 0) = -nu12 / E1;
	S(1, 2) = -nu23 / E3;
	S(2, 0) = -nu13 / E1;
	S(2, 1) = -nu23 / E2;
}

void LamUtils::CalcStiffnessMat(const Matrix <double, 6, 6>& S, Matrix <double, 6, 6>& C)
{
	// First initialize the stiffness matrix
	C.setZero();

	C = S.inverse();

}

void LamUtils::AddLamTransfMat(Matrix <double, 6, 6>& stiff, Matrix <double, 6, 6>& stiffkth, const double t, const double tk)
{
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			stiff(i, j) = stiff(i, j) + (tk / t) * stiffkth(i, j);
		}
	}
}