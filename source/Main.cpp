#include<iostream>
#include<Eigen/Dense>
#include<math.h>
#include<vector>

#include "LamUtils.h"


using namespace std;
using namespace Eigen;

int main() {

	vector<double> thetaLaminaDegrees{0.0,90.0};

	double tk = 1.5;
	double t = tk * size(thetaLaminaDegrees);
	LamUtils M;
	Matrix3d A;
	Matrix <double, 6, 6> T;
	Matrix <double, 6, 6> SDash, S, SLam, CDash, CLam;

	//Initialize C matrix
	CLam.setZero();

	// Calculate S'
	M.CalcComplianceMat(145880, 13312, 13312, 4386, 4386, 4528, 0.263, 0.263, 0.470, SDash);

	for (int k = 0; k < thetaLaminaDegrees.size(); k++)
	{
		// Calculate T
		M.CalcLamTransfMat(thetaLaminaDegrees[k], A);

		M.CalcRotMat(A, T);

		// Calculate rotated S from Sdash and T
		S = (T.transpose() * SDash) * T;

		// Calculate the Stiffness
		M.CalcStiffnessMat(S, CDash);

		CLam = CLam + (tk / t) * CDash;
		//M.AddLamTransfMat(CLam, CDash, t, tk);
	}

	SLam = CLam.inverse();

	//cout << endl << CLam ;
	//cout << endl << t << endl << tk << endl;
	//cout << endl << SLam;
	
	cout << "Laminate properties" << endl;
	cout << "===================" << endl;
	cout << "EXX = " << 1.0 / SLam(0, 0) << endl;
	cout << "EYY = " << 1.0 / SLam(1, 1) << endl;
	cout << "EZZ = " << 1.0 / SLam(2, 2) << endl;

	cout << "GXY = " << 1.0 / SLam(5, 5) << endl;
	cout << "GXZ = " << 1.0 / SLam(4, 4) << endl;
	cout << "GYZ = " << 1.0 / SLam(3, 3) << endl;

	cout << "nuXY = " << -SLam(0, 1)*(1.0 / SLam(1, 1)) << endl;
	cout << "nuXZ = " << -SLam(0, 2)*(1.0 / SLam(2, 2)) << endl;
	cout << "nuYZ = " << -SLam(1, 2)*(1.0 / SLam(2, 2)) << endl;
}