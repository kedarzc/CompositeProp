#include<iostream>
#include<Eigen/Dense>
#include<math.h>
#include<vector>

#include "LamUtils.h"


using namespace std;
using namespace Eigen;

int main() {

	vector<double> thetaLaminaDegrees;

	// Add angles
	thetaLaminaDegrees.push_back(0.0);
	thetaLaminaDegrees.push_back(90.0);
	thetaLaminaDegrees.push_back(30.0);
	thetaLaminaDegrees.push_back(-30.0);
	thetaLaminaDegrees.push_back(-30.0);
	thetaLaminaDegrees.push_back(30.0);
	thetaLaminaDegrees.push_back(90.0);
	thetaLaminaDegrees.push_back(0.0);

	double t = 6, tk = t/thetaLaminaDegrees.size();
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
		S = T.transpose() * SDash * T;

		// Calculate the Stiffness
		M.CalcStiffnessMat(S, CDash);

		CLam = CLam + (tk / t) * CDash;
		//M.AddLamTransfMat(CLam, CDash, t, tk);
	}


	SLam = CLam.inverse();

	cout << 1.0/  SLam(0,0) << endl;
	cout << 1.0 / SLam(1, 1) << endl;
	cout << 1.0 / SLam(2, 2) << endl;
	
}