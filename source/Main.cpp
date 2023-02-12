#include<iostream>
#include<Eigen/Dense>
#include<math.h>
#include "LamUtils.h"

using namespace std;
using namespace Eigen;

int main() {

	//double thetaLaminaDegrees = 60.0;
	//Matrix3d A;
	Matrix <double, 6, 6> S, C;
	LamUtils M;

	//M.CalcLamTransfMat(thetaLaminaDegrees, A);
	//M.CalcRotMat(A, T);

	//Composite properties
	M.CalcComplianceMat(145880, 13312, 13312, 4386, 4386, 4528, 0.263, 0.263, 0.470,S);
	
	// Calculate the Stiffness
	M.CalcStiffnessMat(S, C);

	cout << C << endl;
	
}