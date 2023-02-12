#include<math.h>
#include<Eigen/Dense>


using namespace std;
using namespace Eigen;

class LamUtils
{
	public:
		void CalcLamTransfMat(double angle, Matrix3d& a);
		void CalcRotMat(const Matrix3d a, Matrix <double, 6, 6>& T);
		void CalcComplianceMat(const double E1, const double E2, const double E3, const double G12, const double G13, const double G23,
			const double nu12, const double nu13, const double nu23, Matrix <double, 6, 6>& S);
		void CalcStiffnessMat(const Matrix<double, 6, 6>& S, Matrix<double, 6, 6>& C);

		void AddLamTransfMat(Matrix<double, 6, 6>& stiff, Matrix<double, 6, 6>& stiffkth, const double t, const double tk);

};

