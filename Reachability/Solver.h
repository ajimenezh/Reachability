#pragma once

#include <vector>

class Reachability;
class Matrix2D;
class Mesh;

class Solver {
public:
	void Solve();

	std::vector<Reachability> CalcReachability(Mesh mesh, double T,
		double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double u0, double mu);

	void UpdateReachabilityByInput(std::vector<Reachability>& reachability, Mesh mesh, double T,
		double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& F, double u0, double mu);

	double s = 15.0;
	double c1 = 160.0;
	double c2 = -1.6;
	double c3 = 53.0;
	double c4 = -3.5;
	double c5 = 156.0;
	double c6 = 78.0;
};

