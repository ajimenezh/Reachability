#pragma once

#include "Matrix2D.h"

#include <CGAL/Polygon_2.h>
//#include "VPolytope.hh"

#include "Zonotope.h"
#include "Scenario.h"
#include <boost/numeric/ublas/tensor.hpp>

//void assign(Vector& res, std::vector<double> v);

std::vector< Polygon_2 > CalcReachability(const Polygon_2& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu);

std::vector< Polygon_2 > StateDependantReachability(const Polygon_2& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu);

//std::vector< VPolytope > Reachability_2(const VPolytope& initial_value, double t,
//	double delta_t, const Matrix& a, double mu);

std::vector< Zonotope > ReachabilityLinearSystem(const Zonotope& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu);

Zonotope ReachabilityLinearSystemT(const Zonotope& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu);

std::vector< Zonotope > ReachabilityFromInput(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, double u0, double mu);

std::vector< Zonotope > ReachabilityFromState(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& f_star, double mu);

Zonotope ReachabilityFromStateT(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& f_star, double mu);

std::vector< Zonotope > CalcReachability(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double mu);

Mesh ZonotopesToMatrix(const std::vector< Zonotope >& zonotopes, const Mesh& mesh, bool join_polygons);

void PolygonToMesh(const polygon& pol, Mesh& tmp);

class Reachability {
public:
	Reachability(const Zonotope& initial_value, double T,
		double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double u0, double mu);

	Reachability() {}

	void UpdateReachabilityInput(const Zonotope& initial_value, double T,
		double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double u0, double mu);

	std::vector< Zonotope > R_;
	Zonotope R_t_;
	std::vector< Zonotope > R_bar_;
	Zonotope R_bar_t_;
	std::vector< Zonotope > R_up_;
	Zonotope R_up_t_;
	std::vector< Zonotope > R_down_;
	Zonotope R_down_t_;
};
