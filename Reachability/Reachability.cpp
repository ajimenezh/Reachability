#include "Reachability.h"

#include <CGAL/minkowski_sum_2.h>
#include "bops_linear.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Line_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
//#include "VPolytopeList.hh"

#include <numeric>
#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <fstream>

#include "Zonotope.h"

//#include <CGAL/Gmpzf.h>

#include "matrix_exponential.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

const int MAX_VERTICES = 8;

double Sup(const Polygon_2& polygon) {
	double res = 0.0;
	for (const auto& vertex : polygon) {
		const double& x = to_double(vertex.x().exact());
		const double& y = to_double(vertex.y().exact());
		res = std::max(res, sqrt(x * x + y * y));
	}
	return res;
}

//double Sup(const VPolytope& polytope) {
//	double res = 0.0;
//	for (int i = 0; i < polytope.rowdim(); i++) {
//		const double& x = polytope[0][i].get_d();
//		const double& y = polytope[1][i].get_d();
//		res = std::max(res, sqrt(x * x + y * y));
//	}
//	return res;
//}

Polygon_2 Multiply(const Polygon_2& a, double h) {
	double res = 0.0;
	Polygon_2 b;
	for (const auto& vertex : a) {
		const double& x = to_double(vertex.x().exact());
		const double& y = to_double(vertex.y().exact());
		b.push_back(Point_2(x * h, y * h));
	}
	return b;
}

//VPolytope Transform(const Reachability::Matrix& m, const VPolytope& polytope) {
//	::Matrix res(polytope.rowdim(), polytope.coldim());

//	for (int i = 0; i < polytope.coldim(); i++) {
//		for (int j = 0; j < polytope.rowdim(); j++) {
//			res[i][j] = 0.0;
//			for (int k = 0; k < polytope.rowdim(); k++) {
//				res[i][j] += m[j][k] * polytope[i][k].get_d();
//			}
//		}
//	}

//	return res;
//}

//VPolytope Scale(const VPolytope& polytope, double val) {
//	VPolytope res = polytope;

//	for (int i = 0; i < polytope.rowdim(); i++) {
//		for (int j = 0; j < polytope.coldim(); j++) {
//			res[j][i] *= val;
//		}
//	}

//	return res;
//}

Polygon_2 CreateHyperrectangle(double x) {
	Polygon_2   P;
	P.push_back(Point_2(x, -x));
	P.push_back(Point_2(x, x));
	P.push_back(Point_2(-x, x));
	P.push_back(Point_2(-x, -x));
	return P;
}

//Vector createVector(std::vector<double> v) {
//	Vector res(v.size());
//	for (int i = 0; i < v.size(); i++) {
//		res[i] = v[i];
//	}
//	return res;
//}

//void assign(Vector& res, std::vector<double> v) {
//	for (int i = 0; i < v.size(); i++) {
//		res[i] = v[i];
//	}
//}

//VPolytope CreateHyperrectanglePolytope(double x) {
//	::Matrix m(2, 4);
//	assign(m[0], { x, -x });
//	assign(m[1], { x, x });
//	assign(m[2], { -x, x });
//	assign(m[3], { -x, -x });
//	return m;
//}

Polygon_2 MinkowskiSum(const Polygon_2& a, const Polygon_2& b) {
	return CGAL::minkowski_sum_2(a, b).outer_boundary();
}

//VPolytope MinkowskiSum(const VPolytope& a, const VPolytope& b) {
//	VPolytopeList polylist;
//	polylist.append(a);
//	polylist.append(b);
//	polylist.incMinkSumSetup();

//	std::vector< Vector > v;
//	while (polylist.hasMoreVertices()) {
//		Vector tmp = polylist.incExploreStepFast().coord();
//		for (int i = 0; i < tmp.size(); i++) {
//			std::cout << tmp[i].get_d() << std::endl;
//		}
//		v.push_back(tmp);
//	}
//	std::cout << "---" << std::endl;


//	int n = 2;
//	::Matrix m(n, v.size());
//	for (int i = 0; i < v.size(); i++) {
//		for (int j = 0; j < n; j++) {
//			m[i][j] = v[i][j];
//		}
//	}

//	return m;
//}

Polygon_2 ConvexHullOfPolygons(const std::vector< Polygon_2 >& v) {
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Convex_hull_traits_adapter_2<K,
		CGAL::Pointer_property_map<K::Point_2>::type > Convex_hull_traits_2;

	
	std::vector<K::Point_2> points;
	for (const auto& pol : v) {
		for (const auto& vertex : pol) {
			const double& x = to_double(vertex.x().exact());
			const double& y = to_double(vertex.y().exact());
			points.push_back(K::Point_2(x, y));
		}
	}

	std::vector<std::size_t> indices(points.size()), out;
	std::iota(indices.begin(), indices.end(), 0);
	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
		Convex_hull_traits_2(CGAL::make_property_map(points)));

	Polygon_2 res;
	for (std::size_t i : out) {
		const double& x = points[i].x();
		const double& y = points[i].y();
		res.push_back(Point_2(x, y));
	}
	return res;
}

CGAL::Point_2< CGAL::Cartesian<double> > ToPoint(CGAL::Point_2< CGAL::Epeck > p) {
	const double& x = to_double(p.x().exact());
	const double& y = to_double(p.y().exact());
	return CGAL::Point_2< CGAL::Cartesian<double> >(x, y);
}

Polygon_2 MinParalelogram(const Polygon_2& p) {
	Polygon_2 p_m;
	CGAL::minimum_enclosing_parallelogram_2(
		p.vertices_begin(), p.vertices_end(), std::back_inserter(p_m));
	return p_m;
}

Polygon_2 ReduceVertices(const Polygon_2& pol) {
	if (pol.size() < MAX_VERTICES) return pol;

	return MinParalelogram(pol);

	// TODO(): Fix the implementation below to generate something closer to the input polygon.

	Polygon_2 res;
	for (int i = 0; i < pol.size(); i += 4) {
		if (i + 4 > pol.size()) {
			for (int j = i; j < pol.size(); j++) {
				res.push_back(pol[j]);
			}
		}
		else {
			CGAL::Line_2< CGAL::Cartesian<double> > l1(ToPoint(pol[i]), ToPoint(pol[i + 1]));
			CGAL::Line_2< CGAL::Cartesian<double> > l2(ToPoint(pol[i + 3]), ToPoint(pol[(i + 4) % pol.size()]));
			res.push_back(pol[i]);

			CGAL::Object result;
			CGAL::Point_2<CGAL::Cartesian<double>> ipoint;
			result = CGAL::intersection(l1, l2);
			if (CGAL::assign(ipoint, result)) {
				res.push_back(Point_2(ipoint.x(), ipoint.y()));
			}
			else {
			}
		}
	}

	return res;
}

std::vector< Polygon_2 > CalcReachability(const Polygon_2& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu) {
	int n = t / delta_t;

	double a_norm = a.MaxNorm();
	double alpha_r = (exp(delta_t * a_norm) - 1 - delta_t * a_norm) * Sup(initial_value);
	double beta_r = (exp(delta_t * a_norm) - 1) / a_norm * mu;

	Matrix2D a_exp = (a * delta_t).Exponential();

	Polygon_2 P0_1 = a_exp * initial_value;

	Matrix2D a_exp_2 = (a * -delta_t).Exponential();

	Polygon_2 P0_2 = a_exp_2 * initial_value;

	Polygon_2 P0 = MinkowskiSum(Multiply(P0_1, 0.5), Multiply(P0_2, 0.5));

	Polygon_2 Q0 = MinkowskiSum(P0, CreateHyperrectangle(alpha_r + beta_r));

	std::vector< Polygon_2 > res = { Q0 };

	Polygon_2 result = Q0;

	for (int i = 1; i < n; i++) {
		Q0 = ReduceVertices(Q0);

		P0 = a_exp * Q0;
		Q0 = MinkowskiSum(P0, CreateHyperrectangle(beta_r));

		res.push_back(Q0);
	}

	// res.push_back(ConvexHullOfPolygons(res));

	return res;
}

std::vector< Polygon_2 > StateDependantReachability(const Polygon_2& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu) {
	int n = t / delta_t;

	Matrix2D a_exp = (a * t).Exponential();

	Matrix2D a_inv = (a * t).Exponential();

	Matrix2D f = a_inv * (a_exp - Matrix2D::Identity(a.Size()));

	std::vector< Polygon_2 > res = { };

	return res;
}

//std::vector< VPolytope > Reachability_2(const VPolytope& initial_value, double t,
//		double delta_t, const Matrix& a, double mu) {
//	int n = t / delta_t;

//	double a_norm = a.MaxNorm();
//	double alpha_r = (exp(delta_t * a_norm) - 1 - delta_t * a_norm) * Sup(initial_value);
//	double beta_r = (exp(delta_t * a_norm) - 1) / a_norm * mu;

//	Matrix a_exp = (a * delta_t).Exponential();

//	VPolytope P0_1 = Transform(a_exp, initial_value);

//	Matrix a_exp_2 = (a * -delta_t).Exponential();

//	VPolytope P0_2 = Transform(a_exp_2, initial_value);

//	VPolytope P0 = MinkowskiSum(Scale(P0_1, 0.5), Scale(P0_2, 0.5));

//	VPolytope Q0 = MinkowskiSum(P0, CreateHyperrectanglePolytope(alpha_r + beta_r));

//	std::vector< VPolytope > res = { initial_value, Q0 };

//	VPolytope result = Q0;

//	for (int i = 1; i < 10; i++) {
//		//Q0 = ReduceVertices(Q0);

//		P0 = Transform(a_exp, Q0);
//		Q0 = MinkowskiSum(P0, CreateHyperrectanglePolytope(beta_r));

//		res.push_back(Q0);
//	}

//	// res.push_back(ConvexHullOfPolygons(res));

//	return res;
//}

Zonotope ReduceVerticesZonotope(const Zonotope& zonotope) {
	if (zonotope.Generators().size() < MAX_VERTICES) return zonotope;

	std::vector<std::vector<double> > generators;

	std::vector<std::vector<double> > old_generators = zonotope.Generators();

	std::sort(std::begin(old_generators),
		std::end(old_generators),
		[](const std::vector<double>& l, const std::vector<double>& r) {
		return L1Norm(l) - LInfNorm(l) < 
			L1Norm(r) - LInfNorm(r);
	});

	int d = old_generators[0].size();
	int p = old_generators.size();
	int r;
	if (d == 2) {
		r = 4;
	}
	else {
		r = 2;
	}
	int m = p - d * (r - 1);
	for (int i = 0; i < d; i++) {
		std::vector<double> v(d, 0.0);
		for (int j = 0; j < m; j++) {
			v[i] += abs(old_generators[j][i]);
		}
		generators.push_back(v);
	}

	for (int i = m; i < p; i++) {
		generators.push_back(old_generators[i]);
	}
		
	return Zonotope(zonotope.Center(), generators);
}

std::vector< Zonotope > ReachabilityLinearSystem(const Zonotope& initial_value, double t,
		double delta_t, const Matrix2D& a, double mu) {

	int n = t / delta_t;
	int d = initial_value.Center().size();

	double a_norm = a.MaxNorm();
	double alpha_r = (exp(delta_t * a_norm) - 1 - delta_t * a_norm) * Sup(initial_value);
	double beta_r = (exp(delta_t * a_norm) - 1) / a_norm * mu;

	Matrix2D a_exp = (a * delta_t).Exponential();

	Zonotope P0_1 = a_exp * initial_value;

	Matrix2D a_exp_2 = (a * -delta_t).Exponential();

	Zonotope P0_2 = a_exp_2 * initial_value;

	Zonotope P0 = MinkowskiSum(P0_1 * 0.5, P0_2 * 0.5);

	Zonotope Q0 = MinkowskiSum(P0, CreateHyperrectangleZonotope(alpha_r + beta_r, d));

	std::vector< Zonotope > res = { Q0 };

	Zonotope result = Q0;

	for (int i = 1; i <= n; i++) {
		P0 = a_exp * Q0;
		Q0 = MinkowskiSum(P0, CreateHyperrectangleZonotope(beta_r, d));
			
		Q0 = ReduceVerticesZonotope(Q0);

		res.push_back(Q0);
	}

	return res;
}

Zonotope ReachabilityLinearSystemT(const Zonotope& initial_value, double t,
	double delta_t, const Matrix2D& a, double mu) {

	int n = t / delta_t;
	int d = initial_value.Center().size();

	double a_norm = a.MaxNorm();
	double alpha_r = (exp(delta_t * a_norm) - 1 - delta_t * a_norm) * Sup(initial_value);
	double beta_r = (exp(delta_t * a_norm) - 1) / a_norm * mu;

	Matrix2D a_exp = (a * delta_t).Exponential();

	Zonotope P0_1 = a_exp * initial_value;

	Matrix2D a_exp_2 = (a * -delta_t).Exponential();

	Zonotope P0_2 = a_exp_2 * initial_value;

	Zonotope P0 = MinkowskiSum(P0_1 * 0.5, P0_2 * 0.5);

	Zonotope Q0 = MinkowskiSum(P0, CreateHyperrectangleZonotope(alpha_r + beta_r, d));

	Zonotope result = Q0;

	for (int i = 1; i <= n; i++) {
		result = a_exp * result;
	}

	return result;
}

std::vector< Zonotope > ReachabilityFromInput(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, double u_0, double mu) {

	if (abs(A.calculateDeterminant(A.Size())) < 1.0e-6) {
		std::cout << "Error: Non-invertible matrix" << std::endl;
		return {};
	}

	Matrix2D A_inv = A.Inverse();

	Matrix2D I = A.Identity(A.Size());

	Matrix2D A_exp = (A * T).Exponential();

	Matrix2D M = A_inv * (A_exp - I);

	Zonotope u({ u_0 }, { {mu} });

	Zonotope R_u = (M * B) * u;

	int n = T / delta_t;

	std::vector< Zonotope > res = { Zonotope() };

	for (int i = 1; i <= n; i++) {
		Zonotope R = R_u * (1.0 * i / n);
			
		res.push_back(R);
	}

	return res;
}

std::vector< Zonotope > ReachabilityFromState(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& f_star, double mu) {

	if (abs(A.calculateDeterminant(A.Size())) < 1.0e-6) {
		std::cout << "Error: Non-invertible matrix" << std::endl;
		return {};
	}

	Matrix2D A_inv = A.Inverse();

	Matrix2D I = A.Identity(A.Size());

	Matrix2D A_exp = (A * T).Exponential();

	Matrix2D M = A_inv * (A_exp - I);

	Matrix2D R = M * f_star;

	int iter = T / delta_t;

	Matrix2D tmp2 = (A_inv *
		(A_exp - I - A * T - (A * T) * (A * T) * 0.375));

	Matrix2D tmp = (A_inv *
		(A_exp - I - A * T - (A * T) * (A * T) * 0.375) * f_star);

	double eta = tmp.MaxNorm();

	eta = mu;

	int n = T / delta_t;

	std::vector< Zonotope > res = { Zonotope() };

	for (int i = 1; i <= n; i++) {
		//Zonotope R_1 = initial_value + f_star * (i * delta_t);
		Zonotope R_t = MinkowskiSum(Zonotope(R * (1.0 * i / n)), CreateHyperrectangleZonotope(eta, 2));

		//Zonotope R = Merge(Merge(R_1, R_2), R_t);

		res.push_back(R_t);
	}

	return res;
}

Zonotope ReachabilityFromStateT(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& f_star, double mu) {

	if (abs(A.calculateDeterminant(A.Size())) < 1.0e-6) {
		std::cout << "Error: Non-invertible matrix" << std::endl;
		return {};
	}

	Matrix2D A_inv = A.Inverse();

	Matrix2D I = A.Identity(A.Size());

	Matrix2D A_exp = (A * T).Exponential();

	Matrix2D M = A_inv * (A_exp - I);

	Matrix2D R = M * f_star;

	int iter = T / delta_t;

	Matrix2D tmp = (A_inv *
		(A_exp - I - A * T - (A * T) * (A * T) * 0.375) * f_star);

	double eta = tmp.MaxNorm();

	return Zonotope(R);
}

std::vector< Zonotope > CalcReachability(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double mu) {
	const auto& R_1_arr = 
		ReachabilityLinearSystem(initial_value, T, delta_t, A, mu);
	const auto& R_2_arr =
		ReachabilityFromState(initial_value, T, delta_t, A,
			f_star, mu);
	const auto& R_3_arr =
		ReachabilityFromInput(initial_value, T, delta_t, A, B, 0.038, mu);

	std::vector< Zonotope > res = { };
	int d = R_1_arr[0].Dimension();
	for (int i = 0; i < R_1_arr.size(); i++) {
		Zonotope R_2 = R_2_arr[i].ExpandDim(d, true);
		Zonotope R_3 = R_3_arr[i].ExpandDim(d, true);
		Zonotope R = MinkowskiSum(MinkowskiSum(R_1_arr[i], R_2), R_3);

		R = ReduceVerticesZonotope(R);

		res.push_back(R);
	}

	return res;
}

Polygon_2 ZonotopeToPolygon2(const Zonotope& zonotope) {
	polygon hull = zonotope.GetPolygon();
	Polygon_2 p;
	for (int i = 0; i < hull.outer().size() - 1; i++) {
		p.push_back(Point(hull.outer()[i].get<0>(), hull.outer()[i].get<1>()));
	}
	p.reverse_orientation();
	return p;
}

polygon ZonotopeToPolygon(const Zonotope& zonotope) {
	polygon hull = zonotope.GetPolygon();
	return hull;
}

Polygon_2 MergeZonotopes(const std::vector< Zonotope >& zonotopes, int i, int j) {
	if (j <= i) {
		return Polygon_2();
	}
	if (j - i == 1) {
		return ZonotopeToPolygon2(zonotopes[i]);
	}

	int h = (j - i) / 2;
	Polygon_2 p1 = MergeZonotopes(zonotopes, i, i + h);
	Polygon_2 p2 = MergeZonotopes(zonotopes, i + h, j);

	if (p1.size() == 0) {
		return p2;
	}
	else if (p2.size() == 0) {
		return p1;
	}
	else {
		std::cout << "Joining " << i << " " << j << std::endl;
		Polygon_with_holes_2 unionR;
		CGAL::join(p1, p2, unionR);
		return unionR.outer_boundary();
	}
}

bool check_inside(const Polygon_2& pol, const K::Point_2& pt) {
	static std::vector<K::Point_2> points;
	if (points.size() < pol.size()) {
		points.resize(pol.size());
	}
	for (int i = 0; i < pol.size(); i++) {
		points[i] = K::Point_2(CGAL::to_double(pol[i].x()), CGAL::to_double(pol[i].y()));
	}
	switch (CGAL::bounded_side_2(points.begin(), points.begin() + pol.size(), pt, K())) {
	case CGAL::ON_BOUNDED_SIDE:
		return true;;
	case CGAL::ON_BOUNDARY:
		return true;
	case CGAL::ON_UNBOUNDED_SIDE:
		return false;
	}
	//for (int i = 0; i < pol.size(); i++) {
	//	points[i].~K::Point_2();
	//}
	//delete[] points;
}

Mesh ZonotopesToMatrix(const std::vector< Zonotope >& zonotopes, const Mesh& mesh, bool join_polygons) {
	//using namespace boost::numeric::ublas;
	//typedef boost::geometry::model::point<double, 5, boost::geometry::cs::cartesian> point_t;
	//typedef boost::geometry::model::polygon<point_t> polygon;

	//tensor<bool> res{ (size_t) scenario.N(), (size_t) scenario.N() };
	if (join_polygons) {
		Polygon_2 pol1 = MergeZonotopes(zonotopes, 0, zonotopes.size() / 2);
		Polygon_2 pol2 = MergeZonotopes(zonotopes, zonotopes.size() / 2, zonotopes.size());

		int n = mesh.NElems();

		Mesh res(mesh);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				std::pair<double, double> p = mesh.GetPoint(i, j);
				K::Point_2 q = K::Point_2(p.first, p.second);

				if (check_inside(pol1, q) ||
					check_inside(pol2, q)) {
					res[i][j] = 1.0;
				}
				else {
					res[i][j] = 0.0;
				}
			}
		}

		return res;
	}
	else {
		int n = mesh.NElems();

		Mesh res(mesh);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res[i][j] = 0.0;
			}
		}

		for (const auto& zonotope : zonotopes) {

			//Polygon_2 pol = ZonotopeToPolygon2(zonotope);

			//double x_min, x_max, y_min, y_max;
			//x_min = CGAL::to_double(pol[0].x());

			//for (int i = 0; i < n; i++) {
			//	for (int j = 0; j < n; j++) {
			//		std::pair<double, double> p = mesh.GetPoint(i, j);
			//		K::Point_2 q = K::Point_2(p.first, p.second);

			//		if (check_inside(pol, q)) {
			//			res[i][j] = 1.0;
			//		}
			//	}
			//}

			polygon hull = ZonotopeToPolygon(zonotope);

			Mesh tmp(res);

			PolygonToMesh(hull, tmp);

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (abs(tmp[i][j] - 1.0) < 1.0e-4) {
						res[i][j] = 1.0;
					}
				}
			}

			//auto pol = hull.outer();
			//double x_min, x_max, y_min, y_max;
			//if (pol.size() > 0) {
			//	x_min = x_max = pol[0].get<0>();
			//	y_min = y_max = pol[0].get<1>();
			//	for (int i = 0; i < pol.size(); i++) {
			//		x_min = std::min(x_min, pol[i].get<0>());
			//		x_max = std::max(x_max, pol[i].get<0>());

			//		y_min = std::min(y_min, pol[i].get<1>());
			//		y_max = std::max(y_max, pol[i].get<1>());
			//	}

			//	for (int i = 0; i < n; i++) {
			//		for (int j = 0; j < n; j++) {
			//			std::pair<double, double> p = mesh.GetPoint(i, j);

			//			if (p.first < x_min || p.first > x_max) continue;
			//			if (p.second < y_min || p.second > y_max) continue;

			//			int ii, jj, c = 0;
			//			int nvert = pol.size();
			//			for (ii = 0, jj = nvert - 1; ii < nvert; jj = ii++) {
			//				if (((pol[ii].get<1>() > p.first) != (pol[jj].get<1>() > p.first)) &&
			//					(p.first < (pol[jj].get<0>() - pol[ii].get<0>()) * (p.second - pol[ii].get<1>()) / (pol[jj].get<1>() - pol[ii].get<0>()) + pol[ii].get<0>()))
			//					c = !c;
			//			}

			//			if (c) {
			//				res[i][j] = 1.0;
			//			}
			//		}
			//	}
			//}
		}

		return res;
	}

	//std::ofstream myfile;
	//myfile.open("output.txt");
	//myfile << pol1.size() << " ";
	//for (const auto& it : pol1) {
	//	myfile << it << " ";
	//}
	//myfile << pol2.size() << " ";
	//for (const auto& it : pol2) {
	//	myfile << it << " ";
	//}
	//myfile << "\n";
	//myfile.close();
}

Reachability::Reachability(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double u0, double mu) {

	R_up_ =
		ReachabilityLinearSystem(initial_value, T, delta_t, A, mu);
	R_up_t_ =
		ReachabilityLinearSystemT(initial_value, T, delta_t, A, mu);
	R_down_ =
		ReachabilityFromState(initial_value, T, delta_t, A,
			f_star, mu);
	R_down_t_ = ReachabilityFromStateT(initial_value, T, delta_t, A,
		f_star, mu);
	R_bar_ =
		ReachabilityFromInput(initial_value, T, delta_t, A, B, u0, mu);
	R_bar_t_ = R_bar_.back();

	std::vector< Zonotope > res = { };
	for (int i = 0; i < R_up_.size(); i++) {
		Zonotope R_2 = R_up_[i];
		Zonotope R_3 = R_down_[i];
		Zonotope R = MinkowskiSum(MinkowskiSum(R_bar_[i], R_2), R_3);

		R = ReduceVerticesZonotope(R);

		res.push_back(R);
	}

	R_ = res;

	Zonotope R = MinkowskiSum(MinkowskiSum(R_bar_t_, R_up_t_), R_down_t_);

	R_t_ = ReduceVerticesZonotope(R);
}

void PolygonToMesh(const polygon& pol, Mesh& mesh) {
	int  nodes, pixelY, i, j, swap;
	double pixelX;
	std::vector<double> nodeX;

	int height = mesh.NRows();
	auto& hull = pol.outer();
	int polyCorners = hull.size();
	nodeX.resize(polyCorners);

	mesh.Init(0.0);

	//  Loop through the rows of the image.
	for (pixelY = 0; pixelY < height; pixelY++) {

		//  Build a list of nodes.
		nodes = 0; j = polyCorners - 1;
		double y = mesh.GetPoint(0, pixelY).second;
		for (i = 0; i < polyCorners; i++) {
			if (hull[i].get<1>() < y && hull[j].get<1>() >= y
				|| hull[j].get<1>() < y && hull[i].get<1>() >= y) {
				nodeX[nodes++] = hull[i].get<0>() + (y - hull[i].get<1>()) / (hull[j].get<1>() - hull[i].get<1>())
					* (hull[j].get<0>() - hull[i].get<0>());
			}
			j = i;
		}

		if (nodes == 0) continue;

		//  Sort the nodes, via a simple “Bubble” sort.
		sort(nodeX.begin(), nodeX.begin() + nodes);

		//  Fill the pixels between node pairs.
		for (i = 0; i < nodes; i++) {
			if (nodeX[i] >= mesh.XMax()) break;
			if (nodeX[i + 1] > mesh.XMin()) {
				if (nodeX[i] < mesh.XMin()) nodeX[i] = mesh.XMin();
				if (nodeX[i + 1] > mesh.XMax()) nodeX[i + 1] = mesh.XMax();
				for (pixelX = nodeX[i]; pixelX < nodeX[i + 1]; pixelX += mesh.DeltaX()) mesh.set( pixelX, y);
			}
		}
	}
}

void Reachability::UpdateReachabilityInput(const Zonotope& initial_value, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& f_star, double u0, double mu) {
	R_bar_ =
		ReachabilityFromInput(initial_value, T, delta_t, A, B, u0, mu);
	R_bar_t_ = R_bar_.back();

	std::vector< Zonotope > res = { };
	for (int i = 0; i < R_up_.size(); i++) {
		Zonotope R_2 = R_up_[i];
		Zonotope R_3 = R_down_[i];
		Zonotope R = MinkowskiSum(MinkowskiSum(R_bar_[i], R_2), R_3);

		R = ReduceVerticesZonotope(R);

		res.push_back(R);
	}

	R_ = res;

	Zonotope R = MinkowskiSum(MinkowskiSum(R_bar_t_, R_up_t_), R_down_t_);

	R_t_ = ReduceVerticesZonotope(R);
}