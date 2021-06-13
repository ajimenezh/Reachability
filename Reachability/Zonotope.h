#pragma once

#include <vector>
#include <iterator>

#include "Matrix2D.h"
//#include "zonotope_c.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

typedef boost::geometry::model::point<double, 5, boost::geometry::cs::cartesian> point_t;
typedef boost::geometry::model::polygon<point_t> polygon;

#define MAX(a, b) (a) > (b) ? (a) : (b)

class Zonotope {
public:
	Zonotope() {}
	//Zonotope(const std::vector<std::vector<double>>& generators)
	//	: generators_(generators) {};
	Zonotope(std::vector<double> center, std::vector<std::vector<double>> generators)
		: center_(center), generators_(generators) {};

	Zonotope(std::pair<double, double> c)
		: center_({c.first, c.second}) {};

	Zonotope(const Matrix2D& m) {
		std::vector<double> v;
		if (m.NRows() == 1) {
			for (int i = 0; i < m.NColumns(); i++) {
				v.push_back(m[0][i]);
			}
		}
		else if (m.NColumns() == 1) {
			for (int i = 0; i < m.NRows(); i++) {
				v.push_back(m[i][0]);
			}
		}
		else {

		}

		center_ = v;
	}

	const std::vector<double>& operator[](int i) const {
		return generators_[i];
	}

	std::vector<double>& operator[](int i) {
		return generators_[i];
	}

	const std::vector < std::vector<double> >& Generators() const {
		return generators_;
	}

	std::vector<double>& Center() {
		return center_;
	}

	const std::vector<double>& Center() const {
		return center_;
	}

	void append(const Zonotope& other) {
		for (int i = 0; i < other.Generators().size(); i++) {
			generators_.push_back(other[i]);
		}
	}

	//std::vector<double> GetVertices() const {
	//	//double* vertices;

	//	double* generators = new double[generators_.size()* generators_[0].size()];
	//	for (int i = 0; i < generators_.size(); i++) {
	//		for (int j = 0; j < generators_[0].size(); j++) {
	//			generators[j * (int) generators_.size() + i] = 
	//				generators_[i][j];
	//		}
	//	}

	//	double* center = new double[center_.size()];
	//	for (int i = 0; i < center_.size(); i++) {
	//		center[i] = center_[i];
	//	}

	///*	long n_points = zonotope_vertices_double( 
	//		generators_[0].size(), generators_.size(), generators, center, &vertices);
	//	*/

	//	int n_points = 1;
	//	std::vector<double> res(n_points*2);
	//	//for (int i = 0; i < n_points; i++) {
	//	//	res[2 * i] = vertices[2 * i] + center_[0];
	//	//	res[2 * i + 1] = vertices[2 * i + 1] + center_[1];
	//	//}

	//	//for (int i = 0; i < n_points; i++) {
	//	//	res[2 * i] = vertices[2 * i];
	//	//	res[2 * i + 1] = vertices[2 * i + 1];
	//	//}

	//	return res;
	//}

	std::vector<double> GetVertices() const {
		polygon hull = GetPolygon();

		std::vector<double> res;

		for (const auto& it : hull.outer()) {
			res.push_back(it.get<0>());
			res.push_back(it.get<1>());
		}

		return res;
	}

	polygon GetPolygon() const {
		int n = generators_.size();
		int d = generators_[0].size();

		std::vector< point_t> v;
		boost::geometry::model::multi_point<point_t> pts;
		pts.reserve((1 << n));
		for (int i = 0; i < (1 << n); i++) {
			point_t p(center_[0], center_[1]);
			for (int j = 0; j < n; j++) {
				double f = -1;
				if (i & (1 << j)) {
					f = 1;
				}

				p.set<0>(p.get<0>() + f * generators_[j][0]);
				p.set<1>(p.get<1>() + f * generators_[j][1]);
			}
			pts.emplace_back(p);
		}

		//polygon poly(v.data());

		polygon hull;
		boost::geometry::convex_hull(pts, hull);

		return hull;
	}

	int Dimension() const {
		return center_.size();
	}

	Zonotope ExpandDim(int n, bool end_) const {
		if (Dimension() >= n) return (*this);

		std::vector<std::vector<double> > g;
		for (int i = 0; i < generators_.size(); i++) {
			std::vector<double> v(n, 0.0);
			if (end_) {
				int start = n - Dimension();
				for (int j = 0; j < generators_[i].size(); j++) {
					v[start + j] = generators_[i][j];
				}
			}
			else {
				for (int j = 0; j < generators_[i].size(); j++) {
					v[j] = generators_[i][j];
				}
			}
			g.push_back(v);
		}
		std::vector<double> c(n, 0.0);
		int start = n - Dimension();
		for (int j = 0; j < center_.size(); j++) {
			c[start + j] = center_[j];
		}
		return Zonotope(c, g);
	}

	Zonotope subZonotope(int n, int i) const {
		std::vector < std::vector<double> > g;
		for (int k = 0; k < generators_.size(); k++) {
			std::vector<double> v;
			std::copy(generators_[k].begin() + i,
				generators_[k].begin() + i + n, std::back_inserter(v));
			g.push_back(v);
		}
		std::vector<double> c;
		std::copy(center_.begin() + i, center_.begin() + i + n, 
			std::back_inserter(c));
		return Zonotope(c, g);
	}

private:
	std::vector < std::vector<double> > generators_;
	std::vector<double> center_;
};

inline Zonotope operator*(const Matrix2D& a, const Zonotope& zonotope) {
	std::vector<std::vector<double> > generators;
	std::vector<double> center;
	for (int i = 0; i < zonotope.Generators().size(); i++) {
		std::vector<double> v;
		for (int j = 0; j < a.NRows(); j++) {
			double val = 0.0;
			for (int k = 0; k < zonotope.Generators()[0].size(); k++) {
				val += a[j][k] * zonotope[i][k];
			}
			v.push_back(val);
		}
		generators.push_back(v);
	}
	for (int j = 0; j < a.NRows(); j++) {
		double val = 0.0;
		for (int k = 0; k < zonotope.Center().size(); k++) {
			val += a[j][k] * zonotope.Center()[k];
		}
		center.push_back(val);
	}
	return Zonotope(center, generators);
}

inline Zonotope operator+(const Zonotope& zonotope, const Matrix2D& a) {
	std::vector<std::vector<double> > generators;
	std::vector<double> center;
	for (int i = 0; i < zonotope.Generators().size(); i++) {
		std::vector<double> v;
		for (int j = 0; j < a.NRows(); j++) {
			double val = a[j][0] + zonotope.Generators()[i][j];
			v.push_back(val);
		}
		generators.push_back(v);
	}
	for (int j = 0; j < a.NRows(); j++) {
		double val = a[j][0] * zonotope.Center()[j];
		center.push_back(val);
	}
	return Zonotope(center, generators);
}

inline Zonotope operator*(const Zonotope& zonotope, double v) {
	Zonotope other = zonotope;
	for (int i = 0; i < zonotope.Generators().size(); i++) {
		for (int j = 0; j < zonotope.Generators()[0].size(); j++) {
			other[i][j] = zonotope[i][j] * v;
		}
	}
	for (int i = 0; i < other.Center().size(); i++) {
		other.Center()[i] *= v;
	}
	return other;
}

inline Zonotope MinkowskiSum(const Zonotope& a, const Zonotope& b) {
	if (a.Center().size() == 0) {
		return b;
	} 
	if (b.Center().size() == 0) {
		return a;
	}
	Zonotope res = a;
	res.append(b);
	for (int i = 0; i < b.Center().size(); i++) {
		res.Center()[i] += b.Center()[i];
	}
	return res;
}

inline Zonotope CreateHyperrectangleZonotope(double a, int d) {
	std::vector<double> c(d, 0.0);
	std::vector<double> g1(d, 0.0); g1[0] = a;
	std::vector<double> g2(d, 0.0); g2[1] = a;

	return Zonotope(c, { g1, g2 });
}

inline Zonotope ReduceVertices(const Zonotope& a) {
	return a;
}

inline double L1Norm(const std::vector<double>& v) {
	double res = 0;
	for (double val : v) {
		res += abs(val);
	}
	return res;
}

inline double LInfNorm(const std::vector<double>& v) {
	double res = v[0];
	for (double val : v) {
		res = MAX(res, abs(val));
	}
	return res;
}

inline Zonotope Merge(const Zonotope& a, Zonotope& b) {
	std::vector < std::vector<double> > g;
	std::vector<double> c;

	Zonotope a_ = a.ExpandDim(a.Dimension() + b.Dimension(), false);
	Zonotope b_ = b.ExpandDim(a.Dimension() + b.Dimension(), true);

	for (int i = 0; i < a_.Generators().size(); i++) {
		std::vector<double> v;
		std::copy(a_.Generators()[i].begin(),
			a_.Generators()[i].end(), std::back_inserter(v));
		g.push_back(v);
	}

	for (int i = 0; i < b_.Generators().size(); i++) {
		std::vector<double> v;
		std::copy(b_.Generators()[i].begin(),
			b_.Generators()[i].end(), std::back_inserter(v));
		g.push_back(v);
	}

	std::copy(a.Center().begin(),
		a.Center().end(), std::back_inserter(c));
	std::copy(b.Center().begin(),
		b.Center().end(), std::back_inserter(c));

	return Zonotope(c, g);
}

double Sup(const Zonotope& z);