#include "Zonotope.h"

double Sup(const Zonotope& z) {
	std::vector<double> v = z.GetVertices();
	double res = 0.0;
	for (int i = 0; i < v.size() / 2; i += 2) {
		res = std::max(res, sqrt(v[2 * i] * v[2 * i] + v[2 * i + 1] * v[2 * i + 1]));
	}
	return res;
}
