#include "MathUtils.h"

namespace {
	const double EPS = 1.0e-6;
	const double PI = acos(-1.0);

	std::optional<Line> CalcTangentLine(const Point& p, double s1, double s2) {
		double z = p.x * p.x + p.y * p.y;
		double delta = z - (s2 - s1)* (s2 - s1);
		
		if (delta < -EPS) {
			return std::nullopt;
		}

		double a = (p.x * (s2 - s1) + p.y * sqrt(delta)) / (z);
		double b = (p.y * (s2 - s1) - p.x * sqrt(delta)) / (z);
		double c = s1;

		return Line(a, b, c);
	}

}

std::vector<std::optional<Line>> FindTangentLines(const Circle& c1, const Circle& c2) {
	double r1 = c1.Radius();
	double r2 = c2.Radius();

	std::vector<std::optional<Line>> result;

	for (int i = -1; i <= 1; i += 2) {
		double s1 = i * r1;
		for (int j = -1; j <= 1; j += 2) {
			double s2 = j * r2;
			std::optional<Line> res = CalcTangentLine(c2.Center() - c1.Center(), s1, s2);

			if (res.has_value()) {
				Line line = res.value();
				line.c -= line.a* c1.Center().x + line.b * c1.Center().y;
				result.push_back(line);
			}
			else {
				result.push_back(std::nullopt);
			}
		}
	}

	return result;
}

double Distance(const Circle& c, const Line& line) {
	return abs(line.a * c.Center().x + line.b * c.Center().y + line.c) / sqrt(line.a * line.a + line.b * line.b);
}

double CalcDistance(const Point& a, const Point& b) {
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

std::optional<Point> GetTangentPoint(const Line& tangent, const Circle& circle) {
	Line normal;
	normal.a = -tangent.b;
	normal.b = tangent.a;
	normal.c = -(normal.a * circle.Center().x + normal.b * circle.Center().y);

	return LineIntersection(tangent, normal);
}

std::optional<Point> LineIntersection(const Line& a, const Line& b) {
	if (abs(a.a * b.b - a.b * b.a) < EPS) {
		return std::nullopt;
	}

	if (abs(a.a) < EPS) {
		double y = -a.c / a.b;
		return Point(-(b.c + b.b * y) / b.a, y );
	} else if (abs(a.b) < EPS) {
		double x = -a.c / a.a;
		return Point( x, -(b.c + b.a * x) / b.b );
	}
	else {
		double x = -(b.c + b.b*(-a.c / a.b)) / (b.a + b.b * (-a.a / a.b));
		return Point(x, -(a.c + a.a * x) / a.b);
	}
}

double CrossProduct(const Point& a, const Point& b) {
	return a.x * b.y - a.y * b.x;
}

double CalcAngle(const Point& from, const Point& to, const Point& orientation, const Point& circle_center) {
	double cross_prod = CrossProduct(to - circle_center, orientation);

	double angle_diff = (atan2(to.y - circle_center.y, to.x - circle_center.x) - 
		atan2(from.y - circle_center.y, from.x - circle_center.x));

	if (angle_diff < 0.0) {
		angle_diff = 2 * PI + angle_diff;
	}

	if (cross_prod < 0.0) {
		angle_diff = 2*PI - angle_diff;
	}

	return angle_diff;
}