#pragma once

#include <vector>
#include <optional>

struct Point {
	Point(double xx, double yy) : x(xx), y (yy) {}
	Point() {}

	double x;
	double y;

	Point operator-(const Point& p) const {
		return { x - p.x, y - p.y };
	}

	Point operator+(const Point& p) const  {
		return { x + p.x, y + p.y };
	}

	Point operator*(double f) const  {
		return { x * f, y * f };
	}
};

class Circle {
public:
	Circle(const Point& p, double _r) : c(p), r(_r) {}

	double Radius() const {
		return r;
	}

	Point Center() const {
		return c;
	}

private:
	Point c;
	double r;
};

class Line {
public:
	Line(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {};

	Line() {}

	double a;
	double b;
	double c;
};

std::vector<std::optional<Line>> FindTangentLines(const Circle& c1, const Circle& c2);

double Distance(const Circle& c, const Line& line);

double CalcDistance(const Point& a, const Point& b);

std::optional<Point> GetTangentPoint(const Line& tangent, const Circle& circle);

std::optional<Point> LineIntersection(const Line& a, const Line& b);

double CrossProduct(const Point& a, const Point& b);

double CalcAngle(const Point& from, const Point& to, const Point& orientation, const Point& circle_center);
