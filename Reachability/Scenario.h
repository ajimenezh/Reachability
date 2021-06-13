#pragma once

#include <vector>


struct Interval {
	double l_;
	double r_;

	Interval(double l, double r) : l_(l), r_(r) {}
	Interval() {}
};

class Scenario {
public:
	Scenario(Interval x, Interval y, Interval x3, Interval x4, Interval x5) : 
		x_(x), y_(y), x3_(x3), x4_(x4), x5_(x5) {
	}

	Scenario() {
		x_ = Interval(0, 80);
		y_ = Interval(-4, 4);
		x3_ = Interval(-0.5, 0.5);
		x4_ = Interval(-0.5, 0.5);
		x5_ = Interval(-0.5, 0.5);
	}

	int N() const {
		return n_;
	}

	int M() const {
		return m_;
	}

	std::vector<double> GetPoint(int i, int j) const {
		double delta_x = (x_.r_ - x_.l_) / n_;
		double delta_y = (y_.r_ - y_.l_) / n_;
		return { x_.l_ + i * delta_x, y_.l_ + j * delta_y };
	}

private:
	Interval x_;
	Interval y_;
	Interval x3_;
	Interval x4_;
	Interval x5_;

	int n_ = 50;
	int m_ = 40;
};

