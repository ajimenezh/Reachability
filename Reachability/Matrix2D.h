#pragma once

#include <vector>

#include "bops_linear.h"
#include "matrix_exponential.hpp"
#include <CGAL/Polygon_2.h>

class Matrix2D {
public:
	Matrix2D(int n, int m) {
		data_.resize(n);
		for (int i = 0; i < n; i++) data_[i].resize(m, 0.0);
		n_ = n;
		m_ = m;
	}

	int NRows() const {
		return n_;
	}

	int NColumns() const {
		return m_;
	}

	std::vector<double>& operator[](int i) { return data_[i]; }

	const std::vector<double>& operator[](int i) const { return data_[i]; }

	Matrix2D operator*(double x) const {
		Matrix2D a = *this;
		for (int i = 0; i < n_; i++) {
			for (int j = 0; j < m_; j++) {
				a[i][j] *= x;
			}
		}
		return a;
	}

	Polygon_2 operator*(const Polygon_2& a) const {
		Polygon_2 b;
		if (n_ != 2) return b;
		for (const auto& vertex : a) {
			const double& x = to_double(vertex.x().exact());
			const double& y = to_double(vertex.y().exact());
			b.push_back(Point_2(data_[0][0] * x + data_[0][1] * y,
				data_[1][0] * x + data_[1][1] * y));
		}
		return b;
	}

	Matrix2D operator-(const Matrix2D& other) const {
		Matrix2D a = *this;
		for (int i = 0; i < n_; i++) {
			for (int j = 0; j < m_; j++) {
				a[i][j] -= other[i][j];
			}
		}
		return a;
	}

	Matrix2D operator*(const Matrix2D& other) const {
		Matrix2D a(n_, other.NColumns());
		for (int i = 0; i < n_; i++) {
			for (int j = 0; j < other.NColumns(); j++) {
				for (int k = 0; k < n_; k++) {
					a[i][j] += data_[i][k] * other[k][j];
				}
			}
		}
		return a;
	}

	double MaxNorm() const {
		if (data_.empty() || data_[0].empty()) return 0.0;
		int n = data_.size();
		int m = data_[0].size();
		double res = data_[0][0];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				res = std::max(res, abs(data_[i][j]));
			}
		}

		//double res = calculateDeterminant(n);
		return res;
	}

	Matrix2D Exponential() const  {
		int n = data_.size();

		double* a = new double[n * n];

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				a[i * n + j] = data_[i][j];
			}
		}

		double* res = r8mat_expm1(n, a);

		Matrix2D res_matrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res_matrix[i][j] = res[i * n + j];
			}
		}

		delete[] res;
		delete[] a;

		return res_matrix;
	}

	Matrix2D Inverse() const {
		// Calculate the inverse of the determinant of m.
		double det = calculateDeterminant(n_);
		double inverseDet = 1.0f / det;

		Matrix2D result(n_, n_);

		for (int j = 0; j < n_; j++) {
			for (int i = 0; i < n_; i++) {
				// Get minor of element (j, i) - not (i, j) because
				// this is where the transpose happens.
				double minor = calculateMinor(n_, j, i);

				// Multiply by (−1)^{i+j}
				double factor = ((i + j) % 2 == 1) ? -1.0 : 1.0;
				double cofactor = minor * factor;

				result[i][j] = inverseDet * cofactor;
			}
		}

		return result;
	}

	double calculateDeterminant(int n) const {
		double det = 0.0;

		if (n == 1) return data_[0][0];

		if (n == 2) {
			return data_[0][0] * data_[1][1] -
				data_[0][1] * data_[1][0];
		}

		for (int i = 0; i < n; i++) {
			// Get minor of element (0, i)
			double minor = calculateMinor(n, 0, i);

			// If this is an odd-numbered row, negate the value.
			double factor = (i % 2 == 1) ? -1.0 : 1.0;

			det += factor * data_[0][i] * minor;
		}

		return det;
	}

	double calculateMinor(int n, int row, int col) const {
		auto minorSubmatrix = getMinor(n, row, col);
		return minorSubmatrix.calculateDeterminant(n - 1);
	}

	Matrix2D getMinor(int n, int row, int col) const {
		int colCount = 0, rowCount = 0;

		Matrix2D dest(n - 1, n - 1);
		for (int i = 0; i < n; i++) {
			if (i != row) {
				colCount = 0;
				for (int j = 0; j < n; j++) {
					if (j != col) {
						dest[rowCount][colCount] = data_[i][j];
						colCount++;
					}
				}
				rowCount++;
			}
		}

		return dest;
	}

	Matrix2D getSubmatrix(int n, int row, int col) const {
		Matrix2D dest(std::min(n, n_), std::min(n, m_));
		for (int i = row; i < std::min(row + n, n_); i++) {
			for (int j = col; j < std::min(col + n, m_); j++) {
				dest[i-row][j-col] = data_[i][j];
			}
		}

		return dest;
	}

	static Matrix2D Identity(int n) {
		Matrix2D res(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res[i][j] = 0.0;
			}
			res[i][i] = 1.0;
		}
		return res;
	}

	int Size() const { return n_; }

private:
	std::vector<std::vector<double> > data_;
	int n_;
	int m_;
};

class Mesh : public Matrix2D {
public:
	Mesh(int n, double x_min, double x_max, double y_min, double y_max) : 
		Matrix2D(n, n), x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max) {
		n_elems_ = n;
		delta_x_ = (x_max_ - x_min_) / (n);
		delta_y_ = (y_max_ - y_min_) / (n);
	}

	Mesh(int n, double x_min, double x_max) :
		Matrix2D(n, 1), x_min_(x_min), x_max_(x_max) {
		n_elems_ = n;
		delta_x_ = (x_max_ - x_min_) / (n);
		one_dimensional_ = true;
	}

	Mesh() : Matrix2D(0, 0) {}

	void Init(double val) {
		for (int i = 0; i < n_elems_; i++) {
			if (!one_dimensional_) {
				for (int j = 0; j < n_elems_; j++) {
					(*this)[i][j] = val;
				}
			}
			else {
				(*this)[i][0] = val;
			}
		}
	}

	int NElems() const {
		return n_elems_;
	}

	std::pair<double, double> GetPoint(int i, int j) const {
		return std::make_pair(x_min_ + i * delta_x_, y_min_ + j * delta_y_);
	}

	std::pair<double, double> GetPoint(int idx) const {
		if (!one_dimensional_) {
			std::pair<int, int> p = IdToIndices(idx);
			return std::make_pair(x_min_ + p.first * delta_x_, y_min_ + p.second * delta_y_);
		}
		else {
			return std::make_pair(x_min_ + idx * delta_x_, 0.0);
		}
	}

	std::pair<double, double> GetMidPoint(int i, int j) const {
		return std::make_pair(x_min_ + i * delta_x_ + delta_x_ / 2.0, y_min_ + j * delta_y_ + delta_y_ / 2.0);
	}

	double DeltaX() const {
		return delta_x_;
	}

	double DeltaY() const {
		return delta_y_;
	}

	std::pair<int, int> IdToIndices(int id) const {
		return std::pair(id / n_elems_, id % n_elems_);
	}

	double XMin() const {
		return x_min_;
	}

	double XMax() const {
		return x_max_;
	}

	void set(double x, double y) {
		int ii = (x - x_min_) / (x_max_ - x_min_) * n_elems_;
		int jj = (y - y_min_) / (y_max_ - y_min_) * n_elems_;
		if ((ii < 0 || ii >= n_elems_) || (jj < 0 || jj >= n_elems_)) {
			std::cout << "Out of bounds: " << x << " " << y << std::endl;
			return;
		}
		(*this)[ii][jj] = 1.0;
	}

	void set(double x) {
		if (!one_dimensional_) return;
		int ii = (x - x_min_) / (x_max_ - x_min_) * n_elems_;
		if (ii < 0 || ii >= n_elems_) {
			std::cout << "Out of bounds: " << x << std::endl;
			return;
		}
		(*this)[ii][0] = 1.0;
	}

	void add(double x, double y, double val) {
		int ii = (x - x_min_) / (x_max_ - x_min_) * n_elems_;
		int jj = (y - y_min_) / (y_max_ - y_min_) * n_elems_;
		if ((ii < 0 || ii >= n_elems_) || (jj < 0 || jj >= n_elems_)) {
			std::cout << "Out of bounds in add: " << x << " " << y << std::endl;
			return;
		}
		(*this)[ii][jj] += val;
	}

	void add(double x, double val) {
		if (!one_dimensional_) return;
		int ii = (x - x_min_) / (x_max_ - x_min_) * n_elems_;
		if (ii < 0 || ii >= n_elems_) {
			std::cout << "Out of bounds in add: " << x << std::endl;
			return;
		}
		(*this)[ii][0] += val;
	}

	void Normalize() {
		double sum = 0.0;
		for (int i = 0; i < NRows(); i++) {
			for (int j = 0; j < NColumns(); j++) {
				sum += (*this)[i][j];
			}
		}
		if (sum > 0.0) {
			for (int i = 0; i < NRows(); i++) {
				for (int j = 0; j < NColumns(); j++) {
					(*this)[i][j] /= sum;
				}
			}
		}
	}

	double Sum() {
		double sum = 0.0;
		for (int i = 0; i < NRows(); i++) {
			for (int j = 0; j < NColumns(); j++) {
				sum += (*this)[i][j];
			}
		}
		return sum;
	}

private:
	double x_min_;
	double x_max_;
	double y_min_;
	double y_max_;
	int n_elems_;
	double delta_x_;
	double delta_y_;
	bool one_dimensional_ = false;
};

void PlotMeshData(const std::vector<Mesh>& data);
void PlotMeshData(const Mesh& data, std::string filename);