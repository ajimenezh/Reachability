#include "stdafx.h"

#include "Reachability.h"
#include "Matrix2D.h"
#include "../LibPngUtils/PngImage.h"
#include "Zonotope.h"
#include "Solver.h"

#include <CGAL/Polygon_2.h>

#include <numeric>
#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <fstream>

#include <assert.h>
//#include "CommandlineOptions.hh"
//#include "Array.hh"
//#include "VPolytope.hh"
//#include "VPolytopeList.hh"
//#include "IntegerMatrix.hh"
//#include "LPinter.hh"
//#include "RSinter.hh"
//
//#include "zonotope_c.h"

void SolveCGAL() {
	Polygon_2   P;
	P.push_back(Point_2(1.1, -0.1));
	P.push_back(Point_2(1.1, 0.1));
	P.push_back(Point_2(0.9, 0.1));
	P.push_back(Point_2(0.9, -0.1));

	Matrix2D A(2, 2);
	A[0][0] = -1;
	A[0][1] = -4;
	A[1][0] = 4;
	A[1][1] = -1;

	const auto res = CalcReachability(P, 1, 0.005, A, 0.01);

	std::ofstream myfile;
	myfile.open("output.txt");
	myfile << P << "\n";
	//myfile << Q << "\n";
	for (const auto& it : res) {
		myfile << it << "\n";
	}
	myfile.close();
}

//void SolveMinkSum() {
//	dd_set_global_constants();  // First, this must be called to use cdd
//
//	Matrix m(2, 4);
//	Reachability::assign(m[0], { 1.1, -0.1 });
//	Reachability::assign(m[1], { 1.1, 0.1 });
//	Reachability::assign(m[2], { 0.9, 0.1 });
//	Reachability::assign(m[3], { 0.9, -0.1 });
//
//	Reachability::Matrix A(2, 2);
//	A[0][0] = -1;
//	A[0][1] = -4;
//	A[1][0] = 4;
//	A[1][1] = -1;
//
//	const auto& res = Reachability::Reachability_2(m, 1, 0.01, A, 0.05);
//
//	std::ofstream myfile;
//	myfile.open("output.txt");
//	//myfile << P << "\n";
//	//myfile << Q << "\n";
//	for (const auto& it : res) {
//		myfile << it.coldim() << " ";
//		for (int j = 0; j < it.coldim(); j++) {
//			for (int i = 0; i < it.rowdim(); i++) {
//				myfile << it[j][i].get_d() << " ";
//			}
//		}
//		myfile << "\n";
//	}
//	myfile.close();
//
//	dd_free_global_constants();  // Nice
//}

//void SolveZonotope() {
//	Reachability::Matrix A(2, 2);
//	A[0][0] = -1;
//	A[0][1] = -4;
//	A[1][0] = 4;
//	A[1][1] = -1;
//
//	Zonotope zonotope({1.0, 0.0}, { {0.1, 0.0}, { 0.0, 0.1 } });
//
//	const auto& res = Reachability::ReachabilityZonotope(zonotope, 1, 0.01, A, 0.05);
//
//	std::ofstream myfile;
//	myfile.open("output.txt");
//	for (const auto& it : res) {
//		std::vector<double> vertices = it.GetVertices();
//		myfile << vertices.size()/2 << " ";
//		for (int i = 0; i < vertices.size(); i++) {
//			myfile << vertices[i] << " ";
//		}
//		myfile.flush();
//		myfile << "\n";
//	}
//	myfile.close();
//}

void Solve2() {
	int n = 5;
	Matrix2D A(n, n);

	double s = 15.0;
	double c1 = 160.0;
	double c2 = -1.6;
	double c3 = 53.0;
	double c4 = -3.5;
	double c5 = 156.0;
	double c6 = 78.0;

	double x_1 = 4.5;
	double x_2 = 2.0;
	double x_3 = 0.0;
	double x_4 = 0.0;
	double x_5 = 0.0;

	A[0][2] = -s * sin(x_3);
	A[1][2] = s * cos(x_3);
	A[2][3] = 1;
	A[3][3] = -c1 / s;
	A[3][4] = -c2;
	A[4][3] = -1 - c4 / (s * s);
	A[4][4] = -c5 / s;

	Matrix2D B(5, 1);

	B[0][0] = c3;
	B[1][0] = c6 / s;

	Matrix2D f_star(n, 1);

	f_star[0][0] = s * cos(x_3);
	f_star[1][0] = s * sin(x_3);
	f_star[2][0] = x_4;
	f_star[3][0] = -c1 / s * x_4 - c2 * x_5;
	f_star[4][0] = (-1 - c4 / (s * s)) * x_4 - c5 / s * x_5;

	Zonotope zonotope({ x_1, x_2, x_3,x_4, x_5 },
		{ {1.5, 0.0, 0.0, 0.0, 0.0},
		  {0.0, 0.2, 0.0, 0.0, 0.0},
		  {0.0, 0.0, 0.01, 0.0, 0.0},
		  {0.0, 0.0, 0.0, 0.01, 0.0}, 
		  {0.0, 0.0, 0.0, 0.0, 0.01} });

	double T = 0.8;
	double delta_t = 0.05;

	auto res = 
		Reachability(zonotope, T, delta_t, A.getSubmatrix(2, 3, 3), B.getSubmatrix(2, 3, 0), f_star.getSubmatrix(2, 3, 0), 0.038, 0.005);

	//ZonotopesToMatrix(res, Scenario());

	//std::ofstream myfile;
	//myfile.open("output2.txt");
	//for (const auto& it : res) {
	//	std::vector<double> vertices = it.GetVertices();
	//	myfile << vertices.size() / 2 << " ";
	//	for (int i = 0; i < vertices.size(); i++) {
	//		myfile << vertices[i] << " ";
	//	}
	//	myfile.flush();
	//	myfile << "\n";
	//}
	//myfile.close();

}

int main() {

	//SolveCGAL();
	//SolveMinkSum();

	//SolveZonotope();
	Solver().Solve();

	//jl_init();

	//system("\"C:\\Users\\Alex-User\\AppData\\Local\\Programs\\Julia 1.5.3\\bin\\julia\" main.jl");

	//jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 2);
	//jl_array_t* x = jl_alloc_array_2d(array_type, 10, 5);

	//// Get array pointer
	//double* p = (double*)jl_array_data(x);
	//// Get number of dimensions
	//int ndims = jl_array_ndims(x);
	//// Get the size of the i-th dim
	//size_t size0 = jl_array_dim(x, 0);
	//size_t size1 = jl_array_dim(x, 1);

	//// Fill array with data
	//for (size_t i = 0; i < size1; i++)
	//	for (size_t j = 0; j < size0; j++)
	//		p[j + size0 * i] = i + j;

	/* strongly recommended: notify Julia that the
		 program is about to terminate. this allows
		 Julia time to cleanup pending write requests
		 and run all finalizers
	*/
	//jl_atexit_hook(0);

	// Construct the triangle.
	//Polygon_2   P;
	//P.push_back(Point_2(1.1, -0.1));
	//P.push_back(Point_2(1.1, 0.1));
	//P.push_back(Point_2(0.9, 0.1));
	//P.push_back(Point_2(0.9, -0.1));

	//Matrix A(2, 2);
	//A[0][0] = -1;
	//A[0][1] = -4;
	//A[1][0] = 4;
	//A[1][1] = -1;

	//const auto& res = Reachability(P, 1, 0.01, A, 0.05);

	//std::ofstream myfile;
	//myfile.open("output.txt");
	//myfile << P << "\n";
	////myfile << Q << "\n";
	//for (const auto& it : res) {
	//	myfile << it << "\n";
	//}
	//myfile.close();

	//CommandlineOptions::init(argc, argv);

	//dd_set_global_constants();  // First, this must be called to use cdd

	//// Read the polytopes

	//VPolytopeList polylist;

	//Matrix m(3, 8);
	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 8; j++) {
	//		m[j][i] = -1 + 2 * ((j / (1 << i)) % 2);
	//	}
	//}
	//VPolytope(m).pretty_print(std::cout);
	//polylist.append(VPolytope(m));
	//polylist.append(VPolytope(m));

	//polylist.incMinkSumSetup();

	//std::cout << "[";
	//bool first = true;
	//while (polylist.hasMoreVertices()) {
	//	if (!first)
	//		std::cout << "," << std::endl;
	//	std::cout << polylist.incExploreStepFast().coord();
	//	first = false;
	//}
	//std::cout << "]" << std::endl;

	//dd_free_global_constants();  // Nice

	return 0;
}
