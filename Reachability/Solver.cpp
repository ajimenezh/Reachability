#include "Solver.h"

#include "Matrix2D.h"
#include "Zonotope.h"
#include "Reachability.h"
#include "MarkovChainUtils.h"
#include <sstream>

std::vector<Reachability> Solver::CalcReachability(Mesh mesh, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& F, double u0, double mu) {
	std::vector<Reachability> reachability;

	int n = mesh.NElems();
	for (int i = 0; i < n * n; i++) {
		std::pair<int, int> p = mesh.IdToIndices(i);
		// std::pair<double, double> q = mesh.GetMidPoint(p.first, p.second);
		std::pair<double, double> q = mesh.GetPoint(p.first, p.second);

		double x_4 = q.first;
		double x_5 = q.second;

		Matrix2D f_star = F;
		f_star[0][0] = -c1 / s * x_4 - c2 * x_5 + c3 * u0;
		f_star[1][0] = (-1 - c4 / (s * s)) * x_4 - c5 / s * x_5 + c6 / s * u0;

		Zonotope zonotope({ q.first, q.second },
			{ {mesh.DeltaX(), 0.0},
			  {0.0, mesh.DeltaY()} });

		auto res =
			Reachability(zonotope, T, delta_t, A, B, f_star, u0, mu);

		reachability.push_back(res);
	}
	return reachability;
}

void Solver::UpdateReachabilityByInput(std::vector<Reachability>& reachability, Mesh mesh, double T,
	double delta_t, const Matrix2D& A, const Matrix2D& B, const Matrix2D& F, double u0, double mu) {
	int n = mesh.NElems();
	for (int i = 0; i < n * n; i++) {
		std::pair<int, int> p = mesh.IdToIndices(i);
		// std::pair<double, double> q = mesh.GetMidPoint(p.first, p.second);
		std::pair<double, double> q = mesh.GetPoint(p.first, p.second);

		double x_4 = q.first;
		double x_5 = q.second;

		Matrix2D f_star = F;
		f_star[0][0] = -c1 / s * x_4 - c2 * x_5 + c3 * u0;
		f_star[1][0] = (-1 - c4 / (s * s)) * x_4 - c5 / s * x_5 + c6 / s * u0;

		Zonotope zonotope({ q.first, q.second },
			{ {mesh.DeltaX(), 0.0},
			  {0.0, mesh.DeltaY()} });

		reachability[i].UpdateReachabilityInput(zonotope, T, delta_t, A, B, f_star, u0, mu);
	}
}


void Solver::Solve() {
	Matrix2D A(2, 2);

	A[0][0] = -c1 / s;
	A[0][1] = -c2;
	A[1][0] = -1 - c4 / (s * s);
	A[1][1] = -c5 / s;

	Matrix2D B(2, 1);

	B[0][0] = c3;
	B[1][0] = c6 / s;

	Matrix2D f_star(2, 1);

	//Zonotope zonotope({ x_1, x_2, x_3,x_4, x_5 },
	//	{ {1.5, 0.0, 0.0, 0.0, 0.0},
	//	  {0.0, 0.2, 0.0, 0.0, 0.0},
	//	  {0.0, 0.0, 0.01, 0.0, 0.0},
	//	  {0.0, 0.0, 0.0, 0.01, 0.0},
	//	  {0.0, 0.0, 0.0, 0.0, 0.01} });

	double T = 0.2;
	double delta_t = T / 10;

	double total_t = 3.0;
	int N = int(total_t / T);

	int n = 40;
	Mesh mesh(n, -0.4, 0.6, -0.4, 0.6);

	double u_max = 0.038;

	std::vector<Reachability> reachability = CalcReachability(mesh, T, delta_t, A, B, f_star, u_max, u_max / 10);

	int cnt = 0;
	std::ofstream myfile;
	myfile.open("output2.txt");
	for (const auto& it : reachability[0].R_) if (it.Center().size() > 0) {
		//if (cnt++ > 5) break;
		std::vector<double> vertices = it.GetVertices();
		myfile << vertices.size() / 2 << " ";
		for (int i = 0; i < vertices.size(); i++) {
			myfile << vertices[i] << " ";
		}
		myfile.flush();
		myfile << "\n";
	}
	//for (const auto& it : reachability[0].R_up_) if (it.Center().size() > 0) {
	//	//if (cnt++ > 5) break;
	//	std::vector<double> vertices = it.GetVertices();
	//	myfile << vertices.size() / 2 << " ";
	//	for (int i = 0; i < vertices.size(); i++) {
	//		myfile << vertices[i] << " ";
	//	}
	//	myfile.flush();
	//	myfile << "\n";
	//}
	myfile.close();

	// Initial conditions
	mesh.set(0.0, 0.0);

	Mesh x3_mesh(40, -0.4, 0.4);
	x3_mesh.set(0.0);

	Mesh position_mesh(80, 0, 80, -4, 4);
	for (double y0 = 3; y0 <= 6; y0 += position_mesh.DeltaY()) {
		for (double x0 = -1.8; x0 >= -2.2; x0 -= position_mesh.DeltaX()) {
			position_mesh.set(y0, x0);
		}
	}
	position_mesh.Normalize();

	{
		std::stringstream filename3;
		filename3 << "test_pos_start_" << ".png";
		PlotMeshData(position_mesh, filename3.str());
	}

	Mesh cur_state = mesh;
	std::vector<Mesh> states = { cur_state };
	double u0 = u_max;
	for (int i = 0; i < N + 2; i++) {
		//UpdateReachabilityByInput(reachability, mesh, T, delta_t, A, B, f_star, -u_max, u_max / 10);
		//if (i < N / 4 || (i >= 3 * N / 4 && i < N)) {
		//	if (u0 != u_max) {
		//		UpdateReachabilityByInput(reachability, mesh, T, delta_t, A, B, f_star, u_max, u_max / 10);
		//	}
		//	u0 = u_max;
		//}
		//if (i >= N / 4 && i < 3 * N / 4) {
		//	if (u0 != u_max) {
		//		UpdateReachabilityByInput(reachability, mesh, T, delta_t, A, B, f_star, -u_max, u_max / 10);
		//	}
		//	u0 = -u_max;
		//}
		//if (i >= N) {
		//	if (u0 != 0.0) {
		//		UpdateReachabilityByInput(reachability, mesh, T, delta_t, A, B, f_star, 0.0, u_max / 10);
		//	}
		//	u0 = 0.0;
		//}

		MarkovChain markov_chain = CalculateMarkovChainProbabilities(
			reachability, Scenario(), cur_state);

		std::stringstream filename;
		filename << "test_p_" << i << ".png";
		PlotMeshData(markov_chain.p_, filename.str());

		std::stringstream filename2;
		filename2 << "test_p_t_" << i << ".png";
		PlotMeshData(markov_chain.p_t_, filename2.str());

		Mesh new_x3(x3_mesh);
		new_x3.Init(0.0);
		Mesh new_position(position_mesh);
		new_position.Init(0.0);

		for (int x4 = 0; x4 < markov_chain.p_.NElems(); x4++) {
			for (int x5 = 0; x5 < markov_chain.p_.NElems(); x5++) if (markov_chain.p_[x4][x5] > 0.0) {
				double x4_double = cur_state.GetPoint(x4, x5).first;
				for (int x3 = 0; x3 < x3_mesh.NElems(); x3++) if (x3_mesh[x3][0] > 0.0) {
					double x3_double = x3_mesh.GetPoint(x3).first;
					double tmp = x3_double + T * x4_double;
					new_x3.add(tmp, x3_mesh[x3][0] * markov_chain.p_[x4][x5]);
				}
			}
		}
		new_x3.Normalize();
		x3_mesh = new_x3;
		{
			std::stringstream filename3;
			filename3 << "test_ang_" << i << ".png";
			PlotMeshData(x3_mesh, filename3.str());
		}

		double norm_x3 = 0.0;
		for (int x3 = 0; x3 < x3_mesh.NElems(); x3++) if (x3_mesh[x3][0] > 0.0) {
			norm_x3 += x3_mesh[x3][0];
		}
		std::cout << "Norm x3: " << norm_x3 << std::endl;

		for (int x1 = 0; x1 < position_mesh.NElems(); x1++) {
			for (int x2 = 0; x2 < position_mesh.NElems(); x2++) if (position_mesh[x1][x2] > 0.0) {
				double x1_double = position_mesh.GetPoint(x1, x2).first;
				double x2_double = position_mesh.GetPoint(x1, x2).second;
				for (int x3 = 0; x3 < x3_mesh.NElems(); x3++) if (x3_mesh[x3][0] > 0.0) {
					double x3_double = x3_mesh.GetPoint(x3).first;
					double new_x = x1_double + s * cos(x3_double) * T;
					double new_y = x2_double + s * sin(x3_double) * T;
					new_position.add(new_x, new_y, position_mesh[x1][x2] * x3_mesh[x3][0]);
				}
			}
		}
		std::cout << "Norm pos: " << position_mesh.Sum() << std::endl;

		position_mesh = new_position;

		double norm_pos = position_mesh.Sum();
		std::cout << "Norm pos: " << norm_pos << std::endl;

		std::stringstream filename3;
		filename3 << "test_pos_" << i << ".png";
		PlotMeshData(new_position, filename3.str());

		std::cout << "Iteration " << i << std::endl;

		cur_state = markov_chain.p_t_;
		states.push_back(cur_state);
	}

	//PlotMeshData(states);
}
