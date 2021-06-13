#include "MarkovChainUtils.h"
#include "assert.h"

double CalcVolume(const Mesh& mesh) {
	int n = mesh.NElems();

	double res = 0.0;
	// double cell_vol = mesh.DeltaX() * mesh.DeltaY();
	double cell_vol = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res += cell_vol * mesh[i][j];
		}
	}

	return res;
}

MarkovChain CalculateMarkovChainProbabilities(
	const std::vector<Reachability>& reachability, const Scenario& scenario, const Mesh& previous_state) {

	int n = reachability.size();

	assert(n == scenario.M()* scenario.M());

	// Transition probabilities of state dependant Markov Chain. (can be calculated offline)
	Matrix2D phi_star(n, n);

	for (int j = 0; j < n; j++) {
		Zonotope R = MinkowskiSum(reachability[j].R_up_t_, reachability[j].R_down_t_);

		Mesh mesh = ZonotopesToMatrix({ R }, previous_state, /*join_polygons=*/false);

		double vol = CalcVolume(mesh);

		if (j == 0) {
			std::stringstream filename;
			filename << "test_phi_bar_" << ".png";
			PlotMeshData(mesh, filename.str());

			int cnt = 0;
			std::ofstream myfile;
			myfile.open("output2.txt");
			for (const auto& it : { R }) if (it.Center().size() > 0) {
				//if (cnt++ > 5) break;
				std::vector<double> vertices = it.GetVertices();
				myfile << vertices.size() / 2 << " ";
				for (int i = 0; i < vertices.size(); i++) {
					myfile << vertices[i] << " ";
				}
				myfile.flush();
				myfile << "\n";
			}
		}

		if (vol == 0.0) continue;

		for (int i = 0; i < n; i++) {
			std::pair<int, int> p = previous_state.IdToIndices(i);
			phi_star[i][j] = mesh[p.first][p.second] / vol;
		}
	}

	Matrix2D phi(n, n);

	for (int j = 0; j < n; j++) {
		Mesh mesh = ZonotopesToMatrix(reachability[j].R_, previous_state, /*join_polygons=*/false);

		if (j == 0) {
			std::stringstream filename;
			filename << "test_phi_bar_" << ".png";
			PlotMeshData(mesh, filename.str());

			int cnt = 0;
			std::ofstream myfile;
			myfile.open("output2.txt");
			for (const auto& it : reachability[j].R_) if (it.Center().size() > 0) {
				//if (cnt++ > 5) break;
				std::vector<double> vertices = it.GetVertices();
				myfile << vertices.size() / 2 << " ";
				for (int i = 0; i < vertices.size(); i++) {
					myfile << vertices[i] << " ";
				}
				myfile.flush();
				myfile << "\n";
			}
		}

		double vol = CalcVolume(mesh);

		if (vol == 0.0) continue;

		for (int i = 0; i < n; i++) {
			std::pair<int, int> p = previous_state.IdToIndices(i);
			phi[i][j] = mesh[p.first][p.second] / vol;
		}
	}

	// Transition probabilities of input dependant Markov Chain.
	Matrix2D phi_bar(n, n);

	for (int j = 0; j < n; j++) {
		std::pair<int, int> idx = previous_state.IdToIndices(j);
		std::pair<double, double> p = previous_state.GetMidPoint(idx.first, idx.second);
		Zonotope R = MinkowskiSum(reachability[j].R_bar_t_,
			Zonotope({ p.first, p.second }, { {previous_state.DeltaX(), 0.0}, {0.0, previous_state.DeltaY()} }));
		Mesh mesh = ZonotopesToMatrix({ R }, previous_state, /*join_polygons=*/false);

		if (j == 0) {
			std::stringstream filename;
			filename << "test_phi_bar_" << ".png";
			PlotMeshData(mesh, filename.str());

			int cnt = 0;
			std::ofstream myfile;
			myfile.open("output2.txt");
			for (const auto& it : { R }) if (it.Center().size() > 0) {
				//if (cnt++ > 5) break;
				std::vector<double> vertices = it.GetVertices();
				myfile << vertices.size() / 2 << " ";
				for (int i = 0; i < vertices.size(); i++) {
					myfile << vertices[i] << " ";
				}
				myfile.flush();
				myfile << "\n";
			}
		}

		double vol = CalcVolume(mesh);

		if (vol == 0.0) continue;

		for (int i = 0; i < n; i++) {
			std::pair<int, int> p = previous_state.IdToIndices(i);
			phi_bar[i][j] = mesh[p.first][p.second] / vol;
		}
	}

	MarkovChain res;

	Mesh tmp1 = Mesh(previous_state);
	for (int i = 0; i < tmp1.NElems(); i++) {
		for (int j = 0; j < tmp1.NElems(); j++) {
			tmp1[i][j] = phi_bar[i * tmp1.NElems() + j][0];
		}
	}
	{
		std::stringstream filename;
		filename << "test_phi_bar_" << ".png";
		PlotMeshData(tmp1, filename.str());
	}
	Mesh tmp2 = Mesh(previous_state);
	for (int i = 0; i < tmp1.NElems(); i++) {
		for (int j = 0; j < tmp1.NElems(); j++) {
			tmp2[i][j] = phi_star[i * tmp1.NElems() + j][0];
		}
	}
	{
		std::stringstream filename;
		filename << "test_phi_star_" << ".png";
		PlotMeshData(tmp2, filename.str());
	}


	res.p_t_ = Mesh(previous_state);
	Matrix2D sum(n, 1);
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < n; j++) {
			std::pair<int, int> p = previous_state.IdToIndices(j);
			sum[k][0] += phi_bar[k][j] * previous_state[p.first][p.second];
		}
	}
	for (int i = 0; i < n; i++) {
		// p[i] = phi_star[i][k] * phi_bar[k][j] * previus_state[j];
		std::pair<int, int> q = previous_state.IdToIndices(i);
		res.p_t_[q.first][q.second] = 0.0;
		for (int k = 0; k < n; k++) {
			res.p_t_[q.first][q.second] += phi_star[i][k] * sum[k][0];
		}
	}

	res.p_ = Mesh(previous_state);
	for (int i = 0; i < n; i++) {
		// p[i] = phi_[i][j] * previus_state[j];
		std::pair<int, int> q = previous_state.IdToIndices(i);
		res.p_[q.first][q.second] = 0.0;
		for (int j = 0; j < n; j++) {
			std::pair<int, int> p = previous_state.IdToIndices(j);
			res.p_[q.first][q.second] += phi[i][j] * previous_state[p.first][p.second];
		}
	}

	return res;
}