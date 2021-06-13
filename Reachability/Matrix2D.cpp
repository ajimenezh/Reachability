#include "Matrix2D.h"
#include "PngImage.h"
#include <string>

void PlotMeshData(const std::vector<Mesh>& data) {
	int width = 500;
	int height = 500;
	PngImage image = createWhite(width, height);

	int result = image.writeImage(_strdup("test_1.png"), _strdup("This is my test image"));
}

void PlotMeshData(const Mesh& data, std::string filename) {
	int K = 20;
	int width = data.NRows() * K;
	int height = data.NColumns() * K;
	PngImage image = createWhite(width, height);

	double max_val = 1.0e-4;
	for (int i = 0; i < data.NRows(); i++) {
		for (int j = 0; j < data.NColumns(); j++) {
			max_val = std::max(max_val, data[i][j]);
		}
	}

	for (int i = 0; i < data.NRows(); i++) {
		for (int j = 0; j < data.NColumns(); j++) {
			for (int ii = 0; ii < K; ii++) {
				for (int jj = 0; jj < K; jj++) {
					image.set(i * K + ii, j * K + jj, 255, 255 * (max_val - data[i][j]) / max_val, 255);
				}
			}
		}
	}

	int result = image.writeImage(&filename[0], _strdup(""));
}