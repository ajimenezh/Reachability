#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <png.h>
#include <string.h>
#include <vector>
#include "PngImage.h"

// Creates a test image for saving. Creates a Mandelbrot Set fractal of size width x height
PngImage createMandelbrotImage(int width, int height, float xS, float yS, float rad, int maxIteration);

PngImage createRed(int width, int height);

// This takes the float value 'val', converts it to red, green & blue values, then 
// sets those values into the image memory buffer location pointed to by 'ptr'
inline void setRGB(png_byte* ptr, float val);

// This function actually writes out the PNG image file. The string 'title' is
// also written into the image file
int writeImage(char* filename, int width, int height, const PngImage& buffer, char* title);


int main(int argc, char* argv[])
{
	// Specify an output image size
	int width = 500;
	int height = 300;

	// Create a test image - in this case a Mandelbrot Set fractal
	// The output is a 1D array of floats, length: width * height
	printf("Creating Image\n");
	//PngImage image = createMandelbrotImage(width, height, -0.802, -0.177, 0.011, 110);
	PngImage image = createRed(width, height);

	// Save the image to a PNG file
	// The 'title' string is stored as part of the PNG file
	printf("Saving PNG\n");
	int result = image.writeImage(_strdup("test.png"), _strdup("This is my test image"));

	return result;
}

PngImage createMandelbrotImage(int width, int height, float xS, float yS, float rad, int maxIteration) {
	//float* buffer = (float*)malloc(width * height * sizeof(float));
	//if (buffer == NULL) {
	//	fprintf(stderr, "Could not create image buffer\n");
	//	return NULL;
	//}
	PngImage png_image(width, height);

	// Create Mandelbrot set image

	int xPos, yPos;
	float minMu = maxIteration;
	float maxMu = 0;

	for (yPos = 0; yPos < height; yPos++)
	{
		float yP = (yS - rad) + (2.0f * rad / height) * yPos;

		for (xPos = 0; xPos < width; xPos++)
		{
			float xP = (xS - rad) + (2.0f * rad / width) * xPos;

			int iteration = 0;
			float x = 0;
			float y = 0;

			while (x * x + y * y <= 4 && iteration < maxIteration)
			{
				float tmp = x * x - y * y + xP;
				y = 2 * x * y + yP;
				x = tmp;
				iteration++;
			}

			if (iteration < maxIteration) {
				float modZ = sqrt(x * x + y * y);
				float mu = iteration - (log(log(modZ))) / log(2);
				if (mu > maxMu) maxMu = mu;
				if (mu < minMu) minMu = mu;
				//buffer[yPos * width + xPos] = mu;
				png_image.set(xPos, yPos, mu);
			}
			else {
				//buffer[yPos * width + xPos] = 0;
				png_image.set(xPos, yPos, 0);
			}
		}
	}

	png_image.scale(minMu, maxMu);

	return png_image;
}

PngImage createRed(int width, int height) {
	PngImage png_image(width, height);

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			png_image.set(i, j, 255, 0, 0);
		}
	}

	return png_image;
}