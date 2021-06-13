#include "lib_png_utils.h"

PngImage createRed(int width, int height) {
	PngImage png_image(width, height);

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			png_image.set(i, j, 255, 0, 0);
		}
	}

	return png_image;
}

PngImage CreateWhiteImage(int width, int height) {
	PngImage im(width, height);
	im.Fill(255, 255, 255);
	return im;
}