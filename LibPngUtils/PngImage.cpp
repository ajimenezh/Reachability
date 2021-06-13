#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "PngImage.h"

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <png.h>
#include <string.h>
#include <vector>

int PngImage::writeImage(char* filename, char* title) {
	int code = 0;
	FILE* fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		finalise(fp, png_ptr, info_ptr, row);
		return 1;
	}

	//// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
		finalise(fp, png_ptr, info_ptr, row);
		return 1;
	}

	//// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
		finalise(fp, png_ptr, info_ptr, row);
		return 1;
	}

	//// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
		finalise(fp, png_ptr, info_ptr, row);
		return 1;
	}

	png_init_io(png_ptr, fp);

	//// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width_, height_,
		8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	//// Set title
	if (title != NULL) {
		png_text title_text;
		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		title_text.key = _strdup("Title");
		title_text.text = title;
		png_set_text(png_ptr, info_ptr, &title_text, 1);
	}

	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = (png_bytep)malloc(3 * width_ * sizeof(png_byte));

	// Write image data
	int x, y;
	for (y = 0; y < height_; y++) {
		for (x = 0; x < width_; x++) {
			setRGB(&(row[x * 3]), get(x, y));
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	finalise(fp, png_ptr, info_ptr, row);

	return 0;
}

void PngImage::finalise(FILE* fp, png_structp png_ptr, png_infop info_ptr, png_bytep row) {
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);
}

void PngImage::FillPoints(const std::vector<Point>& points, int size, int r, int g, int b) {
	for (const auto& p : points) {
		int xx = GetX(p.x);
		int yy = GetY(p.y);

		for (int j = -(size + 1) / 2; j <= (size + 1) / 2; j++) if (xx + j >= 0 && xx + j < width_)  {
			for (int i = -((size + 1) / 2 - abs(j)); i <= ((size + 1) / 2 - abs(j)); i++) {
				if (yy + i >= 0 && yy + i < height_) {
					set(xx + j, yy + i, r, g, b);
				}
			}
		}
	}
}

void PngImage::DrawLine(const Point& point, double m, int size, int r, int g, int b) {
	double lx = (xl_ - xr_) / (width_);
	double ly = (yl_ - yr_) / (height_);
	for (int i = 0; i < size; i++) {
		double x = point.x + cos(m) * lx * i;
		double y = point.y + sin(m) * ly * i;
		FillPoints({ { x, y } }, 1, r, g, b);
	}
}

int PngImage::GetY(double y) {
	return height_ - (int)((y - yl_) / (yr_ - yl_) * (height_ - 1)) - 1;
}

int PngImage::GetX(double x) {
	return (int)((x - xl_) / (xr_ - xl_) * (width_ - 1));
}

void PngImage::DrawPoint(const Point& point, int size, int r, int g, int b) {
	for (int i = -size; i <= size; i++) {
		int xc = GetX(point.x) + i;
		int yc = GetY(point.y);
		if (xc >= 0 && xc < width_ && yc >= 0 && yc < height_) {
			set(xc, yc, r, g, b);
		}

		xc = GetX(point.x);
		yc = GetY(point.y) + i;
		if (xc >= 0 && xc < width_ && yc >= 0 && yc < height_) {
			set(xc, yc, r, g, b);
		}
	}
}

int PngImage::ReadPngInit(FILE* infile) {
	unsigned char sig[8];

	fread(sig, 1, 8, infile);
	if (!png_check_sig(sig, 8))
		return 1;   /* bad signature */

	png_struct* png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,
		NULL);
	if (!png_ptr)
		return 4;   /* out of memory */

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_read_struct(&png_ptr, NULL, NULL);
		return 4;   /* out of memory */
	}

	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 2;
	}

	png_init_io(png_ptr, infile);
	png_set_sig_bytes(png_ptr, 8);
	png_read_info(png_ptr, info_ptr);

	unsigned int width, height;
	int bit_depth, color_type;

	width = png_get_image_width(png_ptr, info_ptr);
	height = png_get_image_height(png_ptr, info_ptr);
	color_type = png_get_color_type(png_ptr, info_ptr);
	bit_depth = png_get_bit_depth(png_ptr, info_ptr);

	width_ = width;
	height_ = height;

	int number_of_passes = png_set_interlace_handling(png_ptr);
	
	if (bit_depth == 16)
		png_set_strip_16(png_ptr);

	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png_ptr);

	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand_gray_1_2_4_to_8(png_ptr);

	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png_ptr);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if (color_type == PNG_COLOR_TYPE_RGB ||
		color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);

	if (color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png_ptr);

	png_read_update_info(png_ptr, info_ptr);

	/* read file */
	//if (setjmp(png_jmpbuf(png_ptr)))
	//	return 4;

	png_bytep* row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
	for (int y = 0; y < height_; y++)
		row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png_ptr, info_ptr));

	png_read_image(png_ptr, row_pointers);

	data_.resize(height);
	int sum = 0;
	for (int y = 0; y < height_; y++) {
		data_[y].resize(width_);
		for (int x = 0; x < width_; x++) {
			data_[y][x] = row_pointers[y][4 * x + 2] * 65536 + row_pointers[y][4 * x + 1] * 256 + row_pointers[y][4 * x];
		}

		free(row_pointers[y]);
	}

	free(row_pointers);

	return 0;
}

PngImage ReadPngImage(const char* filename) {
	FILE* pFile = ReadFile(_strdup(filename));

	PngImage image;
	image.ReadPngInit(pFile);

	fclose(pFile);

	return image;
}

FILE* ReadFile(const char* filename) {
	FILE* pFile;

	pFile = fopen(filename, "rb");
	if (pFile == NULL) perror("Error opening file");

	return pFile;
}

bool IsWhite(float val) {
	return val >= ((1 << 24) - 1);
}

PngImage createWhite(int width, int height) {
	PngImage png_image(width, height);

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			png_image.set(i, j, 255, 255, 255);
		}
	}

	return png_image;
}