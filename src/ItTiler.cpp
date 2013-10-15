//
//  ItTiler.cpp
//  pixelbridge
//
//  Created by Dave Estes on 10/13/13.
//  Copyright (c) 2013 Dave Estes. All rights reserved.
//

#include <cmath>
#include <iostream>

#include "ItTiler.h"
#include "PixelBridgeFeatures.h"

#define DEFAULT_QUALITY        1

/**
 * Using a #define because it's used for floats and ints.
 */
#define CLAMP(x, min, max)   (x < min ? min : x >= max ? max : x)

ItTiler::ItTiler(BaseNddiDisplay* display, bool quiet)
: display_(display),
  quiet_(quiet)
{
    initZigZag();
    setQuality(DEFAULT_QUALITY);
}

void ItTiler::initZigZag() {
    // TODO(CDE): Right now this is a simple one-to-one mapping. Sort out a proper ordering.
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        zigZag_[i] = i;
    }
}

void ItTiler::setQuality(uint32_t quality) {
    static uint16_t v[6][3] = {
        10, 16, 13,
        11, 18, 14,
        13, 20, 16,
        14, 23, 18,
        16, 25, 20,
        18, 29, 23
    };
    static uint16_t m[6][3] = {
        13107, 5243, 8066,
        11916, 4660, 7490,
        10082, 4194, 6554,
        9362, 3547, 5825,
        8192, 3355, 5243,
        7282, 2893, 4559
    };
    
    // Don't adjust qp if it wrapped below zero or above 100
    if (quality > 100) return;
    
    // Set qp (the quality)
    qp = quality;
    
    // Initialize Mf4 and Vi4
    Mf4[0] = m[qp % 6][0];
    Mf4[1] = m[qp % 6][2];
    Mf4[2] = m[qp % 6][0];
    Mf4[3] = m[qp % 6][2];
    
    Mf4[4] = m[qp % 6][2];
    Mf4[5] = m[qp % 6][1];
    Mf4[6] = m[qp % 6][2];
    Mf4[7] = m[qp % 6][1];
    
    Mf4[8] = m[qp % 6][0];
    Mf4[9] = m[qp % 6][2];
    Mf4[10] = m[qp % 6][0];
    Mf4[11] = m[qp % 6][2];
    
    Mf4[12] = m[qp % 6][2];
    Mf4[13] = m[qp % 6][1];
    Mf4[14] = m[qp % 6][2];
    Mf4[15] = m[qp % 6][1];
    
    Vi4[0] = v[qp % 6][0];
    Vi4[1] = v[qp % 6][2];
    Vi4[2] = v[qp % 6][0];
    Vi4[3] = v[qp % 6][2];
    
    Vi4[4] = v[qp % 6][2];
    Vi4[5] = v[qp % 6][1];
    Vi4[6] = v[qp % 6][2];
    Vi4[7] = v[qp % 6][1];
    
    Vi4[8] = v[qp % 6][0];
    Vi4[9] = v[qp % 6][2];
    Vi4[10] = v[qp % 6][0];
    Vi4[11] = v[qp % 6][2];
    
    Vi4[12] = v[qp % 6][2];
    Vi4[13] = v[qp % 6][1];
    Vi4[14] = v[qp % 6][2];
    Vi4[15] = v[qp % 6][1];
}

void ItTiler::matrixMultiply(double *D, double *S1, double *S2) {
    double T[BLOCK_SIZE];
    for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
        for (size_t i = 0; i < BLOCK_WIDTH; i++) {
            T[j * BLOCK_WIDTH + i] =
            S1[j * BLOCK_WIDTH + 0] * S2[0 * BLOCK_WIDTH + i] +
            S1[j * BLOCK_WIDTH + 1] * S2[1 * BLOCK_WIDTH + i] +
            S1[j * BLOCK_WIDTH + 2] * S2[2 * BLOCK_WIDTH + i] +
            S1[j * BLOCK_WIDTH + 3] * S2[3 * BLOCK_WIDTH + i];
        }
    }
    memcpy(D, T, sizeof(T));
}

void ItTiler::hadamardMultiply(double *D, double *S1, double *S2) {
    double T[BLOCK_SIZE];
    for (size_t i = 0; i < BLOCK_WIDTH; i++) {
        for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
            T[j * BLOCK_WIDTH + i] =
            S1[j * BLOCK_WIDTH + i] * S2[j * BLOCK_WIDTH + i];
        }
    }
    memcpy(D, T, sizeof(T));
}

void ItTiler::scalarMultiply(double *D, double s1, double *S2) {
    double T[BLOCK_SIZE];
    for (size_t i = 0; i < BLOCK_WIDTH; i++) {
        for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
            T[j * BLOCK_WIDTH + i] =
            s1 * S2[j * BLOCK_WIDTH + i];
        }
    }
    memcpy(D, T, sizeof(T));
}

/*
 * Y = round(Cf4 . X . Cf4T o Mf4 . 1/2^(15 + floor(qp/6)))
 */
void ItTiler::forwardIntegerTransform(int *Y, int *X) {
    double YY[BLOCK_SIZE];
    double XX[BLOCK_SIZE];
    int shifter;
    
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        XX[i] = (double)X[i];
    }
    
    // 1) Y = Cf4 . X
    matrixMultiply(YY, Cf4, XX);
    // 2) Y = Y . Cf4T
    matrixMultiply(YY, YY, Cf4T);
    // 3) Y = Y o Mf4
    hadamardMultiply(YY, YY, Mf4);
    // 4) Y = Y . 1/2^(15 + floor(qp/6))
    shifter = 15 + qp / 6;
    scalarMultiply(YY, 1.0 / (double)(1 << shifter), YY);
    // 5) round
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        Y[i] = (int)floor(YY[i] + .5); // TODO(CDE): Sort out this rounding...likely by converting all of the math to integer math
    }
}

/*
 * Z = round(Ci4T . [Y o Vi4 . 2^floor(qp/6)] . Ci4 . 1/2^6
 */
void ItTiler::inverseIntegerTransform(int *Z, int *Y) {
    double ZZ[BLOCK_SIZE];
    double YY[BLOCK_SIZE];
    int shifter;
    
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        YY[i] = (double)Y[i];
    }
    
    // 1) Z = Y o Vi4
    hadamardMultiply(ZZ, YY, Vi4);
    // 2) Z = Z . 2^floor(qp/6)
    shifter = qp/6;
    scalarMultiply(ZZ, (double)(1 << shifter), ZZ);
    // 3) Z = Ci4T . Z
    matrixMultiply(ZZ, Ci4T, ZZ);
    // 4) Z = Z . Ci4
    matrixMultiply(ZZ, ZZ, Ci4);
    // 5) Z = Z . 1/2^6
    shifter = 6;
    scalarMultiply(ZZ, 1 / (double)(1 << shifter), ZZ);
    // 6) round
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        Z[i] = (int)floor(ZZ[i] + .5); // TODO(CDE): Sort out this rounding...likely by converting all of the math to integer math
    }
}

/**
 * Initializes the Coefficient Planes for this tiler. The coefficient matrices
 * for each plane will pick a plane from the frame volume.
 */
void ItTiler::InitializeCoefficientPlanes() {
    
	// Setup the coefficient matrix to complete 3x3 identity initially
	vector< vector<int> > coeffs;
    coeffs.resize(3);
    coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
    coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
    coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(0);
    
	// Setup start and end points to (0,0,0) initially
	vector<unsigned int> start, end;
	start.push_back(0); start.push_back(0); start.push_back(0);
	end.push_back(0); end.push_back(0); end.push_back(0);
    
	// Break the display into macroblocks and initialize each 4x4x48 cube of coefficients to pick out the proper block from the frame volume
	for (int k = 0; k < FRAMEVOLUME_DEPTH; k++) {
		for (int j = 0; j < (display_->DisplayHeight() / BLOCK_HEIGHT); j++) {
			for (int i = 0; i < (display_->DisplayWidth() / BLOCK_WIDTH); i++) {
				coeffs[2][0] = -i * BLOCK_WIDTH;
				coeffs[2][1] = -j * BLOCK_HEIGHT;
				coeffs[2][2] = k;
				start[0] = i * BLOCK_WIDTH; start[1] = j * BLOCK_HEIGHT; start[2] = k;
				end[0] = (i + 1) * BLOCK_WIDTH - 1; end[1] = (j + 1) * BLOCK_HEIGHT - 1; end[2] = k;
				if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
				if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }
				display_->FillCoefficientMatrix(coeffs, start, end);
			}
		}
	}
    
	// Fill each scaler in every plane with 0
	start[0] = 0; start[1] = 0; start[2] = 0;
	end[0] = display_->DisplayWidth() - 1;
	end[1] = display_->DisplayHeight() - 1;
    end[2] = NUM_COEFFICIENT_PLANES - 1;
    display_->FillScaler(0, start, end);
}

/**
 * Initializes the Frame Volume for this tiler by pre-rendering each
 * of the 16 basis functions into 4x4 planes in the Frame Volume. They're
 * rendered for each color channel and stored in those groups of three in
 * zig-zag order.
 */
void ItTiler::InitializeFrameVolume() {

	Pixel  *pixels;
	size_t  pixels_size = sizeof(Pixel)            // pixel size
                         * BLOCK_SIZE              // size of each x,y plane
                         * FRAMEVOLUME_DEPTH;      // number of basis functions for each color channel
    
	pixels = (Pixel *)malloc(pixels_size);
	memset(pixels, 0x00, pixels_size);
    
    int Y[BLOCK_SIZE];
    int Z[BLOCK_SIZE];
    
	// Pre-render each basis function
#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (int j = 0; j < BASIS_BLOCKS_TALL; j++) {
        for (int i = 0; i < BASIS_BLOCKS_WIDE; i++) {
            // Initialize the coefficients to all off except the special one
            memset(Y, 0x00, sizeof(Y));
            Y[j * BLOCK_WIDTH + i] = MAX_IT_COEFF;
            
            // Perform the Inverse Transform
            inverseIntegerTransform(Z, Y);
            
        	size_t p = zigZag_[j * BASIS_BLOCKS_WIDE + i] * BLOCK_SIZE * 3;

            for (int y = 0; y < BLOCK_HEIGHT; y++) {
                for (int x = 0; x < BLOCK_WIDTH; x++) {
                    
                    // Set the color channels
                    uint32_t m = Z[y * BLOCK_WIDTH + x];

                    // Set the color channels with the magnitude clamped to 127
                    pixels[p + BLOCK_SIZE * 0].r = m;
                    pixels[p + BLOCK_SIZE * 0].a = 0xff;
                    pixels[p + BLOCK_SIZE * 1].g = m;
                    pixels[p + BLOCK_SIZE * 1].a = 0xff;
                    pixels[p + BLOCK_SIZE * 2].b = m;
                    pixels[p + BLOCK_SIZE * 2].a = 0xff;
                    
                    // Move to the next pixel
                    p++;
                }
            }
        }
    }
    
    // Update the frame volume with the basis function renderings and gray block in bulk.
	vector<unsigned int> start, end;
	start.push_back(0); start.push_back(0); start.push_back(0);
	end.push_back(BLOCK_WIDTH - 1); end.push_back(BLOCK_HEIGHT - 1); end.push_back(FRAMEVOLUME_DEPTH - 1);
	display_->CopyPixels(pixels, start, end);
    
    // Free the pixel memory
    free(pixels);
}

void ItTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height) {
#define HACKADACK
#ifdef HACKADACK
	///////////////////DEBUG/////////////////////
	cout << "Update" << endl;
    
	vector<unsigned int> start(3, 0), end(3, 0);
    
	// Clear all scalers
	start[0] = 0; start[1] = 0; start[2] = 0;
	end[0] = display_->DisplayWidth() - 1;
	end[1] = display_->DisplayHeight() - 1;
    end[2] = NUM_COEFFICIENT_PLANES - 1;
    display_->FillScaler(0, start, end);
    
	// Then select the proper planes to just render the 64 basis functions
	for (int j = 0; j < BASIS_BLOCKS_TALL; j++) {
    	start[1] = 4 * j * BLOCK_HEIGHT; end[1] = 4 * ((j + 1) * BLOCK_HEIGHT - 1);
        for (int i = 0; i < BASIS_BLOCKS_WIDE; i++) {
        	start[0] = 4 * i * BLOCK_WIDTH; end[0] = 4 * ((i + 1) * BLOCK_WIDTH - 1);
        	size_t p = zigZag_[j * BASIS_BLOCKS_WIDE + i];
        	start[2] = p * 3; end[2] = (p + 1) * 3 - 1;
            display_->FillScaler(NUM_COEFFICIENT_PLANES, start, end);
        }
    }
	///////////////////DEBUG/////////////////////
#else // HACKADACK
    // TODO(CDE): Implement
#endif
}
