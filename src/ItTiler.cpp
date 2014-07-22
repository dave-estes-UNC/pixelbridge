//
//  ItTiler.cpp
//  pixelbridge
//
//  Created by Dave Estes on 10/13/13.
//  Copyright (c) 2013 Dave Estes. All rights reserved.
//

#include <cmath>
#include <iostream>

#include "PixelBridgeFeatures.h"
#include "Configuration.h"
#include "ItTiler.h"

/**
 * Using a #define because it's used for floats and ints.
 */
#define CLAMP(x, min, max)   (x < min ? min : x >= max ? max : x)
#define MAX(x, y)            (x > y ? x : y)
#define MIN(x, y)            (x < y ? x : y)
#define CEIL(x, y)           (1 + ((x - 1) / y))

/* Static Initialization */
int ItTiler::Cf4[] = {
        2,   2,   2,   2,
        4,   2,  -2,  -4,
        2,  -2,  -2,   2,
        2,  -4,   4,  -2
};
int ItTiler::Cf4T[] = {
        2,   4,   2,   2,
        2,   2,  -2,  -4,
        2,  -2,  -2,   4,
        2,  -4,   2,  -2
};
int ItTiler::Ci4[] = {
        2,   2,   2,   2,
        2,   1,  -1,  -2,
        2,  -2,  -2,   2,
        1,  -2,  2,  -1
};
int ItTiler::Ci4T[] = {
        2,   2,   2,   1,
        2,   1,  -2,  -2,
        2,  -1,  -2,   2,
        2,  -2,   2,  -1
};


ItTiler::ItTiler(size_t display_width, size_t display_height,
                 size_t quality)
: qp(0)
{
    quiet_ = globalConfiguration.headless || !globalConfiguration.verbose;

    // 3 dimensional matching the Macroblock Width x Height x 64+3+1
    vector<unsigned int> fvDimensions;
    fvDimensions.push_back(BLOCK_WIDTH);
    fvDimensions.push_back(BLOCK_HEIGHT);
    fvDimensions.push_back(FRAMEVOLUME_DEPTH);

#ifndef NO_CL
    display_ = new ClNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 FRAMEVOLUME_DEPTH,             // Number of coefficient planes
                                 3);                            // input vector size (x, y, 1)
#else
    display_ = new GlNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 FRAMEVOLUME_DEPTH,             // Number of coefficient planes
                                 3);                            // input vector size (x, y, 1)
#endif


    // Set the full scaler value
    display_->SetFullScaler(MAX_IT_COEFF);

    initZigZag();
    setQuality(quality);
    display_->SetPixelByteSignMode(SIGNED_MODE);

    // Initialize Input Vector
    vector<int> iv;
    iv.push_back(1);
    display_->UpdateInputVector(iv);

    // Initialize Coefficient Planes
    InitializeCoefficientPlanes();

    // Initialize Frame Volume
    InitializeFrameVolume();
}

/**
 * Returns the Display created and initialized by the tiler.
 */
GlNddiDisplay* ItTiler::GetDisplay() {
    return display_;
}

/**
 * See big-honking comment for DctTiler::initZigZag().
 */
void ItTiler::initZigZag() {
    size_t  x = 0, y = 0;
    bool    up = true;

    // Setup the zigZag_ table
    for (size_t i = 0; i < BLOCK_SIZE; i++) {
        zigZag_[y * BLOCK_WIDTH + x] = i;
        if (up) {
            if (x < (BLOCK_WIDTH - 1)) {
                x++;
                if (y > 0) {
                    y--;
                } else {
                    up = false;
                }
            } else {
                y++;
                up = false;
            }
        } else {
            if (y < (BLOCK_HEIGHT - 1)) {
                y++;
                if (x > 0) {
                    x--;
                } else {
                    up = true;
                }
            } else {
                x++;
                up = true;
            }
        }
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
    qp6 = quality / 6;

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

/*
 * Y = round(Cf4 . X . Cf4T o Mf4 . 1/2^(15 + floor(qp/6)))
 */
void ItTiler::forwardIntegerTransform(int *Y, int *X) {
    int32_t T[BLOCK_SIZE];

    // T = Cf4 . X
    for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
        for (size_t i = 0; i < BLOCK_WIDTH; i++) {
            T[j * BLOCK_WIDTH + i] =
            Cf4[j * BLOCK_WIDTH + 0] * X[0 * BLOCK_WIDTH + i] +
            Cf4[j * BLOCK_WIDTH + 1] * X[1 * BLOCK_WIDTH + i] +
            Cf4[j * BLOCK_WIDTH + 2] * X[2 * BLOCK_WIDTH + i] +
            Cf4[j * BLOCK_WIDTH + 3] * X[3 * BLOCK_WIDTH + i];
        }
    }
    // Y = T . Cf4T o Mf4 . 1/2^(15 + floor(qp/6))
    for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
        for (size_t i = 0; i < BLOCK_WIDTH; i++) {
            Y[j * BLOCK_WIDTH + i] =
            T[j * BLOCK_WIDTH + 0] * Cf4T[0 * BLOCK_WIDTH + i] +
            T[j * BLOCK_WIDTH + 1] * Cf4T[1 * BLOCK_WIDTH + i] +
            T[j * BLOCK_WIDTH + 2] * Cf4T[2 * BLOCK_WIDTH + i] +
            T[j * BLOCK_WIDTH + 3] * Cf4T[3 * BLOCK_WIDTH + i];
            Y[j * BLOCK_WIDTH + i] *= Mf4[j * BLOCK_WIDTH + i];
            // Shift right by 2 more to compensate for the x2 performed on Cf4 and Cf4T
            Y[j * BLOCK_WIDTH + i] >>= 15 + qp6 + 2;
        }
    }
}

/*
 * Z = round(Ci4T . [Y o Vi4 . 2^floor(qp/6)] . Ci4 . 1/2^6
 */
void ItTiler::inverseIntegerTransform(int *Z, int *Y) {
    int32_t T1[BLOCK_SIZE];
    int32_t T2[BLOCK_SIZE];

    // T1 = Y o Vi4 . 2^floor(qp/6)
    for (size_t i = 0; i < BLOCK_WIDTH; i++) {
        for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
            T1[j * BLOCK_WIDTH + i] =
            Y[j * BLOCK_WIDTH + i] * Vi4[j * BLOCK_WIDTH + i];
            T1[j * BLOCK_WIDTH + i] <<= qp6;
        }
    }
    // T2 = Ci4T . T1
    for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
        for (size_t i = 0; i < BLOCK_WIDTH; i++) {
            T2[j * BLOCK_WIDTH + i] =
            Ci4T[j * BLOCK_WIDTH + 0] * T1[0 * BLOCK_WIDTH + i] +
            Ci4T[j * BLOCK_WIDTH + 1] * T1[1 * BLOCK_WIDTH + i] +
            Ci4T[j * BLOCK_WIDTH + 2] * T1[2 * BLOCK_WIDTH + i] +
            Ci4T[j * BLOCK_WIDTH + 3] * T1[3 * BLOCK_WIDTH + i];
        }
    }
    // Z = T2 . Ci4 . 1/2^6
    for (size_t j = 0; j < BLOCK_HEIGHT; j++) {
        for (size_t i = 0; i < BLOCK_WIDTH; i++) {
            Z[j * BLOCK_WIDTH + i] =
            T2[j * BLOCK_WIDTH + 0] * Ci4[0 * BLOCK_WIDTH + i] +
            T2[j * BLOCK_WIDTH + 1] * Ci4[1 * BLOCK_WIDTH + i] +
            T2[j * BLOCK_WIDTH + 2] * Ci4[2 * BLOCK_WIDTH + i] +
            T2[j * BLOCK_WIDTH + 3] * Ci4[3 * BLOCK_WIDTH + i];
            // Shift right by 2 more to compensate for the x2 performed on Ci4 and Ci4T
            Z[j * BLOCK_WIDTH + i] >>= 6 + 2;
        }
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
    end[2] = display_->NumCoefficientPlanes() - 1;
    Scaler s;
    s.r = s.g = s.b = 0;
    display_->FillScaler(s, start, end);
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

    // Pre-render each basis function
#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (int j = 0; j < BASIS_BLOCKS_TALL; j++) {
        for (int i = 0; i < BASIS_BLOCKS_WIDE; i++) {

            int Y[BLOCK_SIZE];
            int Z[BLOCK_SIZE];

            // Initialize the coefficients to all off except the special one
            memset(Y, 0x00, sizeof(Y));
            Y[j * BLOCK_WIDTH + i] = MAX_IT_COEFF;

            // Perform the Inverse Transform
            inverseIntegerTransform(Z, Y);

            size_t p = zigZag_[j * BASIS_BLOCKS_WIDE + i] * BLOCK_SIZE;

            for (int y = 0; y < BLOCK_HEIGHT; y++) {
                for (int x = 0; x < BLOCK_WIDTH; x++) {

                    // Set the color channels
                    uint32_t m = CLAMP(Z[y * BLOCK_WIDTH + x], -127, 127);

                    // Set the color channels with the magnitude clamped to 127
                    pixels[p].r = m;
                    pixels[p].g = m;
                    pixels[p].b = m;
                    pixels[p].a = 0xff;

                    // Move to the next pixel
                    p++;
                }
            }
        }
    }

    // Update the frame volume with the basis function renderings in bulk.
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(BLOCK_WIDTH - 1); end.push_back(BLOCK_HEIGHT - 1); end.push_back(FRAMEVOLUME_DEPTH - 1);
    display_->CopyPixels(pixels, start, end);

    // Free the pixel memory
    free(pixels);
}

void ItTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height) {
    vector<unsigned int> start(3, 0), end(3, 0);
    vector<unsigned int> size(2, 0);
    size_t lastNonZeroPlane;
    static size_t largestNonZeroPlaneSeen = 0;
    Scaler s;

    size[0] = BLOCK_WIDTH;
    size[1] = BLOCK_HEIGHT;

    /* Clear all scalers */
    // TODO(CDE): Don't do this in the future. Just write enough planes to overwrite the last known non-zero plane.
    start[0] = 0; start[1] = 0; start[2] = 0;
    end[0] = display_->DisplayWidth() - 1;
    end[1] = display_->DisplayHeight() - 1;
    end[2] = FRAMEVOLUME_DEPTH - 1;
    s.packed = 0;
    display_->FillScaler(s, start, end);

    /*
     * Produces the coefficients for the input buffer using a forward 4x4 integer transform
     */
    for (size_t j = 0; j < CEIL(display_->DisplayHeight(), BLOCK_HEIGHT); j++) {
        for (size_t i = 0; i < CEIL(display_->DisplayWidth(), BLOCK_WIDTH); i++) {

            /* The coefficients are stored in this array in zig-zag order */
            vector<uint64_t> coefficients(BLOCK_WIDTH * BLOCK_HEIGHT, 0);
            lastNonZeroPlane = 0;

            /* Scratch blocks used for forwardIntegerTransform() */
            int redImgBlock[BLOCK_SIZE];
            int redCoeffsBlock[BLOCK_SIZE];
            int greenImgBlock[BLOCK_SIZE];
            int greenCoeffsBlock[BLOCK_SIZE];
            int blueImgBlock[BLOCK_SIZE];
            int blueCoeffsBlock[BLOCK_SIZE];

            /* Copy color channels for this block into their temporary single channel blocks. */
            size_t offset;
            for (size_t v = 0; v < BLOCK_HEIGHT; v++) {
                for (size_t u = 0; u < BLOCK_WIDTH; u++) {
                    offset = ((j * BLOCK_HEIGHT + v) * width + (i * BLOCK_WIDTH + u));
                    if (offset >= width * height) {
                        redImgBlock[v * BLOCK_WIDTH + u] = 0;
                        greenImgBlock[v * BLOCK_WIDTH + u] = 0;
                        blueImgBlock[v * BLOCK_WIDTH + u] = 0;
                    } else {
                        offset *= 3;
                        redImgBlock[v * BLOCK_WIDTH + u] = buffer[offset++];
                        greenImgBlock[v * BLOCK_WIDTH + u] = buffer[offset++];
                        blueImgBlock[v * BLOCK_WIDTH + u] = buffer[offset++];
                    }
                }
            }

            /* Perform Forward Integer Transform on each of the temporary single channel blocks. */
            forwardIntegerTransform(redCoeffsBlock, redImgBlock);
            forwardIntegerTransform(greenCoeffsBlock, greenImgBlock);
            forwardIntegerTransform(blueCoeffsBlock, blueImgBlock);

            /* Then update the coefficients */
            for (int v = 0; v < BLOCK_HEIGHT; v++) {
                for (int u = 0; u < BLOCK_WIDTH; u++) {
                    size_t p = zigZag_[v * BLOCK_WIDTH + u];
                    s.r = redCoeffsBlock[v * BLOCK_WIDTH + u];
                    s.g = greenCoeffsBlock[v * BLOCK_WIDTH + u];
                    s.b = blueCoeffsBlock[v * BLOCK_WIDTH + u];
                    coefficients[p] = s.packed;
                    if (s.packed != 0) lastNonZeroPlane = MAX(lastNonZeroPlane, p);
                }
            }

            /* Resize the coefficients vector */
            coefficients.resize(lastNonZeroPlane + 1);
            if (lastNonZeroPlane > largestNonZeroPlaneSeen) {
                largestNonZeroPlaneSeen = lastNonZeroPlane;
            }

            /* Send the NDDI command to update this macroblock's coefficients, one plane at a time. */
            start[0] = i * BLOCK_WIDTH;
            start[1] = j * BLOCK_HEIGHT;
            display_->FillScalerTileStack(coefficients, start, size);
        }
    }
}
