#include <cmath>
#include <iostream>

#include "PixelBridgeFeatures.h"
#include "DctTiler.h"

#define PI 3.14159265

/*
 * Frame Volume
 *
 * Dimensions are 8 x 8 x 193 and holds the 64 basis function for the red, green, and then blue channels and then
 * one extra plane for medium gray. The pre-rendered basis functions are arranged in the frame volume
 * using zig-zag ordering, with the colors interleaved.
 *
 *      0      1      2      3      4      5      6      7      8      9     10     11
 * R(0,0) G(0,0) B(0,0) R(1,0) G(1,0) B(1,0) R(0,1) G(0,1) B(0,1) R(0,2) G(0,2) B(0,2)
 *     12     13     14     15     16     17     18     19     20     21     22     23
 * R(1,1) G(1,1) B(1,1) R(2,0) G(2,0) B(2,0) R(3,0) G(3,0) B(3,0) R(2,1) G(2,1) B(2,1)
 *   ...
 *    189    190    191
 * R(7,7) G(7,7) B(7,7)
 *
 *
 * Coefficient Plane
 *
 * The Coefficient Planes are then arranged to match, where each coefficient matrix in planes 0 - 191 correspond
 * to the same plane in the Frame Volume with the proper translation values. The last plane (192) picks the
 * extra medium medium gray plane from the frame volume.
 *
 * The scalers for each 8x8 block of coefficient matrices are set per macro block using the DCT coefficients
 * calculated for that macro block. The coefficients are of course in zig-zag order, and so an efficient fill
 * mechanism can be use to pass the non-zero coefficients down the "stack" of 8x8 coefficient planes regions.
 * The scaler for the final medium gray plane is set to full (256) as a way of adding back the 128 to each
 * channel as the final step of this NDDI-based IDCT.
 */

/**
 * Using a #define because it's used for floats and ints.
 */
#define CLAMP(x, min, max)   (x < min ? min : x >= max ? max : x)
#define MAX(x, y)            (x > y ? x : y)
#define MIN(x, y)            (x < y ? x : y)
#define CEIL(x, y)           (1 + ((x - 1) / y))

/**
 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the FlatTiler should be destroyed and re-created.
 */
DctTiler::DctTiler (size_t display_width, size_t display_height,
                    size_t quality, bool quiet)
: quiet_(quiet)
{
    // 3 dimensional matching the Macroblock Width x Height x 64+3+1
    vector<unsigned int> fvDimensions;
    fvDimensions.push_back(BLOCK_WIDTH);
    fvDimensions.push_back(BLOCK_HEIGHT);
    fvDimensions.push_back(FRAMEVOLUME_DEPTH);
    
#ifndef NO_CL
#error CL not supported for DCT Tiler yet.
#else
    display_ = new GlNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 256,                           // Number of coefficient planes
                                 3); 						    // input vector size (x, y, 1)
#endif

    initZigZag();
    initQuantizationMatrix(quality);
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
GlNddiDisplay* DctTiler::GetDisplay() {
    return display_;
}

/*
 * A particular 8x8 macroblock has pixels index by x,y with the origin at the top left:
 *
 *   0,0 1,0 2,0 3,0 4,0 5,0 6,0 7,0
 *   0,1 1,1 2,1 3,1 4,1 5,1 6,1 7,1
 *   0,2 1,2 2,2 3,2 4,2 5,2 6,2 7,2
 *   0,3 1,3 2,3 3,3 4,3 5,3 6,3 7,3
 *   0,4 1,4 2,4 3,4 4,4 5,4 6,4 7,4
 *   0,5 1,5 2,5 3,5 4,5 5,5 6,5 7,5
 *   0,6 1,6 2,6 3,6 4,6 5,6 6,6 7,6
 *   0,7 1,7 2,7 3,7 4,7 5,7 6,7 7,7
 *
 * Zig-Zag order is therefore:
 *
 *                 0,0
 *               1,0 0,1
 *             0,2 1,1 2,0
 *           3,0 2,1 1,2 0,3
 *         0,4 1,3 2,2 3,1 4,0
 *       5,0 4,1 3,2 2,3 1,4 0,5
 *     0,6 1,5 2,4 3,3 4,2 5,1 6,0
 *   7,0 6,1 5,2 4,3 3,4 2,5 1,6 0,7
 *     1,7 2,6 3,5 4,4 5,3 6,2 7,1
 *       7,2 6,3 5,4 4,5 3,6 2,7
 *         3,7 4,6 5,5 6,4 7,3
 *           7,4 6,5 5,6 4,7
 *             5,7 6,6 7,5
 *               7,6 6,7
 *                 7,7
 *
 * This zig-zag ordering is then stored in the zigZag_ array. Each (x,y) pair
 * is used as a unique index into the array by multiplying the pair in a row-major
 * manner (y*8+x). Then the number of the corresponding frame volume frame is stored.
 * Specifically, each plane in the frame volume holds a rendering of a basis function
 * and each of the color channels has its own rendering. So the value stored at the
 * particular entry in the zigZag_ array is actually a logical value that must be
 * multiplied by 3 and then the proper frames for that color are chosen independently.
 */
void DctTiler::initZigZag() {
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

/*
 * Uses simple algorithm from M. Nelson, "The Data Compression Book," San Mateo, CA, M&T Books, 1992.
 */
void DctTiler::initQuantizationMatrix(size_t quality) {
	for (size_t v = 0; v < BLOCK_HEIGHT; v++) {
		for (size_t u = 0; u < BLOCK_WIDTH; u++) {
			quantizationMatrix_[v * BLOCK_WIDTH + u] = 1 + (1 + u + v) * quality;
		}
	}
}

/**
 * Initializes the Coefficient Planes for this tiler. The coefficient matrices
 * for each plane will pick a plane from the frame volume.
 */
void DctTiler::InitializeCoefficientPlanes() {

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

	// Break the display into macroblocks and initialize each 8x8x193 cube of coefficients to pick out the proper block from the frame volume
	for (int j = 0; j < CEIL(display_->DisplayHeight(), BLOCK_HEIGHT); j++) {
		for (int i = 0; i < CEIL(display_->DisplayWidth(), BLOCK_WIDTH); i++) {
			coeffs[2][0] = -i * BLOCK_WIDTH;
			coeffs[2][1] = -j * BLOCK_HEIGHT;
			start[0] = i * BLOCK_WIDTH; start[1] = j * BLOCK_HEIGHT; start[2] = 0;
			end[0] = (i + 1) * BLOCK_WIDTH - 1; end[1] = (j + 1) * BLOCK_HEIGHT - 1; end[2] = FRAMEVOLUME_DEPTH - 1;
			if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
			if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }
			display_->FillCoefficientMatrix(coeffs, start, end);
		}
	}
	// Finish up by setting the proper k for every plane
	start[0] = 0; start[1] = 0;
	end[0] = display_->DisplayWidth() - 1; end[1] = display_->DisplayHeight() - 1;
	for (int k = 0; k < FRAMEVOLUME_DEPTH; k++) {
		start[2] = k; end[2] = k;
		display_->FillCoefficient(k, 2, 2, start, end);
	}

	// Fill each scaler in every plane with 0
	start[0] = 0; start[1] = 0; start[2] = 0;
	end[0] = display_->DisplayWidth() - 1;
	end[1] = display_->DisplayHeight() - 1;
    end[2] = display_->NumCoefficientPlanes() - 1;
    display_->FillScaler(0, start, end);

	// Fill the scalers for the medium gray plane to full on
    start[2] = end[2] = FRAMEVOLUME_DEPTH - 1;
    display_->FillScaler(display_->NumCoefficientPlanes(), start, end);
}

/**
 * Initializes the Frame Volume for this tiler by pre-rendering each
 * of the 16 basis functions into 4x4 planes in the Frame Volume. They're
 * rendered for each color channel and stored in those groups of three in
 * zig-zag order.
 */
void DctTiler::InitializeFrameVolume() {

	Pixel  *pixels;
	size_t  pixels_size = sizeof(Pixel)            // pixel size
						 * BLOCK_SIZE              // size of each x,y plane
						 * FRAMEVOLUME_DEPTH;      // number of basis functions for each color channel
	                                               //   plus one medium gray plane


	pixels = (Pixel *)malloc(pixels_size);
	memset(pixels, 0x00, pixels_size);

	// Pre-render each basis function
#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (int j = 0; j < BASIS_BLOCKS_TALL; j++) {
        for (int i = 0; i < BASIS_BLOCKS_WIDE; i++) {
            
        	size_t p = zigZag_[j * BASIS_BLOCKS_WIDE + i] * BLOCK_SIZE * 3;
            
            for (int y = 0; y < BLOCK_HEIGHT; y++) {
                for (int x = 0; x < BLOCK_WIDTH; x++) {
                    double m = 0.0f;
                    bool neg = false;

                    for (int v = 0; v < BLOCK_HEIGHT; v++) {
                        for (int u = 0; u < BLOCK_WIDTH; u++) {
                            double p = 1.0f;
                            p *= (u == 0) ? sqrt((double)0.125) : sqrt((double)0.25);   // alpha(u)
                            p *= (v == 0) ? sqrt((double)0.125) : sqrt((double)0.25);   // alpha(v)
                            p *= ((u == i) && (v == j)) ? double(MAX_DCT_COEFF) : 0.0;  // DCT coefficient (maximum on or off)
                            p *= cos(PI / 8.0 * ((double)x + 0.5) * (double)u);         // cos with x, u
                            p *= cos(PI / 8.0 * ((double)y + 0.5) * (double)v);         // cos with y, v
                            m += p;
                        }
                    }

                    // Set the neg flag if the magnitude is negative and then remove the sign from the magnitude
                    if (m < 0.0f) {
                        neg = true;
                        m *= -1.0;
                    }

                    // Clamp to 127 and cast to byte and convert to two's complement if necessary
                    unsigned int c = (unsigned char)CLAMP(m, 0.0, 127.0);
                    if (neg) {
                        c = ~c + 1;
                    }

                    // Set the color channels with the magnitude clamped to 127
                    pixels[p + BLOCK_SIZE * 0].r = c;
                    pixels[p + BLOCK_SIZE * 0].a = 0xff;
                    pixels[p + BLOCK_SIZE * 1].g = c;
                    pixels[p + BLOCK_SIZE * 1].a = 0xff;
                    pixels[p + BLOCK_SIZE * 2].b = c;
                    pixels[p + BLOCK_SIZE * 2].a = 0xff;

                    // Move to the next pixel
                    p++;
                }
            }
        }
    }

    // Then render the gray block
    size_t p = BLOCK_SIZE * BASIS_BLOCKS_WIDE * BASIS_BLOCKS_TALL * 3;
    for (int y = 0; y < BLOCK_HEIGHT; y++) {
        for (int x = 0; x < BLOCK_WIDTH; x++) {
        	unsigned int c = 0x7f;
        	pixels[p].r = c;
        	pixels[p].g = c;
        	pixels[p].b = c;
        	pixels[p].a = 0xff;
        	p++;
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

/**
 * Update the display by calculating the DCT coefficients for each macroblock
 * and updating the coefficient plane scalers.
 *
 * @param buffer Pointer to an RGB buffer
 * @param width The width of the RGB buffer
 * @param height The height of the RGB buffer
 */
void DctTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height)
{
	vector<unsigned int> start(3, 0), end(3, 0);
	vector<unsigned int> size(2, 0);
	size_t lastNonZeroPlane;
	static size_t largestNonZeroPlaneSeen = 0;

	size[0] = BLOCK_WIDTH;
	size[1] = BLOCK_HEIGHT;

	/* Clear all scalers up to but not including the medium gray plane. */
	// TODO(CDE): Don't do this in the future. Just write enough planes to overwrite the last known non-zero plane.
	start[0] = 0; start[1] = 0; start[2] = 0;
	end[0] = display_->DisplayWidth() - 1;
	end[1] = display_->DisplayHeight() - 1;
    end[2] = FRAMEVOLUME_DEPTH - 2;
    display_->FillScaler(0, start, end);

	/*
	 * Produces the de-quantized coefficients for the input buffer using the following steps:
	 *
	 * 1. Shift by subtracting 128
	 * 2. Take the 2D DCT
	 * 3. Quantize
	 * 4. De-quantize
	 */
	for (size_t j = 0; j < CEIL(display_->DisplayHeight(), BLOCK_HEIGHT); j++) {
		for (size_t i = 0; i < CEIL(display_->DisplayWidth(), BLOCK_WIDTH); i++) {

			/* The coefficients are stored in this array in zig-zag order */
			vector<int> coefficients(BLOCK_WIDTH * BLOCK_HEIGHT * 3, 0);
			lastNonZeroPlane = 0;

#ifndef NO_OMP
#pragma omp parallel for ordered
#endif
			for (size_t v = 0; v < BLOCK_HEIGHT; v++)
#ifndef NO_OMP
#pragma omp ordered
#endif
			{
				for (size_t u = 0; u < BLOCK_WIDTH; u++) {

					double c_r = 0.0f, c_g = 0.0f, c_b = 0.0f;

					/* Calculate G for each u, v using g(x,y) shifted (1. and 2.) */
					for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
						for (size_t x = 0; x < BLOCK_WIDTH; x++) {

							double p = 1.0f;

							p *= (u == 0) ? sqrt((double)0.125) : sqrt((double)0.25);                                            // alpha(u)
							p *= (v == 0) ? sqrt((double)0.125) : sqrt((double)0.25);                                            // alpha(v)
							p *= cos(PI / 8.0 * ((double)x + 0.5) * (double)u);                                                  // cos with x, u
							p *= cos(PI / 8.0 * ((double)y + 0.5) * (double)v);                                                  // cos with y, v
							/* Fetch each channel, multiply by product p and then shift */
							if ( ((i * BLOCK_WIDTH + x) < display_->DisplayWidth())
									&& ((j * BLOCK_HEIGHT + y) < display_->DisplayHeight()) ) {
								c_r += p * ((double)buffer[((j * BLOCK_HEIGHT + y) * width + (i * BLOCK_WIDTH + x)) * 3 + 0] - 128.0); // red: g(x,y) - 128
								c_g += p * ((double)buffer[((j * BLOCK_HEIGHT + y) * width + (i * BLOCK_WIDTH + x)) * 3 + 1] - 128.0); // green: g(x,y) - 128
								c_b += p * ((double)buffer[((j * BLOCK_HEIGHT + y) * width + (i * BLOCK_WIDTH + x)) * 3 + 2] - 128.0); // blue: g(x,y) - 128
							}
						}
					}

					int g_r, g_g, g_b;

					/* Quantize G(u,v) (3.) */
					g_r = (int)int(c_r / double(quantizationMatrix_[v * BLOCK_WIDTH + u]) + 0.5);
					g_g = (int)int(c_g / double(quantizationMatrix_[v * BLOCK_WIDTH + u]) + 0.5);
					g_b = (int)int(c_b / double(quantizationMatrix_[v * BLOCK_WIDTH + u]) + 0.5);

					/* De-quantized G(u,v) (4.) */
					g_r *= (int)quantizationMatrix_[v * BLOCK_WIDTH + u];
					g_g *= (int)quantizationMatrix_[v * BLOCK_WIDTH + u];
					g_b *= (int)quantizationMatrix_[v * BLOCK_WIDTH + u];

					/* Set the coefficient in zig-zag order. */
					size_t p = zigZag_[v * BLOCK_WIDTH + u] * 3;
					coefficients[p + 0] = g_r;
					coefficients[p + 1] = g_g;
					coefficients[p + 2] = g_b;
					if (g_r != 0) lastNonZeroPlane = MAX(lastNonZeroPlane, p + 0);
					if (g_g != 0) lastNonZeroPlane = MAX(lastNonZeroPlane, p + 1);
					if (g_b != 0) lastNonZeroPlane = MAX(lastNonZeroPlane, p + 2);
				}
			}
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
