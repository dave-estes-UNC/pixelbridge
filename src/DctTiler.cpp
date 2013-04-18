#include <cmath>

#include "PixelBridgeFeatures.h"
#include "DctTiler.h"

/*
 * Frame Volume
 *
 * Dimensions are 8 x 8 x 193 and holds the 64 basis function for the red, green, and then blue channels and then
 * one extra planes for medium gray. The pre-rendered basis functions are arranged in the frame volume
 * using zig-zag ordering, with the colors interleaved.
 *
 * R(0,0) G(0,0) B(0,0) - R(1,0) G(1,0) B(1,0) - R(0,1) G(0,1) B(0,1) - R(0,2) G(0,2) B(0,2) -
 * R(1,1) G(1,1) B(1,1) - R(2,0) G(2,0) B(2,0) - R(3,0) G(3,0) B(3,0) - R(2,1) G(2,1) B(2,1) -
 * ... -
 * R(7,7) G(7,7) B(7,7)
 *
 *
 * Coefficient Plane
 *
 * The Coefficient Planes are then arranged to match, where each coefficient matrix in planes 0 - 195 correspond
 * to the same plane in the Frame Volume with the proper translation values. The last plane (196) picks the
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

/**
 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the FlatTiler should be destroyed and re-created.
 */
DctTiler::DctTiler (BaseNddiDisplay* display, bool quiet)
: display_(display),
  quiet_(quiet)
{
	size_t  x = 0, y = 0;
	bool    up = true;

	// Setup the zigZag_ table
	for (size_t i = 0; i < 64; i++) {
		zigZag_[i] = y * 8 + x;
		if (up) {
			x++;
			if (y > 0) {
				y--;
			} else {
				up = false;
			}
		} else {
			y++;
			if (x > 0) {
				x--;
			} else {
				up = true;
			}

		}
	}
}

/**
 * Initializes the Coefficient Planes for this tiler. The coefficient matrices
 * for each plane will pick a plane from the frame volume in zig-zag order.
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

	// Break the display into macroblocks and initialize each 8x8x193 cube of coefficients to pick out the proper
	for (int k = 0; k < (BLOCKS_WIDE * BLOCKS_TALL * 3 + 1); k++) {
		for (int j = 0; j < display_->DisplayHeight() / MACROBLOCK_HEIGHT; j++) {
			for (int i = 0; display_->DisplayWidth() / i < MACROBLOCK_WIDTH; i++) {
				coeffs[2][0] = -i * MACROBLOCK_WIDTH;
				coeffs[2][1] = -j * MACROBLOCK_HEIGHT;
				coeffs[2][2] = k;
				start[0] = i * MACROBLOCK_WIDTH; start[1] = j * MACROBLOCK_HEIGHT; start[2] = k;
				end[0] = (i + 1) * MACROBLOCK_WIDTH - 1; end[1] = (j + 1) * MACROBLOCK_HEIGHT - 1; end[2] = k;
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

	// Fill the scalers for the medium gray plane to full on
    start[2] = end[2] = (BLOCKS_WIDE * BLOCKS_TALL * 3 + 1) - 1;
    display_->FillScaler(NUM_COEFFICIENT_PLANES, start, end);
}

/**
 * Initializes the Frame Volume for this tiler by pre-rendering each
 * of the 64 basis functions into 8x8 planes in the Frame Volume.
 */
void DctTiler::InitializeFrameVolume() {

	Pixel  *pixels;
	size_t  p;
	size_t  pixels_size = sizeof(Pixel)                           // pixel size
						 * BLOCK_WIDTH * BLOCK_HEIGHT            // size of each x,y plane
						 * (BLOCKS_WIDE * BLOCKS_TALL * 3 + 1);  // number of basis functions for each color channel
	                                                             //   plus one medium gray plane


	pixels = (Pixel *)malloc(pixels_size);
	memset(pixels, pixels_size, 0x00);

	// Pre-render each basis function
    for (int j = 0; j < BLOCKS_TALL; j++) {
        for (int i = 0; i < BLOCKS_WIDE; i++) {
        	p = 0;
            for (int y = 0; y < BLOCK_HEIGHT; y++) {
                for (int x = 0; x < BLOCK_WIDTH; x++) {
                    double m = 0.0f;
                    bool neg = false;

                    for (int v = 0; v < BLOCK_HEIGHT; v++) {
                        for (int u = 0; u < BLOCK_WIDTH; u++) {
                            double s = 1.0f;
                            s *= (u == 0) ? sqrt((double)0.125) : sqrt((double)0.25);   // alpha(u)
                            s *= (v == 0) ? sqrt((double)0.125) : sqrt((double)0.25);   // alpha(v)
                            s *= ((u == i) && (v == j)) ? double(MAX_DCT_COEFF) : 0.0;  // DCT coefficient (maximum on or off)
                            s *= cos(PI / 8.0 * ((double)x + 0.5) * (double)u);         // cos with x, u
                            s *= cos(PI / 8.0 * ((double)y + 0.5) * (double)v);         // cos with y, v
                            m += s;
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
                    pixels[zigZag_[p]].r                  = c;
                    pixels[zigZag_[p] + BLOCK_SIZE].g     = c;
                    pixels[zigZag_[p] + BLOCK_SIZE * 2].b = c;

                    // Finally set the alpha channel and move to the next pixel
                    pixels[zigZag_[p] + BLOCK_SIZE * 3].a = (unsigned char)255;
                    p++;
                }
            }
        }
    }

    // Update the frame volume in bulk
    //TODO(CDE): Do it or it won't work.

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

}
