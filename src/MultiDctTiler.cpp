#include <cmath>
#include <iostream>

#include "PixelBridgeFeatures.h"
#include "Configuration.h"
#include "MultiDctTiler.h"

#define PI    3.14159265
#define PI_8  0.392699081

// For each scale, determines the optimal amount of planes to zero out.
//#define OPTIMAL_ZEROING

// This will turn on the simple truncation used before that simply takes the top X coefficients.
// The new scheme will instead take the most significant coefficients in a square pattern.
//#define SIMPLE_TRUNCATION

// This uses the scale to determine how often to send the coefficients at scales greater than 1x.
// e.g. The 16x scale coefficients are only sent every 16th frame. The 4x every 4th frame.
//#define SKIP_FRAMES


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

/*
 * Using a #define because it's used for floats and ints.
 */
#define CLAMP(x, min, max)   (x < min ? min : x >= max ? max : x)
#define MAX(x, y)            (x > y ? x : y)
#define MIN(x, y)            (x < y ? x : y)
#define CEIL(x, y)           (1 + ((x - 1) / y))
#define SQRT_125             0.353553391
#define SQRT_250             0.5

MultiDctTiler::MultiDctTiler(size_t display_width, size_t display_height, size_t quality)
    : ScaledDctTiler(display_width, display_height, quality) {}

/**
 * Builds the coefficients for the macroblock specified by (i, j) using the source buffer.
 *
 * @param i The column component of the macroblock
 * @param j The row component of the macroblock
 * @param buffer The source buffer. Note this may or may not be a downscaled image.
 * @param width The width of the source buffer
 * @param height The height of the source buffer
 * @return The vector (in zig-zag order) of 4-channel scalers.
 */
vector<uint64_t> MultiDctTiler::BuildCoefficients(size_t i, size_t j, int16_t* buffer, size_t width, size_t height, bool shift) {

    /*
     * Produces the de-quantized coefficients for the input buffer using the following steps:
     *
     * 1. Shift by subtracting 128
     * 2. Take the 2D DCT
     * 3. Quantize
     * 4. De-quantize
     */

    Scaler s;

    /* The coefficients are stored in this array in zig-zag order */
    vector<uint64_t> coefficients(BLOCK_SIZE, 0);

    for (size_t v = 0; v < BLOCK_HEIGHT; v++) {
#ifndef NO_OMP
#pragma omp parallel for ordered
#endif
        for (size_t u = 0; u < BLOCK_WIDTH; u++) {

            double c_r = 0.0, c_g = 0.0, c_b = 0.0;

            /* Calculate G for each u, v using g(x,y) shifted (1. and 2.) */
            for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
                size_t bufPos = ((j * BLOCK_HEIGHT + y) * width + i * BLOCK_WIDTH + 0) * 3;
                for (size_t x = 0; x < BLOCK_WIDTH; x++) {

                    double p = 1.0;

                    p *= (u == 0) ? SQRT_125 : SQRT_250;                                 // alpha(u)
                    p *= (v == 0) ? SQRT_125 : SQRT_250;                                 // alpha(v)
                    p *= cos(PI_8 * ((double)x + 0.5) * (double)u);                      // cos with x, u
                    p *= cos(PI_8 * ((double)y + 0.5) * (double)v);                      // cos with y, v
                    /* Fetch each channel, multiply by product p and then shift */
                    if ( ((i * BLOCK_WIDTH + x) < width)
                            && ((j * BLOCK_HEIGHT + y) < height) ) {
                        if (shift) {
                            c_r += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // red: g(x,y) - 128
                            c_g += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // green: g(x,y) - 128
                            c_b += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // blue: g(x,y) - 128
                        } else {
                            c_r += p * (double)buffer[bufPos]; bufPos++;           // red: g(x,y)
                            c_g += p * (double)buffer[bufPos]; bufPos++;           // green: g(x,y)
                            c_b += p * (double)buffer[bufPos]; bufPos++;           // blue: g(x,y)
                        }
                    } else {
                        bufPos += 3;
                    }
                }
            }

            int g_r, g_g, g_b;
            size_t matPos = v * BLOCK_WIDTH + u;

            /* Quantize G(u,v) (3.) */
            g_r = (int)int(c_r / double(quantizationMatrix_[matPos]) + 0.5); // 0.5 is for rounding
            g_g = (int)int(c_g / double(quantizationMatrix_[matPos]) + 0.5);
            g_b = (int)int(c_b / double(quantizationMatrix_[matPos]) + 0.5);

            /* De-quantized G(u,v) (4.) */
            g_r *= (int)quantizationMatrix_[matPos];
            g_g *= (int)quantizationMatrix_[matPos];
            g_b *= (int)quantizationMatrix_[matPos];

            /* Set the coefficient in zig-zag order. */
            size_t p = zigZag_[matPos];

            /* Skip the last block, because we've used it for the medium gray block */
            if (p == BLOCK_SIZE - 1) continue;

            /* Build the scaler from the three coefficients */
            s.packed = 0;
            s.r = g_r;
            s.g = g_g;
            s.b = g_b;
            coefficients[p] = s.packed;
        }
    }

    return coefficients;
}


/**
 * Takes the coefficients computed for a macroblock and fills them to a region of the display.
 * If a factor greater than one is provided in the config, then this implies that the coefficients were built
 * for a particular macroblock size, but that they are being filled to a larger region of the
 * the display. e.g. Given an 8x8 macroblock, a macroblock location (i, j) of (2, 1) and a factor of
 * 4: the region of the display that will be updated is from (64, 32) to (95, 63) inclusive.
 *
 * @param coefficients The vector (in zig-zag order) of coefficients for the macroblock
 * @param i The column component of the macroblock
 * @param j The row component of the macroblock
 * @param c Index into globalConfiguration.dctScales for information about the factor by
 *          which we're scaling and which planes to fill.
 */
void MultiDctTiler::FillCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, size_t first) {

    vector<unsigned int> start(3, 0);
    vector<unsigned int> size(2, 0);

    scale_config_t config = globalConfiguration.dctScales[c];

    start[0] = i * BLOCK_WIDTH * config.scale_multiplier;
    start[1] = j * BLOCK_HEIGHT * config.scale_multiplier;
    start[2] = first;

    size[0] = BLOCK_WIDTH * config.scale_multiplier;
    size[1] = BLOCK_HEIGHT * config.scale_multiplier;

    /* If any any coefficients have changed, send the NDDI command to update them */
    if (start[2] < display_->NumCoefficientPlanes()) {
        display_->FillScalerTileStack(coefficients, start, size);
    }
}

/**
 * Given the set of coefficients, this routine will simulate a rendering using the specified set of planes.
 * The result is rendered back to the buffer provided into the macroblock specified by (i, j).
 *
 * @param coefficients The vector (in zig-zag order) of coefficients for the macroblock
 * @param i The column component of the macroblock
 * @param j The row component of the macroblock
 * @param c The current scale from the globalConfiguration.dctScales global.
 * @param buffer The destination buffer
 * @param width The width of the destination buffer
 * @param width The height of the destination buffer
 */
void MultiDctTiler::PrerenderCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, int16_t* buffer, size_t width, size_t height, bool shift) {

    scale_config_t config = globalConfiguration.dctScales[c];

#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
        for (size_t x = 0; x < BLOCK_WIDTH; x++) {

            if (i * BLOCK_WIDTH + x >= width || j * BLOCK_HEIGHT + y >= height)
                continue;

            int rAccumulator = 0, gAccumulator = 0, bAccumulator = 0;

#ifdef SIMPLE_TRUNCATION
            for (size_t p = 0; p < coefficients.size(); p++) {
                Scaler s;
                s.packed = coefficients[p];

                size_t bfo = p * BLOCK_WIDTH * BLOCK_HEIGHT + y * BLOCK_WIDTH + x;
                rAccumulator += (int8_t)(basisFunctions_[bfo].r) * s.r;
                gAccumulator += (int8_t)(basisFunctions_[bfo].g) * s.g;
                bAccumulator += (int8_t)(basisFunctions_[bfo].b) * s.b;
            }
#else
            // For an 8x8 region of coefficients, just use simple truncation
            if (config.edge_length == 8) {
                for (size_t p = 0; p < coefficients.size(); p++) {
                    Scaler s;
                    s.packed = coefficients[p];

                    size_t bfo = p * BLOCK_WIDTH * BLOCK_HEIGHT + y * BLOCK_WIDTH + x;
                    rAccumulator += (int8_t)(basisFunctions_[bfo].r) * s.r;
                    gAccumulator += (int8_t)(basisFunctions_[bfo].g) * s.g;
                    bAccumulator += (int8_t)(basisFunctions_[bfo].b) * s.b;
                }
            // Else consider only the most significant coefficients within with the square of edge_length
            } else {
                for (size_t p = 0, yy = 0; yy < config.edge_length && p < coefficients.size(); yy++) {
                    for (size_t xx = 0; xx < config.edge_length && p < coefficients.size(); xx++) {
                        Scaler s;
                        s.packed = coefficients[p];

                        size_t pp = zigZag_[yy * BLOCK_WIDTH + xx];

                        size_t bfo = pp * BLOCK_WIDTH * BLOCK_HEIGHT + y * BLOCK_WIDTH + x;
                        rAccumulator += (int8_t)(basisFunctions_[bfo].r) * s.r;
                        gAccumulator += (int8_t)(basisFunctions_[bfo].g) * s.g;
                        bAccumulator += (int8_t)(basisFunctions_[bfo].b) * s.b;

                        p++;
                    }
                }
            }
#endif

            size_t offset = ((j * BLOCK_HEIGHT + y) * width + i * BLOCK_WIDTH + x) * 3;

            if (shift) {
                buffer[offset + 0] = rAccumulator / display_->GetFullScaler() + 128;
                buffer[offset + 1] = gAccumulator / display_->GetFullScaler() + 128;
                buffer[offset + 2] = bAccumulator / display_->GetFullScaler() + 128;
            } else {
                buffer[offset + 0] = rAccumulator / display_->GetFullScaler();
                buffer[offset + 1] = gAccumulator / display_->GetFullScaler();
                buffer[offset + 2] = bAccumulator / display_->GetFullScaler();
            }
        }
    }
}

/**
 * Returns the Display created and initialized by the tiler.
 */
GlNddiDisplay* MultiDctTiler::GetDisplay() {
    return display_;
}

/**
 * Update the display by calculating the DCT coefficients for each macroblock
 * and updating the coefficient plane scalers.
 *
 * 0. Convert image to signed pixels if we don't already have
 * For each scale level
 *   1. For each NxN macroblock in the image
 *     a. Build coefficients - BuildCoefficients()
 *   2. For this current scale, snap coefficients to zero if configured and then figure out the optimal planes
 *      to zero out in bulk.
 *     a. Snap to zero if configured to do so - SnapCoefficientsToZero()
 *     b. Zero out optimal planes if we've got coefficients cached already - ZeroPlanes()
 *   3. For each NxN macroblock in the image
 *     a. Trim a copy of the coefficients - TrimCoefficients()
 *     b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
 *     c. Perform simulated blending of basis functions and store back to the image - PrerenderCoefficients()
 *   4. Subtract the results from the original image - AdjustFrame()
 *   5. Cleanup
 *
 * @param buffer Pointer to an RGB buffer
 * @param width The width of the RGB buffer
 * @param height The height of the RGB buffer
 */
void MultiDctTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height) {
    assert(width >= display_->DisplayWidth());
    assert(height >= display_->DisplayHeight());

    /* Set the initial delta and planes for snapping and trimming below. */
    size_t delta = 0, planes = 0;

#ifdef SKIP_FRAMES
    static size_t frame = 0;
#endif

    /*
     * 0. Convert image to signed pixels if we don't already have
     */
    int16_t* signedBuf = ConvertToSignedPixels(buffer, width, height);

    // For each scale level
    for (int c = 0; c < globalConfiguration.dctScales.size(); c++) {

        scale_config_t config = globalConfiguration.dctScales[c];

        /*
         * 1. Downsample the image - DownSample()
         */
        int16_t* downBuf;
        size_t downW = CEIL(width, config.scale_multiplier);
        size_t downH = CEIL(height, config.scale_multiplier);
#if 0 // CDE: Fix
        if (config.scale_multiplier == 1)
            downBuf = signedBuf;
        else
            downBuf = DownSample(config.scale_multiplier, signedBuf, width, height);
#else
        downBuf = signedBuf;
#endif

        /*
         * 2. For each macroblock in downsampled image
         *    a. Build coefficients - BuildCoefficients()
         */
        int16_t* rendBuf = (int16_t*)calloc(downW * downH * 3, sizeof(int16_t));

        size_t tilesWide = CEIL(downW, BLOCK_WIDTH);
        size_t tilesHigh = CEIL(downH, BLOCK_HEIGHT);

        vector< vector< vector<uint64_t> > > coefficientsForCurrentScale;

        for (size_t i = 0; i < tilesWide; i++) {
            coefficientsForCurrentScale.push_back(vector< vector<uint64_t> >());
            for (size_t j = 0; j < tilesHigh; j++) {
                // a. Build coefficients - BuildCoefficients()
                vector<uint64_t> coefficients = BuildCoefficients(i, j, downBuf, downW, downH, c == 0);
                SelectCoefficientsForScale(coefficients, c);
                coefficientsForCurrentScale[i].push_back(coefficients);
            }
        }

        /*
         * 3. For this current scale, snap coefficients to zero if configured and then figure out the optimal planes
         *    to zero out in bulk.
         *   a. Snap to zero if configured to do so - SnapCoefficientsToZero()
         *   b. Zero out optimal planes if we've got coefficients cached already - ZeroPlanes()
         */
        // a. Snap to zero if configured to do so - SnapCoefficientsToZero()
        if (globalConfiguration.dctSnap) {
            CalculateSnapCoefficientsToZero(coefficientsForCurrentScale, c, delta, planes);
            SnapCoefficientsToZero(coefficientsForCurrentScale, c, delta, planes);
        }
#ifdef OPTIMAL_ZEROING
        // b. Zero out optimal planes if we've got coefficients cached already - ZeroPlanes()
        if (cachedCoefficients_.size() > c)
            ZeroPlanes(coefficientsForCurrentScale, c);
#endif

        /*
         * 4. For each macroblock in downsampled image
         *   a. Trim a copy of the coefficients - TrimCoefficients()
         *   b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
         *   c. Perform simulated blending of basis functions and store back to downsampled image - PrerenderCoefficients()
         */
        CalculateTrimCoefficients(coefficientsForCurrentScale, c, delta, planes);
        for (size_t i = 0; i < tilesWide; i++) {
            for (size_t j = 0; j < tilesHigh; j++) {

#ifdef SKIP_FRAMES
                if (config.scale_multiplier == 1 || frame % config.scale_multiplier == 0) {
                    // a. Trim a copy of the coefficients - TrimCoefficients()
                    vector<uint64_t> coefficients = coefficientsForCurrentScale[i][j];
                    size_t first_plane_idx = TrimCoefficients(coefficients, i, j, c);

                    // b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
                    FillCoefficients(coefficients, i, j, c, first_plane_idx);
                } else {
                    // Otherwise just re-use the cache if we're skipping this frame.
                    coefficientsForCurrentScale[i][j] = cachedCoefficients_[c][i][j];
                }
#else
                // a. Trim a copy of the coefficients - TrimCoefficients()
                vector<uint64_t> coefficients = coefficientsForCurrentScale[i][j];
                size_t first_plane_idx = TrimCoefficients(coefficients, i, j, c, delta, planes);

                // b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
                FillCoefficients(coefficients, i, j, c, first_plane_idx);
#endif

                // c. Perform simulated blending of basis functions - PrerenderCoefficients()
                // Again, only shift on the first plane
                PrerenderCoefficients(coefficientsForCurrentScale[i][j], i, j, c, rendBuf, downW, downH, c == 0);
            }
        }

        /*
         * 5. Upsample the image - UpSample()
         */
        int16_t* upBuf;
#if 0 //CDE: Fix
        if (config.scale_multiplier == 1)
            upBuf = rendBuf;
        else
            upBuf = UpSample(config.scale_multiplier, rendBuf, width, height);
#else
        upBuf = rendBuf;
#endif

        /*
         * 6. Subtract the results from the original image - AdjustFrame()
         */
        AdjustFrame(signedBuf, upBuf, width, height);

        /*
         * 7. Cleanup
         */
        free(rendBuf);
        if (config.scale_multiplier > 1) {
            free(downBuf);
            free(upBuf);
        }
    }

#ifdef SKIP_FRAMES
    frame++;
#endif

    // Finally clean the signedBuf that we've been using throughout
    free(signedBuf);
}
