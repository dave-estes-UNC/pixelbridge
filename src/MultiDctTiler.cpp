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
#define SQRT2                1.414213562

MultiDctTiler::MultiDctTiler(size_t display_width, size_t display_height, size_t quality)
    : ScaledDctTiler(display_width, display_height, quality) {

    // The frame volume is still 64 planes deep, but the width and height is based on the first
    // configuration which is assumed to be the largest.
    fvWidth_ = globalConfiguration.dctScales[0].scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
    fvHeight_ = globalConfiguration.dctScales[0].scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;

}

/**
 * Initializes each of the zig zag patterns for the different scales.
 */
void MultiDctTiler::initZigZag() {

    // For each scale
    for (size_t c = 0; c < globalConfiguration.dctScales.size(); c++) {

        size_t  block_width = globalConfiguration.dctScales[c].scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
        size_t  block_height = globalConfiguration.dctScales[c].scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;
        size_t  block_size = block_width * block_height;

        // Add the new zigZag and set its size
        zigZag_.push_back(vector<int>());
        zigZag_[c].resize(block_size);

        // Then initialize it
        size_t  x = 0, y = 0;
        bool    up = true;

        // Setup the zigZag_ table
        for (size_t i = 0; i < block_size; i++) {
            zigZag_[c][y * block_width + x] = i;
            if (up) {
                if (x < (block_width - 1)) {
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
                if (y < (block_height - 1)) {
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
}


/*
 * Uses simple algorithm from M. Nelson, "The Data Compression Book," San Mateo, CA, M&T Books, 1992.
 */
void MultiDctTiler::initQuantizationMatrix(size_t quality) {

    // For each scale
    for (size_t c = 0; c < globalConfiguration.dctScales.size(); c++) {

        size_t  block_width = globalConfiguration.dctScales[c].scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
        size_t  block_height = globalConfiguration.dctScales[c].scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;
        size_t  block_size = block_width * block_height;

        // Add the new quantizationMatrix and set its size
        quantizationMatrix_.push_back(vector<uint8_t>());
        quantizationMatrix_[c].resize(block_size);

        for (size_t v = 0; v < block_height; v++) {
            for (size_t u = 0; u < block_width; u++) {
                quantizationMatrix_[c][v * block_width + u] = 1 + (1 + u + v) * quality;
            }
        }
    }
}

/**
 * Initializes the Coefficient Planes for this tiler. The coefficient matrices
 * for each plane will pick a plane from the frame volume.
 */
void MultiDctTiler::InitializeCoefficientPlanes() {

    // Perform the basic initialization first and then do the specific
    // multiscale initialization
    DctTiler::InitializeCoefficientPlanes();

    // Setup start and end points to (0,0,0) initially
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(0); end.push_back(0); end.push_back(0);

    //
    // Now go through the multiscale configuration and adjust the coefficients.
    //

    // For each of the configurations
    for (size_t c = 0; c < globalConfiguration.dctScales.size(); c++) {

        scale_config_t  config = globalConfiguration.dctScales[c];
        size_t          sm = config.scale_multiplier;

        // Adjust the tx and ty coefficients for each supermacroblock
        if (sm > 1) {

            size_t scaledBlockWidth = BLOCK_WIDTH * sm;
            size_t scaledBlockHeight = BLOCK_HEIGHT * sm;
            size_t scaledTilesWide = CEIL(display_->DisplayWidth(), scaledBlockWidth);
            size_t scaledTilesHigh = CEIL(display_->DisplayHeight(), scaledBlockHeight);

#ifndef NO_OMP
#pragma omp parallel for
#endif
            // Break up the display into supermacroblocks at this current scale
            for (int j = 0; j < scaledTilesHigh; j++) {
                for (int i = 0; i < scaledTilesWide; i++) {
                    for (int y = 0; y < scaledBlockHeight && j * scaledBlockHeight + y < display_->DisplayHeight(); y++) {
                        for (int x = 0; x < scaledBlockWidth && i * scaledBlockWidth + x < display_->DisplayWidth(); x++) {
                            vector<unsigned int> start(3), end(3);

                            start[0] = i * scaledBlockWidth + x;
                            start[1] = j * scaledBlockHeight + y;
                            start[2] = config.first_plane_idx;

                            end[0] = i * scaledBlockWidth + x;
                            end[1] = j * scaledBlockHeight + y;
                            end[2] = config.first_plane_idx + config.plane_count - 1;

                            int tx = -i * scaledBlockWidth;
                            int ty = -j * scaledBlockHeight;

                            display_->FillCoefficient(tx, 0, 2, start, end);
                            display_->FillCoefficient(ty, 1, 2, start, end);
                        }
                    }
                }
            }
        }

        // Adjust the k coefficient for each supermacroblock. This will be done one plane
        // at a time starting with the topmost plane (0) and the first configuration. The
        // k used won't necessarily match the plane, because of the configuration.
        start[0] = 0; start[1] = 0;
        end[0] = display_->DisplayWidth() - 1; end[1] = display_->DisplayHeight() - 1;
#ifdef SIMPLE_TRUNCATION
        for (size_t p = 0; p < config.plane_count; p++) {
            size_t k = config.first_plane_idx + p;
            start[2] = k; end[2] = k;
            display_->FillCoefficient(p, 2, 2, start, end);
        }
#else
        // If this is an 8x8 region of macroblocks, then just use simple truncation
        if (config.edge_length == 8) {
            for (size_t p = 0; p < config.plane_count; p++) {
                size_t k = config.first_plane_idx + p;
                start[2] = k; end[2] = k;
                display_->FillCoefficient(p, 2, 2, start, end);
            }
        // Else consider only the most significant coefficients within with the square of edge_length
        } else {
            for (size_t p = 0, y = 0; y < config.edge_length && p < config.plane_count; y++) {
                for (size_t x = 0; x < config.edge_length && p < config.plane_count; x++) {
                    size_t k = config.first_plane_idx + p;
                    start[2] = k; end[2] = k;
                    display_->FillCoefficient(zigZag_[c][y * BLOCK_WIDTH + x], 2, 2, start, end);
                    p++;
                }
            }
        }
#endif
    }
}


/**
 * Initializes the Frame Volume for this tiler by pre-rendering each
 * of the 16 basis functions into 4x4 planes in the Frame Volume. They're
 * rendered for each color channel and stored in those groups of three in
 * zig-zag order.
 */
void MultiDctTiler::InitializeFrameVolume() {

    // Allocate the basisFunctions_. Assumes that the first configuration is the largest
    // scale. The basisFunctions_ (and thus the frame volume) will be filled with 64 planes
    // of that size. Therefore, we're using the simple 8x8 definitions for BASIS_BLOCKS_WIDE/TALL.
    basisFunctions_ = (Pixel *)calloc(fvWidth_ * fvHeight_ * BASIS_BLOCKS_WIDE * BASIS_BLOCKS_TALL,
                                      sizeof(Pixel));

    // For each scale level
    for (int c = 0; c < globalConfiguration.dctScales.size(); c++) {

        // Now we'll need to use the correct block_width/height, block_size, and basis_blocks_wide/tall
        // for this particular configuration
        scale_config_t  config = globalConfiguration.dctScales[c];
        size_t          sm = config.scale_multiplier;
        size_t          block_width = sm * UNSCALED_BASIC_BLOCK_WIDTH;
        size_t          block_height = sm * UNSCALED_BASIC_BLOCK_HEIGHT;
        size_t          block_size = block_width * block_height;
        size_t          basis_blocks_wide = block_width;
        size_t          basis_blocks_tall = block_height;
        double          alpha0 = 1 / (SQRT2 * 2.0);    // Note: The alpjas here are not adjusted. The are simply the alpha from
        double          alphaX = 1 / (2.0);            //       the JPEG DCT along with a 1/2 (which is 1/2 * 1/2 = 1/4 for both).
        double          scaledPi = PI / (8.0 * sm);

        // Pre-render each basis function
#ifndef NO_OMP
#pragma omp parallel for
#endif
        for (int j = 0; j < basis_blocks_tall; j++) {
            for (int i = 0; i < basis_blocks_wide; i++) {

                // Don't process the final basis block.
                if (i == basis_blocks_wide - 1 && j == basis_blocks_tall) continue;

                // Only process the number of planes in this configuration
                if (zigZag_[c][j * basis_blocks_wide + i] >= config.plane_count) continue;

                size_t p = zigZag_[c][j * basis_blocks_wide + i] * block_size;
                for (int y = 0; y < block_height; y++) {
                    for (int x = 0; x < block_width; x++) {
                        double m = 0.0;
                        bool neg = false;

                        for (int v = 0; v < block_height; v++) {
                            for (int u = 0; u < block_width; u++) {
                                double px = 1.0;
                                px *= (u == 0) ? alpha0 : alphaX;   // alpha(u)
                                px *= (v == 0) ? alpha0 : alphaX;   // alpha(v)
                                px *= ((u == i) && (v == j)) ? double(MAX_DCT_COEFF) : 0.0;  // DCT coefficient (maximum on or off)
                                px *= cos(scaledPi * ((double)x + 0.5) * (double)u);         // cos with x, u
                                px *= cos(scaledPi * ((double)y + 0.5) * (double)v);         // cos with y, v
                                m += px;
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
                        basisFunctions_[p].r = c;
                        basisFunctions_[p].g = c;
                        basisFunctions_[p].b = c;
                        basisFunctions_[p].a = 0xff;

                        // Move to the next pixel
                        p++;
                    }
                }
            }
        }
    }

    // Then render the gray block as the actual last basis block
    size_t p = fvWidth_ * fvHeight_ * (FRAMEVOLUME_DEPTH - 1);
    for (int y = 0; y < BLOCK_HEIGHT; y++) {
        for (int x = 0; x < BLOCK_WIDTH; x++) {
            unsigned int c = 0x7f;
            basisFunctions_[p].r = c;
            basisFunctions_[p].g = c;
            basisFunctions_[p].b = c;
            basisFunctions_[p].a = 0xff;
            p++;
        }
    }

    // Update the frame volume with the basis function renderings and gray block in bulk.
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(fvWidth_ - 1); end.push_back(fvHeight_ - 1); end.push_back(FRAMEVOLUME_DEPTH - 1);
    display_->CopyPixels(basisFunctions_, start, end);
}

/**
 * Builds the coefficients for the macroblock specified by (i, j) using the source buffer.
 *
 * @param i The column component of the macroblock
 * @param j The row component of the macroblock
 * @param buffer The source buffer. Note this may or may not be a downscaled image.
 * @param width The width of the source buffer
 * @param height The height of the source buffer
 * @param c Index into globalConfiguration.dctScales for information about the factor by
 *          which we're scaling and which planes to fill.
 * @return The vector (in zig-zag order) of 4-channel scalers.
 */
vector<uint64_t> MultiDctTiler::BuildCoefficients(size_t i, size_t j, int16_t* buffer, size_t width, size_t height, size_t c, bool shift) {

    /*
     * Produces the de-quantized coefficients for the input buffer using the following steps:
     *
     * 1. Shift by subtracting 128
     * 2. Take the 2D DCT
     * 3. Quantize
     * 4. De-quantize
     */

    Scaler  s;
    size_t  sm = globalConfiguration.dctScales[c].scale_multiplier;
    size_t  block_width = sm * UNSCALED_BASIC_BLOCK_WIDTH;
    size_t  block_height = sm * UNSCALED_BASIC_BLOCK_HEIGHT;
    size_t  block_size = block_width * block_height;
    double  alpha0 = 1 / (SQRT2 * 2.0 * sm);    // Note: The alphas here don't just include the 1/2 for each. The 1/2 is actually scaled when
    double  alphaX = 1 / (2.0 * sm);            //       building the coefficients so that the accumulated sums divide by MAX_DCT_COEFF cleanly.
    double  scaledPi = PI / (8.0 * sm);

    /* The coefficients are stored in this array in zig-zag order */
    vector<uint64_t> coefficients(block_size, 0);

    for (size_t v = 0; v < block_height; v++) {
#ifndef NO_OMP
#pragma omp parallel for ordered
#endif
        for (size_t u = 0; u < block_width; u++) {

            double c_r = 0.0, c_g = 0.0, c_b = 0.0;

            /* Calculate G for each u, v using g(x,y) shifted (1. and 2.) */
            for (size_t y = 0; y < block_height; y++) {
                size_t bufPos = ((j * block_height + y) * width + i * block_width + 0) * 3;
                for (size_t x = 0; x < block_width; x++) {

                    double p = 1.0;

                    p *= (u == 0) ? alpha0 : alphaX;                                     // alpha(u)
                    p *= (v == 0) ? alpha0 : alphaX;                                     // alpha(v)
                    p *= cos(scaledPi * ((double)x + 0.5) * (double)u);                  // cos with x, u
                    p *= cos(scaledPi * ((double)y + 0.5) * (double)v);                  // cos with y, v
                    /* Fetch each channel, multiply by product p and then shift */
                    if ( ((i * block_width + x) < width)
                            && ((j * block_height + y) < height) ) {
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
            size_t matPos = v * block_width + u;

            /* Quantize G(u,v) (3.) */
            g_r = (int)int(c_r / double(quantizationMatrix_[c][matPos]) + 0.5); // 0.5 is for rounding
            g_g = (int)int(c_g / double(quantizationMatrix_[c][matPos]) + 0.5);
            g_b = (int)int(c_b / double(quantizationMatrix_[c][matPos]) + 0.5);

            /* De-quantized G(u,v) (4.) */
            g_r *= (int)quantizationMatrix_[c][matPos];
            g_g *= (int)quantizationMatrix_[c][matPos];
            g_b *= (int)quantizationMatrix_[c][matPos];

            /* Set the coefficient in zig-zag order. */
            size_t p = zigZag_[c][matPos];

            /* Skip the last block, because we've used it for the medium gray block */
            if (p == block_size - 1) continue;

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

    scale_config_t  config = globalConfiguration.dctScales[c];
    size_t          block_width = config.scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
    size_t          block_height = config.scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;

    start[0] = i * block_width;
    start[1] = j * block_height;
    start[2] = first;

    size[0] = block_width;
    size[1] = block_height;

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

    scale_config_t  config = globalConfiguration.dctScales[c];
    size_t          block_width = config.scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
    size_t          block_height = config.scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;

    // For this configuration, determine the offset into the basisFunctions_ array since there are
    // likely other planes before it from prior configurations. This offset is in terms of pixels and
    // not bytes.
    size_t ibfo = 0;
    for (size_t d = 0; d < c; d++) {
        ibfo += fvWidth_ * fvHeight_ * globalConfiguration.dctScales[d].plane_count;
    }

#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (size_t y = 0; y < block_height; y++) {
        for (size_t x = 0; x < block_width; x++) {

            if (i * block_width + x >= width || j * block_height + y >= height)
                continue;

            int rAccumulator = 0, gAccumulator = 0, bAccumulator = 0;

#ifdef SIMPLE_TRUNCATION
            for (size_t p = 0; p < coefficients.size(); p++) {
                Scaler s;
                s.packed = coefficients[p];

                size_t bfo = ibfo + p * block_width * block_height + y * block_width + x;
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

                    size_t bfo = ibfo + p * block_width * block_height + y * block_width + x;
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

                        size_t pp = zigZag_[c][yy * block_width + xx];

                        size_t bfo = ibfo + pp * block_width * block_height + y * block_width + x;
                        rAccumulator += (int8_t)(basisFunctions_[bfo].r) * s.r;
                        gAccumulator += (int8_t)(basisFunctions_[bfo].g) * s.g;
                        bAccumulator += (int8_t)(basisFunctions_[bfo].b) * s.b;

                        p++;
                    }
                }
            }
#endif

            size_t offset = ((j * block_height + y) * width + i * block_width + x) * 3;

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

        scale_config_t  config = globalConfiguration.dctScales[c];
        size_t          block_width = config.scale_multiplier * UNSCALED_BASIC_BLOCK_WIDTH;
        size_t          block_height = config.scale_multiplier * UNSCALED_BASIC_BLOCK_HEIGHT;

        /*
         * 2. For each macroblock in downsampled image
         *    a. Build coefficients - BuildCoefficients()
         */
        int16_t* rendBuf = (int16_t*)calloc(width * height * 3, sizeof(int16_t));

        size_t tilesWide = CEIL(width, block_width);
        size_t tilesHigh = CEIL(height, block_height);

        vector< vector< vector<uint64_t> > > coefficientsForCurrentScale;

        for (size_t i = 0; i < tilesWide; i++) {
            coefficientsForCurrentScale.push_back(vector< vector<uint64_t> >());
            for (size_t j = 0; j < tilesHigh; j++) {
                // a. Build coefficients - BuildCoefficients()
                vector<uint64_t> coefficients = BuildCoefficients(i, j, signedBuf, width, height, c, c == 0);
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
                PrerenderCoefficients(coefficientsForCurrentScale[i][j], i, j, c, rendBuf, width, height, c == 0);
            }
        }

        /*
         * 6. Subtract the results from the original image - AdjustFrame()
         */
        AdjustFrame(signedBuf, rendBuf, width, height);

        /*
         * 7. Cleanup
         */
        free(rendBuf);
    }

#ifdef SKIP_FRAMES
    frame++;
#endif

    // Finally clean the signedBuf that we've been using throughout
    free(signedBuf);
}
