#include <cmath>
#include <iostream>

#include "PixelBridgeFeatures.h"
#include "Configuration.h"
#include "DctTiler.h"

#define PI    3.14159265
#define PI_8  0.392699081

// This will turn on the simple truncation used before that simply takes the top X coefficients.
// The new scheme will instead take the most significant coefficients in a square pattern.
//#define SIMPLE_TRUNCATION

// Snapping uses a DCT Delta to snap calculated coefficients to zero if they're within the delta.
// Snapping only occurs in the last 1x scale.
//   SNAP_TO_ZERO is the top level flag that turns this on.
//   FORCE_SNAP_DC snaps the DC regarless of the delta.
//   SNAP_DC_ONLY restricts to the snapping with delta to just the DC
//#define SNAP_TO_ZERO
//#define FORCE_SNAP_DC
//#define SNAP_DC_ONLY

// Trimming with delta is seperate from snapping. This changes the regular trim routine
// That doesn't resend the unchanged coefficients at the top and bottom of a stack to
// consider coefficients within a certain delta to be unchanged.
//#define TRIM_WITH_DELTA

// This uses the scale to determine how often to send the coefficients at scales greater than 1x.
// e.g. The 16x scale coefficients are only sent every 16th frame. The 4x every 4th frame.
#define SKIP_FRAMES


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


/**
 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the FlatTiler should be destroyed and re-created.
 */
DctTiler::DctTiler (size_t display_width, size_t display_height,
                    size_t quality)
{
    quiet_ = globalConfiguration.headless || !globalConfiguration.verbose;

    /* 3 dimensional matching the Macroblock Width x Height x 64 */
    vector<unsigned int> fvDimensions;
    fvDimensions.push_back(BLOCK_WIDTH);
    fvDimensions.push_back(BLOCK_HEIGHT);
    fvDimensions.push_back(FRAMEVOLUME_DEPTH);

    /*
     * Pre-calculate the number of tile used for the display and initialize
     * the tile stack heights that are used to efficiently update the scalers
     * making sure to not update more than needed.
     */
    displayTilesWide_ = CEIL(display_width, BLOCK_WIDTH);
    displayTilesHigh_ = CEIL(display_height, BLOCK_HEIGHT);
    tileStackHeights_ = (uint8_t*)calloc(displayTilesWide_ * displayTilesHigh_, sizeof(uint8_t));

#ifndef NO_CL
    display_ = new ClNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 FRAMEVOLUME_DEPTH,             // Number of coefficient planes
                                 3);                            // Input vector size (x, y, 1)
#else
    display_ = new GlNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 FRAMEVOLUME_DEPTH,             // Number of coefficient planes
                                 3);                            // Input vector size (x, y, 1)
#endif

    /* Set the full scaler value and the sign mode */
    display_->SetFullScaler(MAX_DCT_COEFF);
    display_->SetPixelByteSignMode(SIGNED_MODE);

    /*
     * Initialize the zig-zag order used throughout and the calculate
     * the quantization matrix
     */
    initZigZag();
    initQuantizationMatrix(quality);

    /* Initialize Input Vector */
    vector<int> iv;
    iv.push_back(1);
    display_->UpdateInputVector(iv);

    /* Initialize Coefficient Planes */
    InitializeCoefficientPlanes();

    /* Initialize Frame Volume */
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

    Scaler s;

    /*
     * Perform the basic initialization, which ignores the scaling. Will
     * adjust those coefficients afterwards.
     */

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

    // Break the display into macroblocks and initialize each cube of coefficients to pick out the proper block from the frame volume
    for (int j = 0; j < displayTilesHigh_; j++) {
        for (int i = 0; i < displayTilesWide_; i++) {
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
    s.packed = 0;
    display_->FillScaler(s, start, end);

    // Fill the scalers for the medium gray plane to full on
    start[2] = display_->NumCoefficientPlanes() - 1;
    s.r = s.g = s.b = display_->GetFullScaler();
    display_->FillScaler(s, start, end);


    /*
     * Now go through the multiscale configuration and adjust the coefficients.
     */

    // For each of the configurations
    for (size_t c = 0; c < globalConfiguration.dctScales.size(); c++) {

        scale_config_t config = globalConfiguration.dctScales[c];

        // Adjust the tx and ty coefficients for each supermacroblock
        if (config.scale_multiplier > 1) {

            size_t scaledBlockWidth = BLOCK_WIDTH * config.scale_multiplier;
            size_t scaledBlockHeight = BLOCK_HEIGHT * config.scale_multiplier;
            size_t scaledTilesWide = CEIL(display_->DisplayWidth(), scaledBlockWidth);
            size_t scaledTilesHigh = CEIL(display_->DisplayHeight(), scaledBlockHeight);

            // Break up the display into supermacroblocks at this current scale
            for (int j = 0; j < scaledTilesHigh; j++) {
                for (int i = 0; i < scaledTilesWide; i++) {
                    for (int y = 0; y < scaledBlockHeight && j * scaledBlockHeight + y < display_->DisplayHeight(); y++) {
                        for (int x = 0; x < scaledBlockWidth && i * scaledBlockWidth + x < display_->DisplayWidth(); x++) {

                            start[0] = i * scaledBlockWidth + x;
                            start[1] = j * scaledBlockHeight + y;
                            start[2] = config.first_plane_idx;

                            end[0] = i * scaledBlockWidth + x;
                            end[1] = j * scaledBlockHeight + y;
                            end[2] = config.first_plane_idx + config.plane_count - 1;

                            int tx = -i * scaledBlockWidth - x + x / config.scale_multiplier;
                            int ty = -j * scaledBlockHeight - y + y / config.scale_multiplier;

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
                    display_->FillCoefficient(zigZag_[y * BLOCK_WIDTH + x], 2, 2, start, end);
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
void DctTiler::InitializeFrameVolume() {

    basisFunctions_ = (Pixel *)calloc(BLOCK_SIZE * FRAMEVOLUME_DEPTH, sizeof(Pixel));

    // Pre-render each basis function
#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (int j = 0; j < BASIS_BLOCKS_TALL; j++) {
        for (int i = 0; i < BASIS_BLOCKS_WIDE; i++) {

            // Don't process the final basis block.
            if (i == BASIS_BLOCKS_WIDE - 1 && j == BASIS_BLOCKS_TALL) continue;

            size_t p = zigZag_[j * BASIS_BLOCKS_WIDE + i] * BLOCK_SIZE;

            for (int y = 0; y < BLOCK_HEIGHT; y++) {
                for (int x = 0; x < BLOCK_WIDTH; x++) {
                    double m = 0.0;
                    bool neg = false;

                    for (int v = 0; v < BLOCK_HEIGHT; v++) {
                        for (int u = 0; u < BLOCK_WIDTH; u++) {
                            double px = 1.0;
                            px *= (u == 0) ? SQRT_125 : SQRT_250;   // alpha(u)
                            px *= (v == 0) ? SQRT_125 : SQRT_250;   // alpha(v)
                            px *= ((u == i) && (v == j)) ? double(MAX_DCT_COEFF) : 0.0;  // DCT coefficient (maximum on or off)
                            px *= cos(PI_8 * ((double)x + 0.5) * (double)u);         // cos with x, u
                            px *= cos(PI_8 * ((double)y + 0.5) * (double)v);         // cos with y, v
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

    // Then render the gray block as the actual last basis block
    size_t p = zigZag_[BASIS_BLOCKS_WIDE * BASIS_BLOCKS_TALL - 1] * BLOCK_SIZE;
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
    end.push_back(BLOCK_WIDTH - 1); end.push_back(BLOCK_HEIGHT - 1); end.push_back(FRAMEVOLUME_DEPTH - 1);
    display_->CopyPixels(basisFunctions_, start, end);
}


/**
 * Converts an unsigned byte buffer to a signed byte buffer; which is nothing more than
 * copying each byte to a short.
 *
 * @param buffer The input byte buffer
 * @param width The width of the buffers
 * @param weight The height of the buffers
 * @return mallocs and returns the signed buffer
 */
int16_t* DctTiler::ConvertToSignedPixels(uint8_t* buffer, size_t width, size_t height) {

    int16_t* signedBuf = (int16_t*)calloc(width * height * 3, sizeof(int16_t));

#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < (width * height * 3); i++) {
        signedBuf[i] = (int16_t)buffer[i];
    }

    return signedBuf;
}


/**
 * Takes a 24-bit pixel buffer and downsamples it by a factor with simple averaging. The source
 * buffer is likely not an even multiple of the factor provided, so edge "blocks" will use black
 * for any overscan pixels.
 *
 * @param factor The positive integer to scale by. Should likely be a power of two.
 * @param buffer The 24-bit RGB buffer.
 * @param width The width of the source buffer.
 * @param height The height of the source buffer.
 * @return The newly malloc'd buffer. Callee must free.
 */
int16_t* DctTiler::DownSample(size_t factor, int16_t* buffer, size_t width, size_t height) {

    size_t scaledWidth = CEIL(width, factor);
    size_t scaledHeight = CEIL(height, factor);

    int16_t* downBuf = (int16_t*)calloc(scaledWidth * scaledHeight * 3, sizeof(int16_t));

#ifndef NO_OMP
#pragma omp parallel for
#endif
   for (size_t j = 0; j < scaledHeight; j++) {
        for (size_t i = 0; i < scaledWidth; i++) {
            int64_t r, g, b;
            r = g = b = 0;

            // Sum the pixels over the larger region
            size_t count = 0;
            for (size_t y = j * factor; y < (j + 1) * factor && y < height; y++) {
                for (size_t x = i * factor; x < (i + 1) * factor && x < width; x++) {
                    size_t p = (y * width + x) * 3;
                    r += buffer[p + 0];
                    g += buffer[p + 1];
                    b += buffer[p + 2];
                    count++;
                }
            }

            // And take the average as the return pixel for (i, j)
            size_t p = (j * scaledWidth + i) * 3;
            downBuf[p + 0] = r / count;
            downBuf[p + 1] = g / count;
            downBuf[p + 2] = b / count;
        }
    }

    return downBuf;
}

/**
 * Takes a previously downsampled 24-bit pixel buffer and upsamples it by a factor. The destination
 * buffer is likely not an even multiple of the factor provided, so any overscanned pixels will
 * be clamped.
 *
 * @param factor The positive integer to scale by. Should likely be a power of two.
 * @param buffer The 24-bit RGB source buffer.
 * @param width The width of the *destination* buffer.
 * @param height The height of the *destination* buffer.
 * @return The newly malloc'd destination buffer. Callee must free.
 */
int16_t* DctTiler::UpSample(size_t factor, int16_t* buffer, size_t width, size_t height) {

    size_t scaledHeight = CEIL(height, factor);
    size_t scaledWidth = CEIL(width, factor);

    int16_t* upBuf = (int16_t*)calloc(width * height * 3, sizeof(int16_t));

#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (size_t j = 0; j < scaledHeight; j++) {
        for (size_t i = 0; i < scaledWidth; i++) {
            int16_t r, g, b;

            size_t p = (j * scaledWidth + i) * 3;
            r = buffer[p + 0];
            g = buffer[p + 1];
            b = buffer[p + 2];

            for (size_t y = j * factor; y < (j + 1) * factor && y < height; y++) {
                for (size_t x = i * factor; x < (i + 1) * factor && x < width; x++) {
                    size_t p = (y * width + x) * 3;
                    upBuf[p + 0] = r;
                    upBuf[p + 1] = g;
                    upBuf[p + 2] = b ;
                }
            }
       }
    }

    return upBuf;
}

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
vector<uint64_t> DctTiler::BuildCoefficients(size_t i, size_t j, int16_t* buffer, size_t width, size_t height, bool shift) {

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
 * Selects which coefficients to use at this scale. This function will keep track of our budget of planes
 * and if we run out it will either truncate the extra coefficients in a simple manner or will use the
 * region of most significant coefficients defined by the edge_length.
 *
 * @param coefficients The coefficients to be reduced to fit within the budget.
 * @param c The current scale from the globalConfiguration.dctScales global.
 */
void DctTiler::SelectCoefficientsForScale(vector<uint64_t> &coefficients, size_t c) {
    scale_config_t config = globalConfiguration.dctScales[c];
#ifdef SIMPLE_TRUNCATION
    coefficients.resize(config.plane_count);
#else
    // For an 8x8 region of coefficients, just use simple truncation
    if (config.edge_length == 8) {
        coefficients.resize(config.plane_count);
    // Else consider only the most significant coefficients within with the square of edge_length
    } else {
        size_t p = 0;
        vector<uint64_t> coeffs;
        for (size_t y = 0; y < config.edge_length && p < config.plane_count; y++) {
            for (size_t x = 0; x < config.edge_length && p < config.plane_count; x++) {
                coeffs.push_back(coefficients[zigZag_[y * BLOCK_WIDTH + x]]);
                p++;
            }
        }
        coefficients = coeffs;
    }
#endif
}


/**
 * Considers snapping to zero any coefficients that are within a delta of zero. While it will work on
 * any scale, this is only used in practice on the 1x scale.
 *
 * @param coefficientsForScale All of the coefficient stacks for this scale.
 */
void DctTiler::SnapCoefficientsToZero(vector< vector< vector<uint64_t> > > &coefficientsForScale) {

    for (size_t i = 0; i < coefficientsForScale.size(); i++) {
        for (size_t j = 0; j < coefficientsForScale[i].size(); j++) {
#if defined FORCE_SNAP_DC
            coefficientsForScale[i][j][0] = 0;
#elif defined SNAP_DC_ONLY
            Scaler s;
             s.packed = coefficientsForScale[i][j][0];
             if (abs(s.r) <= globalConfiguration.dctDelta &&
                 abs(s.g) <= globalConfiguration.dctDelta &&
                 abs(s.b) <= globalConfiguration.dctDelta)
             {
                 coefficientsForScale[i][j][0] = 0;
             }
#else
            for (int k = 0; k < coefficientsForScale[i][j].size(); k++) {
                Scaler s;
                s.packed = coefficientsForScale[i][j][k];
                if (abs(s.r) <= globalConfiguration.dctDelta &&
                    abs(s.g) <= globalConfiguration.dctDelta &&
                    abs(s.b) <= globalConfiguration.dctDelta)
                {
                    coefficientsForScale[i][j][k] = 0;
                }
            }
#endif
        }
    }
}


/**
 * Estimates what the cost of updating the coefficients will be when all but a defined range (front to last inclusive) are
 * zeroed out.
 *
 * @param coefficientsForScale The calculated coefficients for this scale
 * @param c The current scale from the globalConfiguration.dctScales global.
 * @param first The first non-zeroed plane.
 * @param last The last non-zeroed plane.
 * @return The estimated cost in terms of bytes sent over the wire.
 */
size_t DctTiler::EstimateCost(vector< vector< vector<uint64_t> > > &coefficientsForScale, size_t c, size_t first, size_t last) {

    scale_config_t config = globalConfiguration.dctScales[c];
    size_t scaleFirst = config.first_plane_idx;
    size_t scaleLast = config.first_plane_idx + config.plane_count - 1;

#if defined FIRST_NON_ZEROED_PLANE && defined LAST_NON_ZEROED_PLANE
    if (first == CLAMP(FIRST_NON_ZEROED_PLANE, scaleFirst, scaleLast + 1) && last == LAST_NON_ZEROED_PLANE)
        return 0;
    else
        return 1000;
#elif defined FIRST_NON_ZEROED_PLANE
    if (first == CLAMP(FIRST_NON_ZEROED_PLANE, scaleFirst, scaleLast + 1))
        return 0;
    else
        return 1000;
#elif defined LAST_NON_ZEROED_PLANE
    if (last == CLAMP(FIRST_NON_ZEROED_PLANE, scaleFirst, scaleLast + 1))
        return 0;
    else
        return 1000;
#endif

    /*
     * Use the dimensions of the cachedCoefficients for this scale to set the
     * dimensions of the tops and bottoms vector.
     */
    vector< vector<int> > firsts(cachedCoefficients_[c].size(), vector<int>(cachedCoefficients_[c][0].size(), -1));
    vector< vector<int> >  lasts(cachedCoefficients_[c].size(), vector<int>(cachedCoefficients_[c][0].size(), -1));

    /*
     * Work through the planes to and setup the tops and bottoms.
     */
    for (size_t i = 0; i < cachedCoefficients_[c].size(); i++) {
        for (size_t j = 0; j < cachedCoefficients_[c][i].size(); j++) {
            /* Work from the front (top) and track the tops */
            for (int k = config.first_plane_idx; k < config.first_plane_idx + config.plane_count; k++) {
                if (firsts[i][j] == -1) {
                    if (k < first && coefficientsForScale[i][j][k] != 0)
                        firsts[i][j] = k;
                    else if (k > last && coefficientsForScale[i][j][k] != 0)
                        firsts[i][j] = k;
                    else if (coefficientsForScale[i][j][k] != cachedCoefficients_[c][i][j][k])
                        firsts[i][j] = k;
                }
            }
            /* Then work from the last (bottom) and track the bottoms */
            for (int k = config.first_plane_idx + config.plane_count - 1; k+1 > config.first_plane_idx; k--) {
                if (lasts[i][j] == -1) {
                    if (k > last && coefficientsForScale[i][j][k] != 0)
                        lasts[i][j] = k;
                    else if (k < first && coefficientsForScale[i][j][k] != 0)
                        lasts[i][j] = k;
                    else if (coefficientsForScale[i][j][k] != cachedCoefficients_[c][i][j][k])
                        lasts[i][j] = k;
                }
            }
        }
    }

    /*
     * Then tally the score
     */
    size_t score = 0;

    // Add the cost of zeroing sections (if any)
    if (first > config.first_plane_idx)
        score += BYTES_PER_SCALER * 1 + CALC_BYTES_FOR_CP_COORD_TRIPLES(2);

    if (last < config.first_plane_idx + config.plane_count - 1)
        score += BYTES_PER_SCALER * 1 + CALC_BYTES_FOR_CP_COORD_TRIPLES(2);

    for (size_t i = 0; i < cachedCoefficients_[c].size(); i++) {
        for (size_t j = 0; j < cachedCoefficients_[c][i].size(); j++) {
            // Add cost of filling a stack of tiles
            if (firsts[i][j] != -1 && lasts[i][j] != -1) {
                size_t stack_height = (lasts[i][j] - firsts[i][j] + 1);
                score += BYTES_PER_SCALER * stack_height +
                         CALC_BYTES_FOR_CP_COORD_TRIPLES(1) +
                         CALC_BYTES_FOR_TILE_COORD_DOUBLES(1);
            }
        }
    }

    return score;
}


/**
 * Compares all of the coefficients computed for this current scale to what's in the cache for this scale.
 * Will then choose which planes to zero. Two directions are considered, which may result in the first X
 * planes and the last Y planes cleared. Alternatively, all planes might be cleared.
 *
 * After the operation, the fill scaler commands are sent to the NDDI display and the cached coefficients
 * are updated.
 */
void DctTiler::ZeroPlanes(vector< vector< vector<uint64_t> > > &coefficientsForScale, size_t c) {

    scale_config_t config = globalConfiguration.dctScales[c];
    size_t scaleFirst = config.first_plane_idx;
    size_t scaleLast = config.first_plane_idx + config.plane_count - 1;

    vector<unsigned int> start(3, 0);
    vector<unsigned int> end(3, 0);

    /*
     * Determine the best regions to clear. Start from the front, and work from the last.
     */
    int first = scaleFirst;
    int last  = scaleLast;

    // Get a baseline cost over the whole set of planes
    size_t minCost = EstimateCost(coefficientsForScale, c, first, last);
    size_t cost;

    // Then adjust first until we find a minimal cost
    if (first != last)
        for (int i = first + 1; i <= last; i++) {
            cost = EstimateCost(coefficientsForScale, c, i, last);
            if (cost < minCost) {
                minCost = cost;
                first = i;
            }
        }

    // Then adjust last until we find a minimal cost
    if (first != last)
        for (int i = last - 1; i + 1 > first; i--) {
            cost = EstimateCost(coefficientsForScale, c, first, i);
            if (cost < minCost) {
                minCost = cost;
                last = i;
            }
        }

    // Lastly check if it's cheaper to zero everything out by setting first past the end
    cost = EstimateCost(coefficientsForScale, c, config.first_plane_idx + config.plane_count, config.first_plane_idx + config.plane_count);
    if (cost < minCost) {
        minCost = cost;
        first = last = config.first_plane_idx + config.plane_count;
    }

    if (globalConfiguration.verbose)
        cout << "Scale [" << scaleFirst << "," << scaleLast << "] -- Not Zeroed [" << first << "," << last << "]" << endl;

    /*
     * Configure the first region to clear, then clear the scalers and then the cached coefficients.
     */
    start[0] = start[1] = 0;
    start[2] = scaleFirst;

    end[0] = display_->DisplayWidth() - 1;
    end[1] = display_->DisplayHeight() - 1;
    end[2] = first - 1;

    if (first > scaleFirst) {
        /* Zero the scalers for these planes first. */
        Scaler s;
        s.packed = 0;
        display_->FillScaler(s, start, end);

        /* Then update the cached coefficients. */
        if (cachedCoefficients_.size() > c) {
            for (size_t i = 0; i < cachedCoefficients_[c].size(); i++) {
                for (size_t j = 0; j < cachedCoefficients_[c][i].size(); j++) {
                    for (size_t k = start[2]; k <= end[2]; k++) {
                        cachedCoefficients_[c][i][j][k - scaleFirst] = 0;
                    }
                }
            }
        }
    }

    /*
     * Configure the last region to clear, then clear the scalers and then the cached coefficients.
     */
    start[0] = start[1] = 0;
    start[2] = last + 1;

    end[0] = display_->DisplayWidth() - 1;
    end[1] = display_->DisplayHeight() - 1;
    end[2] = scaleLast;

    if (last < scaleLast) {
        /* Zero the scalers for these planes first. */
        Scaler s;
        s.packed = 0;
        display_->FillScaler(s, start, end);

        /* Then update the cached coefficients. */
        if (cachedCoefficients_.size() > c) {
            for (size_t i = 0; i < cachedCoefficients_[c].size(); i++) {
                for (size_t j = 0; j < cachedCoefficients_[c][i].size(); j++) {
                    for (size_t k = start[2]; k <= end[2]; k++) {
                        cachedCoefficients_[c][i][j][k - scaleFirst] = 0;
                    }
                }
            }
        }
    }
}


/**
 * Will trim the leading and trailing coefficients that don't need to be updated.
 *
 * @param coefficients The vector (in zig-zag order) of coefficients for the macroblock
 * @param c Index into globalConfiguration.dctScales for information about the factor by
 *          which we're scaling and which planes to fill.
 * @return The first plane to be updated.
 */
size_t DctTiler::TrimCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c) {

    /* Get the set of cached coefficients for (c, i, j), building out the cache the first time. */
    if (cachedCoefficients_.size() <= c)
        cachedCoefficients_.push_back(vector< vector< vector<uint64_t> > >());
    if (cachedCoefficients_[c].size() <= i)
        cachedCoefficients_[c].push_back(vector< vector<uint64_t> >());
    if (cachedCoefficients_[c][i].size() <= j) {
        /* Initialize the coefficients to 0 */
        cachedCoefficients_[c][i].push_back(vector<uint64_t>(globalConfiguration.dctScales[c].plane_count, 0));
    }

    /* first_plane_idx will be the return value. Start it at the configured first plane. */
    size_t first_plane_idx = globalConfiguration.dctScales[c].first_plane_idx;

    /* Figure out the first and last unchanged coefficient */
    size_t first = 0, last = coefficients.size() - 1;
    bool foundFirst = false;
    for (size_t k = 0; k < coefficients.size(); k++) {
#ifdef TRIM_WITH_DELTA
        Scaler coeff, cachedCoeff;
        coeff.packed = coefficients[k];
        cachedCoeff.packed = cachedCoefficients_[c][i][j][k];
        if (abs(coeff.r - cachedCoeff.r) > globalConfiguration.dctDelta ||
            abs(coeff.g - cachedCoeff.g) > globalConfiguration.dctDelta ||
            abs(coeff.b - cachedCoeff.b) > globalConfiguration.dctDelta)
        {
#else
        if (coefficients[k] != cachedCoefficients_[c][i][j][k]) {
#endif
            last = k;
            if (!foundFirst) {
                first = k;
                foundFirst = true;
            }
        }
    }

    /* Store this set of coefficients in the cache. */
    cachedCoefficients_[c][i][j] = vector<uint64_t>(coefficients);

    /* Trim the zeros from the end and then the beginning */
    if (foundFirst) {
        if (last < coefficients.size() - 1)
            coefficients.erase(coefficients.begin() + last + 1, coefficients.end());
        coefficients.erase(coefficients.begin(), coefficients.begin() + first);
        return first_plane_idx + first;
    } else {
        // If we didn't find a first, then that means that we've trimmed everything away.
        // As in there are no changes to make. So return a first plan that's out of bounds.
        return display_->NumCoefficientPlanes();
    }
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
void DctTiler::FillCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, size_t first) {

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
void DctTiler::PrerenderCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, int16_t* buffer, size_t width, size_t height, bool shift) {

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
 * Subtracts the pre-rendered buffer from the original buffer.
 *
 * @param buffer The original buffer and also the destination buffer
 * @param renderedBuffer The pre-rendered buffer from the previous scale
 * @param width The width of both buffers
 * @param height The height of both buffers
 */
void DctTiler::AdjustFrame(int16_t* buffer, int16_t* renderedBuffer, size_t width, size_t height) {
#ifndef NO_OMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < (width * height * 3); i++) {
        int16_t adjustedChannel = buffer[i] - renderedBuffer[i];

        buffer[i] = adjustedChannel;
    }
}

void DctTiler::UpdateScaledDisplay(uint8_t* buffer, size_t width, size_t height) {

    assert(width >= display_->DisplayWidth());
    assert(height >= display_->DisplayHeight());

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
        if (config.scale_multiplier == 1)
            downBuf = signedBuf;
        else
            downBuf = DownSample(config.scale_multiplier, signedBuf, width, height);

        /*
         * 2. For each macroblock in downsampled image
         *    a. Build coefficients - BuildCoefficients()
         */
        int16_t* rendBuf = (int16_t*)calloc(downW * downH * 3, sizeof(int16_t));

        size_t tilesWide = CEIL(downW, BLOCK_WIDTH);
        size_t tilesHigh = CEIL(downH, BLOCK_HEIGHT);

        vector< vector< vector<uint64_t> > > coefficientsForCurrentScale;

        // Only shift pixels by 128 for this first configuration. Subsequent configurations
        // will only build only address the differences.
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
         * 4. For this current scale, calculate optimal planes to zero and trim coefficients
         *   If we've initialized the coefficient cache
         *   a. Zero out optimal planes - ZeroPlanes()
         */
#ifdef SNAP_TO_ZERO
        if (c > 0 && config.scale_multiplier == 1)
            SnapCoefficientsToZero(coefficientsForCurrentScale);
#endif
#ifdef USE_DCT_OPTIMAL_ZEROING
        if (cachedCoefficients_.size() > c)
            ZeroPlanes(coefficientsForCurrentScale, c);
#endif

        /*
         * 5. For each macroblock in downsampled image
         *   a. Trim a copy of the coefficients - TrimCoefficients()
         *   b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
         *   c. Perform simulated blending of basis functions and store back to downsampled image - PrerenderCoefficients()
         */
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
                size_t first_plane_idx = TrimCoefficients(coefficients, i, j, c);

                // b. Fill trimmed coefficients to super-macroblock's coefficient scalers - FillCoefficients()
                FillCoefficients(coefficients, i, j, c, first_plane_idx);
#endif

                // c. Perform simulated blending of basis functions - PrerenderCoefficients()
                // Again, only shift on the first plane
                PrerenderCoefficients(coefficientsForCurrentScale[i][j], i, j, c, rendBuf, downW, downH, c == 0);
            }
        }

        /*
         * 6. Upsample the image - UpSample()
         */
        int16_t* upBuf;
        if (config.scale_multiplier == 1)
            upBuf = rendBuf;
        else
            upBuf = UpSample(config.scale_multiplier, rendBuf, width, height);

        /*
         * 7. Subtract the results from the original image - AdjustFrame()
         */
        AdjustFrame(signedBuf, upBuf, width, height);

        /*
         * 8. Cleanup
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
#ifdef USE_MULTISCALE_DCT
    UpdateScaledDisplay(buffer, width, height);
#else
    vector<unsigned int> start(3, 0);
    vector<unsigned int> size(2, 0);
    Scaler s;

    assert(width >= display_->DisplayWidth());
    assert(height >= display_->DisplayHeight());

    start[0] = 0; start[1] = 0; start[2] = 0;
    size[0] = BLOCK_WIDTH;
    size[1] = BLOCK_HEIGHT;

    /*
     * Produces the de-quantized coefficients for the input buffer using the following steps:
     *
     * 1. Shift by subtracting 128
     * 2. Take the 2D DCT
     * 3. Quantize
     * 4. De-quantize
     */
    for (size_t j = 0; j < displayTilesHigh_; j++) {
        for (size_t i = 0; i < displayTilesWide_; i++) {

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
                            if ( ((i * BLOCK_WIDTH + x) < display_->DisplayWidth())
                                    && ((j * BLOCK_HEIGHT + y) < display_->DisplayHeight()) ) {
                                c_r += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // red: g(x,y) - 128
                                c_g += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // green: g(x,y) - 128
                                c_b += p * ((double)buffer[bufPos] - 128.0); bufPos++;           // blue: g(x,y) - 128
                            } else {
                                bufPos += 3;
                            }
                        }
                    }

                    int g_r, g_g, g_b;
                    size_t matPos = v * BLOCK_WIDTH + u;

                    /* Quantize G(u,v) (3.) */
                    g_r = (int)int(c_r / double(quantizationMatrix_[matPos]) + 0.5);
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

            /*
             * Determine the minimum height of the stack that needs to be sent and
             * update the coefficients vector size. The minimum stack height is
             * determined by the non-zero planes in this update, but it might be
             * slightly larger if we need to overwrite the larger non-zero planes from
             * the last update. So the stack height is the largest of the old stack height
             * and the current stack height.
             */
            size_t h = BLOCK_SIZE - 2;
            while (h >= 0 && coefficients[h] == 0)
                 h--;
            if (h < tileStackHeights_[j * displayTilesWide_ + i])
                coefficients.resize(tileStackHeights_[j * displayTilesWide_ + i] + 1);
            else
                coefficients.resize(h + 1);
            tileStackHeights_[j * displayTilesWide_ + i] = h;

            /* Send the NDDI command to update this macroblock's coefficients, one plane at a time. */
            start[0] = i * BLOCK_WIDTH;
            start[1] = j * BLOCK_HEIGHT;
            display_->FillScalerTileStack(coefficients, start, size);
        }
    }
#endif
}
