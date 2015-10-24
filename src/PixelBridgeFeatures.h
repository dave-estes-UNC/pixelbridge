/*
 * PixelBridgeFeatures.h
 *
 *  Created on: Sep 24, 2012
 *      Author: cdestes
 */

#ifndef PIXELBRIDGEFEATURES_H_
#define PIXELBRIDGEFEATURES_H_

/*
 * These features are all on by default. Uncomment a line (and cross fingers) to turn one off.
 */
#ifdef DEBUG
#define NO_OMP
#endif

//#define NO_OMP
//#define NO_GL
#define NO_CL


/*
 * Running make with NO_HACKS=1 should turn these off if using the Makefile for Linux,
 * but uncommenting these will accomplish the same thing. And in OS X it's the only way
 * since these hacks are hardcoded in the project for lack of recent Makefile support
 * for OS X.
 */
//#undef SUPRESS_EXCESS_RENDERING
//#undef SKIP_COMPUTE_WHEN_SCALER_ZERO

/*
 * This turns on the OpenCL profiling code which updates the cost model with time
 * NOTE: This triggers a lot of clFinish()
 */
//#define CL_PROFILING_ENABLED

/*
 * Sets the checksum calculating mode.
 */
#define CRC 1
#define ADLER 2
#define TRIVIAL 3

#define CHECKSUM_CALCULATOR CRC
//#define CHECKSUM_CALCULATOR ADLER
//#define CHECKSUM_CALCULATOR TRIVIAL

/*
 * Use the narrowed data fields instead of 4-byte words for EVERYTHING.
 */
#define USE_NARROW_DATA_FIELDS

/*
 * Uses the nDDI extension to update groups of tiles
 */
#define USE_COPY_PIXEL_TILES

/*
 * Uses multiple scaled 8x8 macroblocks for DCT or
 * uses true multi-DCT with 8x8, 16x16, etc. macroblocks.
 */
//#define USE_SCALED_DCT
#define USE_MULTI_DCT

/*
 * Resets the cost model immediately after the tiler setup is complete.
 * Useful for considering the performance of a tiler with said setup.
 */
#define CLEAR_COST_MODEL_AFTER_SETUP

/*
 * Creates random frame(s) instead of decoding video.
 */
//#define USE_RANDOM_PLAYER

/*
 * Turns on the VERY BROKEN asynchronous decoder.
 */
// TODO(CDE): Fix this junk
//#define USE_ASYNC_DECODER

/*
 * Used to dramatically narrow the various data stores. Can lead to bugs, so proceed carefully.
 */
#define NARROW_DATA_STORES

/*
 * Used to reduce the rendered frame to a much smaller region for testing on memory constrained devices.
 */
//#define USE_SMALLER_WINDOW

/*
 * Divides the average optical flow by the diagonal.
 */
#define NORMALIZE_FLOW_FOR_RESOLUTION

#endif /* PIXELBRIDGEFEATURES_H_ */
