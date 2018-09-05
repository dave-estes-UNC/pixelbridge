/*
 * PixelBridgeFeatures.h
 *
 *  Created on: Sep 24, 2012
 *      Author: cdestes
 */

#ifndef PIXELBRIDGEFEATURES_H_
#define PIXELBRIDGEFEATURES_H_

#include "nddi/Features.h"

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
 * When defined, then the RAM saving features for the coefficient planes are used by
 * any tiler that supports them.
 */

#define USE_RAM_SAVING_COEFFICIENT_PLANE_FEATURES

/*
 * Divides the average optical flow by the diagonal.
 */
#define NORMALIZE_FLOW_FOR_RESOLUTION

/*
 * Chooses which DCT tiler to use when running with --dctscales.
 * The ScaledDctTiler is the original implementation which always used
 * 8x8 macroblocks for the DCT transform. The results were then scaled.
 * The MultiDctTiler is a variation that actually uses macroblocks of each
 * size and then prerenders all of those basis functions into the frame volume.
 */
//#define USE_SCALED_DCT
#define USE_MULTI_DCT

#endif /* PIXELBRIDGEFEATURES_H_ */
