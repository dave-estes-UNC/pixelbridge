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
//#define NO_OMP
//#define NO_GL
#define NO_CL

/*
 * Running make with NO_HACKS=1 still turn this off, but uncommenting this will accomplish the same thing.
 */
//#undef SUPRESS_EXCESS_RENDERING

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
 * Uses the nDDI extension to update groups of tiles
 */
#define USE_COPY_PIXEL_TILES

/*
 * Creates random frame(s) instead of decoding video.
 */
//#define USE_RANDOM_PLAYER

#endif /* PIXELBRIDGEFEATURES_H_ */
