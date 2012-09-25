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
//#define NO_CL

/*
 * Running make with NO_HACKS=1 till turn this off, but uncommenting this will accomplish the same thing.
 */
//#undef SUPRESS_EXCESS_RENDERING

/*
 * This turns on the OpenCL profiling code which updates the cost model with time
 */
#define CL_PROFILING_ENABLED

#endif /* PIXELBRIDGEFEATURES_H_ */
