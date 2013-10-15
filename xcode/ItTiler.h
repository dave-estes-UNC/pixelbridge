/*
 * ItTiler.h
 *  pixelbridge
 *
 *  Implements the Integer Transform Tiler
 *
 *  Created by Dave Estes on 10/13/13.
 *  Copyright (c) 2013 Dave Estes. All rights reserved.
 */

#ifndef __pixelbridge__ItTiler__
#define __pixelbridge__ItTiler__

#include "Tiler.h"
#include "BaseNddiDisplay.h"

using namespace nddi;
using namespace std;

/**
 * This tiler will split provided frames into 4x4 macroblocks and will perform the forward
 * Integer Transform from AVC and and create coefficients that will be used as the scalers
 * for the coefficient planes.
 */
class ItTiler : public Tiler {
    
public:
	/**
	 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
	 * dimensions change, then the DctTiler should be destroyed and re-created.
	 *
	 * @param display A pointer to the NDDI display
	 */
	ItTiler(BaseNddiDisplay* display,
            bool quiet);
    
	~ItTiler() {
	}
    
    /**
     * Initializes the Coefficient Planes for this tiler.
     */
    void InitializeCoefficientPlanes();
    
    /**
     * Initializes the Frame Volume for this tiler.
     */
    void InitializeFrameVolume();
    
	/**
	 * Update the scalers and then the NDDI display based on the frame that's passed in.
	 *
	 * @param buffer Pointer to the return frame buffer
	 * @param width The width of that frame buffer
	 * @param height The height of that frame buffer
	 */
	void UpdateDisplay(uint8_t* buffer, size_t width, size_t height);


private:
    void initZigZag();
    void setQuality(uint32_t quality);
    void matrixMultiply(double *D, double *S1, double *S2);
    void hadamardMultiply(double *D, double *S1, double *S2);
    void scalarMultiply(double *D, double s1, double *S2);
    void forwardIntegerTransform(int *Y, int *X);
    void inverseIntegerTransform(int *Z, int *Y); // Added for completeness, but not used.

public:
	static const size_t  BLOCK_WIDTH = 4;
	static const size_t  BLOCK_HEIGHT = 4;
    static const size_t  FRAMEVOLUME_DEPTH  = BLOCK_WIDTH * BLOCK_HEIGHT * 3;
    
private:
	static const size_t  BLOCK_SIZE = BLOCK_WIDTH * BLOCK_HEIGHT;
	static const size_t  BASIS_BLOCKS_WIDE = 4;
	static const size_t  BASIS_BLOCKS_TALL = 4;

	static const size_t  MAX_IT_COEFF = 1 << 8;
    
    double Cf4[BLOCK_SIZE] = {
        1,   1,   1,   1,
        2,   1,  -1,  -2,
        1,  -1,  -1,   1,
        1,  -2,   2,  -1
    };
    double Cf4T[BLOCK_SIZE] = {
        1,   2,   1,   1,
        1,   1,  -1,  -2,
        1,  -1,  -1,   2,
        1,  -2,   1,  -1
    };
    double Ci4[BLOCK_SIZE] = {
        1,   1,   1,   1,
        1,  .5, -.5,  -1,
        1,  -1,  -1,   1,
        .5,  -1,   1, -.5
    };
    double Ci4T[BLOCK_SIZE] = {
        1,   1,   1,  .5,
        1,  .5,  -1,  -1,
        1, -.5,  -1,   1,
        1,  -1,   1, -.5
    };
    double Mf4[BLOCK_SIZE];
    double Vi4[BLOCK_SIZE];
    
    uint32_t qp = 0;
    
	BaseNddiDisplay*  display_;
	bool              quiet_;
	int               zigZag_[BLOCK_SIZE];
};

#endif /* defined(__pixelbridge__ItTiler__) */
