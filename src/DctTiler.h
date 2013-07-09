#ifndef DCT_TILER_H
#define DCT_TILER_H
/*
 *  DctTiler.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 11/26/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include "Tiler.h"
#include "BaseNddiDisplay.h"

using namespace nddi;
using namespace std;

/**
 * Macro Block dimensions and frame volume depth
 */
#define MACROBLOCK_WIDTH       8
#define MACROBLOCK_HEIGHT      8
#define DCT_FRAMEVOLUME_DEPTH  (64 * 3 + 1)

/**
 * This tiler will split provided frames into macroblocks and will perform the forward DCT
 * and create coefficients that will be used as the scalers for the coefficient planes.
 */
class DctTiler : public Tiler {

public:
	/**
	 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
	 * dimensions change, then the DctTiler should be destroyed and re-created.
	 *
	 * @param display A pointer to the NDDI display
	 */
	DctTiler(BaseNddiDisplay* display,
			 bool quiet);

	~DctTiler() {
	}

	void initZigZag();
	void initQuantizationMatrix(unsigned char quality);

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
	static const size_t  BLOCK_WIDTH = MACROBLOCK_WIDTH;
	static const size_t  BLOCK_HEIGHT = MACROBLOCK_HEIGHT;
	static const size_t  BLOCK_SIZE = BLOCK_WIDTH * BLOCK_HEIGHT;
	static const size_t  BASIS_BLOCKS_WIDE = 8;
	static const size_t  BASIS_BLOCKS_TALL = 8;
	static const size_t  MAX_DCT_COEFF = 256;
	static const double  PI = 3.14159265;

	BaseNddiDisplay*  display_;
	bool              quiet_;
	int               zigZag_[BLOCK_WIDTH * BLOCK_HEIGHT];
	unsigned char     quantizationMatrix_[BLOCK_WIDTH * BLOCK_HEIGHT];
};
#endif // DCT_TILER_H
