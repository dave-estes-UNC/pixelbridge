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
#include "nddi/BaseNddiDisplay.h"

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
     * @param display_width The width of the display
     * @param display_height The height of the display
     * @param quality The quality factor used for DCT.
     */
    ItTiler(size_t display_width, size_t display_height, size_t quality);

    ~ItTiler() {
    }

    /**
     * Returns the Display created and initialized by the tiler.
     */
    GlNddiDisplay* GetDisplay();

    /**
     * Update the scalers and then the NDDI display based on the frame that's passed in.
     *
     * @param buffer Pointer to the return frame buffer
     * @param width The width of that frame buffer
     * @param height The height of that frame buffer
     */
    void UpdateDisplay(uint8_t* buffer, size_t width, size_t height);

    /**
     * Calculates the costs for rendering without actually rendering.
     */
    void SimulateRenderCosts(bool force = false);

private:
    void InitializeCoefficientPlanes();
    void InitializeFrameVolume();
    void initZigZag();
    void setQuality(uint32_t quality);
    void forwardIntegerTransform(int *Y, int *X);
    void inverseIntegerTransform(int *Z, int *Y); // Added for completeness, but not used.

private:
    static const size_t  BLOCK_WIDTH = 4;
    static const size_t  BLOCK_HEIGHT = 4;
    static const size_t  BLOCK_SIZE = BLOCK_WIDTH * BLOCK_HEIGHT;
    static const size_t  FRAMEVOLUME_DEPTH  = BLOCK_SIZE;
    static const size_t  BASIS_BLOCKS_WIDE = 4;
    static const size_t  BASIS_BLOCKS_TALL = 4;

    static const size_t  MAX_IT_COEFF = 64;

    static int Cf4[BLOCK_SIZE];
    static int Cf4T[BLOCK_SIZE];
    static int Ci4[BLOCK_SIZE];
    static int Ci4T[BLOCK_SIZE];

    int Mf4[BLOCK_SIZE];
    int Vi4[BLOCK_SIZE];

    uint32_t qp;
    uint32_t qp6;

    GlNddiDisplay*    display_;
    bool              quiet_;
    int               zigZag_[BLOCK_SIZE];
};

#endif /* defined(__pixelbridge__ItTiler__) */
