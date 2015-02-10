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
 * This tiler will split provided frames into macroblocks and will perform the forward DCT
 * and create coefficients that will be used as the scalers for the coefficient planes.
 */
class DctTiler : public Tiler {

public:
/**
 * The DctTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the DctTiler should be destroyed and re-created.
 *
 * @param display_width The width of the display
 * @param display_height The height of the display
 * @param quality The quality factor used for DCT.
 * @param quiet Used to squelch extra information output.
 */
DctTiler(size_t display_width, size_t display_height, size_t quality);

~DctTiler() {
    if (tileStackHeights_)
        free(tileStackHeights_);
    if (basisFunctions_)
        free(basisFunctions_);
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

private:
    void InitializeCoefficientPlanes();
    void InitializeFrameVolume();
    void initZigZag();
    void initQuantizationMatrix(size_t quality);
    int16_t* ConvertToSignedPixels(uint8_t* buffer, size_t width, size_t height);
    int16_t* DownSample(size_t factor, int16_t* buffer, size_t width, size_t height);
    int16_t* UpSample(size_t factor, int16_t* buffer, size_t width, size_t height);
    vector<uint64_t> BuildCoefficients(size_t i, size_t j, int16_t* buffer, size_t width, size_t height, bool adjustPixels);
    size_t EstimateCost(vector< vector< vector<uint64_t> > > &coefficientsForScale, size_t c, size_t first, size_t last);
    void ZeroPlanes(vector< vector< vector<uint64_t> > > &coefficientsForScale, size_t c);
    size_t TrimCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c);
    void FillCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, size_t first);
    void PrerenderCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, int16_t* renderedBuffer, size_t width, size_t height, bool shift);
    void AdjustFrame(int16_t* buffer, int16_t* renderedBuffer, size_t width, size_t height);
    void UpdateScaledDisplay(uint8_t* buffer, size_t width, size_t height);


private:
    static const size_t  BLOCK_WIDTH = 8;
    static const size_t  BLOCK_HEIGHT = 8;
    static const size_t  BLOCK_SIZE = BLOCK_WIDTH * BLOCK_HEIGHT;
    static const size_t  FRAMEVOLUME_DEPTH = BLOCK_SIZE;
    static const size_t  BASIS_BLOCKS_WIDE = 8;
    static const size_t  BASIS_BLOCKS_TALL = 8;
    static const size_t  MAX_DCT_COEFF = 256;

    GlNddiDisplay  *display_;
    bool            quiet_;
    int             zigZag_[BLOCK_WIDTH * BLOCK_HEIGHT];
    uint8_t         quantizationMatrix_[BLOCK_WIDTH * BLOCK_HEIGHT];
    Pixel          *basisFunctions_;

    uint32_t        displayTilesWide_, displayTilesHigh_;
    uint8_t        *tileStackHeights_;

    vector< vector< vector< vector<uint64_t> > > >  cachedCoefficients_;
};
#endif // DCT_TILER_H
