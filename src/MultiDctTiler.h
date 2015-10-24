#ifndef MULTI_DCT_TILER_H
#define MULTI_DCT_TILER_H
/*
 *  ScaledDctTiler.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 11/26/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include "ScaledDctTiler.h"
#include "BaseNddiDisplay.h"

using namespace nddi;
using namespace std;

/**
 * This tiler will split provided frames into macroblocks and will perform the forward DCT
 * and create coefficients that will be used as the scalers for the coefficient planes.
 */
class MultiDctTiler : public ScaledDctTiler {

public:
    /**
     * The MultiDctTiler is created based on the dimensions of the NDDI display that's passed in. If those
     * dimensions change, then the MultiDctTiler should be destroyed and re-created.
     *
     * @param display_width The width of the display
     * @param display_height The height of the display
     * @param quality The quality factor used for DCT.
     * @param quiet Used to squelch extra information output.
     */
        MultiDctTiler(size_t display_width, size_t display_height, size_t quality);

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
    void initZigZag();
    void initQuantizationMatrix(size_t quality);
    void InitializeCoefficientPlanes();
    void InitializeFrameVolume();
    vector<uint64_t> BuildCoefficients(size_t i, size_t j, int16_t* buffer, size_t width, size_t height, size_t c, bool adjustPixels);
    void SelectCoefficientsForScale(vector<uint64_t> &coefficients, size_t c);
    void FillCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, size_t first);
    void PrerenderCoefficients(vector<uint64_t> &coefficients, size_t i, size_t j, size_t c, int16_t* renderedBuffer, size_t width, size_t height, bool shift);

private:
    vector< vector< vector< vector<uint64_t> > > >  cachedCoefficients_;
    vector< vector <int> >                          zigZag_;
    vector< vector<uint8_t> >                       quantizationMatrix_;
    size_t                                          fvWidth_, fvHeight_;

};
#endif // MULTI_DCT_TILER_H
