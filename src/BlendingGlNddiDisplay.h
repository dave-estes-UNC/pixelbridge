//
//  BlendingGlNddiDisplay.h
//  pixelbridge
//
//  Created by Dave Estes on 1/31/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_BlendingGlNddiDisplay_h
#define pixelbridge_BlendingGlNddiDisplay_h

#include "GlNddiDisplay.h"
#include "NDimensionalDisplayInterfaceExtended.h"

using namespace nddi;

/**
 * Blending version of a GL NDDI Display.
 */
class BlendingGlNddiDisplay : public GlNddiDisplay, public NDimensionalDisplayInterfaceExtended {
    
public:
    BlendingGlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                          int inputVectorSize);
    BlendingGlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                          int displayWidth, int displayHeight,
                          int inputVectorSize);
    BlendingGlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                          int displayWidth, int displayHeight,
                          int inputVectorSize,
                          unsigned int planes);
    ~BlendingGlNddiDisplay();
    
    // Overridden to handle multiple Coefficient Planes
    void PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> location);
    void FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> start, std::vector<unsigned int> end);
    void FillCoefficient(int coefficient, int row, int col, std::vector<unsigned int> start, std::vector<unsigned int> end);
    
    // To satisfy the NDimensionalDisplayInterfaceExtended interface
    void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest, bool blend);
    void CopyPixelTiles(Pixel* p, std::vector<std::vector<unsigned int> > starts, std::vector<unsigned int> size) {};

    nddi::Pixel* GetFrameBuffer();
    
private:
    nddi::Pixel BlendPixel(nddi::Pixel pTo, nddi::Pixel pFrom);
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y);
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y, int* iv, nddi::Pixel* fv);
    void Render();
   
    int numPlanes_;
    CoefficientPlane** coefficientPlanes_;

};

#endif
