//
//  BlendingGlNddiDisplay.cpp
//  pixelbridge
//
//  Created by Dave Estes on 1/31/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>

#include "PixelBridgeFeatures.h"
#include "BlendingGlNddiDisplay.h"

#define CLAMPTOBYTE(color)                  \
if((color) & (~255)) {                      \
color = (unsigned char)((-(color)) >> 31);  \
} else {                                    \
color = (unsigned char)(color);             \
}

using namespace nddi;

// public

BlendingGlNddiDisplay::BlendingGlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                             int inputVectorSize)
: GlNddiDisplay(frameVolumeDimensionalSizes, inputVectorSize)
{
    numPlanes_ = 1;
    coefficientPlanes_ = &coefficientPlane_;
    costModel = new CostModel();
}

BlendingGlNddiDisplay::BlendingGlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                             int displayWidth, int displayHeight,
                                             int inputVectorSize)
: GlNddiDisplay(frameVolumeDimensionalSizes, displayWidth, displayHeight, inputVectorSize)
{
    numPlanes_ = 1;
    coefficientPlanes_ = &coefficientPlane_;
    costModel = new CostModel();
}

BlendingGlNddiDisplay::BlendingGlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                             int displayWidth, int displayHeight,
                                             int inputVectorSize,
                                             unsigned int planes)
: GlNddiDisplay(frameVolumeDimensionalSizes, displayWidth, displayHeight, inputVectorSize)
{
    numPlanes_ = planes;
    
    // Allocate the coefficient planes and copy over the initial one created by parent contructor
    coefficientPlanes_ = (CoefficientPlane**)malloc(sizeof(CoefficientPlane*) * planes);
    coefficientPlanes_[0] = coefficientPlane_;

	// Expand coefficient planes with zeroed coefficient matrices
    for (int i = 1; i < numPlanes_; i++) {
        coefficientPlanes_[i] = new CoefficientPlane(costModel, displayWidth, displayHeight, CM_WIDTH, CM_HEIGHT);
    }

    // Create the CostModel
    costModel = new CostModel();
}

void BlendingGlNddiDisplay::PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) {
	
    if (numPlanes_ > 1) {
        assert(location[2] <= numPlanes_);
    } else {
        if (location.size() == 2) {
            location.push_back(0);
        }
    }
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1), // One Coefficient Plane Coordinate triple
                                          0);
    
    // Update the coefficient matrix
    coefficientPlanes_[location[2]]->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BlendingGlNddiDisplay::FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &start, vector<unsigned int> &end) {

    if (numPlanes_ > 1) {
        assert(start[2] <= numPlanes_);
        assert(end[2] <= numPlanes_);
    } else {
        if (start.size() == 2) {
            start.push_back(0);
        }
        if (end.size() == 2) {
            end.push_back(0);
        }
    }
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);
    
    // For each plane in the range, fill the coefficient matrices
    for (int i = start[2]; i <= end[2]; i++) {
        coefficientPlanes_[i]->FillCoefficientMatrix(coefficientMatrix, start, end);
    }
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BlendingGlNddiDisplay::FillCoefficient(int coefficient, int row, int col, vector<unsigned int> &start, vector<unsigned int> &end) {

    if (numPlanes_ > 1) {
        assert(start[2] <= numPlanes_);
        assert(end[2] <= numPlanes_);
    } else {
        if (start.size() == 2) {
            start.push_back(0);
        }
        if (end.size() == 2) {
            end.push_back(0);
        }
    }
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * 1 +                // One coefficient
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(1) + // One Coefficient Matrix Coordinate double
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2),  // Two Coefficient Plane Coordinate triples
                                          0);
    
    // For each plane in the range, fill the coefficient matrices
    for (int i = start[2]; i <= end[2]; i++) {
        coefficientPlanes_[i]->FillCoefficient(coefficient, row, col, start, end);
    }
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

// Implements the frame volume function directly using getPixel and setPixel since frame volume does not yet
// support blending.
void BlendingGlNddiDisplay::CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest, bool blend) {

	vector<unsigned int> positionFrom = start;
	vector<unsigned int> positionTo = dest;
	bool copyFinished = false;
	int pixelsCopied = 0;
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_FV_COORD_TUPLES(3), // Three Coordinate Tuples
                                          0);
    
	// Move from start to end, filling in each location with the provided pixel
	do {
		// Set pixel in frame volume at position
        if (!blend) {
            frameVolume_->setPixel(positionTo, frameVolume_->getPixel(positionFrom));
        } else {
            // Register the cost of the blending
            costModel->registerPixelBlendCharge(1);
            
            frameVolume_->setPixel(positionTo, BlendPixel(frameVolume_->getPixel(positionTo), frameVolume_->getPixel(positionFrom)));
        }
		pixelsCopied++;
		
		// Move to the next position
		int fvDim = 0;
		bool overflow;
		do {
			overflow = false;
			positionFrom[fvDim]++;
			positionTo[fvDim]++;
			if ( (positionFrom[fvDim] >= frameVolumeDimensionalSizes_[fvDim])
				|| (positionFrom[fvDim] > end[fvDim]) ) {
				overflow = true;
				positionFrom[fvDim] = start[fvDim];
				positionTo[fvDim] = dest[fvDim];
				if (++fvDim >= frameVolumeDimensionalSizes_.size())
					copyFinished = true;
			}
		} while (overflow && !copyFinished);
		
	} while (!copyFinished);
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

// private

void BlendingGlNddiDisplay::Render() {
    
    timeval startTime, endTime; // Used for timing data
    if (!quiet_)
    	gettimeofday(&startTime, NULL);
    
    // Even though the InputVector can be invoked concurrently, it's really slow. So
    // we'll use a local copy instead and update the cost model in bulk later.
#ifndef NO_OMP
    int   * iv = inputVector_->data();
    Pixel * fv = frameVolume_->data();
#pragma omp parallel for
#endif // !NO_OMP
    for (int j = 0; j < displayHeight_; j++) {
        for (int i = 0; i < displayWidth_; i++) {
#ifndef NO_OMP
			frameBuffer_[j * displayWidth_ + i].packed = ComputePixel(i, j, iv, fv).packed;
#else
            frameBuffer_[j * displayWidth_ + i].packed = this->ComputePixel(i, j).packed;
#endif
        }   
    }
    
    // Update the cost model for the in bulk now if we are using OpenMP since we bypassed the traditional
    // getters for input vector, frame volume, and coefficient matrix.
#ifndef NO_OMP
    costModel->registerBulkMemoryCharge(INPUT_VECTOR_COMPONENT,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * (CM_WIDTH - 2) * numPlanes_,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * (CM_WIDTH - 2) * numPlanes_ * BYTES_PER_IV_VALUE,
                                        0);
    costModel->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * CM_WIDTH * numPlanes_,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * CM_WIDTH * numPlanes_ * BYTES_PER_COEFF,
                                        0);
    costModel->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                        displayWidth_ * displayHeight_ * numPlanes_,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * numPlanes_ * BYTES_PER_PIXEL,
                                        0);
    costModel->registerPixelMappingCharge(displayWidth_ * displayHeight_ * numPlanes_);
    if (numPlanes_ > 1) {
        costModel->registerPixelBlendCharge(displayWidth_ * displayHeight_ * (numPlanes_ - 1));
    }
#endif
    
    if (!quiet_) {
		gettimeofday(&endTime, NULL);
		printf("Render Statistics:\n  Size: %dx%d - FPS: %f\n",
			   displayWidth_,
			   displayHeight_,
			   1.0f / ((double)(endTime.tv_sec * 1000000
								+ endTime.tv_usec
								- startTime.tv_sec * 1000000
								- startTime.tv_usec) / 1000000.0f)
			   );
    }
}


// Don't register the blending cost here, since we want to bypass it sometimes for speed.
nddi::Pixel BlendingGlNddiDisplay::BlendPixel(nddi::Pixel pTo, nddi::Pixel pFrom) {
    int r, g, b;

    r = ((int)pFrom.r * (int)pFrom.a + (int)pTo.r * (255 - (int)pFrom.a)) / 255;
    g = ((int)pFrom.g * (int)pFrom.a + (int)pTo.g * (255 - (int)pFrom.a)) / 255;
    b = ((int)pFrom.b * (int)pFrom.a + (int)pTo.b * (255 - (int)pFrom.a)) / 255;

    CLAMPTOBYTE(r);
    CLAMPTOBYTE(g);
    CLAMPTOBYTE(b);

    pTo.r = (unsigned char)r;
    pTo.g = (unsigned char)g;
    pTo.b = (unsigned char)b;

    return pTo;
}


nddi::Pixel BlendingGlNddiDisplay::ComputePixel(unsigned int x, unsigned int y) {
    
    nddi::Pixel computedPixel;
    
    // For each coefficient plan, calculate the pixel value by blending the other layers on top
    for (int p = 0; p < numPlanes_; p++) {
    
        // Grab the coefficient matrix
        CoefficientMatrix * matrix = coefficientPlanes_[p]->getCoefficientMatrix(x, y, 0);
        
        // Compute the position vector for the proper pixel in the frame volume.
        vector<unsigned int> fvPosition;
        // Matrix multiply the input vector by the coefficient matrix
        for (int j = 0; j < CM_HEIGHT; j++) {
            fvPosition.push_back(0);
            // No need to read the x and y from the input vector, just multiply directly.
            fvPosition[j] += matrix->getCoefficient(0, j) * x;
            fvPosition[j] += matrix->getCoefficient(1, j) * y;
            // Then multiply the remainder of the input vector
            for (int i = 2; i < CM_WIDTH; i++) {
                fvPosition[j] += matrix->getCoefficient(i, j) * inputVector_->getValue(i);
            }
        }
        
        // Register the cost of the mapping
        costModel->registerPixelMappingCharge(1);

        // Calculate the index into the frameVolume_ int array
        if (p == 0) {
            computedPixel = frameVolume_->getPixel(fvPosition);
        } else {
            // Register the cost of the blending
            costModel->registerPixelBlendCharge(1);
            
            computedPixel = BlendPixel(computedPixel, frameVolume_->getPixel(fvPosition));
        }
    }
    
    return computedPixel;
}

Pixel BlendingGlNddiDisplay::ComputePixel(unsigned int x, unsigned int y, int* iv, Pixel* fv) {
	
    nddi::Pixel computedPixel;
    
    // For each coefficient plan, calculate the pixel value by blending the other layers on top
    for (int p = 0; p < numPlanes_; p++) {
        
        // Grab the coefficient matrix
        int * cm = coefficientPlanes_[p]->getCoefficientMatrix(x, y, 0)->data();
        
        // Compute the position vector for the proper pixel in the frame volume.
        vector<unsigned int> fvPosition;
        // Matrix multiply the input vector by the coefficient matrix
        for (int j = 0; j < CM_HEIGHT; j++) {
            // Initialize to zero
            fvPosition.push_back(0);
            // No need to read the x and y from the input vector, just multiply directly.
            fvPosition[j] += cm[j * CM_WIDTH + 0] * x;
            fvPosition[j] += cm[j * CM_WIDTH + 1] * y;
            // Then multiply the remainder of the input vector
            for (int i = 2; i < CM_WIDTH; i++) {
                fvPosition[j] += cm[j * CM_WIDTH + i] * iv[i];
            }
        }
        
        // Compute the offset and grab the pixel directly from the frame volume
        unsigned int offset = 0;
        unsigned int multiplier = 1;
        
        for (int i = 0; i < fvPosition.size(); i++) {
            offset += fvPosition[i] * multiplier;
            multiplier *= frameVolumeDimensionalSizes_[i];
        }
        
        // Either grab the pixel or blend
        if (p == 0) {
            computedPixel = fv[offset];
        } else {
            computedPixel = BlendPixel(computedPixel, fv[offset]);
        }
    }
   
    return computedPixel;
}

nddi::Pixel* BlendingGlNddiDisplay::GetFrameBuffer() {
#ifdef SUPRESS_EXCESS_RENDERING
	Render();
#endif
	return frameBuffer_;
}
