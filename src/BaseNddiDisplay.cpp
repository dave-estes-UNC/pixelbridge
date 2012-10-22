#include <iostream>

#include "PixelBridgeFeatures.h"
#include "BaseNddiDisplay.h"

using namespace nddi;

// public

BaseNddiDisplay::BaseNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                                 int inputVectorSize) {}
BaseNddiDisplay::BaseNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                                 int displayWidth, int displayHeight,
                                 int inputVectorSize) {}
BaseNddiDisplay::~BaseNddiDisplay() {}

int BaseNddiDisplay::DisplayWidth() {
	return displayWidth_;
}

int BaseNddiDisplay::DisplayHeight() {
	return displayHeight_;
}

void BaseNddiDisplay::PutPixel(Pixel p, std::vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + frameVolumeDimensionalSizes_.size()), 0);
    
    // Set the single pixel
	frameVolume_->PutPixel(p, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
	int dimensionToCopyAlong;
	bool dimensionFound = false;
	
    // Find the dimension to copy along
	for (int i = 0; !dimensionFound && (i < frameVolumeDimensionalSizes_.size()); i++) {
		if (start[i] != end[i]) {
			dimensionToCopyAlong = i;
			dimensionFound = true;
		}
	}
	
    // Register transmission cost now that we know the length of the strip sent
    costModel->registerTransmissionCharge(4 * ((end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1) + 2 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Copy the pixels
    frameVolume_->CopyPixelStrip(p, start, end);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < start.size(); i++) {
        pixelsToCopy *= end[i] - start[i] + 1;
    }
    costModel->registerTransmissionCharge(4 * (pixelsToCopy + 2 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Copy pixels
    frameVolume_->CopyPixels(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillPixel(Pixel p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + 2 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Fill pixels
    frameVolume_->FillPixel(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Copy pixels
    frameVolume_->CopyFrameVolume(start, end, dest);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::UpdateInputVector(std::vector<int> input) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * input.size(), 0);
	
    // Update the input vector
	inputVector_->UpdateInputVector(input);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix,
                                           std::vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + frameVolumeDimensionalSizes_.size()), 0);
    
    // Update the coefficient matrix
    coefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix,
                                            std::vector<unsigned int> start,
                                            std::vector<unsigned int> end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + 2 * 2), 0);
    
    // Fill the coefficient matrices
    coefficientPlane_->FillCoefficientMatrix(coefficientMatrix, start, end);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficient(int coefficient,
                                      int row, int col,
                                      std::vector<unsigned int> start,
                                      std::vector<unsigned int> end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 + 2 * 2), 0);
	
    // Fill the coefficient matrices
    coefficientPlane_->FillCoefficient(coefficient, row, col, start, end);

#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

CostModel* BaseNddiDisplay::GetCostModel() {
    return costModel;
}
