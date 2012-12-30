#include <iostream>

#include "PixelBridgeFeatures.h"
#include "BaseNddiDisplay.h"

using namespace nddi;

// public

BaseNddiDisplay::BaseNddiDisplay() :
		displayWidth_(0),
		displayHeight_(0),
		inputVector_(NULL),
		frameVolume_(NULL),
		coefficientPlane_(NULL),
		frameBuffer_(NULL),
		costModel(NULL),
		quiet_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> frameVolumeDimensionalSizes,
                                 int inputVectorSize) :
		displayWidth_(0),
		displayHeight_(0),
		inputVector_(NULL),
		frameVolumeDimensionalSizes_(frameVolumeDimensionalSizes),
		frameVolume_(NULL),
		coefficientPlane_(NULL),
		frameBuffer_(NULL),
		costModel(NULL),
		quiet_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> frameVolumeDimensionalSizes,
                                 int displayWidth, int displayHeight,
                                 int inputVectorSize) :
		displayWidth_(displayWidth),
		displayHeight_(displayHeight),
		inputVector_(NULL),
		frameVolumeDimensionalSizes_(frameVolumeDimensionalSizes),
		frameVolume_(NULL),
		coefficientPlane_(NULL),
		frameBuffer_(NULL),
		costModel(NULL),
		quiet_(false)
{}

BaseNddiDisplay::~BaseNddiDisplay() {}

int BaseNddiDisplay::DisplayWidth() {
	return displayWidth_;
}

int BaseNddiDisplay::DisplayHeight() {
	return displayHeight_;
}

void BaseNddiDisplay::PutPixel(Pixel p, vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + frameVolumeDimensionalSizes_.size()), 0);
    
    // Set the single pixel
	frameVolume_->PutPixel(p, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyPixelStrip(Pixel* p, vector<unsigned int> start, vector<unsigned int> end) {
	
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

void BaseNddiDisplay::CopyPixels(Pixel* p, vector<unsigned int> start, vector<unsigned int> end) {
	
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

void BaseNddiDisplay::FillPixel(Pixel p, vector<unsigned int> start, vector<unsigned int> end) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + 2 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Fill pixels
    frameVolume_->FillPixel(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyFrameVolume(vector<unsigned int> start, vector<unsigned int> end, vector<unsigned int> dest) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Copy pixels
    frameVolume_->CopyFrameVolume(start, end, dest);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::UpdateInputVector(vector<int> input) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * input.size(), 0);
	
    // Update the input vector
	inputVector_->UpdateInputVector(input);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::PutCoefficientMatrix(vector< vector<int> > coefficientMatrix,
                                           vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + frameVolumeDimensionalSizes_.size()), 0);
    
    // Update the coefficient matrix
    coefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficientMatrix(vector< vector<int> > coefficientMatrix,
                                            vector<unsigned int> start,
                                            vector<unsigned int> end) {
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
                                      vector<unsigned int> start,
                                      vector<unsigned int> end) {
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
