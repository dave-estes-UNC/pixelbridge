#include <iostream>
#include <cassert>

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

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
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

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
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

void BaseNddiDisplay::PutPixel(Pixel p, vector<unsigned int> &location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + frameVolumeDimensionalSizes_.size()), 0);
    
    // Set the single pixel
	frameVolume_->PutPixel(p, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) {
	
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

void BaseNddiDisplay::CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) {
	
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

void BaseNddiDisplay::CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) {

	size_t tile_count = p.size();

	// Ensure parameter vectors' sizes match
	assert(starts.size() == tile_count);
	assert(size.size() == 2);

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < size.size(); i++) {
        pixelsToCopy *= size[i];
    }
    pixelsToCopy *= starts.size();
    costModel->registerTransmissionCharge(4 * (pixelsToCopy + (tile_count + 1) * frameVolumeDimensionalSizes_.size()), 0);

    // Copy pixels
    vector<unsigned int> end;
	end.resize(size.size());
    for (int i = 0; i < starts.size(); i++) {
    	assert(starts[i].size() == frameVolumeDimensionalSizes_.size());
		end[0] = starts[i][0] + size[0]- 1; if (end[0] >= frameVolumeDimensionalSizes_[0]) end[0] = frameVolumeDimensionalSizes_[0] - 1;
		end[1] = starts[i][1] + size[1] - 1; if (end[1] >= frameVolumeDimensionalSizes_[1]) end[1] = frameVolumeDimensionalSizes_[1] - 1;
    	for (int j = 2; j < frameVolumeDimensionalSizes_.size(); j++) {
    		end[j] = starts[i][j];
    	}
    	frameVolume_->CopyPixels(p[i], starts[i], end);
    }

    #ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + 2 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Fill pixels
    frameVolume_->FillPixel(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 * frameVolumeDimensionalSizes_.size()), 0);
    
    // Copy pixels
    frameVolume_->CopyFrameVolume(start, end, dest);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::UpdateInputVector(vector<int> &input) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * input.size(), 0);
	
    // Update the input vector
	inputVector_->UpdateInputVector(input);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                           vector<unsigned int> &location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + frameVolumeDimensionalSizes_.size()), 0);
    
    // Update the coefficient matrix
    coefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                            vector<unsigned int> &start,
                                            vector<unsigned int> &end) {
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
                                      vector<unsigned int> &start,
                                      vector<unsigned int> &end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 + 2 * 2), 0);
	
    // Fill the coefficient matrices
    coefficientPlane_->FillCoefficient(coefficient, row, col, start, end);

#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficientTiles(vector<int> &coefficients,
		                                   vector<vector<unsigned int> > &positions,
		                                   vector<vector<unsigned int> > &starts,
		                                   vector<unsigned int> &size) {

	size_t tile_count = coefficients.size();

	// Ensure parameter vectors' sizes match
	assert(positions.size() == tile_count);
	assert(starts.size() == tile_count);
	assert(size.size() == 2);

    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * ((1 + 3 + 2) * tile_count + 2), 0);

    // Fill the coefficient matrices
    vector<unsigned int> end;
    end.push_back(0); end.push_back(0); end.push_back(0);
    for (size_t i = 0; i < tile_count; i++) {
    	assert(positions[i].size() == 2);
    	assert(starts[i].size() == 3);
    	end[0] = starts[i][0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
    	end[1] = starts[i][1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
    	coefficientPlane_->FillCoefficient(coefficients[i], positions[i][0], positions[i][1], starts[i], end);
    }

#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillScaler(int scaler,
                                 vector<unsigned int> &start,
                                 vector<unsigned int> &end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + 3 * 2), 0);

    // Fill the coefficient matrices
    coefficientPlane_->FillScaler(scaler, start, end);

#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void BaseNddiDisplay::FillScalerTiles(vector<int> &scalers,
                                      vector<vector<unsigned int> > &starts,
                                      vector<unsigned int> &size) {
	// TODO(CDE): Implement
}

CostModel* BaseNddiDisplay::GetCostModel() {
    return costModel;
}
