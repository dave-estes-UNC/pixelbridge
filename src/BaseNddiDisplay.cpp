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
        pixelSignMode_(UNSIGNED_MODE),
		costModel(NULL),
		quiet_(false),
		changed_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                 int inputVectorSize) :
		displayWidth_(0),
		displayHeight_(0),
		inputVector_(NULL),
		frameVolumeDimensionalSizes_(frameVolumeDimensionalSizes),
		frameVolume_(NULL),
		coefficientPlane_(NULL),
        pixelSignMode_(UNSIGNED_MODE),
		costModel(NULL),
		quiet_(false),
		changed_(false)
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
        pixelSignMode_(UNSIGNED_MODE),
		costModel(NULL),
		quiet_(false),
		changed_(false)
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
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(1), // One Coordinate Tuple
                                          0);
    
    // Set the single pixel
	frameVolume_->PutPixel(p, location);
	
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
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
    int pixelsToCopy = end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1;
	
    // Register transmission cost now that we know the length of the strip sent
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // A strip of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);
    
    // Copy the pixels
    frameVolume_->CopyPixelStrip(p, start, end);
	
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) {
	
    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < start.size(); i++) {
        pixelsToCopy *= end[i] - start[i] + 1;
    }
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // Range of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);
    
    // Copy pixels
    frameVolume_->CopyPixels(p, start, end);
    
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) {

	size_t tile_count = p.size();

	// Ensure parameter vectors' sizes match
	assert(starts.size() == tile_count);
	assert(starts[0].size() == frameVolumeDimensionalSizes_.size());
	assert(size.size() == 2);

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < size.size(); i++) {
        pixelsToCopy *= size[i];
    }
    pixelsToCopy *= tile_count;
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +            // t tiles of x by y tiles of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(tile_count + 1) + // t start coordinate tuples + 1 tuple for tile size dimensions
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),            // 1 X by Y tile dimension double
                                          0);

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

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2), // Two Coordinate Tuples
                                          0);
    
    // Fill pixels
    frameVolume_->FillPixel(p, start, end);
    
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_FV_COORD_TUPLES(3), // Three Coordinate Tuples
                                          0);
    
    // Copy pixels
    frameVolume_->CopyFrameVolume(start, end, dest);
    
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::UpdateInputVector(vector<int> &input) {
    
    assert(input.size() == inputVector_->getSize() - 2);
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_IV_UPDATE(), // Input Vector
                                          0);
	
    // Update the input vector
	inputVector_->UpdateInputVector(input);
	
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                           vector<unsigned int> &location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1), // One Coefficient Plane Coordinate triple
                                          0);
    
    // Update the coefficient matrix
    coefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                            vector<unsigned int> &start,
                                            vector<unsigned int> &end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);
    
    // Fill the coefficient matrices
    coefficientPlane_->FillCoefficientMatrix(coefficientMatrix, start, end);
	
#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillCoefficient(int coefficient,
                                      int row, int col,
                                      vector<unsigned int> &start,
                                      vector<unsigned int> &end) {
    assert(row >= 0 && row < CM_HEIGHT);
    assert(col >= 0 && col < CM_WIDTH);
    assert(start.size() == 3);
    assert(end.size() == 3);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * 1 +                // One coefficient
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(1) + // One Coefficient Matrix Coordinate double
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2),  // Two Coefficient Plane Coordinate triples
                                          0);
	
    // Fill the coefficient matrices
    coefficientPlane_->FillCoefficient(coefficient, row, col, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
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
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * tile_count +                 // t coefficients
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(tile_count) +  // t Coefficient Matrix Coordinate doubles
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(tile_count) +  // t Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),          // 1 X by Y tile dimension double
                                          0);

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

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillScaler(int scaler,
                                 vector<unsigned int> &start,
                                 vector<unsigned int> &end) {
    assert(start.size() == 3);
    assert(end.size() == 3);
    
    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALAR * 1 +              // One Scalar
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    coefficientPlane_->FillScaler(scaler, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillScalerTiles(vector<int> &scalers,
                                      vector<vector<unsigned int> > &starts,
                                      vector<unsigned int> &size) {
	size_t tile_count = scalers.size();

	// Ensure parameter vectors' sizes match
	assert(starts.size() == tile_count);
	assert(size.size() == 2);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALAR * tile_count +                // t scalers
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(tile_count) +  // t Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),          // One X by Y tile dimension double
                                          0);

    vector<unsigned int> end;
    end.push_back(0); end.push_back(0); end.push_back(0);
    for (size_t i = 0; i < tile_count; i++) {
    	assert(starts[i].size() == 3);
    	end[0] = starts[i][0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
    	end[1] = starts[i][1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
    	end[2] = starts[i][2];
    	coefficientPlane_->FillScaler(scalers[i], starts[i], end);
    }

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::FillScalerTileStack(vector<int> &scalers,
                                          vector<unsigned int> &start,
                                          vector<unsigned int> &size) {
	vector<unsigned int> st = start;

	size_t tile_count = scalers.size();

	assert(start.size() == 3);

    // Register transmission cost first
	/*
    costModel->registerTransmissionCharge(BYTES_PER_SCALAR * tile_count +       // t scalers
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1) +  // One Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1), // One X by Y tile dimension double
                                          0);
                                          */

    vector<unsigned int> end;
    end.push_back(0); end.push_back(0); end.push_back(0);
    for (size_t i = 0; i < tile_count; i++) {
    	end[0] = start[0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
    	end[1] = start[1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
    	end[2] = start[2];
    	coefficientPlane_->FillScaler(scalers[i], st, end);
    	st[2]++;
    }

#ifdef SUPRESS_EXCESS_RENDERING
	changed_ = true;
#else
	Render();
#endif
}

void BaseNddiDisplay::SetPixelByteSignMode(SignMode mode) {
    assert(mode == UNSIGNED_MODE || mode == SIGNED_MODE);
    pixelSignMode_ = mode;
}

CostModel* BaseNddiDisplay::GetCostModel() {
    return costModel;
}
