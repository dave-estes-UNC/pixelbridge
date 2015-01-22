#include <iostream>
#include <cassert>
#include <cmath>

#include "PixelBridgeFeatures.h"
#include "BaseNddiDisplay.h"

using namespace nddi;

// public

BaseNddiDisplay::BaseNddiDisplay() :
        displayWidth_(0),
        displayHeight_(0),
        inputVector_(NULL),
        frameVolume_(NULL),
        coefficientPlanes_(NULL),
        costModel(NULL),
        quiet_(false),
        changed_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                 int numCoefficientPlanes, int inputVectorSize) :
        displayWidth_(0),
        displayHeight_(0),
        numPlanes_(numCoefficientPlanes),
        inputVector_(NULL),
        frameVolumeDimensionalSizes_(frameVolumeDimensionalSizes),
        frameVolume_(NULL),
        coefficientPlanes_(NULL),
        costModel(NULL),
        quiet_(false),
        changed_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                                 int displayWidth, int displayHeight,
                                 int numCoefficientPlanes, int inputVectorSize) :
        displayWidth_(displayWidth),
        displayHeight_(displayHeight),
        numPlanes_(numCoefficientPlanes),
        inputVector_(NULL),
        frameVolumeDimensionalSizes_(frameVolumeDimensionalSizes),
        frameVolume_(NULL),
        coefficientPlanes_(NULL),
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

int BaseNddiDisplay::NumCoefficientPlanes() {
    return numPlanes_;
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
    coefficientPlanes_->PutCoefficientMatrix(coefficientMatrix, location);

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
    coefficientPlanes_->FillCoefficientMatrix(coefficientMatrix, start, end);

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
    coefficientPlanes_->FillCoefficient(coefficient, row, col, start, end);

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
        coefficientPlanes_->FillCoefficient(coefficients[i], positions[i][0], positions[i][1], starts[i], end);
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScaler(Scaler scaler,
                                 vector<unsigned int> &start,
                                 vector<unsigned int> &end) {
    assert(start.size() == 3);
    assert(end.size() == 3);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * 1 +              // One Scaler
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    coefficientPlanes_->FillScaler(scaler, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScalerTiles(vector<uint64_t> &scalers,
                                      vector<vector<unsigned int> > &starts,
                                      vector<unsigned int> &size) {
    size_t tile_count = scalers.size();
    Scaler s;

    // Ensure parameter vectors' sizes match
    assert(starts.size() == tile_count);
    assert(size.size() == 2);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * tile_count +                // t scalers
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
        s.packed = scalers[i];
        coefficientPlanes_->FillScaler(s, starts[i], end);
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScalerTileStack(vector<uint64_t> &scalers,
                                          vector<unsigned int> &start,
                                          vector<unsigned int> &size) {
    vector<unsigned int> st = start;

    size_t stack_height = scalers.size();
    Scaler s;

    assert(start.size() == 3);
    assert(size.size() == 2);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * stack_height +     // h scalers
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1) +  // One Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1), // One X by Y tile dimension double
                                          0);

    vector<unsigned int> end;
    end.push_back(0); end.push_back(0); end.push_back(0);
    for (size_t i = 0; i < stack_height; i++) {
        end[0] = start[0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
        end[1] = start[1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
        end[2] = start[2];
        s.packed = scalers[i];
        coefficientPlanes_->FillScaler(s, st, end);
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

void BaseNddiDisplay::SetFullScaler(uint16_t scaler) {
    if (scaler & (scaler - 1)) {
        cout << "ERROR: THE FULL_SCALER specified is not a power of two." << endl;
    } else {
        // Set the full scaler
        fullScaler_ = scaler;

        // Initialize the shifter used during accumulation
        double s = log2((double)fullScaler_);
        accumulatorShifter_ = int(s);
    }
}

CostModel* BaseNddiDisplay::GetCostModel() {
    return costModel;
}
