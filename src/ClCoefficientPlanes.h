//
//  ClCoefficientPlanes.h
//  pixelbridge
//
//  Created by Dave Estes on 4/1/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_ClCoefficientPlanes_h
#define pixelbridge_ClCoefficientPlanes_h

#include "CoefficientPlanes.h"

using namespace nddi;

class ClCoefficientPlanes : public CoefficientPlanes {

private:
    cl_mem            clCoeffBuffer_;
    cl_mem            clScalerBuffer_;
    cl_context        clContext_;
    cl_command_queue  clQueue_;

    unsigned int      matrixWidth_;
    unsigned int      matrixHeight_;
    unsigned int      matrixSize_;      // In bytes
    unsigned int      coefficientSize_; // In bytes
    unsigned int      scalerSize_;      // In bytes
    unsigned int      numPlanes_;

#ifdef NARROW_DATA_STORES
    int16_t             * coefficients_;
    int16_t             * scalers_;
#else
    int                 * coefficients_;
    int                 * scalers_;
#endif

    inline unsigned int calcOffset(vector<unsigned int> &location) {
        assert(location.size() == 3);
        return (location[2] * height_ * width_ + location[1] * width_ + location[0]) * matrixWidth_ * matrixHeight_;
    }

public:

    ClCoefficientPlanes(CostModel* costModel,
                       unsigned int displayWidth, unsigned int displayHeight,
                       unsigned int numPlanes,
                       unsigned int matrixWidth, unsigned int matrixHeight) {

        costModel_ = costModel;
        width_ = displayWidth;
        height_ = displayHeight;
        numPlanes_ = numPlanes;
        coefficientMatrices_ = NULL;

        matrixWidth_ = matrixWidth;
        matrixHeight_ = matrixHeight;
#ifdef NARROW_DATA_STORES
        coefficientSize_ = sizeof(int16_t);
        scalerSize_ = sizeof(int16_t) * 3;
#else
        coefficientSize_ = sizeof(int);
        scalerSize_ = sizeof(int) * 3;
#endif
        matrixSize_ = matrixWidth_ * matrixHeight_ * coefficientSize_;

#ifdef NARROW_DATA_STORES
        coefficients_ = (int16_t *)malloc(CoefficientMatrix::memoryRequired(matrixWidth, matrixHeight) * displayWidth * displayHeight * numPlanes_);
        scalers_ = (int16_t *)malloc(scalerSize_ * displayWidth * displayHeight * numPlanes_);
#else
        coefficients_ = (int *)malloc(CoefficientMatrix::memoryRequired(matrixWidth, matrixHeight) * displayWidth * displayHeight * numPlanes_);
        scalers_ = (int *)malloc(scalerSize_ * displayWidth * displayHeight * numPlanes_);
#endif
    }

    ~ClCoefficientPlanes() {

        if (coefficients_)
            free(coefficients_);

        // Release CL mem buffers
        if (clCoeffBuffer_ != 0)
            clReleaseMemObject(clCoeffBuffer_);
        if (clScalerBuffer_ != 0)
            clReleaseMemObject(clScalerBuffer_);
    }

    void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) {

        assert(location.size() == 3);
        assert(location[0] < width_);
        assert(location[1] < height_);
        assert(location[2] < numPlanes_);

        unsigned int offset = calcOffset(location);
#ifdef NARROW_DATA_STORES
        int16_t * coefficientPtr = coefficients_ + offset;
#else
        int * coefficientPtr = coefficients_ + offset;
#endif

        for (int y = 0; y < matrixHeight_; y++) {
            for (int x = 0; x < matrixWidth_; x++) {
                if (coefficientMatrix[x][y] != COFFICIENT_UNCHANGED) {
                    *coefficientPtr = coefficientMatrix[x][y];
                }
                assert(0); // TODO(CDE): Shouldn't coefficientPtr be moved to the next coefficient?
            }
        }

        // Enqueue CL commands
        int err = clEnqueueWriteBuffer(clQueue_, clCoeffBuffer_, CL_TRUE,
                                       offset * coefficientSize_,
                                       matrixSize_,
                                       coefficients_,
                                       0, NULL, NULL);
        if (err != CL_SUCCESS) {
            cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << endl;
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + offset, matrixSize_, 0);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + offset, matrixSize_, 0);
    }

    void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                               vector<unsigned int> &start,
                               vector<unsigned int> &end) {

        assert(start.size() == 3);
        assert(start[0] < width_);
        assert(start[1] < height_);
        assert(start[2] < numPlanes_);
        assert(end.size() == 3);
        assert(end[0] < width_);
        assert(end[1] < height_);
        assert(end[2] < numPlanes_);

        vector<unsigned int> position = start;
        unsigned int stride = end[0] - start[0] + 1;
        unsigned int rowCount = end[1] - start[1] + 1;
        unsigned int planeCount = end[2] - start[2] + 1;

        // Copy matrices into coefficients_ and enqueue as one large write
        unsigned int offset;
        for (position[2] = start[2]; position[2] <= end[2]; position[2]++) {
            for (position[1] = start[1]; position[1] <= end[1]; position[1]++) {
                for (position[0] = start[0]; position[0] <= end[0]; position[0]++) {
                    offset = calcOffset(position);
                    for (int y = 0; y < matrixHeight_; y++) {
                        for (int x = 0; x < matrixWidth_; x++) {
                            if (coefficientMatrix[x][y] != COFFICIENT_UNCHANGED) {
                                coefficients_[offset] = coefficientMatrix[x][y];
                            }
                            offset++;
                        }
                    }
                }
            }
        }

        // Copy the rectangular region of coefficients to the compute device
        size_t row_pitch = matrixSize_ * width_;
        size_t slice_pitch = matrixSize_ * width_ * height_;
        const size_t origin[3] = {start[0] * matrixSize_, start[1], start[2]};
        const size_t region[3] = {(end[0] - start[0] + 1) * matrixSize_, end[1] - start[1] + 1, end[2] - start[2] + 1};

        int err = clEnqueueWriteBufferRect(clQueue_, clCoeffBuffer_, CL_FALSE,
                                           origin, origin, region,
                                           row_pitch, slice_pitch, row_pitch, slice_pitch,
                                           coefficients_, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            cout << __FUNCTION__ << " - Err: (" << err << ") - Failed to create enqueue write buffer rect command." << endl;
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + calcOffset(start), matrixSize_* stride * rowCount * planeCount, 0);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + calcOffset(start), matrixSize_* stride * rowCount * planeCount, 0);
    }

    void FillCoefficient(int coefficient,
                         int row, int col,
                         vector<unsigned int> &start,
                         vector<unsigned int> &end) {

        assert(start.size() == 3);
        assert(start[0] < width_);
        assert(start[1] < height_);
        assert(start[2] < numPlanes_);
        assert(end.size() == 3);
        assert(end[0] < width_);
        assert(end[1] < height_);
        assert(end[2] < numPlanes_);

        vector<unsigned int> position = start;
        unsigned int stride = end[0] - start[0] + 1;
        unsigned int rowCount = end[1] - start[1] + 1;
        unsigned int planeCount = end[2] - start[2] + 1;

        // Copy matrices into coefficients_ and enqueue as one large write
        unsigned int offset;
        for (position[2] = start[2]; position[2] <= end[2]; position[2]++) {
            for (position[1] = start[1]; position[1] <= end[1]; position[1]++) {
                for (position[0] = start[0]; position[0] <= end[0]; position[0]++) {
                    offset = calcOffset(position);
                    offset += row * matrixWidth_ + col;
                    coefficients_[offset] = coefficient;
                }
            }
        }

        // Copy the rectangular region of coefficients to the compute device
        size_t row_pitch = matrixSize_ * width_;
        size_t slice_pitch = matrixSize_ * width_ * height_;
        const size_t origin[3] = {start[0] * matrixSize_, start[1], start[2]};
        const size_t region[3] = {(end[0] - start[0] + 1) * matrixSize_, end[1] - start[1] + 1, end[2] - start[2] + 1};

        int err = clEnqueueWriteBufferRect(clQueue_, clCoeffBuffer_, CL_FALSE,
                                           origin, origin, region,
                                           row_pitch, slice_pitch, row_pitch, slice_pitch,
                                           coefficients_, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            cout << __FUNCTION__ <<  " - Err: (" << err << ") - Failed to create enqueue write buffer rect command." << endl;
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + calcOffset(start), coefficientSize_* stride * rowCount * planeCount, 0);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + calcOffset(start), coefficientSize_* stride * rowCount * planeCount, 0);
    }

    void FillScaler(Scaler scaler,
                    vector<unsigned int> &start,
                    vector<unsigned int> &end) {

        assert(start.size() == 3);
        assert(start[0] < width_);
        assert(start[1] < height_);
        assert(start[2] < numPlanes_);
        assert(end.size() == 3);
        assert(end[0] < width_);
        assert(end[1] < height_);
        assert(end[2] < numPlanes_);

        vector<unsigned int> position = start;
        unsigned int stride = end[0] - start[0] + 1;
        unsigned int rowCount = end[1] - start[1] + 1;
        unsigned int planeCount = end[2] - start[2] + 1;

        // Copy matrices into coefficients_ and enqueue as one large write
        unsigned int offset;
        for (size_t p = start[2]; p <= end[2]; p++) {
            for (size_t y = start[1]; y <= end[1]; y++) {
                for (size_t x = start[0]; x <= end[0]; x++) {
                    offset = p * height_ * width_ + y * width_ + x;
                    scalers_[offset] = scaler.r; offset++;
                    scalers_[offset] = scaler.g; offset++;
                    scalers_[offset] = scaler.b; offset++;
                }
            }
        }

        // Copy the rectangular region of scalers to the compute device
        size_t row_pitch = matrixSize_ * width_;
        size_t slice_pitch = matrixSize_ * width_ * height_;
        const size_t origin[3] = {start[0] * scalerSize_, start[1], start[2]};
        const size_t region[3] = {(end[0] - start[0] + 1) * scalerSize_, end[1] - start[1] + 1, end[2] - start[2] + 1};

        int err = clEnqueueWriteBufferRect(clQueue_, clScalerBuffer_, CL_FALSE,
                                           origin, origin, region,
                                           row_pitch, slice_pitch, row_pitch, slice_pitch,
                                           scalers_, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            cout << __FUNCTION__ <<  " - Err: (" << err << ") - Failed to create enqueue write buffer rect command." << endl;
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, scalers_, scalerSize_* stride * rowCount * planeCount, 0);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, scalers_, scalerSize_* stride * rowCount * planeCount, 0);
    }

    cl_mem initializeCl(cl_context context, cl_command_queue queue) {

        // Set CL variables
        clContext_ = context;
        clQueue_ = queue;

        // Create CL mem buffers
        clCoeffBuffer_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY,
                                   matrixSize_ * width_ * height_ * numPlanes_, NULL, NULL);
        clScalerBuffer_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY,
                                   scalerSize_ * width_ * height_ * numPlanes_, NULL, NULL);

        return clCoeffBuffer_ ;
    }

    // TODO(CDE): Remove this garbage
    cl_mem getClCoeffBuffer() {
        return clCoeffBuffer_;
    }

    // TODO(CDE): Remove this garbage
    cl_mem getClScalerBuffer() {
        return clScalerBuffer_;
    }

};

#endif
