//
//  ClCoefficientPlane.h
//  pixelbridge
//
//  Created by Dave Estes on 4/1/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_ClCoefficientPlane_h
#define pixelbridge_ClCoefficientPlane_h

#include "CoefficientPlane.h"

using namespace nddi;

class ClCoefficientPlane : public CoefficientPlane {
    
private:
    cl_mem            clBuffer_;
    cl_context        clContext_;
    cl_command_queue  clQueue_;
    
    unsigned int      matrixWidth_;
    unsigned int      matrixHeight_;
    unsigned int      matrixSize_;   // In bytes
    
    int             * coefficients_;
    
    inline unsigned int calcOffset(std::vector<unsigned int> location) {
        return (location[1] * width_ + location[0]) * matrixWidth_ * matrixHeight_;
    }

public:
    
    ClCoefficientPlane(CostModel* costModel,
                       unsigned int displayWidth, unsigned int displayHeight,
                       unsigned int matrixWidth, unsigned int matrixHeight) {
        
        costModel_ = costModel;
        width_ = displayWidth;
        height_ = displayHeight;
        coefficientMatrices_ = NULL;
            
        matrixWidth_ = matrixWidth;
        matrixHeight_ = matrixHeight;
        matrixSize_ = matrixWidth_ * matrixHeight_ * sizeof(int);
        
        coefficients_ = (int *)malloc(matrixSize_ * width_ * height_);
    }
    
    ~ClCoefficientPlane() {
        
        if (coefficients_)
            free(coefficients_);

        // Release CL mem buffer
        if (clBuffer_ != 0)
            clReleaseMemObject(clBuffer_);
    }
    
    void PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> location) {
        
        // TODO(CDE): Add future support for multiple planes
        assert(location.size() == 2);
        assert(location[0] < width_);
        assert(location[1] < height_);
        
        unsigned int offset = calcOffset(location);
        int * coefficientPtr = coefficients_ + offset;

        for (int y = 0; y < matrixHeight_; y++) {
            for (int x = 0; x < matrixWidth_; x++) {
                if (coefficientMatrix[x][y] != COFFICIENT_UNCHANGED) {
                    *coefficientPtr = coefficientMatrix[x][y];
                }
            }
        }

        // Enqueue CL commands
        int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, offset, matrixSize_, coefficients_, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + offset, matrixSize_);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + offset, matrixSize_);
    }
    
    void FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix,
                               std::vector<unsigned int> start,
                               std::vector<unsigned int> end) {
        
        // TODO(CDE): Add future support for multiple planes
        assert(start.size() == 2);
        assert(start[0] < width_);
        assert(start[1] < height_);
        assert(end.size() == 2);
        assert(end[0] < width_);
        assert(end[1] < height_);
        
        std::vector<unsigned int> position = start;
        unsigned int stride = end[0] - start[0] + 1;
        unsigned int rowCount = end[1] - start[1] + 1;
        
        // Copy matrices into coefficients_ and enqueue as one large write
        unsigned int offset;
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

        // TODO(CDE): DIRTY DIRTY DIRTY. Get rid of this ASAP! For now, it saves us from
        //            all those individual writes when initializing the coefficient plane.
        //            This WILL BREAK cached tiled mode.
#define WRITE_ONLY_ON_LAST_TILE
#ifdef WRITE_ONLY_ON_LAST_TILE
        if (end[0] == width_ - 1 && end[1] == height_ - 1) {
            int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, 0, matrixSize_ * width_ * height_, coefficients_, 0, NULL, NULL);
            if (err != CL_SUCCESS) {
                std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
            }
        }
        
        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + calcOffset(start), matrixSize_* width_ * height_);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + calcOffset(start), matrixSize_* width_ * height_);
        
        return;
#endif
        
        // TODO(CDE): Change this to one single copy rect of the entire region of coefficients
        // If the stride matches the entire x dimension of the coefficient plane, then...
        if (stride == width_) {
            // Enqueue CL commands
            int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, calcOffset(start), matrixSize_ * stride * rowCount, coefficients_, 0, NULL, NULL);
            if (err != CL_SUCCESS) {
                std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
            }
        } else {
            // Otherwise enqueue the writes row by row
            position[0] = start[0];
            for (position[1] = start[1]; position[1] <= end[1]; position[1]++) {
                // Enqueue CL commands
                int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, calcOffset(position), matrixSize_ * stride, coefficients_, 0, NULL, NULL);
                if (err != CL_SUCCESS) {
                    std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
                }
            }
        }

        // Update Cost Model
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, coefficients_ + calcOffset(start), matrixSize_* stride * rowCount);
        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, coefficients_ + calcOffset(start), matrixSize_* stride * rowCount);
    }
    
    void FillCoefficient(int coefficient,
                         int row, int col,
                         std::vector<unsigned int> start,
                         std::vector<unsigned int> end) {
        
        // TODO(CDE): Add future support for multiple planes
        assert(start.size() == 2);
        assert(start[0] < width_);
        assert(start[1] < height_);
        assert(end.size() == 2);
        assert(end[0] < width_);
        assert(end[1] < height_);

        std::vector<unsigned int> position = start;
        bool fillFinished = false;
        
        // Move from start to end, filling in each location with the provided pixel
        do {
            // Set coefficient in the coefficient matrix at this position in the coefficient plane
            getCoefficientMatrix(position[0], position[1])->setCoefficient(col, row, coefficient);
            
            // Move to the next position
            if (++position[0] > end[0]) {
                position[0] = start[0];
                if (++position[1] > end[1]) {
                    fillFinished = true;
                }
            }
        } while (!fillFinished);
        
        // Enqueue CL commands
        // TODO(CDE): Add CL support
    }

    cl_mem* initializeCl(cl_context context, cl_command_queue queue) {
        
        // Set CL variables
        clContext_ = context;
        clQueue_ = queue;
        
        // Create CL mem buffer
        clBuffer_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY,
                                   matrixSize_ * width_ * height_, NULL, NULL);
        return &clBuffer_ ;
    }
    
    // TODO(CDE): Remove this garbage
    cl_mem* getClBuffer() {
        return &clBuffer_;
    }
    
};

#endif
