//
//  ClFrameVolume.h
//  pixelbridge
//
//  Created by Dave Estes on 4/1/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_ClFrameVolume_h
#define pixelbridge_ClFrameVolume_h

#include <vector>
#include <cassert>
#include <string.h>

#include "FrameVolume.h"

using namespace nddi;

class ClFrameVolume : public FrameVolume {
    
private:
    cl_mem            clBuffer_;
    cl_context        clContext_;
    cl_command_queue  clQueue_;
    size_t            fv_row_pitch_, fv_slice_pitch_;
    
    inline unsigned int calcOffset(std::vector<unsigned int> location) {
        unsigned int  offset = 0;
        unsigned int  multiplier = 1;
        
        assert(dimensionalSizes_.size() == location.size());
        
        for (int i = 0; i < location.size(); i++) {
            assert(location[i] < dimensionalSizes_[i]);
            
            offset += location[i] * multiplier;
            multiplier *= dimensionalSizes_[i];
        }
        
        return offset;
    }
    
public:
    
    ClFrameVolume(CostModel* costModel,
                vector<unsigned int> frameVolumeDimensionalSizes)
    : FrameVolume(costModel, frameVolumeDimensionalSizes), clBuffer_(NULL), clContext_(NULL), clQueue_(NULL) {
        
        fv_row_pitch_ = fv_slice_pitch_ = 0;
        switch (dimensionalSizes_.size()) {
            case 3:
                fv_slice_pitch_ = dimensionalSizes_[1] * dimensionalSizes_[0] * sizeof(Pixel);
                // No break
            case 2:
                fv_row_pitch_ = dimensionalSizes_[0] * sizeof(Pixel);
                // No break
            case 1:
            default:
                break;
        }
    }
    
    ~ClFrameVolume() {
        
        if (!dimensionalSizes_.empty()) {
            dimensionalSizes_.clear();
        }
        if (pixels_) {
            free(pixels_);
        }
        
        // Release CL mem buffer
        if (clBuffer_ != 0)
            clReleaseMemObject(clBuffer_);
    }
    
    void PutPixel(Pixel p, std::vector<unsigned int> location) {
        
        unsigned int offset = calcOffset(location);
            
        // Enqueue CL commands
        int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, offset, sizeof(unsigned int), &p, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
        }
        
        costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, WRITE_ACCESS, NULL, 4);
    }
    
    // TODO(CDE): move to buffer writes in 3 dimensions. Consider not even supporting strips in any other dimensions
    void CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
        
        int dimensionToCopyAlong;
        bool dimensionFound = false;
        unsigned int stripLength;
        
        // Find the dimension to copy along
        for (int i = 0; !dimensionFound && (i < dimensionalSizes_.size()); i++) {
            if (start[i] != end[i]) {
                dimensionToCopyAlong = i;
                dimensionFound = true;
                stripLength = end[dimensionToCopyAlong] - start[dimensionToCopyAlong];
            }
        }
        
        // If we're copying along the x dimension, then we can copy at at once
        if (dimensionToCopyAlong == 0) {
            
            // Enqueue CL command
            unsigned int offset = calcOffset(start);
            
            // Enqueue CL commands
            int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, offset, sizeof(unsigned int) * stripLength, p, 0, NULL, NULL);
            if (err != CL_SUCCESS) {
                std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
            }

        // Otherwise copy one pixel at a time
        } else {
            std::vector<unsigned int> position = start;
            for (int j = 0; j <= stripLength; j++) {
                
                // Calculate offset
                unsigned int offset = calcOffset(position);
                
                // Enqueue CL command
                int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE, offset, sizeof(unsigned int), p + offset, 0, NULL, NULL);
                if (err != CL_SUCCESS) {
                    std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << std::endl;
                }
                
                // Move to next position
                position[dimensionToCopyAlong]++;
            }
        }

        // Register memory charge
        costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, WRITE_ACCESS, NULL, 4 * stripLength);
    }
    
    void CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
        
        std::vector<unsigned int> position = start;
        bool copyFinished = false;
        int pixelsCopied = 0;

        size_t buffer_origin[3] = { 0, 0, 0 };
        size_t host_origin[3] = { 0, 0, 0 };
        size_t region[3] = { 1, 1, 1 };
        size_t p_row_pitch = 0, p_slice_pitch = 0;
        
        // Setup enqueue parameters that are unchanged within the loop
        switch (position.size()) {
            case 3:
                buffer_origin[2] = position[2];
                region[2] = end[2] - start[2] + 1;
                // slice pitch is effectively region[1] * region[0] * sizeof(Pixel)
                p_slice_pitch = (end[1] - start[1] + 1) * (end[0] - start[0] + 1) * sizeof(Pixel);
                // No break
            case 2:
                buffer_origin[1] = position[1];
                region[1] = end[1] - start[1] + 1;
                // row pitch is effectively region[0] * sizeof(Pixel)
                p_row_pitch = (end[0] - start[0] + 1) * sizeof(Pixel);
                // No break
            case 1:
                buffer_origin[0] = position[0] * sizeof(Pixel);
                region[0] = (end[0] - start[0] + 1) * sizeof(Pixel);
                // No break
            default:
                break;
        }

        // Copy 3D volumes at a time
        do {
            // First copy into pixels array
            void* startPixel = (unsigned char *)pixels_ + buffer_origin[2] * fv_slice_pitch_ + buffer_origin[1] * fv_row_pitch_ + buffer_origin[0];
            memcpy(startPixel, p, region[2] * region[1] * region[0]);
            // Enqueue CL command
            if (clQueue_ && clBuffer_) {
                int err = clEnqueueWriteBufferRect(clQueue_, clBuffer_, CL_TRUE,
                                                   buffer_origin,
                                                   host_origin,
                                                   region,
                                                   fv_row_pitch_, fv_slice_pitch_,
                                                   p_row_pitch, p_slice_pitch,
                                                   p, 0, NULL, NULL);
                if (err != CL_SUCCESS) {
                    std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << err << std::endl;
                }
                
            } else {
                std::cout << __FUNCTION__ << " - Still null?." << std::endl;
            }

            // Update pixel count
            pixelsCopied += region[0] * region[1] * region[2];
            
            if (dimensionalSizes_.size() <= 3) {
                copyFinished = true;
            } else {
                // Move to the next position
                int fvDim = 3; // Starting at 4th dimension since we're copying full xyz volumes at a time
                bool overflow;
                do {
                    overflow = false;
                    position[fvDim]++;
                    if ( (position[fvDim] >= dimensionalSizes_[fvDim])
                        || (position[fvDim] > end[fvDim]) ) {
                        overflow = true;
                        position[fvDim] = start[fvDim];
                        if (++fvDim >= dimensionalSizes_.size())
                            copyFinished = true;
                    }
                } while (overflow && !copyFinished);
            
                // Update host_origin[2] and buffer_origin[2] to copy the next cubic region
                host_origin[2] += p_slice_pitch * region[2];
                buffer_origin[2] += fv_slice_pitch_ * region[2];
            }
        } while (!copyFinished);
        
        // Register memory charge
        costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, WRITE_ACCESS, NULL, 4 * pixelsCopied);
    }
    
    // TODO(CDE): Implement for CL
    void FillPixel(Pixel p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
        
        std::vector<unsigned int> position = start;
        bool fillFinished = false;
        int pixelsFilled = 0;
        
        // Move from start to end, filling in each location with the provided pixel
        do {
            // Set pixel in frame volume at position
            setPixel(position, p);
            pixelsFilled++;
            
            // Move to the next position
            int fvDim = 0;
            bool overflow;
            do {
                overflow = false;
                position[fvDim]++;
                if ( (position[fvDim] >= dimensionalSizes_[fvDim])
                    || (position[fvDim] > end[fvDim]) ) {
                    overflow = true;
                    position[fvDim] = start[fvDim];
                    if (++fvDim >= dimensionalSizes_.size())
                        fillFinished = true;
                }
            } while (overflow && !fillFinished);
            
        } while (!fillFinished);
        
        // Enqueue CL commands
    }
    
    // TODO(CDE): Implement for CL
    void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest) {
        
        std::vector<unsigned int> positionFrom = start;
        std::vector<unsigned int> positionTo = dest;
        bool copyFinished = false;
        int pixelsCopied = 0;
        
        // Move from start to end, filling in each location with the provided pixel
        do {
            // Set pixel in frame volume at position
            setPixel(positionFrom, getPixel(positionTo));
            pixelsCopied++;
            
            // Move to the next position
            int fvDim = 0;
            bool overflow;
            do {
                overflow = false;
                positionFrom[fvDim]++;
                positionTo[fvDim]++;
                if ( (positionFrom[fvDim] >= dimensionalSizes_[fvDim])
                    || (positionFrom[fvDim] > end[fvDim]) ) {
                    overflow = true;
                    positionFrom[fvDim] = start[fvDim];
                    positionTo[fvDim] = dest[fvDim];
                    if (++fvDim >= dimensionalSizes_.size())
                        copyFinished = true;
                }
            } while (overflow && !copyFinished);
            
        } while (!copyFinished);
        
        // Enqueue CL commands
    }
    
    cl_mem* initializeCl(cl_context context, cl_command_queue queue) {
        
        // Set CL variables
        clContext_ = context;
        clQueue_ = queue;
        
        // Create CL mem buffer
        clBuffer_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY,
                                   sizeof(unsigned int) * getSize(), NULL, NULL);

        return &clBuffer_;
    }

    // TODO(CDE): Remove this garbage
    cl_mem* getClBuffer() {
        return &clBuffer_;
    }
    
};

#endif
