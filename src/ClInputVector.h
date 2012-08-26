//
//  ClInputVector.h
//  pixelbridge
//
//  Created by Dave Estes on 3/26/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_ClInputVector_h
#define pixelbridge_ClInputVector_h

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#endif

#include "InputVector.h"

using namespace nddi;

class ClInputVector : public InputVector {

private:
    cl_mem            clBuffer_;
    cl_context        clContext_;
    cl_command_queue  clQueue_;

public:

    // TODO(CDE): Don't use super's ctor
    ClInputVector(CostModel* costModel,
                  unsigned int size)
    : InputVector(costModel, size) {
    }

    ~ClInputVector() {

        if (values_) {
            free(values_);
        }

        // Release CL mem buffer
        if (clBuffer_ != 0)
            clReleaseMemObject(clBuffer_);
    }

    void UpdateInputVector(std::vector<int> input) {

        assert(input.size() + 2 == size_);

        // Update the underlying Input Vector and cost model
        for (int i = 0; (i < input.size()) && ((i + 2) < size_); i++) {
            setValue(i+2, input[i]);
        }

        // Enqueue CL command to write from the underlying Input Vector
        if (clQueue_ && clBuffer_) {
            // TODO(CDE): Change this back to skip the x and y values. There some kind of
            //            bug with it at the moment, so just writing the whole thing.
//            int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE,
//                                           2 * sizeof(int),
//                                           (size_  - 2) * sizeof(int),
//                                           values_, 0, NULL, NULL);
            int err = clEnqueueWriteBuffer(clQueue_, clBuffer_, CL_TRUE,
                                           0,
                                           size_ * sizeof(int),
                                           values_, 0, NULL, NULL);
            if (err != CL_SUCCESS) {
                std::cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << err << std::endl;
            }
        }
    }

    cl_mem* initializeCl(cl_context context, cl_command_queue queue) {

        // Set CL variables
        clContext_ = context;
        clQueue_ = queue;

        // Create CL mem buffer
        clBuffer_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY, sizeof(int) * size_, NULL, NULL);
        return &clBuffer_;
    }

    // TODO(CDE): Remove this garbage
    cl_mem* getClBuffer() {
        return &clBuffer_;
    }

};

#endif
