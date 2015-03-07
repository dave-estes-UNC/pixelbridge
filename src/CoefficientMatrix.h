//
//  CoefficientMatrix.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CoefficientMatrix_h
#define pixelbridge_CoefficientMatrix_h

#include <limits.h>
#include <cassert>

namespace nddi {

    /**
     *  Originally, the CoefficientMatrix class would be instantiated for each and every
     *  coefficient matrix used. This quickly grew out of hand. Now it represents a single
     *  logical coefficient matrix, which acts as a sort of static class that holds just
     *  enough common information to be usefull along with the methods to get and set
     *  coefficients. Those methods are changed to accept a memory location where the coefficients
     *  are logically held with the coefficient planes themselves.
     */
    class CoefficientMatrix {

        CostModel     * costModel_;
        size_t          width_, height_;

    public:

        /**
         * Basic contructor for a coefficent matrix.
         */
        CoefficientMatrix(CostModel* costModel,
                          unsigned int width, unsigned int height)
        : costModel_(costModel), width_(width), height_(height) {
        }

        ~CoefficientMatrix() {
        }

        unsigned int getWidth() {
            return width_;
        }

        unsigned int getHeight() {
            return height_;
        }

        unsigned int getSize() {
            return width_ * height_;
        }

        unsigned int getDataSize() {
            return width_ * height_ * sizeof(Coeff);
        }

        /**
         * When using preallocated memory, this function is first called to
         * determine how much memory the coefficient matrix needs.
         */
        static size_t memoryRequired(unsigned int width, unsigned int height) {
                return sizeof(Coeff) * width * height;
        }
    };
}

#endif
