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
        unsigned int    width_, height_;

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
#ifdef NARROW_DATA_STORES
            return width_ * height_ * sizeof(int16_t);
#else
            return width_ * height_ * sizeof(int);
#endif
        }

        void setCoefficient(unsigned int x, unsigned int y, int value,
#ifdef NARROW_DATA_STORES
                            int16_t * coefficients)
#else
                            int * coefficients)
#endif
        {

            assert(x < width_);
            assert(y < height_);
#ifdef NARROW_DATA_STORES
            assert(value >= SHRT_MIN && value <= SHRT_MAX);
#endif
            assert(!globalConfiguration.headless);

            coefficients[y * width_ + x] = value;
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &coefficients[y * width_ + x], BYTES_PER_COEFF, 0);
        }

        void setCoefficients(vector< vector<int> > &coefficientVector,
#ifdef NARROW_DATA_STORES
                            int16_t * coefficients)
#else
                            int * coefficients)
#endif
        {

            assert(coefficientVector.size() == width_);
            assert(coefficientVector[0].size() == height_);
            assert(!globalConfiguration.headless);

            // Examine each coefficient in the coefficient matrix vector and use it unless it's a COFFICIENT_UNCHANGED
            for (int y = 0; y < height_; y++) {
                for (int x = 0; x < width_; x++) {
                    if (coefficientVector[x][y] != COFFICIENT_UNCHANGED) {
#ifdef NARROW_DATA_STORES
                        assert(coefficientVector[x][y] >= SHRT_MIN && coefficientVector[x][y] <= SHRT_MAX);
#endif
                        coefficients[y * width_ + x] = coefficientVector[x][y];
                        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &coefficients[y * width_ + x], BYTES_PER_COEFF, 0);
                    }
                }
            }
        }

        int getCoefficient(unsigned int x, unsigned int y,
#ifdef NARROW_DATA_STORES
                            int16_t * coefficients)
#else
                            int * coefficients)
#endif
        {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &coefficients[y * width_ + x], BYTES_PER_COEFF, 0);
            return coefficients[y * width_ + x];
        }

        /**
         * When using preallocated memory, this function is first called to
         * determine how much memory the coefficient matrix needs.
         */
        static size_t memoryRequired(unsigned int width, unsigned int height) {
#ifdef NARROW_DATA_STORES
                return sizeof(int16_t) * width * height;
#else
                return sizeof(int) * width * height;
#endif

        }
    };
}

#endif
