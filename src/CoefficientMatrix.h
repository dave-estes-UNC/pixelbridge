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

    class CoefficientMatrix {

        CostModel     * costModel_;
        unsigned int    width_, height_;
#ifdef NARROW_DATA_STORES
        int16_t       * coefficients_;
#else
        int           * coefficients_;
#endif
        bool            preallocated;

    public:

        /**
         * Basic contructor for a coefficent matrix. Will allocate the memory for the
         * coefficients, this setting the preallocated flag to false since there was not
         * memory previously allocated for this coefficient matrix.
         */
        CoefficientMatrix(CostModel* costModel,
                          unsigned int width, unsigned int height)
        : costModel_(costModel), width_(width), height_(height), coefficients_(NULL) {

            if (!globalConfiguration.headless) {
#ifdef NARROW_DATA_STORES
                coefficients_ = (int16_t *)malloc(sizeof(int16_t) * width * height);
                memset(coefficients_, 0x00, sizeof(int16_t) * width * height);
#else
                coefficients_ = (int *)malloc(sizeof(int) * width * height);
                memset(coefficients_, 0x00, sizeof(int) * width * height);
#endif
            }
            preallocated = false;
        }

        /**
         * Alternative constructor that will not allocate its own memory
         * for the coefficients. Using the contructor is ideal if the memory
         * layout for all of the coefficient planes is being carefully controlled.
         */
        CoefficientMatrix(CostModel* costModel,
                          unsigned int width, unsigned int height,
                          void* memory)
        : costModel_(costModel), width_(width), height_(height), coefficients_(NULL) {

#ifdef NARROW_DATA_STORES
            coefficients_ = (int16_t *)memory;
#else
            coefficients_ = (int *)memory;
#endif
            preallocated = true;
        }

        /**
         * Copy constructor
         */
        CoefficientMatrix(CoefficientMatrix* cm)
        : costModel_(cm->costModel_), width_(cm->width_), height_(cm->height_), coefficients_(NULL) {

            if (!globalConfiguration.headless) {
#ifdef NARROW_DATA_STORES
                coefficients_ = (int16_t *)malloc(sizeof(int16_t) * cm->width_ * cm->height_);
                memcpy(coefficients_, cm->coefficients_, sizeof(int16_t) * cm->width_ * cm->height_);
#else
                coefficients_ = (int *)malloc(sizeof(int) * cm->width_ * cm->height_);
                memcpy(coefficients_, cm->coefficients_, sizeof(int) * cm->width_ * cm->height_);
#endif
            }
        }

        ~CoefficientMatrix() {

            if (!preallocated && coefficients_) {
                free(coefficients_);
            }
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

        void setCoefficient(unsigned int x, unsigned int y, int value) {

            assert(x < width_);
            assert(y < height_);
#ifdef NARROW_DATA_STORES
            assert(value >= SHRT_MIN && value <= SHRT_MAX);
#endif
            assert(!globalConfiguration.headless);

            coefficients_[y * width_ + x] = value;
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &coefficients_[y * width_ + x], BYTES_PER_COEFF, 0);
        }

        void setCoefficients(vector< vector<int> > &coefficientVector) {

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
                        coefficients_[y * width_ + x] = coefficientVector[x][y];
                        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &coefficients_[y * width_ + x], BYTES_PER_COEFF, 0);
                    }
                }
            }
        }

        int getCoefficient(unsigned int x, unsigned int y) {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &coefficients_[y * width_ + x], BYTES_PER_COEFF, 0);
            return coefficients_[y * width_ + x];
        }

#ifdef NARROW_DATA_STORES
        int16_t * data() {
#else
        int * data() {
#endif
            assert(!globalConfiguration.headless);

            return coefficients_;
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
