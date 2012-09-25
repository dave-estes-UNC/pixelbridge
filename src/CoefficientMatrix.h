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
        int           * coefficients_;

    public:

        // Basic contructor
        CoefficientMatrix(CostModel* costModel,
                          unsigned int width, unsigned int height)
        : costModel_(costModel), width_(width), height_(height), coefficients_(NULL) {

            coefficients_ = (int *)malloc(sizeof(int) * width * height);
            memset(coefficients_, 0x00, sizeof(int) * width * height);
        }

        // Copy constructor
        CoefficientMatrix(CoefficientMatrix* cm)
        : costModel_(cm->costModel_), width_(cm->width_), height_(cm->height_), coefficients_(NULL) {

            coefficients_ = (int *)malloc(sizeof(int) * cm->width_ * cm->height_);
            memcpy(coefficients_, cm->coefficients_, sizeof(int) * cm->width_ * cm->height_);
        }

        ~CoefficientMatrix() {

            if (coefficients_) {
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

            return width_ * height_ * sizeof(int);
        }

        void setCoefficient(unsigned int x, unsigned int y, int value) {

            assert(x < width_);
            assert(y < height_);
            coefficients_[y * width_ + x] = value;
        }

        void setCoefficients(std::vector< std::vector<int> > coefficientVector) {

            assert(coefficientVector.size() == width_);
            assert(coefficientVector[0].size() == height_);

            // Examine each coefficient in the coefficient matrix vector and use it unless it's a COFFICIENT_UNCHANGED
            for (int y = 0; y < height_; y++) {
                for (int x = 0; x < width_; x++) {
                    if (coefficientVector[x][y] != COFFICIENT_UNCHANGED) {
                        coefficients_[y * width_ + x] = coefficientVector[x][y];
                        costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &coefficients_[y * width_ + x], 4, 0);
                    }
                }
            }
        }

        int getCoefficient(unsigned int x, unsigned int y) {

            assert(x < width_);
            assert(y < height_);

            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &coefficients_[y * width_ + x], 4, 0);

            return coefficients_[y * width_ + x];
        }

        int * data() {

            return coefficients_;
        }
    };
}

#endif
