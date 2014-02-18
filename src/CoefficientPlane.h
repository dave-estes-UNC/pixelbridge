//
//  CoefficientPlane.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CoefficientPlane_h
#define pixelbridge_CoefficientPlane_h

#include <cstdlib>
#include <cassert>

#include "Configuration.h"
#include "NDimensionalDisplayInterface.h"
#include "CoefficientMatrix.h"

using namespace std;

namespace nddi {

    class CoefficientPlane {

    protected:
        CostModel           * costModel_;
        unsigned int          width_, height_, numPlanes_, matrixWidth_, matrixHeight_;
        CoefficientMatrix  ** coefficientMatrices_;
#ifdef NARROW_DATA_STORES
        int16_t             * scalers_;
#else
        int                 * scalers_;
#endif

    public:

        CoefficientPlane() {
        }

        CoefficientPlane(CostModel* costModel,
                         unsigned int displayWidth, unsigned int displayHeight,
                         unsigned int numPlanes,
                         unsigned int matrixWidth, unsigned int matrixHeight)
        : costModel_(costModel),
          width_(displayWidth), height_(displayHeight),
          numPlanes_(numPlanes),
          matrixWidth_(matrixWidth), matrixHeight_(matrixHeight),
          coefficientMatrices_(NULL) {

            coefficientMatrices_ = (CoefficientMatrix **)malloc(sizeof(CoefficientMatrix *) * displayWidth * displayHeight * numPlanes_);
            for (int p = 0; p < numPlanes_; p++) {
                for (int y = 0; y < displayHeight; y++) {
                    for (int x = 0; x < displayWidth; x++) {
                        if (!globalConfiguration.headless)
                            coefficientMatrices_[p * displayWidth * displayHeight + y * displayWidth + x] = new CoefficientMatrix(costModel_, matrixWidth, matrixHeight);
                    }
                }
            }

            if (!globalConfiguration.headless) {
#ifdef NARROW_DATA_STORES
                scalers_ = (int16_t *)malloc(sizeof(int16_t) * 3 * displayWidth * displayHeight * numPlanes_);
#else
                scalers_ = (int *)malloc(sizeof(int) * 3 * displayWidth * displayHeight * numPlanes_);
#endif
            }
        }

        ~CoefficientPlane() {

            CoefficientMatrix * matrix;

            if (coefficientMatrices_) {

                for (int y = 0; y < height_; y++) {
                    for (int x = 0; x < width_; x++) {
                        matrix = coefficientMatrices_[y * width_ + x];
                        if (matrix) {
                            delete(matrix);
                        }
                    }
                }

                free(coefficientMatrices_);
            }
            if (scalers_) free(scalers_);
        }

        unsigned int getWidth() {

            return width_;
        }

        unsigned int getHeight() {

            return height_;
        }

        void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) {

            assert(location.size() == 3);

            if (!globalConfiguration.headless) {
                getCoefficientMatrix(location[0], location[1], location[2])->setCoefficients(coefficientMatrix);
            } else {
                costModel_->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                                     1 * (matrixHeight_ * matrixWidth_),
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     1 * (matrixHeight_ * matrixWidth_) * BYTES_PER_COEFF,
                                                     0);

            }
        }

        void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                   vector<unsigned int> &start,
                                   vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(end.size() == 3);

            vector<unsigned int> position = start;
            bool fillFinished = false;
            int matricesFilled = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Update coefficient matrix in coefficient plane at position
                if (!globalConfiguration.headless)
                    getCoefficientMatrix(position[0], position[1], position[2])->setCoefficients(coefficientMatrix);
                matricesFilled++;

                // Move to the next position
                if (++position[0] > end[0]) {
                    position[0] = start[0];
                    if (++position[1] > end[1]) {
                        position[1] = start[1];
                        if (++position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            if (globalConfiguration.headless)
                costModel_->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                                     matricesFilled * (matrixHeight_ * matrixWidth_),
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     matricesFilled * (matrixHeight_ * matrixWidth_) * BYTES_PER_COEFF,
                                                     0);
        }

        void FillCoefficient(int coefficient,
                             int row, int col,
                             vector<unsigned int> &start,
                             vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(end.size() == 3);

            vector<unsigned int> position = start;
            bool fillFinished = false;
            int coefficientsFilled = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set coefficient in the coefficient matrix at this position in the coefficient plane
                if (!globalConfiguration.headless)
                    getCoefficientMatrix(position[0], position[1], position[2])->setCoefficient(col, row, coefficient);
                coefficientsFilled++;

                // Move to the next position
                if (++position[0] > end[0]) {
                    position[0] = start[0];
                    if (++position[1] > end[1]) {
                        position[1] = start[1];
                        if (++position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            if (globalConfiguration.headless)
                costModel_->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                                     coefficientsFilled,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     coefficientsFilled * BYTES_PER_COEFF,
                                                     0);
        }

        void FillScaler(Scaler scaler,
                        vector<unsigned int> &start,
                        vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(end.size() == 3);

            vector<unsigned int> position = start;
            bool fillFinished = false;
            int scalersFilled = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set scaler at this position in the coefficient plane
                if (!globalConfiguration.headless)
                    putScaler(position[0], position[1], position[2], scaler);
                scalersFilled++;

                // Move to the next position
                if (++position[0] > end[0]) {
                    position[0] = start[0];
                    if (++position[1] > end[1]) {
                        position[1] = start[1];
                        if (++position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            if (globalConfiguration.headless)
                costModel_->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                                     scalersFilled,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     scalersFilled * BYTES_PER_SCALER,
                                                     0);
        }

        CoefficientMatrix* getCoefficientMatrix(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            return coefficientMatrices_[p * width_ * height_ + y * width_ + x];
        }

        void putScaler(unsigned int x, unsigned int y, unsigned int p, Scaler scaler) {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            // TODO(CDE): Get this working properly
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &scalers_[p * width_ * height_ + y * width_ + x], BYTES_PER_SCALER, 0);

            scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 0] = scaler.r;
            scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 1] = scaler.g;
            scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 2] = scaler.b;
        }

        Scaler getScaler(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            // TODO(CDE): Get this working properly
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &scalers_[p * width_ * height_ + y * width_ + x], BYTES_PER_SCALER, 0);

            Scaler s;
            s.packed = 0;
            s.r = scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 0];
            s.g = scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 1];
            s.b = scalers_[(p * width_ * height_ + y * width_ + x) * 3 + 2];

            return s;
        }

#ifdef NARROW_DATA_STORES
        int16_t * data() {
#else
        int * data() {
#endif
            return scalers_;
        }
    };
}

#endif
