//
//  CoefficientPlanes.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CoefficientPlanes_h
#define pixelbridge_CoefficientPlanes_h

#include <cstdlib>
#include <cassert>

#include "Configuration.h"
#include "NDimensionalDisplayInterface.h"
#include "CoefficientMatrix.h"

/**
 * These helper macros will calculate the offset into scaler memory or coefficient memory. This is helpful
 * for tweaking the layout when doing memory profiling experiments.
 */
// R x C x P
#define SC_OFF(x, y, p, c)   ((y * width_ * numPlanes_ + x * numPlanes_ + p) * 3 + c)
//#define CP_OFF(x, y, p)   ((y * width_ * numPlanes_ + x * numPlanes_ + p) * matrixWidth_ * matrixHeight_)

// P x R x C
//#define SC_OFF(x, y, p, c)   ((p * height_ * width_ + y * width_ + x) * 3 + c)
#define CP_OFF(x, y, p)   ((p * height_ * width_ + y * width_ + x) * matrixWidth_ * matrixHeight_)

using namespace std;

namespace nddi {

    class CoefficientPlanes {

    protected:
        CostModel           * costModel_;
        size_t                width_, height_, numPlanes_, matrixWidth_, matrixHeight_;
        CoefficientMatrix   * coefficientMatrix_;
        Coeff               * coefficients_;
        int16_t             * scalers_;

    public:

        CoefficientPlanes() {
        }

        CoefficientPlanes(CostModel* costModel,
                         unsigned int displayWidth, unsigned int displayHeight,
                         unsigned int numPlanes,
                         unsigned int matrixWidth, unsigned int matrixHeight)
        : costModel_(costModel),
          width_(displayWidth), height_(displayHeight),
          numPlanes_(numPlanes),
          matrixWidth_(matrixWidth), matrixHeight_(matrixHeight) {

            // Create the common CoefficientMatrix
            coefficientMatrix_ =  new CoefficientMatrix(costModel_, matrixWidth, matrixHeight);

            // Alloc the actual coefficients and scalers
            if (!globalConfiguration.headless) {
                coefficients_ = (Coeff *)malloc(CoefficientMatrix::memoryRequired(matrixWidth, matrixHeight) * displayWidth * displayHeight * numPlanes_);
                scalers_ = (int16_t *)malloc(sizeof(int16_t) * 3 * displayWidth * displayHeight * numPlanes_);
            }
        }

        ~CoefficientPlanes() {

            if (coefficientMatrix_) delete(coefficientMatrix_);
            if (coefficients_) free(coefficients_);
            if (scalers_) free(scalers_);
        }

        unsigned int getWidth() {

            return width_;
        }

        unsigned int getHeight() {

            return height_;
        }

        Coeff GetCoefficient(vector<unsigned int> &location, int row, int col) {
            assert(location.size() == 3);
            assert(row < matrixWidth_);
            assert(col < matrixHeight_);

            Coeff *cm = dataCoefficient(location[0], location[1], location[2]);
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &cm[col * matrixWidth_ + row], BYTES_PER_COEFF, 0);

            return cm[col * matrixWidth_ + row];
        }

        void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) {

            assert(location.size() == 3);
            assert(coefficientMatrix.size() == matrixWidth_);
            assert(coefficientMatrix[0].size() == matrixHeight_);

            if (!globalConfiguration.headless) {
                // Examine each coefficient in the coefficient matrix vector and use it unless it's a COFFICIENT_UNCHANGED
                for (int col = 0; col < matrixHeight_; col++) {
                    for (int row = 0; row < matrixWidth_; row++) {
                        if (coefficientMatrix[row][col] != COFFICIENT_UNCHANGED) {
    #ifdef NARROW_DATA_STORES
                            assert(coefficientMatrix[row][col] >= SHRT_MIN && coefficientMatrix[row][col] <= SHRT_MAX);
    #endif
                            Coeff * cm = dataCoefficient(location[0], location[1], location[2]);
                            cm[col * matrixWidth_ + row] = coefficientMatrix[row][col];
                            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &cm[col * matrixWidth_ + row], BYTES_PER_COEFF, 0);
                        }
                    }
                }
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
            assert(coefficientMatrix.size() == matrixWidth_);
            assert(coefficientMatrix[0].size() == matrixHeight_);

            vector<unsigned int> position = start;
            bool fillFinished = false;
            int matricesFilled = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Update coefficient matrix in coefficient plane at position
                if (!globalConfiguration.headless) {
                    // Examine each coefficient in the coefficient matrix vector and use it unless it's a COFFICIENT_UNCHANGED
                    for (int col = 0; col < matrixHeight_; col++) {
                        for (int row = 0; row < matrixWidth_; row++) {
                            if (coefficientMatrix[row][col] != COFFICIENT_UNCHANGED) {
        #ifdef NARROW_DATA_STORES
                                assert(coefficientMatrix[row][col] >= SHRT_MIN && coefficientMatrix[row][col] <= SHRT_MAX);
        #endif
                                Coeff * cm = dataCoefficient(position[0], position[1], position[2]);
                                cm[col * matrixWidth_ + row] = coefficientMatrix[row][col];
                                costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &cm[col * matrixWidth_ + row], BYTES_PER_COEFF, 0);
                            }
                        }
                    }

                }
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
            assert(col < matrixWidth_);
            assert(row < matrixHeight_);
#ifdef NARROW_DATA_STORES
            assert(coefficient >= SHRT_MIN && coefficient <= SHRT_MAX);
#endif

            vector<unsigned int> position = start;
            bool fillFinished = false;
            int coefficientsFilled = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set coefficient in the coefficient matrix at this position in the coefficient plane
                if (!globalConfiguration.headless) {
                    Coeff* cm = dataCoefficient(position[0], position[1], position[2]);
                    cm[row * matrixWidth_ + col] = coefficient;
                    costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &cm[row * width_ + col], BYTES_PER_COEFF, 0);

                }
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

        void FillScalerStack(vector<uint64_t> &scalers,
                             vector<unsigned int> &start,
                             vector<unsigned int> &size) {

        }

        CoefficientMatrix* getCoefficientMatrix() {

            return coefficientMatrix_;
        }

        void putScaler(unsigned int x, unsigned int y, unsigned int p, Scaler scaler) {

            assert(x < width_);
            assert(y < height_);
            assert(!globalConfiguration.headless);

            // TODO(CDE): Get this working properly
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &scalers_[SC_OFF(x, y, p, 0)], BYTES_PER_SCALER, 0);

            scalers_[SC_OFF(x, y, p, 0)] = scaler.r;
            scalers_[SC_OFF(x, y, p, 1)] = scaler.g;
            scalers_[SC_OFF(x, y, p, 2)] = scaler.b;
        }

        Scaler getScaler(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);
            assert(p < numPlanes_);
            assert(!globalConfiguration.headless);

            // TODO(CDE): Get this working properly
            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &scalers_[p * width_ * height_ + y * width_ + x], BYTES_PER_SCALER, 0);

            Scaler s;
            s.packed = 0;
            s.r = scalers_[SC_OFF(x, y, p, 0)];
            s.g = scalers_[SC_OFF(x, y, p, 1)];
            s.b = scalers_[SC_OFF(x, y, p, 2)];

            return s;
        }

        int16_t * dataScaler() {
            assert(!globalConfiguration.headless);
            return scalers_;
        }

        Coeff * dataCoefficient(size_t x, size_t y, size_t p) {
            assert(!globalConfiguration.headless);
            return &coefficients_[CP_OFF(x, y, p)];
        }

        inline size_t computeScalerOffset(size_t x, size_t y, size_t p) {
            return SC_OFF(x, y, p, 0);
        }

    };
}

#endif
