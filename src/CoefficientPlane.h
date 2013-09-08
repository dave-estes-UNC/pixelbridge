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

#include "NDimensionalDisplayInterface.h"
#include "CoefficientMatrix.h"

using namespace std;

namespace nddi {

    class CoefficientPlane {

    protected:
        CostModel           * costModel_;
        unsigned int          width_, height_, numPlanes_;
        CoefficientMatrix  ** coefficientMatrices_;
        int                 * scalers_;

    public:

        CoefficientPlane() {
        }

        CoefficientPlane(CostModel* costModel,
                         unsigned int displayWidth, unsigned int displayHeight,
                         unsigned int matrixWidth, unsigned int matrixHeight)
        : costModel_(costModel), width_(displayWidth), height_(displayHeight), numPlanes_(NUM_COEFFICIENT_PLANES), coefficientMatrices_(NULL) {

            coefficientMatrices_ = (CoefficientMatrix **)malloc(sizeof(CoefficientMatrix *) * displayWidth * displayHeight * NUM_COEFFICIENT_PLANES);
            for (int p = 0; p < NUM_COEFFICIENT_PLANES; p++) {
				for (int y = 0; y < displayHeight; y++) {
					for (int x = 0; x < displayWidth; x++) {
						coefficientMatrices_[p * displayWidth * displayHeight + y * displayWidth + x] = new CoefficientMatrix(costModel_, matrixWidth, matrixHeight);
					}
				}
            }

            scalers_ = (int *)malloc(sizeof(int) * displayWidth * displayHeight * NUM_COEFFICIENT_PLANES);
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

            getCoefficientMatrix(location[0], location[1], location[2])->setCoefficients(coefficientMatrix);
        }

        void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                   vector<unsigned int> &start,
                                   vector<unsigned int> &end) {

        	assert(start.size() == 3);
        	assert(end.size() == 3);

            vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Update coefficient matrix in coefficient plane at position
                getCoefficientMatrix(position[0], position[1], position[2])->setCoefficients(coefficientMatrix);

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
        }

        void FillCoefficient(int coefficient,
                             int row, int col,
                             vector<unsigned int> &start,
                             vector<unsigned int> &end) {

        	assert(start.size() == 3);
        	assert(end.size() == 3);

        	vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set coefficient in the coefficient matrix at this position in the coefficient plane
                getCoefficientMatrix(position[0], position[1], position[2])->setCoefficient(col, row, coefficient);

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
        }

        void FillScaler(int scaler,
                        vector<unsigned int> &start,
                        vector<unsigned int> &end) {

        	assert(start.size() == 3);
        	assert(end.size() == 3);

        	vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set scaler at this position in the coefficient plane
                putScaler(position[0], position[1], position[2], scaler);

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
        }

        void putCoefficientMatrix(unsigned int x, unsigned int y, unsigned int p, CoefficientMatrix * matrix) {

            assert(x < width_);
            assert(y < height_);

            CoefficientMatrix* m = coefficientMatrices_[p * width_ * height_ * numPlanes_ + y * width_ + x];
            if (m != NULL) {
                delete(m);
            }
            coefficientMatrices_[p * width_ * height_ + y * width_ + x] = matrix;
        }

        CoefficientMatrix* getCoefficientMatrix(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);
            return coefficientMatrices_[p * width_ * height_ + y * width_ + x];
        }

        void putScaler(unsigned int x, unsigned int y, unsigned int p, int scaler) {

            assert(x < width_);
            assert(y < height_);

            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, WRITE_ACCESS, &scalers_[p * width_ * height_ + y * width_ + x], BYTES_PER_SCALER, 0);

            scalers_[p * width_ * height_ + y * width_ + x] = scaler;
        }

        int getScaler(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);

            costModel_->registerMemoryCharge(COEFFICIENT_PLANE_COMPONENT, READ_ACCESS, &scalers_[p * width_ * height_ + y * width_ + x], BYTES_PER_SCALER, 0);

            return scalers_[p * width_ * height_ + y * width_ + x];
        }
    };
}

#endif
