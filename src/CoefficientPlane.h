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
        unsigned int          width_, height_;
        CoefficientMatrix  ** coefficientMatrices_;

    public:

        CoefficientPlane() {
        }

        CoefficientPlane(CostModel* costModel,
                         unsigned int displayWidth, unsigned int displayHeight,
                         unsigned int matrixWidth, unsigned int matrixHeight)
        : costModel_(costModel), width_(displayWidth), height_(displayHeight), coefficientMatrices_(NULL) {

            coefficientMatrices_ = (CoefficientMatrix **)malloc(sizeof(CoefficientMatrix *) * displayWidth * displayHeight);
            for (int y = 0; y < displayHeight; y++) {
                for (int x = 0; x < displayWidth; x++) {
                    coefficientMatrices_[y * displayWidth + x] = new CoefficientMatrix(costModel_, matrixWidth, matrixHeight);
                }
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
        }

        unsigned int getWidth() {

            return width_;
        }

        unsigned int getHeight() {

            return height_;
        }

        void PutCoefficientMatrix(vector< vector<int> > coefficientMatrix, vector<unsigned int> location) {

            getCoefficientMatrix(location[0], location[1])->setCoefficients(coefficientMatrix);
        }

        void FillCoefficientMatrix(vector< vector<int> > coefficientMatrix,
                                   vector<unsigned int> start,
                                   vector<unsigned int> end) {

            vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Update coefficient matrix in coefficient plane at position
                getCoefficientMatrix(position[0], position[1])->setCoefficients(coefficientMatrix);

                // Move to the next position
                if (++position[0] > end[0]) {
                    position[0] = start[0];
                    if (++position[1] > end[1]) {
                        fillFinished = true;
                    }
                }
            } while (!fillFinished);
        }

        void FillCoefficient(int coefficient,
                             int row, int col,
                             vector<unsigned int> start,
                             vector<unsigned int> end) {

            vector<unsigned int> position = start;
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
        }

        void putCoefficientMatrix(unsigned int x, unsigned int y, CoefficientMatrix * matrix) {

            assert(x < width_);
            assert(y < height_);

            CoefficientMatrix* m = coefficientMatrices_[y * width_ + x];
            if (m != NULL) {
                delete(m);
            }
            coefficientMatrices_[y * width_ + x] = matrix;
        }

        CoefficientMatrix* getCoefficientMatrix(unsigned int x, unsigned int y) {

            assert(x < width_);
            assert(y < height_);
            return coefficientMatrices_[y * width_ + x];
        }
    };
}

#endif
