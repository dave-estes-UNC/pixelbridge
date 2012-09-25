//
//  FrameVolume.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_FrameVolume_h
#define pixelbridge_FrameVolume_h

#include <cstdlib>
#include <vector>
#include <cassert>

#include "NDimensionalDisplayInterface.h"

using namespace std;

namespace nddi {

    class FrameVolume {

    protected:
        CostModel             * costModel_;
        vector<unsigned int>    dimensionalSizes_;
        unsigned int            size_;
        Pixel                 * pixels_;

    public:

        FrameVolume(CostModel* costModel,
                    vector<unsigned int> frameVolumeDimensionalSizes)
        : costModel_(costModel), size_(1), pixels_(NULL) {

            dimensionalSizes_ = frameVolumeDimensionalSizes;
            for (int i = 0; i < dimensionalSizes_.size(); i++) {
                size_ *= dimensionalSizes_[i];
            }
            pixels_ = (Pixel *)malloc(sizeof(Pixel) * size_);
            memset(pixels_, 0x00, sizeof(Pixel) * size_);
        }

        ~FrameVolume() {

            if (!dimensionalSizes_.empty()) {
                dimensionalSizes_.clear();
            }
            if (pixels_) {
                free(pixels_);
            }
        }

        unsigned int getSize() {

            return size_;
        }

        void PutPixel(Pixel p, std::vector<unsigned int> location) {
            setPixel(location, p);
        }

        void CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
            int dimensionToCopyAlong;
            bool dimensionFound = false;

            // Find the dimension to copy along
            for (int i = 0; !dimensionFound && (i < dimensionalSizes_.size()); i++) {
                if (start[i] != end[i]) {
                    dimensionToCopyAlong = i;
                    dimensionFound = true;
                }
            }

            std::vector<unsigned int> position = start;
            for (int j = 0; j <= end[dimensionToCopyAlong] - start[dimensionToCopyAlong]; j++) {
                setPixel(position, p[j]);
                position[dimensionToCopyAlong]++;
            }
        }

        void CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {

            std::vector<unsigned int> position = start;
            bool copyFinished = false;
            int pixelsCopied = 0;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set pixel in frame volume at position
                setPixel(position, p[pixelsCopied]);
                pixelsCopied++;

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
                            copyFinished = true;
                    }
                } while (overflow && !copyFinished);

            } while (!copyFinished);
        }

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
        }

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
        }

        void setPixel(vector<unsigned int> location, Pixel pixel) {

            unsigned int  offset = 0;
            unsigned int  multiplier = 1;

            assert(dimensionalSizes_.size() == location.size());

            for (int i = 0; i < location.size(); i++) {
                assert(location[i] < dimensionalSizes_[i]);

                offset += location[i] * multiplier;
                multiplier *= dimensionalSizes_[i];
            }

            pixels_[offset].packed = pixel.packed;

            costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, WRITE_ACCESS, pixels_ + offset, 4, 0);
        }

        Pixel getPixel(vector<unsigned int> location) {

            Pixel         pixel;
            unsigned int  offset = 0;
            unsigned int  multiplier = 1;

            assert(dimensionalSizes_.size() == location.size());

            for (int i = 0; i < location.size(); i++) {
                assert(location[i] < dimensionalSizes_[i]);

                offset += location[i] * multiplier;
                multiplier *= dimensionalSizes_[i];
            }

            pixel.packed = pixels_[offset].packed;

            costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, READ_ACCESS, pixels_ + offset, 4, 0);

            return pixel;
        }

        Pixel * data() {

            return pixels_;
        }
    };
}

#endif
