#ifndef BASE_NDDI_DISPLAY_H
#define BASE_NDDI_DISPLAY_H

#include "NDimensionalDisplayInterface.h"
#include "InputVector.h"
#include "CoefficientPlanes.h"
#include "FrameVolume.h"

#define CM_WIDTH   (inputVector_->getSize())
#define CM_HEIGHT  (frameVolumeDimensionalSizes_.size())
#define CM_SIZE    (CM_WIDTH * CM_HEIGHT)
#define CM_STORAGE (sizeof(int) * CM_SIZE)

namespace nddi {

    /**
     * This base class implements the NDimensionalDisplayInterface API and is meant,
     * to be extended by various other forms of an NDDI Display, such as a test display
     * and an OpenGL display that will display the framebuffer in a window.
     */
    class BaseNddiDisplay : public NDimensionalDisplayInterface {

    protected:

        /**
         * Renders each pixel of the frame buffer by setting the x, y in the input vector and computing which
         * frame volume pixel to use.
         */
        virtual void Render() = 0;

        /**
         * Holds the Display Width
         */
        size_t displayWidth_;

        /**
         * Holds the Display Height
         */
        size_t displayHeight_;

        /**
         * Holds the number of coefficient planes
         */
        size_t numPlanes_;

        /**
         * The input vector
         */
        InputVector * inputVector_;

        /**
         * Holds the sizes of each dimension of the Frame Volume.
         */
        vector<unsigned int> frameVolumeDimensionalSizes_;

        /**
         * The frameVolume_ is physcially a flat buffer of Pixels, but it logically managed based on the configured
         * dimensions in the constructor.
         */
        FrameVolume * frameVolume_;

        /**
         * The coefficientPlanes_ is a flat buffer of coefficient matrices that are logically
         * arranged into rows and columns maching the dimenions of the display.
         */
        CoefficientPlanes* coefficientPlanes_;

        /**
         * Holds the setting that determines if accumulation math is signed or not.
         */
        SignMode  pixelSignMode_;

        /**
         * The current setting for a Full Scaler. For example, if it was 256, then a value
         * of 256 is full on, a value of 128 is half on, etc. Must be a power of two.
         */
        uint16_t  fullScaler_;

        /**
         * The fullScaler is used to set this shift amount, which is used once an entire stack of
         * pixels is accumulated.
         */
        size_t    accumulatorShifter_;

        /**
         * The CostModel for this display.
         */
        CostModel* costModel;

        /**
         * Tracks whether verbose output should be used or not.
         */
        bool quiet_;

        /**
         * Used to indicate that a render was suppressed.
         */
        bool changed_;

    public:
        BaseNddiDisplay();
        BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                        int numCoefficientPlanes, int inputVectorSize);
        BaseNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                        int displayWidth, int displayHeight,
                        int numCoefficientPlanes, int inputVectorSize);
        ~BaseNddiDisplay();
        int DisplayWidth();
        int DisplayHeight();
        int NumCoefficientPlanes();
        void PutPixel(Pixel p, vector<unsigned int> &location);
        void CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end);
        void CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end);
        void CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size);
        void FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end);
        void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest);
        void UpdateInputVector(vector<int> &input);
        void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location);
        void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &start, vector<unsigned int> &end);
        void FillCoefficient(int coefficient, int row, int col, vector<unsigned int> &start, vector<unsigned int> &end);
        void FillCoefficientTiles(vector<int> &coefficients, vector<vector<unsigned int> > &positions, vector<vector<unsigned int> > &starts, vector<unsigned int> &size);
        void FillScaler(Scaler scaler, vector<unsigned int> &start, vector<unsigned int> &end);
        void FillScalerTiles(vector<uint64_t> &scalers, vector<vector<unsigned int> > &starts, vector<unsigned int> &size);
        void FillScalerTileStack(vector<uint64_t> &scalers, vector<unsigned int> &start, vector<unsigned int> &size);
        void SetPixelByteSignMode(SignMode mode);
        void SetFullScaler(uint16_t scaler);
        uint16_t GetFullScaler() { return fullScaler_; }
        void Mute() { quiet_ = true; }
        void Unmute() { quiet_ = false; }
        CostModel* GetCostModel();
    };

}

#endif // BASE_NDDI_DISPLAY_H
