#ifndef BASE_NDDI_DISPLAY_H
#define BASE_NDDI_DISPLAY_H

#include "NDimensionalDisplayInterface.h"
#include "InputVector.h"
#include "CoefficientPlane.h"
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
        int displayWidth_;
        
        /**
         * Holds the Display Height
         */
        int displayHeight_;
        
        /**
         * The input vector
         */
        InputVector * inputVector_;
        
        /**
         * Holds the sizes of each dimenion of the Frame Volume.
         */
        std::vector<unsigned int> frameVolumeDimensionalSizes_;
        
        /**
         * The frameVolume_ is physcially a flat buffer of Pixels, but it logically managed based on the configured
         * dimensions in the constructor.
         */
        FrameVolume * frameVolume_;
        
        /**
         * The coefficientPlane_ is a flat buffer of coefficient matrices that are logically
         * arranged into rows and columns maching the dimenions of the display.
         */
        CoefficientPlane* coefficientPlane_;
        
        /**
         * The frameBuffer_ holds the rendered pixels.
         */
        Pixel* frameBuffer_;
        
        /**
         * The CostModel for this display.
         */
        CostModel* costModel;
        
    public:
        BaseNddiDisplay() {}
        BaseNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                        int inputVectorSize);
        BaseNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                        int displayWidth, int displayHeight,
                        int inputVectorSize);
        ~BaseNddiDisplay();
        int DisplayWidth();
        int DisplayHeight();
        void PutPixel(Pixel p, std::vector<unsigned int> location);
        void CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end);
        void CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end);
        void FillPixel(Pixel p, std::vector<unsigned int> start, std::vector<unsigned int> end);
        void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest);
        void UpdateInputVector(std::vector<int> input);
        void PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> location);
        void FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> start, std::vector<unsigned int> end);
        void FillCoefficient(int coefficient, int row, int col, std::vector<unsigned int> start, std::vector<unsigned int> end);
        CostModel* GetCostModel();
    };
    
}

#endif // BASE_NDDI_DISPLAY_H
