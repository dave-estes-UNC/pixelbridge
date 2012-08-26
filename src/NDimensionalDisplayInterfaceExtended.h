#ifndef N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
#define N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H

#include <vector>

namespace nddi {
    
    /**
     * This abstract class serves as an extension to the NDDI software interface.
     */
    class NDimensionalDisplayInterfaceExtended  {
        
    public:

        /**
         * Copies from one region of the frame volume to another with blending.
         *
         * @param start The first pixel in the frame volume to be copied from.
         * @param end The last pixel in the frame volume to be copied from.
         * @param dest The first pixel in the frame volume to be copied to.
         * @return The cost of the operation. Can be measured in time, byte-count, or another
         *         measurements based on the display implementation
         */
        virtual void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest, bool blend) = 0;
    };
    
} // namespace nddi {

#endif // N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
