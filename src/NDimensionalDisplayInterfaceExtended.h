#ifndef N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
#define N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H

#include <vector>

using namespace std;

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
    	virtual void CopyFrameVolume(vector<unsigned int> start, vector<unsigned int> end, vector<unsigned int> dest, bool blend) = 0;

        /**
         * Copies the array of pixels into the designated tile regions of the frame volume. The data must be
         * arranged in the array with strides for each dimension of the area. Tiles are not necessarily two-
         * dimensional and can have any dimensionality at most the dimensionality of the frame volume.
         *
         * @param p The pointer to the pixel values to be copied.
         * @param starts Vector holding series of first pixel for each destination tile in the frame volume.
         * @param size The size of each tile. Dimensionality should match that of the individual start vectors.
         * @return The cost of the operation. Can be measured in time, byte-count, or another
         *         measurements based on the display implementation
         */
        virtual void CopyPixelTiles(vector<Pixel*> p, vector<vector<unsigned int> > starts, vector<unsigned int> size) = 0;

        /**
         * For each coefficient, positions, and start; copies the coefficient to the position
         * in the in each coefficient matrix in the tile specified by the start and size.
         *
         * @param coefficients The buffer of coefficients.
         * @param positions The position (row, col) to place the coefficient within the coefficient matrix.
         * @param starts The location (x, y) of the start of the tile in the coefficient plane.
         * @param size The size (w, h) of the tile.
         */
        virtual void FillCoefficientTiles(vector<int> coefficients, vector<vector<unsigned int> > positions, vector<vector<unsigned int> > starts, vector<unsigned int> size) = 0;
};
    
} // namespace nddi {

#endif // N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
