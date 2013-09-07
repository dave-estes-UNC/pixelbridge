#ifndef N_DIMENSIONAL_DISPLAY_INTERFACE_H
#define N_DIMENSIONAL_DISPLAY_INTERFACE_H

#include <vector>
#include "CostModel.h"

using namespace std;

namespace nddi {

    /**
     * When specifying coefficient matrices for purposes of updating the coefficient plane, the NDDI client
     * can use this value for one or more of the elements in the matrix if they would like the element in the
     * same location of the destination coefficient matrix to remain unchanged.
     */
    #define COFFICIENT_UNCHANGED INT_MAX

	/**
	 * The number of coefficient planes is fixed at 256.
	 */
	#define NUM_COEFFICIENT_PLANES 256

    /**
     * Struct representing an RGBA 32-bit pixel.
     */
    typedef union {
        struct {
            unsigned char r;
            unsigned char g;
            unsigned char b;
            unsigned char a;
        };
        unsigned int packed;
    } Pixel;
    
    /**
     * Options for pixel byte sign mode.
     */
    typedef enum {
        UNSIGNED_MODE,
        SIGNED_MODE
    } SignMode;

    /**
     * This abstract class serves as a software interface to an n-Dimensional Display Interface (NDDI) compliant
     * display device. Implementations of this interface may work with a simulated display that writes
     * to a system framebuffer, an embedded display that works with a device driver, or a remote display
     * connected via an IP-based socket and communicating with an NDDI specific protocol.
     *
     */
    class NDimensionalDisplayInterface  {

    public:

        /**
         * Required default constructor for abstract class NDimensionalDisplayInterface
         */
        NDimensionalDisplayInterface() {}

        /**
         * Each NDimensionalDisplayInterface is configured during construction.
         *
         * @param frameVolumeDimensionalSizes This vector is used to configure the frame volume.
         *                                    Each element in the vector represents a dimension and that element's
         *                                    value represents the size of that dimension. e.g. a simple 4x4 2D
         *                                    frame volume will be configured with a two-element vector with 4 and 4 in it.
         * @param inputVectorSize Used to configure the size of the input vector. It must be greater than or equal to two.
         */
        NDimensionalDisplayInterface(vector<unsigned int> &frameVolumeDimensionalSizes,
                                     int inputVectorSize) {}
        /**
         * Each NDimensionalDisplayInterface is configured during contruction. This contructor allows the NDDI client
         * to reduce the size of the displayable area.
         *
         * @param frameVolumeDimensionalSizes This vector is used to configure the frame volume.
         *                                    Each element in the vector represents a dimension and that element's
         *                                    value represents the size of that dimension. e.g. a simple 4x4 2D
         *                                    frame volume will be configured with a two-element vector with 4 and 4 in it.
         * @param displayWidth Used to configure the width of the diplay if it is less than the display device.
         * @param displayHeight Used to configure the width of the diplay if it is less than the display device.
         * @param inputVectorSize Used to configure the size of the input vector. It must be greater than or equal to two.
         */
        NDimensionalDisplayInterface(vector<unsigned int> &frameVolumeDimensionalSizes,
                                     int displayWidth, int displayHeight,
                                     int inputVectorSize) {}

        /**
         * Used to query the display width.
         *
         * @return The width of the display.
         */
        virtual int DisplayWidth() = 0;

        /**
         * Used to query the display height.
         *
         * @return The height of the display.
         */
        virtual int DisplayHeight() = 0;

        /**
         * Copies the provided pixel to the specified location.
         *
         * @param p The pixel value to be copied.
         * @param location The location within the frame volume where the pixel will be copied to.
         */
        virtual void PutPixel(Pixel p, vector<unsigned int> &location) = 0;

        /**
         * Copies the one dimensional array of pixels along a particular dimension in the frame volume. In a
         * two-dimensional frame volume, this can be thought of as a way to copy along a row or along a column, but
         * not both since the input pixels are only one-dimensional.
         *
         * @param p The pointer to the pixel values to be copied.
         * @param start The first pixel in the frame volume to be filled.
         * @param end The last pixel in the frame volume to be filled. All but one of the values in
         *            values in this last pixel should be identical to the start pixel.
         */
        virtual void CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * Copies the array of pixels into the designated region of the frame volume. The data must be
         * arranged in the array with strides for each dimension of the area. So to copy pixels into a
         * 2 x 2 x 2 region in the frame volume, the array must be arranged accordingly:
         * (0,0,0) (1,0,0) (0,1,0) (1,1,0) (0,0,1) (1,0,1) (0,1,1) (1,1,1)
         *
         * @param p The pointer to the pixel values to be copied.
         * @param start The first pixel in the frame volume to be filled.
         * @param end The last pixel in the frame volume to be filled.
         */
        virtual void CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * Copies the array of pixels into the designated tile regions of the frame volume. The data must be
         * arranged in the array with strides for each dimension of the area. Only 2D tiles are supported.
         *
         * @param p The pointer to the pixel values to be copied.
         * @param starts Vector holding series of first pixel for each destination tile in the frame volume.
         * @param size The size of each tile (w, h).
         */
        virtual void CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;

        /**
         * Fills the frame volume with the specified pixel. It can fill in multiple
         * dimensions by starting at the start pixel and filling in each dimension until
         * the end pixel value is reached.
         *
         * @param p The pixel value to be filled.
         * @param start The first pixel in the frame volume to be filled.
         * @param end The last pixel in the frame volume to be filled.
         */
        virtual void FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * Copies pixels from one multi-dimensional region of the frame volume to another region.
         *
         * @param start The starting coordinate of the source region.
         * @param end The ending coordinate of the source region.
         * @param dest The starting coordinate of the destination region.
         */
        virtual void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest) = 0;

        /**
         * Used to update the input vector with the extra values in the input vector.
         *
         * @param input The values to use for the update. The length of this
         *              vector must equal the size of the actual input vector
         *              minus two, since the first two values in the input
         *              vector cannot be changed.
         */
        virtual void UpdateInputVector(vector<int> &input) = 0;

        /**
         * Used to copy the specified coefficientMatrix into the specified location of the coefficient
         * planes.
         *
         * @param coefficientMatrix This two-dimensional vector holds the matrix to be copied.
         *                          It's size must match the configuration of the coefficient matrices
         *                          exactly. Can use COFFICIENT_UNCHANGED for one or more elements.
         * @param location This two-element vector specifies the location in the coefficient plane where the provided
         *                 coefficient matrix will be copied.
         */
        virtual void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) = 0;

        /**
         * Used to copy the specified coefficientMatrix into a range of locations in the coefficient planes.
         *
         * @param coefficientMatrix This two-dimensional vector holds the matrix to be copied.
         *                          It's size must match the configuration of the coefficient matrices
         *                          exactly. Can use COFFICIENT_UNCHANGED for one or more elements.
         * @param start This three-element vector specifies the location in the coefficient planes where the first
         *              coefficient matrix will be copied to.
         * @param end This three-element vector specifies the location in the coefficient planes where the last
         *            coefficient matrix will be copied to.
         */
        virtual void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * Used to copy the specified single coefficient value from a matrix into a range of locations in the coefficient planes.
         *
         * @param coefficient This single value will be placed into each coefficient at the specified location in the coefficient
         *                    matrices of the specified range.
         * @param row The row of the coefficient to be updated in the coefficient matrix.
         * @param col The column of the coefficient to be updated in the coefficient matrix.
         * @param start This three-element vector specifies the location in the coefficient planes where the first
         *              coefficient matrix will be copied to.
         * @param end This three-element vector specifies the location in the coefficient planes where the last
         *            coefficient matrix will be copied to.
         */
        virtual void FillCoefficient(int coefficient, int row, int col, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * For each coefficient, positions, and start; copies the coefficient to the position
         * in the in each coefficient matrix in the 2D tile specified by the start and size.
         *
         * @param coefficients The buffer of coefficients.
         * @param positions The position (row, col) to place the coefficient within the coefficient matrix.
         * @param starts The location (x, y) of the start of the tile in the coefficient planes.
         * @param size The size (w, h) of the tile.
         */
        virtual void FillCoefficientTiles(vector<int> &coefficients, vector<vector<unsigned int> > &positions, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;

        /**
         * Used to copy the specified scaler to a range of locations in the coefficient planes.
         *
         * @param scaler This single scaler will be copied to each location in the range of coefficient planes.
         * @param start This three-element vector specifies the location in the coefficient planes where the
         *              scaler will be copied to.
         * @param end This three-element vector specifies the location in the coefficient planes where the
         *            scalers will be copied to.
         */
        virtual void FillScaler(int scaler, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * Used to copy the specified scalers to a series of 2D ranges of locations (tiles) in the coefficient planes.
         *
         * @param scaler This single scaler will be copied to each location in the range of coefficient planes.
         * @param starts The location (x, y) of the start of the tile in the coefficient planes.
         * @param size The size (w, h) of the tile.
         */
        virtual void FillScalerTiles(vector<int> &scalers, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;
        
        /**
         * Allows the bytes of pixel values to be interpretted as signed values when scaling, accumulating, and clamping
         * in the pixel blending pipeline.
         *
         */
        virtual void SetPixelByteSignMode(SignMode mode) = 0;

        /**
         * Returns the CostModel for this display. The CostModel can be queried by the
         * host application to understand the cost of operations after they complete.
         *
         * @return The CostModel created and maintained by this display.
         */
        virtual CostModel* GetCostModel() = 0;
    };

} // namespace nddi {

#endif // N_DIMENSIONAL_DISPLAY_INTERFACE_H
