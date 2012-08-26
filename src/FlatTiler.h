#ifndef FLAT_TILER_H
#define FLAT_TILER_H
/*
 *  FlatTiler.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 11/26/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include "Tiler.h"
#include "GlNddiDisplay.h"


/**
 * This tiler will split provided frames into tiles and update the NDDI display. It organizes the
 * display into dimensions matching the tile size in the x and y directions. If a tile changes, then it
 * is updated in the frame volume.
 */
class FlatTiler : public Tiler {
	
public:
	/**
	 * The FlatTiler is created based on the dimensions of the NDDI display that's passed in. If those
	 * dimensions change, then the CachedTiler should be destroyed and re-created.
	 *
	 * @param display A pointer to the NDDI display
	 * @param tile_width The width of the tiles
	 * @param tile_height The height of the tiles
	 * @param bits The number of most significant bits to use when computing checksums for a tile match
	 */
	FlatTiler(GlNddiDisplay* display,
			  size_t tile_width, size_t tile_height,
			  size_t bits,
			  bool quiet);
	
	~FlatTiler() {
		tile_map_.clear();
	}
	
    /**
     * Intializes the Coefficient Plane for this tiler.
     *
	 * @return The cost of this operation, including all of the NDDI operations
     */
    void InitializeCoefficientPlane();
    
	/**
	 * Update the tile_map, tilecache, and then the NDDI display based on the frame that's passed in.
	 *
	 * @param buffer Pointer to the return frame buffer
	 * @param width The width of that frame buffer
	 * @param height The height of that frame buffer
	 * @return The cost of this operation, including all of the NDDI operations
	 */
	void UpdateDisplay(uint8_t* buffer, size_t width, size_t height);
	
private:
	void UpdateFrameVolume(nddi::Pixel* pixels, int i_map, int j_map);
	
	GlNddiDisplay* display_;
	size_t tile_width_, tile_height_;
	size_t tile_map_width_, tile_map_height_;
	size_t bits_;
	bool quiet_;
	
	std::vector< std::vector<unsigned long> > tile_map_;

	int unchanged_tiles_, tile_updates_;
};
#endif // FLAT_TILER_H