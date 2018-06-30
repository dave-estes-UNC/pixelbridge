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
#include "nddi/GlNddiDisplay.h"
#ifdef USE_CL
#include "nddi/ClNddiDisplay.h"
#endif

using namespace nddi;
using namespace std;

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
     * @param display_width The width of the display
     * @param display_height The height of the display
	 * @param tile_width The width of the tiles
	 * @param tile_height The height of the tiles
	 * @param bits The number of most significant bits to use when computing checksums for a tile match
	 */
	FlatTiler(size_t display_width, size_t display_height,
			  size_t tile_width, size_t tile_height,
			  size_t bits);

	~FlatTiler() {
		tile_map_.clear();
	}

    /**
     * Returns the Display created and initialized by the tiler.
     */
    virtual GlNddiDisplay* GetDisplay();

	/**
	 * Update the tile_map, tilecache, and then the NDDI display based on the frame that's passed in.
	 *
	 * @param buffer Pointer to the return frame buffer
	 * @param width The width of that frame buffer
	 * @param height The height of that frame buffer
	 */
	void UpdateDisplay(uint8_t* buffer, size_t width, size_t height);

private:
    void InitializeCoefficientPlanes();
#ifndef USE_COPY_PIXEL_TILES
	void UpdateFrameVolume(Pixel* pixels, int i_map, int j_map);
#endif

	GlNddiDisplay*    display_;
	size_t            tile_width_, tile_height_;
	size_t            tile_map_width_, tile_map_height_;
	size_t            bits_;
	bool              quiet_;

	vector< vector<unsigned long> > tile_map_;

	int unchanged_tiles_, tile_updates_;
};
#endif // FLAT_TILER_H
