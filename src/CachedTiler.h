#ifndef CACHED_TILER_H
#define CACHED_TILER_H
/*
 *  CachedTiler.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 10/28/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include "Tiler.h"
#include "GlNddiDisplay.h"
#include "ClNddiDisplay.h"


using namespace std;
using namespace nddi;

/**
 * This struct holds the checksum of a tile as well as the zIndex into the frame volume
 * when the tile is used in a caching configuration.
 */
typedef struct {
	unsigned long checksum;
	size_t zIndex;
} tile_t;


/**
 * This tiler will split provided frames into tiles and update the NDDI display. It organizes the
 * display into dimensions matching the tile size in the x and y directions and then the cache size
 * in the z direction. It will then update the coefficient plane so that a region of the display
 * maps to the appropriate tile in the cache.
 */
class CachedTiler : public Tiler {
	
public:
	/**
	 * The CachedTiler is created based on the dimensions of the NDDI display that's passed in. If those
	 * dimensions change, then the CachedTiler should be destroyed and re-created.
	 *
	 * @param display A pointer to the NDDI display
	 * @param tile_width The width of the tiles
	 * @param tile_height The height of the tiles
	 * @param max_tiles The maximum number of tiles in the cache
	 * @param bits The number of most significant bits to use when computing checksums for a tile match
	 */
	CachedTiler(GlNddiDisplay* display,
				size_t tile_width, size_t tile_height,
				size_t max_tiles, size_t bits,
				bool quiet);
	
	~CachedTiler() {
		tile_cache_.clear();
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
	
	int IsTileInCache(tile_t tile);
	bool IsTileInMap(tile_t tile);
	int GetExpiredCacheTile();
	void InitializeCoefficientMatrices();
#ifndef USE_COPY_PIXEL_TILES
	void UpdateFrameVolume(Pixel* pixels, tile_t tile);
	void UpdateCoefficientMatrices(size_t x, size_t y, tile_t tile);
#else
	void PushTile(tile_t tile, Pixel* pixels, size_t i, size_t j);
	void PushTile(tile_t tile, Pixel* pixels);
	void PushTile(tile_t tile, size_t i, size_t j);
#endif
	
	GlNddiDisplay*                 display_;
	size_t                         tile_width_, tile_height_, max_tiles_;
	size_t                         tile_map_width_, tile_map_height_;
	size_t                         bits_;
	bool                           quiet_;
	
	vector<tile_t>                 tile_cache_;
	vector< vector<tile_t> >       tile_map_;
	
	int                            unchanged_tiles_, cache_hits_, cache_misses_;

#ifdef USE_COPY_PIXEL_TILES
	vector<Pixel *>                tile_pixels_list;
	vector<vector<unsigned int> >  tile_starts_list;

	vector<int>                    coefficients_list;
	vector<vector<unsigned int> >  coefficient_positions_list;
	vector<vector<unsigned int> >  coefficient_plane_starts_list;
#endif
};
#endif // CACHED_TILER_H
