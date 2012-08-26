#ifndef TILER_H
#define TILER_H

/*
 *  Tiler.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 11/24/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */


/**
 * The Tiler base class is used to update an attached NDDI display. It's derivative classes
 * will break up each frame into tiles and update the tiles that are changed.
 */
class Tiler {
	
public:
    
    /**
     * Intializes the Coefficient Plane for this tiler.
     */
    virtual void InitializeCoefficientPlane() = 0;
    
	/**
	 * Update the tile_map, tilecache, and then the NDDI display based on the frame that's passed in.
	 */
	virtual void UpdateDisplay(uint8_t* buffer, size_t width, size_t height) = 0;
	
};
#endif // TILER_H