//
//  Player.h
//  pixelbridge
//
//  Created by Dave Estes on 8/26/13.
//  Copyright (c) 2013 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_Player_h
#define pixelbridge_Player_h

#include <stdlib.h>
#include <stdint.h>

class Player {
public:
    
	/**
	 * Returns the width of the video.
	 *
	 * @returns The width of the video
	 */
	virtual size_t width() = 0;
    
	/**
	 * Returns the height of the video.
	 *
	 * @returns The height of the video
	 */
	virtual size_t height() = 0;
    
	/**
	 * Decodes a frame.
	 *
	 * @returns A pointer to the decoded frame.
	 */
	virtual uint8_t* decodeFrame() = 0;
};


#endif
