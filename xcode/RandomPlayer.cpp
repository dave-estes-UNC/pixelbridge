//
//  RandomPlayer.cpp
//  pixelbridge
//
//  Created by Dave Estes on 8/26/13.
//  Copyright (c) 2013 Dave Estes. All rights reserved.
//

#include "RandomPlayer.h"
#include <time.h>
#include <string.h>


void RandomPlayer::createRandomFrame() {
    srand((unsigned int)time(NULL));

    memset(buffer_, 0x00, numBytes_);
    for (size_t i = 0; i < width_ * height_; i++) {
        // Set the red channel only
        buffer_[i * VIDEO_PIXEL_SIZE] = rand() % 265;
    }
}

RandomPlayer::RandomPlayer() : width_(FIXED_WIDTH), height_(FIXED_HEIGHT), framesRemaining_(FIXED_FRAME_COUNT) {
    
    // Determine required buffer size and allocate buffers
    numBytes_ = VIDEO_PIXEL_SIZE * FIXED_WIDTH * FIXED_HEIGHT;
    buffer_ = new uint8_t[numBytes_];
    
    // Initialize the buffer with random pixels
    createRandomFrame();
}

RandomPlayer::~RandomPlayer() {
    
    // Free the RGB image
    delete [] buffer_;
}

size_t RandomPlayer::width() {
    return width_;
}

size_t RandomPlayer::height() {
    return height_;
}

uint8_t* RandomPlayer::decodeFrame() {
    return buffer_;
}
