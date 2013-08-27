#ifndef FFMPEGPLAYER_H
#define FFMPEGPLAYER_H
/*
 *  FfmpegPlayer.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 10/23/10.
 *  Copyright Dave Estes. All rights reserved.
 *
 */

#include "Player.h"

extern "C" {
#include "libavcodec/avcodec.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
}

#define VIDEO_PIXEL_SIZE 3 // RGB


/**
 * This class will utilize an FFMPEG demuxer, decoder, and color convertor to decode
 * each frame of the provided video.
 */
class FfmpegPlayer : Player {
public:
	/**
	 * The constructor for the player.
	 *
	 * @param fileName The filename of the file to be decoded.
	 */
	FfmpegPlayer(const char* fileName);

	~FfmpegPlayer();

	/**
	 * Returns the width of the video.
	 *
	 * @returns The width of the video
	 */
	size_t width();

	/**
	 * Returns the height of the video.
	 *
	 * @returns The height of the video
	 */
	size_t height();

	/**
	 * Decodes a frame.
	 *
	 * @returns A pointer to the decoded frame.
	 */
	uint8_t* decodeFrame();

private:
	const char* fileName_;
	size_t width_;
	size_t height_;

	AVFormatContext* pFormatCtx_;
	int              videoStream_;
	AVCodecContext*  pCodecCtx_;
	AVCodec*         pCodec_;
	AVFrame*         pFrame_;
	AVFrame*         pFrameRGB_;
	int              numBytes_;
	uint8_t*         buffer_;
	SwsContext*      pSwsCtx_;
};

#endif // FFMPEGPLAYER_H
