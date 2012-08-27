/*
 *  FfmpegPlayer.cpp
 *  pixelbridge
 *
 *  Created by Dave Estes on 10/23/10.
 *  Copyright 2010 Qualcomm. All rights reserved.
 *
 */
#include "FfmpegPlayer.h"

FfmpegPlayer::FfmpegPlayer(const char* fileName) : fileName_(fileName), width_(0), height_(0) {

    // Register all of the codecs
    av_register_all();

    // Allocate a context
    pFormatCtx_ = avformat_alloc_context();

    // Open video file
    if (avformat_open_input(&pFormatCtx_, fileName_, NULL, NULL) != 0)
	return; // Couldn't open file

    // Retrieve stream information
    if (avformat_find_stream_info(pFormatCtx_, NULL) < 0)
	return; // Couldn't find stream information

    // Dump information about file onto standard error
    av_dump_format(pFormatCtx_, 0, fileName_, false);

    // Find the first video stream
    videoStream_ = -1;
    for (int i = 0; i < pFormatCtx_->nb_streams; i++) {
	if (pFormatCtx_->streams[i]->codec->codec_type == AVMEDIA_TYPE_VIDEO) {
	    videoStream_ = i;
	    break;
	}
    }
    if (videoStream_ == -1)
        return; // Didn't find a video stream

    // Get a pointer to the codec context for the video stream
    pCodecCtx_ = pFormatCtx_->streams[videoStream_]->codec;
    width_ = pCodecCtx_->width;
    height_ = pCodecCtx_->height;

    // Find the decoder for the video stream
    pCodec_ = avcodec_find_decoder(pCodecCtx_->codec_id);
    if (pCodec_ == NULL)
        return; // Codec not found

    // Inform the codec that we can handle truncated bitstreams -- i.e.,
    // bitstreams where frame boundaries can fall in the middle of packets
    if (pCodec_->capabilities & CODEC_CAP_TRUNCATED)
        pCodecCtx_->flags |= CODEC_FLAG_TRUNCATED;

    // Open codec
    if (avcodec_open2(pCodecCtx_, pCodec_, NULL) < 0)
        return; // Could not open codec

    // Allocate the video frame that we'll be decoding into
    pFrame_ = avcodec_alloc_frame();

    // Allocate the video frame that we'll be color-converting into
    pFrameRGB_ = avcodec_alloc_frame();
    if (pFrameRGB_ == NULL)
        return;

    // Determine required buffer size and allocate buffers
    numBytes_ = avpicture_get_size(PIX_FMT_RGB24, pCodecCtx_->width, pCodecCtx_->height);
    buffer_ = new uint8_t[numBytes_];

    // Assign appropriate parts of buffer to image planes in pFrameRGB
    avpicture_fill((AVPicture *)pFrameRGB_, buffer_, PIX_FMT_RGB24,
		   pCodecCtx_->width, pCodecCtx_->height);

    // Set up the scaling and color conversion context
    pSwsCtx_ = sws_getContext(pCodecCtx_->width, pCodecCtx_->height, pCodecCtx_->pix_fmt,
			      pCodecCtx_->width, pCodecCtx_->height, PIX_FMT_RGB24,
			      SWS_POINT, NULL, NULL, NULL);
    // TODO(cdestes): Try to use this new method below and get rid of the deprecated call above.
    //pSwsCtx_ = sws_alloc_context();
    //sws_init_context(pSwsCtx_, NULL, NULL);
}

FfmpegPlayer::~FfmpegPlayer() {

    // Free the RGB image
    delete [] buffer_;
    av_free(pFrameRGB_);

    // Free the YUV frame
    av_free(pFrame_);

    // Close the codec
    avcodec_close(pCodecCtx_);

    // Close the video file
    avformat_close_input(&pFormatCtx_);

    // Free the SWS context
    sws_freeContext(pSwsCtx_);

}

size_t FfmpegPlayer::width() {
    return width_;
}

size_t FfmpegPlayer::height() {
    return height_;
}

// TODO(cdestes): Get rid of the nasty gotos from the sample code.
uint8_t* FfmpegPlayer::decodeFrame() {

    static AVPacket packet;
    static int      bytesRemaining = 0;
    static uint8_t  *rawData;
    static bool     fFirstTime = true;
    int             bytesDecoded;
    int             frameFinished = 0;

    // First time we're called, set packet.data to NULL to indicate it
    // doesn't have to be freed
    if (fFirstTime) {
        fFirstTime = false;
        packet.data = NULL;
    }

    // Decode packets until we have decoded a complete frame
    while (true) {
        // Work on the current packet until we have decoded all of it
        while (bytesRemaining > 0) {
	    // Decode the next chunk of data
	    // avcodec_decode_video() has been deprecated
            //bytesDecoded = avcodec_decode_video(pCodecCtx_, pFrame_,
	    //                                    &frameFinished, rawData, bytesRemaining);
            bytesDecoded = avcodec_decode_video2(pCodecCtx_, pFrame_, &frameFinished, &packet);
	    // Was there an error?
            if (bytesDecoded < 0) {
                fprintf(stderr, "Error while decoding frame\n");
		return NULL;
            }

            bytesRemaining -= bytesDecoded;
            rawData += bytesDecoded;

            // Did we finish the current frame? Then we can return
            if (frameFinished)
                goto function_exit;
        }

        // Read the next packet, skipping all packets that aren't for this
        // stream
        do {
            // Free old packet
            if (packet.data != NULL)
                av_free_packet(&packet);

            // Read new packet
            if (av_read_frame(pFormatCtx_, &packet) < 0)
                goto loop_exit;
        } while(packet.stream_index != videoStream_);

        bytesRemaining = packet.size;
        rawData = packet.data;
    }

  loop_exit:

    // Decode the rest of the last frame
    // avcodec_decode_video() has been deprecated
    //bytesDecoded = avcodec_decode_video(pCodecCtx_, pFrame_, &frameFinished,
    //                                    rawData, bytesRemaining);
    bytesDecoded = avcodec_decode_video2(pCodecCtx_, pFrame_, &frameFinished, &packet);

    // Free last packet
    if (packet.data != NULL)
        av_free_packet(&packet);

  function_exit:

    if (frameFinished == 0) {
	return NULL;
    } else {
	// Scale/Color-Convert
	sws_scale(pSwsCtx_, pFrame_->data, pFrame_->linesize,
		  0, pCodecCtx_->height,
		  pFrameRGB_->data, pFrameRGB_->linesize);

	return buffer_;
    }
}
