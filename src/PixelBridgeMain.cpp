#include <iostream>
#include <sys/time.h>
#include <cmath>
#include <assert.h>
#include <queue>
#include <pthread.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <opencv/cv.hpp>

#include "PixelBridgeFeatures.h"
#include "Configuration.h"
#include "CostModel.h"
#include "GlNddiDisplay.h"
#include "ClNddiDisplay.h"
#include "BlendingGlNddiDisplay.h"
#include "CachedTiler.h"
#include "DctTiler.h"
#include "ItTiler.h"
#include "FlatTiler.h"
#include "FfmpegPlayer.h"
#include "RandomPlayer.h"
#include "Rewinder.h"


typedef enum {
    NOT_REWINDING,  // Normal playback, not rewinding
    PLAY_BACKWARDS, // Rewinding frame by frame backwards
    PLAY_FORWARD    // Finished rewinding, playing forward from memory
} rewind_play_t;

// General Globals
size_t displayWidth = 40, displayHeight = 32;
const char* fileName = NULL;
Configuration globalConfiguration = Configuration();

// Helper Objects
Player*  myPlayer;
GlNddiDisplay* myDisplay;
BlendingGlNddiDisplay* myBlendingDisplay;
Tiler* myTiler;
Rewinder* myRewinder = NULL;

// Decoder thread
pthread_t           decoderThread;

// Statistical Instrumentation
CostModel* costModel;
int totalUpdates = 0;
timeval startTime, endTime; // Used for timing data
long totalSE = 0; // Cumulative Square Error used to calculate MSE and then PSNR

// Stores the current and previously decoded frame
uint8_t* videoBuffer = NULL;
uint8_t* lastBuffer = NULL;

#ifdef USE_ASYNC_DECODER
// Buffers decoded frame
queue<uint8_t*> bufferQueue;

// Synchronization variables
pthread_mutex_t      bufferQueueMutex;
#endif

void setupDisplay() {

    Scaler s;

    // Cached-Tiled
    if (globalConfiguration.tiler == CACHE) {

        // Setup Cached Tiler which initializies Coefficient Plane
        myTiler = new CachedTiler(displayWidth, displayHeight,
                globalConfiguration.tileWidth, globalConfiguration.tileHeight,
                globalConfiguration.maxTiles, globalConfiguration.sigBits);

        // Grab the display and cost model
        myDisplay = myTiler->GetDisplay();
        costModel = myDisplay->GetCostModel();

    // DCT-Tiled
    } else if (globalConfiguration.tiler == DCT) {

        // Setup DCT Tiler and initializes Coefficient Plane and Frame Volume
        myTiler = new DctTiler(displayWidth, displayHeight,
                               globalConfiguration.quality);

        // Grab the display and cost model
        myDisplay = myTiler->GetDisplay();
        costModel = myDisplay->GetCostModel();

    // IT-Tiled
    } else if (globalConfiguration.tiler == IT) {

        // Setup IT Tiler and initializes Coefficient Plane and Frame Volume
        myTiler = new ItTiler(displayWidth, displayHeight,
                              globalConfiguration.quality);

        // Grab the display and cost model
        myDisplay = myTiler->GetDisplay();
        costModel = myDisplay->GetCostModel();

    // Frame Volume Blending
    } else if ((globalConfiguration.tiler == SIMPLE) && (globalConfiguration.blend == FRAME_VOLUME)) {

        // 3 dimensional matching the Video Width x Height x 3
        vector<unsigned int> fvDimensions;
        fvDimensions.push_back(displayWidth);
        fvDimensions.push_back(displayHeight);
        fvDimensions.push_back(3);

#ifndef NO_CL
        assert(0 && "FrameVolume blending mode not supported in the OpenCL version of NDDI display.");
#else
        myBlendingDisplay = new BlendingGlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                                      displayWidth, displayHeight, // display size
                                                      1,                           // number of coefficient planes
                                                      3);                            // input vector size (x, y, 1)
#endif
        myDisplay = myBlendingDisplay;

        // Grab the cost model
        costModel = myDisplay->GetCostModel();

        // Initialize Input Vector
        vector<int> iv;
        iv.push_back(0);
        myDisplay->UpdateInputVector(iv);

        // Initialize Frame Volume
        nddi::Pixel p;
        p.r = p.g = p.b = 0x00; p.a = 0x7f;
        vector<unsigned int> start, end;
        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(2);
        myDisplay->FillPixel(p, start, end);

        // Setup Coefficient Plane
        vector< vector<int> > coeffs;
        coeffs.resize(3);
        coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
        coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
        coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(1);

        start.clear(); end.clear();

        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);

        myDisplay->FillCoefficientMatrix(coeffs, start, end);

        // Turn off all planes and then set the 0 plane to full on.
        end[2] = myDisplay->NumCoefficientPlanes() - 1;
        s.packed = 0;
        myDisplay->FillScaler(s, start, end);
        end[2] = 0;
        s.r = s.g = s.b = myDisplay->GetFullScaler();
        myDisplay->FillScaler(s, start, end);

    // Temporal Blending
    } else if ((globalConfiguration.tiler == SIMPLE) && (globalConfiguration.blend == TEMPORAL)) {

        // 3 dimensional matching the Video Width x Height x 2
        vector<unsigned int> fvDimensions;
        fvDimensions.push_back(displayWidth);
        fvDimensions.push_back(displayHeight);
        fvDimensions.push_back(2);

#ifndef NO_CL
        myDisplay = new ClNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      1,                           // number of coefficient planes on the display
                                      3);                          // input vector size (x, y, t)
#else
        myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      1,                           // number of coefficient planes on the display
                                      3);                          // input vector size (x, y, t)
#endif

        // Grab the cost model
        costModel = myDisplay->GetCostModel();

        // Initialize Input Vector
        vector<int> iv;
        iv.push_back(0);
        myDisplay->UpdateInputVector(iv);

        // Initialize Frame Volume
        nddi::Pixel p;
        p.r = p.g = p.b = 0x00; p.a = 0x7f;
        vector<unsigned int> start, end;
        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(1);
        myDisplay->FillPixel(p, start, end);

        // Setup Coefficient Plane
        vector< vector<int> > coeffs;
        coeffs.resize(3);
        coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
        coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
        coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(1);

        start.clear(); end.clear();

        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);

        myDisplay->FillCoefficientMatrix(coeffs, start, end);

        // Turn off all planes and then set the 0 plane to full on.
        end[2] = myDisplay->NumCoefficientPlanes() - 1;
        s.packed = 0;
        myDisplay->FillScaler(s, start, end);
        end[2] = 0;
        s.r = s.g = s.b = myDisplay->GetFullScaler();
        myDisplay->FillScaler(s, start, end);

    // Coefficient Plane Blending
    } else if ((globalConfiguration.tiler == SIMPLE) && (globalConfiguration.blend == COEFFICIENT_PLANE)) {

        // 3 dimensional matching the Video Width x Height x 2
        vector<unsigned int> fvDimensions;
        fvDimensions.push_back(displayWidth);
        fvDimensions.push_back(displayHeight);
        fvDimensions.push_back(2);

#ifndef NO_CL
        myDisplay = new ClNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      2,                           // number of coefficient planes on the display
                                      3);                          // input vector size (x, y, 1)
#else
        myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      2,                           // number of coefficient planes on the display
                                      3);                          // input vector size (x, y, 1)
#endif

        // Grab the cost model
        costModel = myDisplay->GetCostModel();

        // Initialize Input Vector
        vector<int> iv;
        iv.push_back(1);
        myDisplay->UpdateInputVector(iv);

        // Initialize Frame Volume
        nddi::Pixel p;
        p.r = p.g = p.b = 0x00; p.a = 0xff;
        vector<unsigned int> start, end;
        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(1);
        myDisplay->FillPixel(p, start, end);

        // Setup Coefficient Plane 0 to pull plane 0 of the Frame Volume
        vector< vector<int> > coeffs;
        coeffs.resize(3);
        coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
        coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
        coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(0);

        start.clear(); end.clear();

        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);

        myDisplay->FillCoefficientMatrix(coeffs, start, end);

        // Then setup Coefficient Plane 1 to pull plane 1 of the Frame Volume, which will remain black
        coeffs[2][2] = 1;
        start[2] = 1; end[2] = 1;

        myDisplay->FillCoefficientMatrix(coeffs, start, end);

        // Turn off all planes and then set the 0 and 1 planes to half each.
        start[2] = 0;
        end[2] = myDisplay->NumCoefficientPlanes() - 1;
        s.packed = 0;
        myDisplay->FillScaler(s, start, end);
        end[2] = 1;
        s.r = s.g = s.b = myDisplay->GetFullScaler() >> 1;
        myDisplay->FillScaler(s, start, end);

    // Flat-Tiled
    } else if (globalConfiguration.tiler == FLAT) {

        // Set up Flat Tiler and initialize Coefficient Planes
        myTiler = new FlatTiler(displayWidth, displayHeight,
                                globalConfiguration.tileWidth,
                                globalConfiguration.tileHeight,
                                globalConfiguration.sigBits);

        // Grab the display and cost model
        myDisplay = myTiler->GetDisplay();
        costModel = myDisplay->GetCostModel();

    // Simple Framebuffer
    } else {

        // 2 dimensional matching the Video Width x Height
        vector<unsigned int> fvDimensions;
        fvDimensions.push_back(displayWidth);
        fvDimensions.push_back(displayHeight);

#ifndef NO_CL
        myDisplay = new ClNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      1,                           // number of coefficient planes on the display
                                      2);                          // input vector size (x and y only)
#else
        myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      1,                           // number of coefficient planes on the display
                                      2);                          // input vector size (x and y only)
#endif
        // Grab the cost model
        costModel = myDisplay->GetCostModel();

        // Initialize Frame Volume
        nddi::Pixel p;
        p.r = p.g = p.b = p.a = 0xff;
        vector<unsigned int> start, end;
        start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
        myDisplay->FillPixel(p, start, end);

        // Initialize Coefficient Plane
        vector< vector<int> > coeffs;
        coeffs.resize(2);
        coeffs[0].push_back(1); coeffs[0].push_back(0);
        coeffs[1].push_back(0); coeffs[1].push_back(1);

        start.clear(); end.clear();

        start.push_back(0); start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);

        myDisplay->FillCoefficientMatrix(coeffs, start, end);

        // Turn off all planes and then set the 0 plane to full on.
        end[2] = myDisplay->NumCoefficientPlanes() - 1;
        s.packed = 0;
        myDisplay->FillScaler(s, start, end);
        end[2] = 0;
        s.r = s.g = s.b = myDisplay->GetFullScaler();
        myDisplay->FillScaler(s, start, end);

    }

    if (globalConfiguration.verbose)
        myDisplay->Unmute();

#ifdef CLEAR_COST_MODEL_AFTER_SETUP
    costModel->clearCosts();
#endif

    // Renders the initial display
    if (!globalConfiguration.headless)
        myDisplay->GetFrameBufferTex();
    else
        myDisplay->SimulateRender();
}


void updateDisplay(uint8_t* buffer, size_t width, size_t height) {

    // CACHE, DCT, IT, or FLAT
    if ( (globalConfiguration.tiler == CACHE) || (globalConfiguration.tiler == DCT) || (globalConfiguration.tiler == IT) || (globalConfiguration.tiler == FLAT) ) {
        // Update the display with the Tiler
        myTiler->UpdateDisplay(buffer, width, height);
    // SIMPLE
    } else {
        // Update the display like a simple Frame Buffer
        size_t bufferPos = 0;

        // Array of pixels used for framebuffer mode.
        Pixel* frameBuffer = (Pixel*)malloc(sizeof(Pixel) * displayWidth * displayHeight);

        // Transform the buffer into pixels
        for (int j = 0; j < displayHeight; j++) {
#ifdef USE_SMALLER_WINDOW
            // Set the bufferPos to the beginning of the next row
            bufferPos = j * width * 3;
#endif
            for (int i = 0; i < displayWidth; i++) {
                frameBuffer[j * displayWidth + i].r = buffer[bufferPos++];
                frameBuffer[j * displayWidth + i].g = buffer[bufferPos++];
                frameBuffer[j * displayWidth + i].b = buffer[bufferPos++];
                frameBuffer[j * displayWidth + i].a = 0xff;
            }
        }

        // Update the frame volume
        vector<unsigned int> start, end, dest;

        // Not Blending
        if (globalConfiguration.blend == NONE) {
            // Just send the pixels to the single plane
            start.push_back(0); start.push_back(0);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
            myDisplay->CopyPixels(frameBuffer, start, end);

        // Frame Volume Blending
        } else if (globalConfiguration.blend == FRAME_VOLUME) {

            // Set the plane that we're working on. Alternates between 0 and 1
            int plane = totalUpdates & 1;

            // Send the pixels to the plane
            start.push_back(0); start.push_back(0); start.push_back(plane);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(plane);
            myDisplay->CopyPixels(frameBuffer, start, end);

            // Blend the z=2 plane over it
            start[2] = 2; end[2] = 2;
            dest.push_back(0); dest.push_back(0); dest.push_back(plane);
            myBlendingDisplay->CopyFrameVolume(start, end, dest, true);

            // Render
            vector<int> iv;
            iv.push_back(plane);
            myDisplay->UpdateInputVector(iv);
        // Temporal Blending
        } else if (globalConfiguration.blend == TEMPORAL) {

            // Send the pixels to the 0 plane
            start.push_back(0); start.push_back(0); start.push_back(0);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);
            myDisplay->CopyPixels(frameBuffer, start, end);

            // Render the new frame and then render the black frame.
            // Note: This is a loose simulation. Actual temporal blending would carefully
            // drive the input vector, perhaps flipping back and forth several times per frame.
            vector<int> iv(1);
            for (int c = 0; c < globalConfiguration.temporalFlipCountPerFrame; c++) {
                iv[0] = 1;
                myDisplay->UpdateInputVector(iv);
                iv[0] = 0;
                myDisplay->UpdateInputVector(iv);
            }
        // Coefficient Plane Blending
        } else if (globalConfiguration.blend == COEFFICIENT_PLANE) {
            // Send the pixels to the 0 plane
            start.push_back(0); start.push_back(0); start.push_back(0);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);
            myDisplay->CopyPixels(frameBuffer, start, end);
        }

        // Free the temporary frame buffer
        free(frameBuffer);
    }
    totalUpdates++;

}


long calculateSE(uint8_t* videoIn, Pixel* frameOut) {
    long errR = 0, errG = 0, errB = 0;
    size_t i, bufferPos = 0;
    long SE = 0;

    for (i = 0; i < displayWidth * displayHeight; i++) {
        errR = videoIn[bufferPos++] - frameOut[i].r;
        errG = videoIn[bufferPos++] - frameOut[i].g;
        errB = videoIn[bufferPos++] - frameOut[i].b;
        SE += errR * errR + errG * errG + errB * errB;
    }

    return SE;
}


void draw( void ) {

    if (globalConfiguration.tiler > COUNT) {

        if (globalConfiguration.headless) {

            myDisplay->SimulateRender();

        } else {

            // Grab the frame buffer from the NDDI display
            GLuint texture = myDisplay->GetFrameBufferTex();

// TODO(CDE): Temporarily putting this here until GlNddiDisplay and ClNddiDisplay
//            are using the exact same kind of GL textures
#ifndef NO_CL
            glClearColor(0.0f, 0.0f, 1.0f, 1.0f );
            glClear( GL_COLOR_BUFFER_BIT );

            // Draw the texture with ARB Rectangles
            glEnable(GL_TEXTURE_RECTANGLE_ARB);
            glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texture );

            glBegin(GL_QUADS);
            glTexCoord2f(0, 0); glVertex2f(-1.0, 1.0);
            glTexCoord2f(0, displayHeight); glVertex2f(-1.0, -1.0);
            glTexCoord2f(displayWidth, displayHeight); glVertex2f(1.0, -1.0);
            glTexCoord2f(displayWidth, 0); glVertex2f(1.0, 1.0);
            glEnd();

            glDisable(GL_TEXTURE_RECTANGLE_ARB);
#else
            glBindTexture( GL_TEXTURE_2D, texture );

            glBegin( GL_QUADS );
            glTexCoord2d(0.0,0.0); glVertex2d(-1.0,1.0);
            glTexCoord2d(1.0,0.0); glVertex2d(1.0,1.0);
            glTexCoord2d(1.0,1.0); glVertex2d(1.0,-1.0);
            glTexCoord2d(0.0,1.0); glVertex2d(-1.0,-1.0);
            glEnd();

            glBindTexture( GL_TEXTURE_2D, 0 );
#endif

            // Update the window
            glutSwapBuffers();
        }
    }
}

void outputStats(bool exitNow) {

    // Calculate the PSNR regardless of whether or not we were calculating it
    double MSE, PSNR;
    MSE = (double)totalSE / (double)(displayWidth * displayHeight * totalUpdates) / 3.0f;
    PSNR = 10.0f * log10(65025.0f / MSE);

    //
    // Print a detailed, readable report if we're not running configHeadless
    //
    if (!globalConfiguration.headless) {

        cout << endl;

        // General
        //
        cout << "General Information" << endl;
        cout << "  File: " << fileName << endl;
        cout << "  Dimensions: " << displayWidth << " x " << displayHeight << endl;
        cout << "  Coefficient Planes: " << (myDisplay ? myDisplay->NumCoefficientPlanes() : 0) << endl;
        cout << "  Frames Rendered: " << totalUpdates << endl;
        cout << endl;

        // Configuration
        //
        cout << "Configuration Information:" << endl;

        switch (globalConfiguration.tiler) {
            case SIMPLE:
                cout << "  Configuring NDDI as a simple Framebuffer." << endl;
                break;
            case FLAT:
                cout << "  Configuring NDDI as Flat Tiled." << endl;
                cout << "  Tile Dimensions: " << globalConfiguration.tileWidth << "x" << globalConfiguration.tileHeight << endl;
                break;
            case CACHE:
                cout << "  Configuring NDDI as Cached Tiled." << endl;
                cout << "  Tile Dimensions: " << globalConfiguration.tileWidth << "x" << globalConfiguration.tileHeight << endl;
                cout << "  Significant Bits: " << globalConfiguration.sigBits << endl;
                break;
            case DCT:
                cout << "  Configuring NDDI as DCT Tiled." << endl;
                cout << "  Quality: " << globalConfiguration.quality << endl;
                break;
            case IT:
                cout << "  Configuring NDDI as IT Tiled." << endl;
                cout << "  Quality: " << globalConfiguration.quality << endl;
                break;
            case COUNT:
                cout << "  Configuring NDDI to just count." << endl;
                break;
            case FLOW:
                cout << "  Configuring NDDI to compute optical flow." << endl;
                break;
            default:
                break;
        }
        cout << endl;

        // Transmission
        //
        cout << "Transmission Statistics:" << endl;
        // Get total transmission cost
        long totalCost = costModel->getLinkBytesTransmitted();
        cout << "  Total Pixel Data Updated (bytes): " << totalUpdates * displayWidth * displayHeight * BYTES_PER_PIXEL <<
        " Total NDDI Cost (bytes): " << totalCost <<
        " Ratio: " << (double)totalCost / (double)totalUpdates / (double)displayWidth / (double)displayHeight / BYTES_PER_PIXEL << endl;
        cout << endl;


        // Memory
        //
        cout << "Memory Statistics:" << endl;
        cout << "  Input Vector" << endl;
        cout << "    - Num Reads: " << costModel->getReadAccessCount(INPUT_VECTOR_COMPONENT) <<  " - Bytes Read: " << costModel->getBytesRead(INPUT_VECTOR_COMPONENT) << endl;
        cout << "    - Num Writes: " << costModel->getWriteAccessCount(INPUT_VECTOR_COMPONENT) << " - Bytes Written: " << costModel->getBytesWritten(INPUT_VECTOR_COMPONENT) << endl;
        cout << "  Coefficient Plane" << endl;
        cout << "    - Num Reads: " << costModel->getReadAccessCount(COEFFICIENT_PLANE_COMPONENT) <<  " - Bytes Read: " << costModel->getBytesRead(COEFFICIENT_PLANE_COMPONENT) << endl;
        cout << "    - Num Writes: " << costModel->getWriteAccessCount(COEFFICIENT_PLANE_COMPONENT) << " - Bytes Written: " << costModel->getBytesWritten(COEFFICIENT_PLANE_COMPONENT) << endl;
        cout << "  Frame Volume" << endl;
        cout << "    - Num Reads: " << costModel->getReadAccessCount(FRAME_VOLUME_COMPONENT) <<  " - Bytes Read: " << costModel->getBytesRead(FRAME_VOLUME_COMPONENT) << endl;
        cout << "    - Num Writes: " << costModel->getWriteAccessCount(FRAME_VOLUME_COMPONENT) << " - Bytes Written: " << costModel->getBytesWritten(FRAME_VOLUME_COMPONENT) << endl;
        cout << endl;


        // Pixel
        //
        cout << "Pixel Statistics:" << endl;
        cout << "  Pixel Mappings: " << costModel->getPixelsMapped() << endl;
        cout << "  Pixel Blends: " << costModel->getPixelsBlended() << endl;
        cout << endl;

        // Performance
        //
        gettimeofday(&endTime, NULL);
        cout << "Performance Statistics:" << endl;
        cout << "  Average FPS: " << (double)totalUpdates / ((double)(endTime.tv_sec * 1000000
                                                                      + endTime.tv_usec
                                                                      - startTime.tv_sec * 1000000
                                                                      - startTime.tv_usec) / 1000000.0f) << endl;
        cout << endl;

        // Quality
        //
        if (globalConfiguration.PSNR) {
            cout << "Quality Statistics:" << endl << "  MSE: " << MSE << endl << "  Total PSNR: " << PSNR << endl;
        } else {
            cout << "Quality Statistics:" << endl << "  Use --psnr to enable." << endl;
        }
        cout << endl;
    }

    // CSV
    //

    if (globalConfiguration.csv) {
        // Pretty print a heading to stdout, but for headless just spit it to stderr for reference
        if (!globalConfiguration.headless) {
            cout << "CSV Headings:" << endl;
            cout << "Frames,Commands Sent,Bytes Transmitted,IV Num Reads,IV Bytes Read,IV Num Writes,IV Bytes Written,CP Num Reads,CP Bytes Read,CP Num Writes,CP Bytes Written,FV Num Reads,FV Bytes Read,FV Num Writes,FV Bytes Written,FV Time,Pixels Mapped,Pixels Blended,PSNR,File Name" << endl;
        } else {
            cerr << "CSV Headings:" << endl;
            cerr << "Frames,Commands Sent,Bytes Transmitted,IV Num Reads,IV Bytes Read,IV Num Writes,IV Bytes Written,CP Num Reads,CP Bytes Read,CP Num Writes,CP Bytes Written,FV Num Reads,FV Bytes Read,FV Num Writes,FV Bytes Written,FV Time,Pixels Mapped,Pixels Blended,PSNR,File Name" << endl;
        }

        cout
        << totalUpdates << " , "
        << costModel->getLinkCommandsSent() << " , "
        << costModel->getLinkBytesTransmitted() << " , "
        << costModel->getReadAccessCount(INPUT_VECTOR_COMPONENT) << " , "
        << costModel->getBytesRead(INPUT_VECTOR_COMPONENT) << " , "
        << costModel->getWriteAccessCount(INPUT_VECTOR_COMPONENT) << " , "
        << costModel->getBytesWritten(INPUT_VECTOR_COMPONENT) << " , "
        << costModel->getReadAccessCount(COEFFICIENT_PLANE_COMPONENT) << " , "
        << costModel->getBytesRead(COEFFICIENT_PLANE_COMPONENT) << " , "
        << costModel->getWriteAccessCount(COEFFICIENT_PLANE_COMPONENT) << " , "
        << costModel->getBytesWritten(COEFFICIENT_PLANE_COMPONENT) << " , "
        << costModel->getReadAccessCount(FRAME_VOLUME_COMPONENT) << " , "
        << costModel->getBytesRead(FRAME_VOLUME_COMPONENT) << " , "
        << costModel->getWriteAccessCount(FRAME_VOLUME_COMPONENT) << " , "
        << costModel->getBytesWritten(FRAME_VOLUME_COMPONENT) << " , "
        << costModel->getTime(FRAME_VOLUME_COMPONENT) << " , "
        << costModel->getPixelsMapped() << " , "
        << costModel->getPixelsBlended() << " , "
        << PSNR << ","
        << fileName << endl;
    }

    cerr << endl;

    // Warnings about Configuration
#if defined(SUPRESS_EXCESS_RENDERING) || defined(SKIP_COMPUTE_WHEN_SCALER_ZERO) || defined(NO_CL) || defined(NO_GL) || defined(CLEAR_COST_MODEL_AFTER_SETUP)
    cerr << endl << "CONFIGURATION WARNINGS:" << endl;
#ifdef SUPRESS_EXCESS_RENDERING
    cerr << "  - Was compiled with SUPRESS_EXCESS_RENDERING, and so the numbers may be off. Recompile with \"make NO_HACKS=1\"." << endl;
#endif
#ifdef SKIP_COMPUTE_WHEN_SCALER_ZERO
    cerr << "  - Was compiled with SKIP_COMPUTE_WHEN_SCALER_ZERO, and so the numbers may be off when running with NO_OMP." << endl <<
            "    When using OpenMP, the number will be fine regardless because they're register in bulk later. Recompile" << endl <<
            "    with \"make NO_HACKS=1\"." << endl;
#endif
#ifdef NO_CL
    cerr << "  - Was compiled without OpenCL." << endl;
#endif
#ifdef NO_GL
    cerr << "  - Was compiled without OpenGL." << endl;
#endif
#ifdef CLEAR_COST_MODEL_AFTER_SETUP
    cerr << "  - Was compiled with CLEAR_COST_MODEL_AFTER_SETUP, affecting the true cost." << endl;
#endif
#endif

    // Clean up
    if (exitNow) {
        delete myDisplay;
        // TODO(CDE): Commenting out ot avoid a warning. Figure out the warning.
        // if (myPlayer) { delete(myPlayer); }
        if (lastBuffer) { free(lastBuffer); }

        exit(0);
    }
}


void renderFrame() {

    static int framesDecoded = 0;
    static int framesRendered = 0;
    static int currentFrame = 0;
    long transCost = costModel->getLinkBytesTransmitted();
    long SE = 0;

    static rewind_play_t rewindState = NOT_REWINDING;

    if (!globalConfiguration.maxFrames || (totalUpdates < globalConfiguration.maxFrames)) {
        // If we're not in a rewind backwards or forward state...
        if (rewindState == NOT_REWINDING) {
#ifdef USE_ASYNC_DECODER
            // Get a decoded frame
            // TODO(CDE): Don't poll. Use signals instead.
            while (bufferQueue.size() == 0) {
                usleep(5);
             }
            if ((videoBuffer = bufferQueue.front()) != NULL) {
                pthread_mutex_lock(&bufferQueueMutex);
                bufferQueue.pop();
                pthread_mutex_unlock(&bufferQueueMutex);
#else
            // Decode a frame
            if ((videoBuffer = myPlayer->decodeFrame()) != NULL) {
#endif
                framesDecoded++;

                // If we're past the designated start frame, then update the NDDI display
                if (framesDecoded >= globalConfiguration.startFrame) {
                    // If configured to rewind, then check to see if this decoded frame should be stored.
                    if ( myRewinder &&
                        (currentFrame >= (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)) &&
                        (currentFrame < globalConfiguration.rewindStartFrame) ) {
                        myRewinder->CopyFrame(videoBuffer, currentFrame - (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames));
                    }

                    // Update NDDI Display
                    updateDisplay(videoBuffer,
                                  myPlayer->width(), myPlayer->height());

                    framesRendered++;

                    // If we used a lossy mode, then calculate the Square Error for this frame
                    if (globalConfiguration.PSNR) {
                        Pixel* frameBuffer = myDisplay->GetFrameBuffer();
                        SE = calculateSE(videoBuffer, frameBuffer);
                        totalSE += SE;
                    }

                    // Draw
                    glutPostRedisplay();


                    // If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
                    if (myRewinder && currentFrame == globalConfiguration.rewindStartFrame) {
                        rewindState = PLAY_BACKWARDS;
                        currentFrame--;
                    } else {
                        currentFrame++;
                    }
                }
            } else {
                // Output stats and exit
                outputStats(true);
            }

        // Else if we're in the rewind placy backwards state...
        } else if (rewindState == PLAY_BACKWARDS) {
            // Update NDDI Display
            updateDisplay(myRewinder->GetFrame(currentFrame - (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)),
                          myPlayer->width(), myPlayer->height());

            framesRendered++;

            // Draw
            glutPostRedisplay();

            // If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
            if (currentFrame == (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)) {
                rewindState = PLAY_FORWARD;
                currentFrame++;
            } else {
                currentFrame--;
            }
        // Else if we're in the rewind play forward state...
        } else if (rewindState == PLAY_FORWARD) {
            // Update NDDI Display
            updateDisplay(myRewinder->GetFrame(currentFrame - (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)),
                          myPlayer->width(), myPlayer->height());

            framesRendered++;

            // Draw
            glutPostRedisplay();

            // Increment the currentFrame
            currentFrame++;

            // If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
            if (currentFrame == globalConfiguration.rewindStartFrame) {
                rewindState = NOT_REWINDING;
                delete(myRewinder);
            }

        }

        if (!globalConfiguration.headless) {
            if (globalConfiguration.verbose) {
                cout << "PixelBidge Statistics:" << endl;
                cout << "  Decoded Frames: " << framesDecoded << " - Rendered Frames: " << framesRendered << endl;
                cout << "  Transmission: " << costModel->getLinkBytesTransmitted() - transCost << endl;
            }

            double MSE, PSNR = 0.0;
            if (globalConfiguration.PSNR) {
                MSE = (double)SE / (double)(displayWidth * displayHeight) / 3.0f;
                PSNR = 10.0f * log10(65025.0f / MSE);
                if (globalConfiguration.verbose)
                    cout << "  PSNR: " << PSNR << endl;
            }

            if (globalConfiguration.csv)
                cout << "RenderCSV," << framesRendered << "," << costModel->getLinkBytesTransmitted() - transCost << "," << PSNR << endl;
        }

    } else {
        // Output stats and exit
        outputStats(true);
    }

}


/**
 * Operates much like renderFrame(), except it just counts changed pixels instead of
 * rendering the contents. Used exclusively for the COUNT configuration which corresponds
 * to Perfect Pixel Latching AKA the Ideal Pixel Bridge Mode.
 */
void countChangedPixels() {

    static int framesDecoded = 0;
    static int framesCounted = 0;

    if (!globalConfiguration.maxFrames || (totalUpdates < globalConfiguration.maxFrames)) {
        int pixel_count = myPlayer->width() * myPlayer->height();
#ifdef USE_ASYNC_DECODER
        // Get a decoded frame
        // TODO(CDE): Don't poll. Use signals instead.
        while (bufferQueue.size() == 0) {
            usleep(5);
        }
        if ((videoBuffer = bufferQueue.front()) != NULL) {

            pthread_mutex_lock(&bufferQueueMutex);
            bufferQueue.pop();
            pthread_mutex_unlock(&bufferQueueMutex);
#else
            // Decode a frame
            if ((videoBuffer = myPlayer->decodeFrame()) != NULL) {
#endif
            framesDecoded++;

            if (framesDecoded > globalConfiguration.startFrame) {
                // If we're in a rewind period, then we'll be incrementing diffs as well as frame count by 3 instead of 1
                int inc;
                if ( (framesDecoded > (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)) && (framesDecoded <= globalConfiguration.rewindStartFrame) ) {
                    inc = 3;
                } else {
                    inc = 1;
                }

                // If we haven't previously decoded a frame
                if (lastBuffer == NULL) {
                    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(inc * pixel_count), 0);

                    // Allocate the lastBuffer and copy the decode frame into it.
                    lastBuffer = (uint8_t*)malloc(VIDEO_PIXEL_SIZE * pixel_count);
                    memcpy(lastBuffer, videoBuffer, VIDEO_PIXEL_SIZE * pixel_count);

                    // Otherwise compare the newly decode frame to the previous frame, pixel by pixel
                } else {
                    // Count the changed pixels
                    int diffs = 0;

                    for (int offset = 0; offset < (pixel_count * 3);) {
                        if (lastBuffer[offset] != videoBuffer[offset]) {
                            diffs += inc;
                            offset += 3;
                            continue;
                        }
                        offset++;

                        if (lastBuffer[offset] != videoBuffer[offset]) {
                            diffs += inc;
                            offset += 2;
                            continue;
                        }
                        offset++;

                        if (lastBuffer[offset] != videoBuffer[offset]) {
                            diffs += inc;
                        }
                        offset++;
                    }

                    // Update the cost
                    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(diffs), 0);

                    // Then copy to the lastBuffer
                    memcpy(lastBuffer, videoBuffer, VIDEO_PIXEL_SIZE * pixel_count);
                }

                // Update the totalUpdates
                totalUpdates += inc;
                framesCounted += inc;

                // Draw nothing
                glutPostRedisplay();
            }
        } else {
            // Output stats and exit
            outputStats(true);
        }

        if (!globalConfiguration.headless && globalConfiguration.verbose) {
            cout << "PixelBidge Statistics:" << endl << "  Decoded Frames: " << framesDecoded << " - Counted Frames: " << framesCounted << endl;
        }

    } else {
        // Output stats and exit
        outputStats(true);
    }
}

/**
 * Used to compute the optical flow between two frames and to report.
 */
void computeFlow() {
    static int framesDecoded = 0;
    static std::vector<cv::Point2f> features, featuresLast;

    if (!globalConfiguration.maxFrames || (totalUpdates < globalConfiguration.maxFrames)) {
        int pixel_count = myPlayer->width() * myPlayer->height();
#ifdef USE_ASYNC_DECODER
        // Get a decoded frame
        // TODO(CDE): Don't poll. Use signals instead.
        while (bufferQueue.size() == 0) {
            usleep(5);
        }
        if ((videoBuffer = bufferQueue.front()) != NULL) {

            pthread_mutex_lock(&bufferQueueMutex);
            bufferQueue.pop();
            pthread_mutex_unlock(&bufferQueueMutex);
#else
            // Decode a frame
            if ((videoBuffer = myPlayer->decodeFrame()) != NULL) {
#endif
            framesDecoded++;

            if (framesDecoded > globalConfiguration.startFrame) {
                // If we're in a rewind period, then we'll be incrementing diffs as well as frame count by 3 instead of 1
                int inc;
                if ( (framesDecoded > (globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames)) && (framesDecoded <= globalConfiguration.rewindStartFrame) ) {
                    inc = 3;
                } else {
                    inc = 1;
                }

                // If we haven't previously decoded a frame
                if (lastBuffer == NULL) {
#if 0 // Normally you pick points once, and they're tracked throughout. However, they dwindle over time
      // so I'm opting to just find them every frame below instead of here.
                    // Convert the image buffers to OpenCV Mats
                    cv::Mat imgCurrent(displayHeight, displayWidth , CV_8UC3, videoBuffer);

                    // Then convert to grayscale
                    cv::Mat grayCurrent;
                    cv::cvtColor(imgCurrent, grayCurrent, CV_RGB2GRAY);

                    // Then find the features to track
                    cv::goodFeaturesToTrack(grayCurrent, featuresLast,
                                            500,  // the maximum number of features
                                            0.01, // quality level (top 1%)
                                            10);  // min distance between two features
#endif

                    // Allocate the lastBuffer and copy the decode frame into it.
                    lastBuffer = (uint8_t*)malloc(VIDEO_PIXEL_SIZE * pixel_count);
                    memcpy(lastBuffer, videoBuffer, VIDEO_PIXEL_SIZE * pixel_count);

                // Otherwise compute flow between the newly decoded frame to the previous frame
                } else {

                    // Convert the image buffers to OpenCV Mats
                    cv::Mat imgCurrent(displayHeight, displayWidth , CV_8UC3, videoBuffer);
                    cv::Mat imgLast(displayHeight, displayWidth , CV_8UC3, lastBuffer);

                    // Then convert to grayscale
                    cv::Mat grayCurrent, grayLast;
                    cv::cvtColor(imgCurrent, grayCurrent, CV_RGB2GRAY);
                    cv::cvtColor(imgLast, grayLast, CV_RGB2GRAY);

#if 1 // Normally these wouldn't be found each time, they'd just be re-used. See comment in previous if block.
                    // Then find the features to track
                    cv::goodFeaturesToTrack(grayCurrent, featuresLast,
                                            500,  // the maximum number of features
                                            0.01, // quality level (top 1%)
                                            10);  // min distance between two features
#endif

                    vector<uchar> status;
                    vector<float> err;
                    cv::calcOpticalFlowPyrLK(grayLast, grayCurrent,
                                             featuresLast, // input point positions in previous image
                                             features,     // output point positions for the current
                                             status,       // tracking success
                                             err);         // tracking error

                    float dist;
                    int j = 0;
                    for (int i = 0; i < featuresLast.size(); i++) {
                        if (status[i]) {
                            dist += sqrt((features[j].x - featuresLast[i].x) * (features[j].x - featuresLast[i].x) +
                                         (features[j].y - featuresLast[i].y) * (features[j].y - featuresLast[i].y));
                            featuresLast[j] = features[j];
                            j++;
                        }
                    }
                    featuresLast.resize(j);
                    dist = dist / j;

                    if (!globalConfiguration.headless && globalConfiguration.verbose)
                        cout << "FlowCSV," << framesDecoded - 1 << "," << framesDecoded << "," << dist << endl;

                    // Then copy to the lastBuffer
                    memcpy(lastBuffer, videoBuffer, VIDEO_PIXEL_SIZE * pixel_count);
                }

                // Update the totalUpdates
                totalUpdates += inc;

                // Draw nothing
                glutPostRedisplay();
            }
        } else {
            // Output stats and exit
            outputStats(true);
        }

    } else {
        // Output stats and exit
        outputStats(true);
    }
}


void keyboard( unsigned char key, int x, int y ) {

    switch (key) {
        case 27: // Esc
            outputStats(true);
            break;
        case ('u'):
        case('U'):
            outputStats(false);
            break;
        default:
            break;
    }

    glutPostRedisplay();
}


static bool doing_it = false;


void mouse( int button, int state, int x, int y ) {

    if (button == GLUT_LEFT && state == GLUT_DOWN) {
        doing_it = true;
    } else {
        doing_it = false;
    }
    glutPostRedisplay();
}


void motion( int x, int y ) {

    if (doing_it) {
        glutPostRedisplay();
    }
}


void showUsage() {
    cout << "pixelbridge [--mode <fb|flat|cache|dct|count|flow>] [--blend <fv|t|cp|>] [--ts <n> <n>] [--tc <n>] [--bits <1-8>]" << endl <<
            "            [--dctscales x:y[,x:y...]] [--dctdelta <n>] [--dctplanes <n>] [--dctbudget <n>] [--dctsnap] [--dcttrim] [--quality <0/1-100>]" << endl <<
            "            [--start <n>] [--frames <n>] [--rewind <n> <n>] [--psnr] [--verbose] [--csv] [--headless] <filename>" << endl;
    cout << endl;
    cout << "  --mode  Configure NDDI as a framebuffer (fb), as a flat tile array (flat), as a cached tile (cache), using DCT (dct), or using IT (it).\n" <<
            "          Optional the mode can be set to count the number of pixels changed (count) or determine optical flow (flow)." << endl;
    cout << "  --blend  For fb mode, this option alpha-blends using just the frame volume, input vector (temporal), or coefficient plane." << endl;
    cout << "  --ts  Sets the tile size to the width and height provided." << endl;
    cout << "  --tc  Sets the maximum number of tiles in the cache." << endl;
    cout << "  --bits  Sets the number of significant bits per channel when computing checksums." << endl;
    cout << "  --dctscales  For dct mode, this will set a series of scales in the form of comma-separated two tuples holding the scale and then\n" <<
            "               the edge length of a square that determines the number of planes used (i.e. 2 -> 2 x 2 = 4 planes." << endl;
    cout << "  --dctdelta  For dct mode, this is the delta between coefficients that is considered a match which will not be updated.\n" <<
            "              Setting to zero when using a budget will use an optimal setting." << endl;
    cout << "  --dctplanes  For dct mode, this is the number of planes either zeroed or trimmed.\n" <<
            "               Setting to zero when using a budget will use an optimal setting." << endl;
    cout << "  --dctbudget  For dct mode, this is budget in bytes for the transmission of each frame." << endl;
    cout << "  --dctsnap  For dct mode, this turns on snap to zero using either planes or delta." << endl;
    cout << "  --dcttrim  For dct mode, this turns on lossy trimming using either planes or delta." << endl;
    cout << "  --quality  Sets the quality for DCT [1-100] or IT [0-100] mode." << endl;
    cout << "  --start  Will start with this frame, ignoring any decoded frames prior to it." << endl;
    cout << "  --frames  Sets the number of maximum frames that are decoded." << endl;
    cout << "  --rewind  Sets a start point and a number of frames to play in reverse. Once finished, normal playback continues." << endl;
    cout << "  --psnr  Calculates and outputs PSNR. Cannot use with headless." << endl;
    cout << "  --verbose  Outputs frame-by-frame statistics." << endl;
    cout << "  --csv  Outputs CSV data." << endl;
    cout << "  --headless  Removes rendering and excessive data output. Overrides --verbose. Cannot use with psnr." << endl;
}


bool parseArgs(int argc, char *argv[]) {
    argc--;
    argv++;

    while (argc) {
        if (strcmp(*argv, "--mode") == 0) {
            argc--;
            argv++;
            if (strcmp(*argv, "fb") == 0) {
                globalConfiguration.tiler = SIMPLE;
            } else if (strcmp(*argv, "flat") == 0) {
                globalConfiguration.tiler = FLAT;
            } else if (strcmp(*argv, "cache") == 0) {
                globalConfiguration.tiler = CACHE;
            } else if (strcmp(*argv, "dct") == 0) {
                globalConfiguration.tiler = DCT;
            } else if (strcmp(*argv, "it") == 0) {
                globalConfiguration.tiler = IT;
            } else if (strcmp(*argv, "count") == 0) {
                globalConfiguration.tiler = COUNT;
            } else if (strcmp(*argv, "flow") == 0) {
                globalConfiguration.tiler = FLOW;
            } else {
                showUsage();
                return false;
            }
            argc--;
            argv++;
        } else if (strcmp(*argv, "--dctscales") == 0) {
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            argc--;
            argv++;
            // Clear the default dct scales
            globalConfiguration.clearDctScales();
            // Then pull out each tuple x:y
            char *p = strtok(*argv, ",");
            size_t first = 0;
            while (p != NULL) {
                // Scan the scale and count
                int scale, edge;
                sscanf(p, "%d:%d", &scale, &edge);
                if (edge < 0 || edge > 8) {
                    showUsage();
                    return false;
                }
                // Add the dct scale
                globalConfiguration.addDctScale(scale, first, edge);
                first += edge * edge;
                // Move to next tuple
                p = strtok(NULL, ",");
            }
            argc--;
            argv++;
        } else if (strcmp(*argv, "--dctdelta") == 0) {
            argc--;
            argv++;
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            globalConfiguration.dctDelta = atoi(*argv);
            argc--;
            argv++;
        } else if (strcmp(*argv, "--dctplanes") == 0) {
            argc--;
            argv++;
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            globalConfiguration.dctPlanes = atoi(*argv);
            argc--;
            argv++;
        } else if (strcmp(*argv, "--dctbudget") == 0) {
            argc--;
            argv++;
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            globalConfiguration.dctBudget = atoi(*argv);
            argc--;
            argv++;
        } else if (strcmp(*argv, "--dctsnap") == 0) {
            argc--;
            argv++;
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            globalConfiguration.dctSnap = true;
        } else if (strcmp(*argv, "--dcttrim") == 0) {
            argc--;
            argv++;
            if (globalConfiguration.tiler != DCT) {
                showUsage();
                return false;
            }
            globalConfiguration.dctTrim = true;
        } else if (strcmp(*argv, "--blend") == 0) {
            if (globalConfiguration.tiler != SIMPLE) {
                showUsage();
                return false;
            }
            argc--;
            argv++;
            if (strcmp(*argv, "fv") == 0) {
                globalConfiguration.blend = FRAME_VOLUME;
            } else if (strcmp(*argv, "t") == 0) {
                globalConfiguration.blend = TEMPORAL;
            } else if (strcmp(*argv, "cp") == 0) {
                globalConfiguration.blend = COEFFICIENT_PLANE;
            } else {
                showUsage();
                return false;
            }
            argc--;
            argv++;
        } else if (strcmp(*argv, "--ts") == 0) {
            globalConfiguration.tileWidth = atoi(argv[1]);
            globalConfiguration.tileHeight = atoi(argv[2]);
            if ((globalConfiguration.tileWidth == 0) || (globalConfiguration.tileHeight == 0)) {
                showUsage();
                return false;
            }
            argc -= 3;
            argv += 3;
        } else if (strcmp(*argv, "--tc") == 0) {
            globalConfiguration.maxTiles = atoi(argv[1]);
            if (globalConfiguration.maxTiles == 0) {
                showUsage();
                return false;
            }
            argc -= 2;
            argv += 2;
        } else if (strcmp(*argv, "--bits") == 0) {
            globalConfiguration.sigBits = atoi(argv[1]);
            if ( (globalConfiguration.sigBits == 0) || (globalConfiguration.sigBits > 8) ) {
                showUsage();
                return false;
            }
            argc -= 2;
            argv += 2;
        } else if (strcmp(*argv, "--quality") == 0) {
            globalConfiguration.quality = atoi(argv[1]);
            if ( (globalConfiguration.quality == 0) || (globalConfiguration.quality > 100) ) {
                showUsage();
                return false;
            }
            argc -= 2;
            argv += 2;
        } else if (strcmp(*argv, "--start") == 0) {
            globalConfiguration.startFrame = atoi(argv[1]);
            if (globalConfiguration.startFrame == 0) {
                showUsage();
                return false;
            }
            argc -= 2;
            argv += 2;
        } else if (strcmp(*argv, "--frames") == 0) {
            globalConfiguration.maxFrames = atoi(argv[1]);
            if (globalConfiguration.maxFrames == 0) {
                showUsage();
                return false;
            }
            argc -= 2;
            argv += 2;
        } else if (strcmp(*argv, "--rewind") == 0) {
            globalConfiguration.rewindStartFrame = atoi(argv[1]);
            globalConfiguration.rewindFrames = atoi(argv[2]);
            if ((globalConfiguration.rewindStartFrame == 0) || (globalConfiguration.rewindFrames == 0) || (globalConfiguration.rewindFrames > globalConfiguration.rewindStartFrame)) {
                showUsage();
                return false;
            }
            argc -= 3;
            argv += 3;
        } else if (strcmp(*argv, "--psnr") == 0) {
            if (globalConfiguration.headless) {
                showUsage();
                return false;
            }
            globalConfiguration.PSNR = true;
            argc--;
            argv++;
        } else if (strcmp(*argv, "--verbose") == 0) {
            globalConfiguration.verbose = true;
            argc--;
            argv++;
        } else if (strcmp(*argv, "--csv") == 0) {
            globalConfiguration.csv = true;
            argc--;
            argv++;
        } else if (strcmp(*argv, "--headless") == 0) {
            if (globalConfiguration.PSNR) {
                showUsage();
                return false;
            }
            globalConfiguration.headless = true;
            argc--;
            argv++;
        } else {
            fileName = *argv;
            argc--;
            argv++;
        }
    }

    if (!fileName) {
        showUsage();
        return false;
    } else {
        return true;
    }
}


#ifdef USE_ASYNC_DECODER
void* decoderRun(void* none) {
    uint8_t* buf;

    while ((buf = myPlayer->decodeFrame()) != NULL) {
        // TODO(CDE): Don't poll. Use signals instead.
        while (bufferQueue.size() >= 20) {
            usleep(50);
        }
        pthread_mutex_lock(&bufferQueueMutex);
        bufferQueue.push(buf);
        pthread_mutex_unlock(&bufferQueueMutex);
    }
    // Push a final NULL to trigger the end of stream
    bufferQueue.push(NULL);

    return NULL;
}
#endif


int main(int argc, char *argv[]) {

    // Parse command line arguments
    if (!parseArgs(argc, argv)) {
        return -1;
    }

    // Initialize ffmpeg and set dimensions
#ifndef USE_RANDOM_PLAYER
    myPlayer = (Player*)new FfmpegPlayer(fileName);
#else
    myPlayer = (Player*)new RandomPlayer();
#endif
#ifndef USE_SMALLER_WINDOW
    displayWidth = myPlayer->width();
    displayHeight = myPlayer->height();
#endif
    if (globalConfiguration.rewindStartFrame > 0) {
        myRewinder = new Rewinder(globalConfiguration.rewindStartFrame - globalConfiguration.rewindFrames, displayWidth, displayHeight);
    }
    if ( (displayWidth == 0) || (displayHeight == 0) ) {
        cerr << "Error: Could not get video dimensions." << endl;
        return -1;
    }

#ifdef USE_ASYNC_DECODER
    // Start decoded thread
    pthread_create(&decoderThread, NULL, decoderRun, NULL);
#endif

    // Initialize GLUT
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
    glutInitWindowSize( displayWidth, displayHeight );

    glutCreateWindow(fileName);

    glutDisplayFunc( draw );
    glutKeyboardFunc( keyboard );
    glutMouseFunc( mouse );
    glutMotionFunc( motion );

    glEnable( GL_TEXTURE_2D );

    if (globalConfiguration.tiler > COUNT) {
        glutIdleFunc(renderFrame);
    } else if (globalConfiguration.tiler == COUNT){
        glutIdleFunc(countChangedPixels);
    } else if (globalConfiguration.tiler == FLOW){
        glutIdleFunc(computeFlow);
    }

    // Setup NDDI Display
    if (globalConfiguration.tiler > COUNT) {

        // Force tilesize to be 8x8 for DCT since we're using 8x8 macroblocks
        if (globalConfiguration.tiler == DCT) {
            globalConfiguration.tileWidth = globalConfiguration.tileHeight = 8;
        // Else if the tile size wasn't specified, then dynamically calculate it
        } else if ((globalConfiguration.tileWidth == 0) || (globalConfiguration.tileHeight == 0)) {
            // Calculate the tile as a square who's edge is 1/40 of the longest display dimension
            size_t edge = ((displayWidth > displayHeight) ? displayWidth : displayHeight) / 40;
            // But don't go smaller than 8x8 when automatically calculating the tile size
            if (edge < 8)
                edge = 8;
            globalConfiguration.tileWidth = globalConfiguration.tileHeight = edge;
        }

        // Setup the GlNddiDisplay and Tiler if required
        setupDisplay();
    } else {
        costModel = new CostModel();
    }

    // Take the start time stamp
    gettimeofday(&startTime, NULL);

    // Run main loop
    glutMainLoop();

    return 0;
}
