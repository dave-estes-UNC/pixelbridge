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

#include "PixelBridgeFeatures.h"
#include "CostModel.h"
#include "GlNddiDisplay.h"
#include "ClNddiDisplay.h"
#include "BlendingGlNddiDisplay.h"
#include "CachedTiler.h"
#include "FlatTiler.h"
#include "FfmpegPlayer.h"
#include "Rewinder.h"


using namespace std;

typedef enum {
	COUNT,  // Just count (aka Perfect Pixel Latching)
	SIMPLE, // Simple Framebuffer
	FLAT,   // Tiled, but not cached
	CACHE   // Tiled and cached
} config_t;

typedef enum {
    NONE = 0,         // No blending
    FRAME_VOLUME,     // Uses frame volume blending operations
    TEMPORAL,         // Uses the input vector to switch rapidly between planes
    COEFFICIENT_PLANE // Uses multiple coefficient planes
} blend_t;

#define TEMPORAL_FLIP_PER_FRAME_COUNT 4

typedef enum {
	NOT_REWINDING,  // Normal playback, not rewinding
	PLAY_BACKWARDS, // Rewinding frame by frame backwards
	PLAY_FORWARD    // Finished rewinding, playing forward from memory
} rewind_play_t;


// General Globals
int displayWidth = 320, displayHeight = 240;
const char* fileName = NULL;

// Helper Objects
FfmpegPlayer*  myPlayer;
GlNddiDisplay* myDisplay;
BlendingGlNddiDisplay* myBlendingDisplay;
Tiler* myTiler;
Rewinder* myRewinder = NULL;

// Decoder thread
pthread_t           decoderThread;

// Configuration Options
config_t config = CACHE;
blend_t configBlend = NONE;
size_t configTileWidth = 0;
size_t configTileHeight = 0;
size_t configMaxTiles = 1000;
size_t configSigBits = 8;
size_t configStartFrame = 0;
size_t configMaxFrames = 0;
size_t configRewindStartFrame = 0;
size_t configRewindFrames = 0;
bool configHeadless = false;

// Statistical Instrumentation
CostModel* costModel;
int totalUpdates = 0;
timeval startTime, endTime; // Used for timing data
long totalSE = 0; // Cumulative Square Error used to calculate MSE and then PSNR

// Stores the current and previously decoded frame
uint8_t* videoBuffer = NULL;
uint8_t* lastBuffer = NULL;

// Buffers decoded frame
std::queue<uint8_t*> bufferQueue;

// Synchronization variables
pthread_mutex_t      bufferQueueMutex;


void setupDisplay() {
    
	// Cached-Tiled
	if (config == CACHE) {

		// 3 dimensional matching the Tile Width x Height x max tiles
		std::vector<unsigned int> fvDimensions;
		fvDimensions.push_back(configTileWidth);
		fvDimensions.push_back(configTileHeight);
		fvDimensions.push_back(configMaxTiles);
		
		myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
									  displayWidth, displayHeight, // display size
									  3); 						   // input vector size (x, y, and z)
        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Input Vector
		std::vector<int> iv;
		iv.push_back(1);
		myDisplay->UpdateInputVector(iv);
		
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = p.a = 0xff;
		std::vector<unsigned int> start, end;
		start.push_back(0); start.push_back(0); start.push_back(0);
		end.push_back(configTileWidth); end.push_back(configTileHeight); end.push_back(configMaxTiles);
		myDisplay->FillPixel(p, start, end);
		
		// Setup Cached Tiler which initializies Coefficient Plane
		myTiler = new CachedTiler(myDisplay,
								  configTileWidth, configTileHeight,
								  configMaxTiles, configSigBits,
								  configHeadless);

        myTiler->InitializeCoefficientPlane();

    // Frame Volume Blending
    } else if ((config == SIMPLE) && (configBlend == FRAME_VOLUME)) {
		
		// 3 dimensional matching the Video Width x Height x 3
		std::vector<unsigned int> fvDimensions;
		fvDimensions.push_back(displayWidth);
		fvDimensions.push_back(displayHeight);
		fvDimensions.push_back(3);
		
		myBlendingDisplay = new BlendingGlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                                      displayWidth, displayHeight, // display size
                                                      3); 						   // input vector size (x, y, 1)
        myDisplay = myBlendingDisplay;
		
        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Input Vector
		std::vector<int> iv;
		iv.push_back(0);
		myDisplay->UpdateInputVector(iv);
		
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = 0x00; p.a = 0x7f;
		std::vector<unsigned int> start, end;
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
		
		start.push_back(0); start.push_back(0); 
		end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
		
        myDisplay->FillCoefficientMatrix(coeffs, start, end);

    // Temporal Blending
    } else if ((config == SIMPLE) && (configBlend == TEMPORAL)) {
		
		// 3 dimensional matching the Video Width x Height x 2
		std::vector<unsigned int> fvDimensions;
		fvDimensions.push_back(displayWidth);
		fvDimensions.push_back(displayHeight);
		fvDimensions.push_back(2);
		
		myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
									  displayWidth, displayHeight, // display size
									  3); 						   // input vector size (x, y, t)
		
        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Input Vector
		std::vector<int> iv;
		iv.push_back(0);
		myDisplay->UpdateInputVector(iv);
		
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = 0x00; p.a = 0x7f;
		std::vector<unsigned int> start, end;
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
		
		start.push_back(0); start.push_back(0); 
		end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
		
		myDisplay->FillCoefficientMatrix(coeffs, start, end);
        
    // Coefficient Plane Blending
    } else if ((config == SIMPLE) && (configBlend == COEFFICIENT_PLANE)) {
		
		// 3 dimensional matching the Video Width x Height x 2
		std::vector<unsigned int> fvDimensions;
		fvDimensions.push_back(displayWidth);
		fvDimensions.push_back(displayHeight);
		fvDimensions.push_back(2);
		
		myBlendingDisplay = new BlendingGlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                                      displayWidth, displayHeight, // display size
                                                      3, 						   // input vector size (x, y, 1)
                                                      2);                          // 2 coefficient planes
        myDisplay = myBlendingDisplay;
		
        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Input Vector
		std::vector<int> iv;
		iv.push_back(1);
		myDisplay->UpdateInputVector(iv);
		
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = 0x00; p.a = 0x7f;
		std::vector<unsigned int> start, end;
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
        
    // Flat-Tiled
	} else if (config == FLAT) {

		// 2 dimensional matching the Video Width x Height
		std::vector<unsigned int> fvDimensions;
		
		fvDimensions.push_back(displayWidth);
		fvDimensions.push_back(displayHeight);
#ifndef NO_CL
		myDisplay = new ClNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      2); 						   // input vector size (x and y only)
#else
		myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
                                      displayWidth, displayHeight, // display size
                                      2); 						   // input vector size (x and y only)
#endif

        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = p.a = 0xff;
		std::vector<unsigned int> start, end;
		start.push_back(0); start.push_back(0);
		end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
		myDisplay->FillPixel(p, start, end);
		
        // Set up Flat Tiler and initialize Coefficient plane`
        myTiler = new FlatTiler(myDisplay,
                                configTileWidth, configTileHeight, configSigBits,
                                configHeadless);
        
        myTiler->InitializeCoefficientPlane();
    // Simple Framebuffer
	} else {
		
		
		// 2 dimensional matching the Video Width x Height
		std::vector<unsigned int> fvDimensions;
		fvDimensions.push_back(displayWidth);
		fvDimensions.push_back(displayHeight);
		
#ifndef NO_CL
		myDisplay = new ClNddiDisplay(fvDimensions,                // framevolume dimensional sizes
									  displayWidth, displayHeight, // display size
									  2); 						   // input vector size (x and y only)
#else
		myDisplay = new GlNddiDisplay(fvDimensions,                // framevolume dimensional sizes
									  displayWidth, displayHeight, // display size
									  2); 						   // input vector size (x and y only)
#endif		
        // Grab the cost model
        costModel = myDisplay->GetCostModel();
        
		// Initialize Frame Volume
		nddi::Pixel p;
		p.r = p.g = p.b = p.a = 0xff;
		std::vector<unsigned int> start, end;
		start.push_back(0); start.push_back(0);
		end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
		myDisplay->FillPixel(p, start, end);
		
        // Initialize Coefficient Plane
        vector< vector<int> > coeffs;
        coeffs.resize(2);
        coeffs[0].push_back(1); coeffs[0].push_back(0);
        coeffs[1].push_back(0); coeffs[1].push_back(1);
        
        start.clear(); end.clear();
        
        start.push_back(0); start.push_back(0);
        end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
        
        myDisplay->FillCoefficientMatrix(coeffs, start, end);
    }
	totalUpdates++;
    
    // Renders the initial display
    myDisplay->GetFrameBuffer();
}


void updateDisplay(uint8_t* buffer, size_t width, size_t height) {

    // CACHE or FLAT
	if ( (config == CACHE) || (config == FLAT) ) {
		// Update the display with the Tiler
		myTiler->UpdateDisplay(buffer, width, height);
    // SIMPLE
	} else {
		// Update the display like a simple Frame Buffer
		size_t bufferPos = 0;
		
        // Array of pixels used for framebuffer mode.
        nddi::Pixel* frameBuffer = (nddi::Pixel*)malloc(sizeof(nddi::Pixel) * width * height);

		// Transform the buffer into pixels
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				frameBuffer[j * width + i].r = buffer[bufferPos++];
				frameBuffer[j * width + i].g = buffer[bufferPos++];
				frameBuffer[j * width + i].b = buffer[bufferPos++];
				frameBuffer[j * width + i].a = 0xff;
			}
		}
		
		// Update the frame volume
		std::vector<unsigned int> start, end, dest;

        // Not Blending
        if (configBlend == NONE) {
            // Just send the pixels to the single plane
            start.push_back(0); start.push_back(0);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1);
            myDisplay->CopyPixels(frameBuffer, start, end);
            
        // Frame Volume Blending
        } else if (configBlend == FRAME_VOLUME) {
            
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
            std::vector<int> iv;
            iv.push_back(plane);
            myDisplay->UpdateInputVector(iv);
        // Temporal Blending
        } else if (configBlend == TEMPORAL) {
            
            // Send the pixels to the 0 plane
            start.push_back(0); start.push_back(0); start.push_back(0);
            end.push_back(displayWidth - 1); end.push_back(displayHeight - 1); end.push_back(0);
            myDisplay->CopyPixels(frameBuffer, start, end);
            
            // Render the new frame and then render the black frame.
            // Note: This is a loose simulation. Actual temporal blending would carefully
            // drive the input vector, perhaps flipping back and forth several times per frame.
            std::vector<int> iv;
            for (int c = 0; c < TEMPORAL_FLIP_PER_FRAME_COUNT; c++) {
                iv.push_back(0);
                myDisplay->UpdateInputVector(iv);
                iv.push_back(1);
                myDisplay->UpdateInputVector(iv);
            }
        // Coefficient Plane Blending
        } else if (configBlend == COEFFICIENT_PLANE) {
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


long calculateSE(uint8_t* videoIn, nddi::Pixel* frameOut) {
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
	
	if (config > COUNT) {
		
		// Grab the frame buffer from the NDDI display
        GLuint texture = myDisplay->GetFrameBuffer();
		
		// If we used a lossy mode, then calculate the Square Error for this frame
#if 0   // TODO(CDE): Add this back! Will likely be done in the CachedTiler
		if (videoBuffer && (config == CACHE) && (configSigBits < 8)) {
			totalSE += calculateSE(videoBuffer, frameBuffer);
		}
#endif
        
        if (!configHeadless) {
            
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
    
    double MSE, PSNR;
    long totalCost = costModel->getLinkBytesTransmitted();
    
    MSE = (double)totalSE / (double)(displayWidth * displayHeight * totalUpdates) / 3.0f;
    PSNR = 10.0f * log10(65025.0f / MSE);
    
    gettimeofday(&endTime, NULL);
    
    // Configuration
    //
    cout << "Configuration Information:" << endl;
    
    switch (config) {
        case SIMPLE:
            cout << "  Configuring NDDI as a simple Framebuffer." << endl;
            break;
        case FLAT:
            cout << "  Configuring NDDI as Flat Tiled." << endl;
            cout << "  Tile Dimensions: " << configTileWidth << "x" << configTileHeight << endl;
            break;
        case CACHE:
            cout << "  Configuring NDDI as Cached Tiled." << endl;
            cout << "  Tile Dimensions: " << configTileWidth << "x" << configTileHeight << endl;
            break;
        case COUNT:
            cout << "  Configuring NDDI to just count." << endl;
            break;
        default:
            break;
    }
    
    cout << "  File: " << fileName << endl;
    
    // Transmission
    //
    cout << "Transmission Statistics:" << endl;
    
    cout << "  Total Pixel Data Updated (bytes): " << totalUpdates * displayWidth * displayHeight * 4 <<
    " Total NDDI Cost (bytes): " << totalCost <<  
    " Ratio: " << (double)totalCost / (double)totalUpdates / (double)displayWidth / (double)displayHeight / 4.0f << endl;
    
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
    
    
    // Pixel
    //
    cout << "Pixel Statistics:" << endl;
    cout << "  Pixel Mappings: " << costModel->getPixelsMapped() << endl;
    cout << "  Pixel Blends: " << costModel->getPixelsBlended() << endl;
    
    // Performance
    //
    cout << "Performance Statistics:" << endl;
    cout << "  Average FPS: " << (double)totalUpdates / ((double)(endTime.tv_sec * 1000000
                                                                  + endTime.tv_usec
                                                                  - startTime.tv_sec * 1000000
                                                                  - startTime.tv_usec) / 1000000.0f) << endl;
    
    // Quality
    //
    cout << "Quality Statistics:" << endl << "  Total PSNR: " << PSNR << endl;
    
    // CSV
    //
    cout << "CSV" << endl;
    cout << "Commands Sent\tBytes Transmitted\tIV Num Reads\tIV Bytes Read\tIV Num Writes\tIV Bytes Written\tCP Num Reads\tCP Bytes Read\tCP Num Writes\tCP Bytes Written\tFV Num Reads\tFV Bytes Read\tFV Num Writes\tFV Bytes Written\tFV Time\tPixels Mapped\tPixels Blended" << endl;
    cout
    << costModel->getLinkCommandsSent() << "\t" << costModel->getLinkBytesTransmitted() << "\t"
    << costModel->getReadAccessCount(INPUT_VECTOR_COMPONENT) << "\t" << costModel->getBytesRead(INPUT_VECTOR_COMPONENT) << "\t"
    << costModel->getWriteAccessCount(INPUT_VECTOR_COMPONENT) << "\t" << costModel->getBytesWritten(INPUT_VECTOR_COMPONENT) << "\t"
    << costModel->getReadAccessCount(COEFFICIENT_PLANE_COMPONENT) << "\t" << costModel->getBytesRead(COEFFICIENT_PLANE_COMPONENT) << "\t"
    << costModel->getWriteAccessCount(COEFFICIENT_PLANE_COMPONENT) << "\t" << costModel->getBytesWritten(COEFFICIENT_PLANE_COMPONENT) << "\t"
    << costModel->getReadAccessCount(FRAME_VOLUME_COMPONENT) << "\t" << costModel->getBytesRead(FRAME_VOLUME_COMPONENT) << "\t"
    << costModel->getWriteAccessCount(FRAME_VOLUME_COMPONENT) << "\t" << costModel->getBytesWritten(FRAME_VOLUME_COMPONENT) << "\t"
    << costModel->getTime(FRAME_VOLUME_COMPONENT) << "\t"
    << costModel->getPixelsMapped() << "\t" << costModel->getPixelsBlended() << endl;
    
    // Warnings about Configuration
#if defined(SUPRESS_EXCESS_RENDERING) || defined(NO_CL) || defined(NO_GL)
    cout << endl << "CONFIGURATION WARNINGS:" << endl;
#ifdef SUPRESS_EXCESS_RENDERING
    cout << "- Was compiled with SUPRESS_EXCESS_RENDERING, and so the numbers may be off. Recompile with \"make NO_HACKS=1\"." << endl;
#endif
#ifdef NO_CL
    cout << - "Was compiled without OpenCL." << endl;
#endif
#ifdef NO_GL
    cout << - "Was compiled without OpenGL." << endl;
#endif
#endif

    // Clean up
    if (exitNow) {
        delete myDisplay;
        if (myPlayer) { delete(myPlayer); }
        if (lastBuffer) { free(lastBuffer); }
        
        exit(0);
    }
}


void renderFrame() {
	
	static int framesDecoded = 0;
	static int framesRendered = 0;
	static int currentFrame = 0;
	static rewind_play_t rewindState = NOT_REWINDING;
	
	if (!configMaxFrames || (totalUpdates <= configMaxFrames)) {
		// If we're not in a rewind backwards or forward state...
		if (rewindState == NOT_REWINDING) {
			// Get a decoded frame
            // TODO(CDE): Don't poll. Use signals instead.
            while (bufferQueue.size() == 0) {
                usleep(5);
            }
            if ((videoBuffer = bufferQueue.front()) != NULL) {
                
                pthread_mutex_lock(&bufferQueueMutex);
                bufferQueue.pop();
                pthread_mutex_unlock(&bufferQueueMutex);
				framesDecoded++;

				// If we're past the designated start frame, then update the NDDI display
				if (framesDecoded >= configStartFrame) {
					// If configured to rewind, then check to see if this decoded frame should be stored.
					if ( myRewinder &&
						(currentFrame >= (configRewindStartFrame - configRewindFrames)) &&
						(currentFrame < configRewindStartFrame) ) {
						myRewinder->CopyFrame(videoBuffer, currentFrame - (configRewindStartFrame - configRewindFrames));
					}
					
					// Update NDDI Display
					updateDisplay(videoBuffer,
								  myPlayer->width(), myPlayer->height());
					
					framesRendered++;
					
					// Draw
                    glutPostRedisplay();
					
					
					// If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
					if (myRewinder && currentFrame == configRewindStartFrame) {
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
			updateDisplay(myRewinder->GetFrame(currentFrame - (configRewindStartFrame - configRewindFrames)),
						  myPlayer->width(), myPlayer->height());
			
			framesRendered++;
			
			// Draw
            glutPostRedisplay();
			
			// If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
			if (currentFrame == (configRewindStartFrame - configRewindFrames)) {
				rewindState = PLAY_FORWARD;
				currentFrame++;
			} else {
				currentFrame--;
			}
		// Else if we're in the rewind play forward state...
		} else if (rewindState == PLAY_FORWARD) {
			// Update NDDI Display
			updateDisplay(myRewinder->GetFrame(currentFrame - (configRewindStartFrame - configRewindFrames)),
						  myPlayer->width(), myPlayer->height());
			
			framesRendered++;
			
			// Draw
            glutPostRedisplay();
            
			// Increment the currentFrame
			currentFrame++;
			
			// If this is the final frame stored, then set the rewindState to PLAY_BACKWARDS
			if (currentFrame == configRewindStartFrame) {
				rewindState = NOT_REWINDING;
				delete(myRewinder);
			}

		}

		if (!configHeadless) {
			cout << "PixelBidge Statistics:" << endl << "  Decoded Frames: " << framesDecoded << " - Rendered Frames: " << framesRendered << endl;
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

	if (!configMaxFrames || (totalUpdates < configMaxFrames)) {
		int pixel_count = myPlayer->width() * myPlayer->height();
		
		// Get a decoded frame
        // TODO(CDE): Don't poll. Use signals instead.
        while (bufferQueue.size() == 0) {
            usleep(5);
        }
		if ((videoBuffer = bufferQueue.front()) != NULL) {
			
            pthread_mutex_lock(&bufferQueueMutex);
            bufferQueue.pop();
            pthread_mutex_unlock(&bufferQueueMutex);
			framesDecoded++;
			
			if (framesDecoded > configStartFrame) {
				// If we're in a rewind period, then we'll be incrementing diffs as well as frame count by 3 instead of 1
				int inc;
				if ( (framesDecoded > (configRewindStartFrame - configRewindFrames)) && (framesDecoded <= configRewindStartFrame) ) {
					inc = 3;
				} else {
					inc = 1;
				}
				
				// If we haven't previously decoded a frame
				if (lastBuffer == NULL) {
                    costModel->registerTransmissionCharge(4 * inc * pixel_count, 0);

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
                    costModel->registerTransmissionCharge(4 * diffs, 0);

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
		
		if (!configHeadless) {
			cout << "PixelBidge Statistics:" << endl << "  Decoded Frames: " << framesDecoded << " - Counted Frames: " << framesCounted << endl;
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
	cout << "pixelbridge [--mode <fb|flat|cache|count>] [--blend <fv|t|cp|>] [--ts <n> <n>] [--tc <n>] [--bits <1-8>] [--start <n>] [--frames <n>] [--rewind <n> <n>] <filename>" << endl;
	cout << endl;
	cout << "  --mode  Configure NDDI as a framebuffer (fb), as a flat tile array (flat), or as a cached tile (cache).\n" <<
	        "          Optional the mode can be set to count the number of pixels changed (count)." << endl;
	cout << "  --blend  For fb mode, this option alpha-blends using just the frame volume, input vector (temporal), or coefficient plane." << endl;
	cout << "  --ts  Sets the tile size to the width and height provided." << endl;
	cout << "  --tc  Sets the maximum number of tiles in the cache." << endl;
	cout << "  --bits  Sets the number of significant bits per channel when computing checksums." << endl;
	cout << "  --start  Will start with this frame, ignoring any decoded frames prior to it." << endl;
	cout << "  --frames  Sets the number of maximum frames that are decoded." << endl;
	cout << "  --rewind  Sets a start point and a number of frames to play in reverse. Once finished, normal playback continues." << endl;
	cout << "  --headless  Removes rendering and excessive data output. Great for batch processing." << endl;
}


bool parseArgs(int argc, char *argv[]) {
	argc--;
	argv++;
	
	while (argc) {
		if (strcmp(*argv, "--mode") == 0) {
			argc--;
			argv++;
			if (strcmp(*argv, "fb") == 0) {
				config = SIMPLE;
			} else if (strcmp(*argv, "flat") == 0) {
				config = FLAT;
			} else if (strcmp(*argv, "cache") == 0) {
				config = CACHE;
			} else if (strcmp(*argv, "count") == 0) {
				config = COUNT;
			}
			argc--;
			argv++;
		} else if (strcmp(*argv, "--blend") == 0) {
			if (config != SIMPLE) {
				showUsage();
				return false;
			}
			argc--;
			argv++;
            if (strcmp(*argv, "fv") == 0) {
                configBlend = FRAME_VOLUME;
            } else if (strcmp(*argv, "t") == 0) {
                configBlend = TEMPORAL;
            } else if (strcmp(*argv, "cp") == 0) {
                configBlend = COEFFICIENT_PLANE;
            } else {
				showUsage();
				return false;
            }
			argc--;
			argv++;
		} else if (strcmp(*argv, "--ts") == 0) {
			configTileWidth = atoi(argv[1]);
			configTileHeight = atoi(argv[2]);
			if ((configTileWidth == 0) || (configTileHeight == 0)) {
				showUsage();
				return false;
			}
			argc -= 3;
			argv += 3;
		} else if (strcmp(*argv, "--tc") == 0) {
			configMaxTiles = atoi(argv[1]);
			if (configMaxTiles == 0) {
				showUsage();
				return false;
			}
			argc -= 2;
			argv += 2;
		} else if (strcmp(*argv, "--bits") == 0) {
			configSigBits = atoi(argv[1]);
			if ( (configSigBits == 0) || (configSigBits > 8) ) {
				showUsage();
				return false;
			}
			argc -= 2;
			argv += 2;
		} else if (strcmp(*argv, "--start") == 0) {
			configStartFrame = atoi(argv[1]);
			if (configStartFrame == 0) {
				showUsage();
				return false;
			}
			argc -= 2;
			argv += 2;
		} else if (strcmp(*argv, "--frames") == 0) {
			configMaxFrames = atoi(argv[1]);
			if (configMaxFrames == 0) {
				showUsage();
				return false;
			}
			argc -= 2;
			argv += 2;
		} else if (strcmp(*argv, "--rewind") == 0) {
			configRewindStartFrame = atoi(argv[1]);
			configRewindFrames = atoi(argv[2]);
			if ((configRewindStartFrame == 0) || (configRewindFrames == 0) || (configRewindFrames > configRewindStartFrame)) {
				showUsage();
				return false;
			}
			argc -= 3;
			argv += 3;
		} else if (strcmp(*argv, "--headless") == 0) {
			configHeadless = true;
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


int main(int argc, char *argv[]) {

	// Parse command line arguments
	if (!parseArgs(argc, argv)) {
		return -1;
	}
	
	// Initialize ffmpeg and set dimensions
	myPlayer = new FfmpegPlayer(fileName);
	displayWidth = myPlayer->width();
	displayHeight = myPlayer->height();
	if (configRewindStartFrame > 0) {
		myRewinder = new Rewinder(configRewindStartFrame - configRewindFrames, displayWidth, displayHeight);
	}
	if ( (displayWidth == 0) || (displayHeight == 0) ) {
		cout << "Error: Could not get video dimensions." << endl;
		return -1;
	}
    
    // Start decoded thread
    pthread_create(&decoderThread, NULL, decoderRun, NULL);
    
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
    
	if (config > COUNT) {
		glutIdleFunc(renderFrame);
	} else {
		glutIdleFunc(countChangedPixels);
	}
    
	// Setup NDDI Display
    if (config > COUNT) {
		// If the tile size wasn't specified, then dynamically calculate it
		if ((configTileWidth == 0) || (configTileHeight == 0)) {
			configTileWidth = configTileHeight = ((displayWidth > displayHeight) ? displayWidth : displayHeight) / 40;
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
