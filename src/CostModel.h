//
//  CostModel.h
//  pixelbridge
//
//  Created by Dave Estes on 10/18/11.
//  Copyright 2011 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CostModel_h
#define pixelbridge_CostModel_h

/*
 * Definitions for widths for the various
 * data in the NDDI display. Modify to perform
 * various cost experiments.
 */
#define BYTES_PER_PIXEL     3
#define BYTES_PER_COORD     4
#define BYTES_PER_IV_VALUE  4
#define BYTES_PER_COEFF     4
#define BYTES_PER_SCALER    4

/*
 * Helper macros to be used when registering
 * cost charges.
 */
#define CALC_BYTES_FOR_PIXELS(c)              (BYTES_PER_PIXEL * c)
#define CALC_BYTES_FOR_FV_COORD_TUPLES(c)     (BYTES_PER_COORD * frameVolumeDimensionalSizes_.size() * c)
#define CALC_BYTES_FOR_TILE_COORD_DOUBLES(c)  (BYTES_PER_COORD * 2 * c)
#define CALC_BYTES_FOR_IV_UPDATE()            (BYTES_PER_IV_VALUE * (inputVector_->getSize() - 2))
#define CALC_BYTES_FOR_CMS(c)                 (BYTES_PER_COEFF * inputVector_->getSize() * frameVolumeDimensionalSizes_.size() * c)
#define CALC_BYTES_FOR_CM_COORD_DOUBLES(c)    (BYTES_PER_COORD * 2 * c)
#define CALC_BYTES_FOR_CP_COORD_TRIPLES(c)    (BYTES_PER_COORD * 3 * c)

namespace nddi {

    typedef enum {
    	NDDI_LINK_COMPONENT,
        INPUT_VECTOR_COMPONENT,
        COEFFICIENT_PLANE_COMPONENT,
        FRAME_VOLUME_COMPONENT
    } component_t;
    
    /**
     * Types of charges. Each NDDI Command is first assessed a link charge,
     * The the memory, pixel blends, and pixel mappings resulting from that command
     * are then charged.
     */
    typedef enum {
        LINK_CHARGE,
        MEMORY_CHARGE,
        PIXEL_BLEND_CHARGE,
        PIXEL_MAPPING_CHARGE
    } charge_type_t;
    
    typedef enum {
        READ_ACCESS,
        WRITE_ACCESS
    } memory_access_t;
    
    /**
     * Represents number of bytes read from or written to a component.
     */
    typedef struct {
        component_t      component;
        memory_access_t  access;
        void*            address;
        unsigned long    numBytes;
        unsigned long    time;
    } memory_charge_t;

    /**
     * Represents the number of bytes sent over the NDDI Link.
     */
    typedef struct {
        unsigned long    numBytes;
    } link_charge_t;

    /**
     * Represents the number of pixel blend operations where each operation alpha
     * blends two pixels.
     */
    typedef struct {
        unsigned long     numBlends;
    } pixel_blend_charge_t;

    /**
     * Represents the number of pixel mapping operations, where each operation is
     * a matrix multiplication.
     */
    typedef struct {
        unsigned long     numMappings;
    } pixel_mapping_charge_t;

    /**
     * Holds an instance of a type of charge.
     */
    typedef struct {
        unsigned int      sequenceNumber;
        charge_type_t     chargeType;
        union {
            memory_charge_t         memory;
            link_charge_t           link;
            pixel_blend_charge_t    pixel_blend;
            pixel_mapping_charge_t  pixel_mapping;
        } u;
    } charge_t;

    
    /**
     * The CostModel allows different types of charges to be made and will run reports
     * later.
     */
    class CostModel {
        
    private:

        unsigned long linkCommandsSent;
        unsigned long linkBytesSent;
        unsigned long linkTime;

        unsigned long pixelsBlended;
        unsigned long pixelBlendingTime;
        
        unsigned long pixelsMapped;
        unsigned long pixelMappingTime;
        
        unsigned long inputVectorReads;
        unsigned long inputVectorWrites;
        unsigned long inputVectorBytesRead;
        unsigned long inputVectorBytesWritten;
        unsigned long inputVectorTime;

        unsigned long coefficientPlaneReads;
        unsigned long coefficientPlaneWrites;
        unsigned long coefficientPlaneBytesRead;
        unsigned long coefficientPlaneBytesWritten;
        unsigned long coefficientPlaneTime;

        unsigned long frameVolumeReads;
        unsigned long frameVolumeWrites;
        unsigned long frameVolumeBytesRead;
        unsigned long frameVolumeBytesWritten;
        unsigned long frameVolumeTime;

    public:
        
        CostModel()
        : linkCommandsSent(0),
          linkBytesSent(0),
          linkTime(0),
          pixelsBlended(0),
          pixelBlendingTime(0),
          pixelsMapped(0),
          pixelMappingTime(0),
          inputVectorReads(0),
          inputVectorWrites(0),
          inputVectorBytesRead(0),
          inputVectorBytesWritten(0),
          inputVectorTime(0),
          coefficientPlaneReads(0),
          coefficientPlaneWrites(0),
          coefficientPlaneBytesRead(0),
          coefficientPlaneBytesWritten(0),
          coefficientPlaneTime(0),
          frameVolumeReads(0),
          frameVolumeWrites(0),
          frameVolumeBytesRead(0),
          frameVolumeBytesWritten(0),
          frameVolumeTime(0) {
        }
        
        void registerMemoryCharge(component_t component,
                                  memory_access_t access,
                                  void* address,
                          	      unsigned long numBytes,
                                  unsigned long time) {

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        inputVectorReads++;
#pragma omp atomic
                        inputVectorBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        inputVectorWrites++;
#pragma omp atomic
                        inputVectorBytesWritten += numBytes;
#pragma omp atomic
                        inputVectorTime += time;
                    }
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        coefficientPlaneReads++;
#pragma omp atomic
                        coefficientPlaneBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        coefficientPlaneWrites++;
#pragma omp atomic
                        coefficientPlaneBytesWritten += numBytes;
                    }
#pragma omp atomic
                    coefficientPlaneTime += time;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        frameVolumeReads++;
#pragma omp atomic
                        frameVolumeBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        frameVolumeWrites++;
#pragma omp atomic
                        frameVolumeBytesWritten += numBytes;
                    }
#pragma omp atomic
                    frameVolumeTime += time;
                    break;
                default:
                    break;
            }
        }
        
        void registerBulkMemoryCharge(component_t component,
                                      unsigned long accessCount,
                                      memory_access_t access,
                                      void* address,
                                      unsigned long numBytes,
                                      unsigned long time) {
            
            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        inputVectorReads += accessCount;
#pragma omp atomic
                        inputVectorBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        inputVectorWrites += accessCount;
#pragma omp atomic
                        inputVectorBytesWritten += numBytes;
                    }
                    inputVectorTime += time;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        coefficientPlaneReads += accessCount;
#pragma omp atomic
                        coefficientPlaneBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        coefficientPlaneWrites += accessCount;
#pragma omp atomic
                        coefficientPlaneBytesWritten += numBytes;
                    }
                    coefficientPlaneTime += time;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        frameVolumeReads += accessCount;
#pragma omp atomic
                        frameVolumeBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        frameVolumeWrites += accessCount;
#pragma omp atomic
                        frameVolumeBytesWritten += numBytes;
                    }
                    frameVolumeTime += time;
                    break;
                default:
                    break;
            }
        }
        
        void registerTransmissionCharge(unsigned long numBytes, unsigned long time) {
#pragma omp atomic
            linkCommandsSent++;
            
#pragma omp atomic
            linkBytesSent += numBytes;

#pragma omp atomic
            linkTime += time;
        }
        
        void registerPixelBlendCharge(unsigned long numBlends) {
#pragma omp atomic
            pixelsBlended += numBlends;
        }
        
        void registerPixelMappingCharge(unsigned long numMappings) {
#pragma omp atomic
            pixelsMapped += numMappings;
        }
        
        unsigned long getLinkCommandsSent() {
            return linkCommandsSent;
        }

        unsigned long getLinkBytesTransmitted() {
            return linkBytesSent;
        }

        unsigned long getReadAccessCount(component_t component) {
            
            unsigned long count = 0;
            
            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorReads;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    count = coefficientPlaneReads;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeReads;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getWriteAccessCount(component_t component) {
            
            unsigned long count = 0;
            
            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorWrites;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    count = coefficientPlaneWrites;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeWrites;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getBytesRead(component_t component) {
            
        	unsigned long count = 0;
            
            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorBytesRead;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    count = coefficientPlaneBytesRead;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeBytesRead;
                    break;
                default:
                    break;
            }
            return count;
        }
        
        unsigned long getBytesWritten(component_t component) {
            
        	unsigned long count = 0;
            
            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorBytesWritten;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    count = coefficientPlaneBytesWritten;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeBytesWritten;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getTime(component_t component) {

        	unsigned long time = 0;

            switch (component) {
            	case NDDI_LINK_COMPONENT:
            		time = linkTime;
            		break;
                case INPUT_VECTOR_COMPONENT:
                    time = inputVectorTime;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    time = coefficientPlaneTime;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    time = frameVolumeTime;
                    break;
                default:
                    break;
            }
            return time;
        }

        unsigned long getPixelsBlended() {
            return pixelsBlended;
        }
        
        unsigned long getPixelsMapped() {
            return pixelsMapped;
        }
    };
}
#endif
