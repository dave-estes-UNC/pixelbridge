//
//  CostModel.h
//  pixelbridge
//
//  Created by Dave Estes on 10/18/11.
//  Copyright 2011 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CostModel_h
#define pixelbridge_CostModel_h

#include <vector>

using namespace std;

namespace nddi {
    
    typedef enum {
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
        unsigned int           sequenceNumber;
        charge_type_t  chargeType;
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
        
        void registerBulkMemoryCharge(component_t component, long accessCount, memory_access_t access, void* address, long numBytes) {
            
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
                    break;
                default:
                    break;
            }
        }
        
        void registerTransmissionCharge(long numBytes) {
#pragma omp atomic
            linkCommandsSent++;
            
#pragma omp atomic
            linkBytesSent += numBytes;

        }
        
        void registerPixelBlendCharge(long numBlends) {
#pragma omp atomic
            pixelsBlended += numBlends;
        }
        
        void registerPixelMappingCharge(long numMappings) {
#pragma omp atomic
            pixelsMapped += numMappings;
        }
        
        long getLinkCommandsSent() {
            return linkCommandsSent;
        }

        long getLinkBytesTransmitted() {
            return linkBytesSent;
        }

        long getReadAccessCount(component_t component) {
            
            long count = 0;
            
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

        long getWriteAccessCount(component_t component) {
            
            long count = 0;
            
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

        long getBytesRead(component_t component) {
            
            long count = 0;
            
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
        
        long getBytesWritten(component_t component) {
            
            long count = 0;
            
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

        long getTime(component_t component) {

            long count = 0;

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorTime;
                    break;
                case COEFFICIENT_PLANE_COMPONENT:
                    count = coefficientPlaneTime;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeTime;
                    break;
                default:
                    break;
            }
            return count;
        }

        long getPixelsBlended() {
            return pixelsBlended;
        }
        
        long getPixelsMapped() {
            return pixelsMapped;
        }
    };
}
#endif
