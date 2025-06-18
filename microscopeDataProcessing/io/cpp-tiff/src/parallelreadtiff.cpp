#include <cstdint>
#include <cmath>
#include <cstring>
#include <vector>
#include <limits.h>
#include <omp.h>
#include "tiffio.h"
#include "../src/helperfunctions.h"


// Backup method in case there are errors reading strips
uint8_t readTiffParallelBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
	uint8_t err = 0;
	char errString[10000];
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){
		if(err) continue;
        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif){
			sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
			err = 1;
		}
        void* buffer = malloc(x*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice || err) break;

            int counter = 0;
            while(!TIFFSetDirectory(tif, (uint64_t)dir) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,dir,counter+1);
                counter++;
            }

            for (int64_t i = 0; i < y; i++)
            {
                TIFFReadScanline(tif, buffer, i, 0);
                if(!flipXY){
                    // Probably need to fix this
                    memcpy(tiff+(((i*y)+((dir-startSlice)*(x*y)))*bytes),buffer,x*bytes);
                    continue;
                }
                //loading the data into a buffer
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j];
                        }
                            break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j];
                        }
                            break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((float*)buffer)[j];
                        }
                            break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((double*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((double*)buffer)[j];
                        }
                            break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
	if(err){
		printf(errString);
	}
	return err;
}

uint8_t readTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    uint8_t errBak = 0;
    uint8_t err = 0;
    char errString[10000];
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        uint8_t outCounter = 0;
        TIFF* tif = TIFFOpen(fileName, "r");
        while(!tif){
            tif = TIFFOpen(fileName, "r");
            if(outCounter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                }
                continue;
            }
            outCounter++;
        }

        void* buffer = malloc(x*stripSize*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice || err) break;

            uint8_t counter = 0;
            while(!TIFFSetDirectory(tif, (uint64_t)dir) && counter<3){
                counter++;
                if(counter == 3){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                    }
                }
            }
            if(err) break;
            for (int64_t i = 0; i*stripSize < y; i++)
            {

                //loading the data into a buffer
                int64_t cBytes = TIFFReadEncodedStrip(tif, i, buffer, stripSize*x*bytes);
                if(cBytes < 0){
                    #pragma omp critical
                    {
                        errBak = 1;
                        err = 1;
                        sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,i);
                    }
                    break;
                }
                if(!flipXY){
                    memcpy(tiff+(((i*stripSize*x)+((dir-startSlice)*(x*y)))*bytes),buffer,cBytes);
                    continue;
                }
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((float*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((float*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((double*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((double*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
    if(err){
        if(errBak) return readTiffParallelBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else {
			printf(errString);
		}
    }
	return err;
}

// Backup method in case there are errors reading strips
uint8_t readTiffParallel2DBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (y-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
	uint8_t err = 0;
	char errString[10000];
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){
		if(err) continue;
        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif) {
			sprintf(errString,"tiff:threadError","Thread %d: File \"%s\" cannot be opened\n",w,fileName);
			err = 1;
		}
        void* buffer = malloc(x*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice || err) break;

            int counter = 0;
            while(!TIFFSetDirectory(tif, startSlice) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,dir,counter+1);
                counter++;
            }

            for (int64_t i = (w*batchSize); i < ((w+1)*batchSize); i++)
            {
                if(i >= y) break;
                TIFFReadScanline(tif, buffer, i, 0);
                if(!flipXY){
                    memcpy(tiff+((i*x)*bytes),buffer,x*bytes);
                    continue;
                }
                //loading the data into a buffer
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j];
                        }
                            break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j];
                        }
                            break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((float*)buffer)[j];
                        }
                            break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((double*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((double*)buffer)[j];
                        }
                            break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
	if(err){
		printf(errString);
	}
	return err;
}


uint8_t readTiffParallel2D(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    int32_t batchSize = (stripsPerDir-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    uint8_t err = 0;
    uint8_t errBak = 0;
    char errString[10000];
    uint16_t compressed = 1;
    TIFF* tif = TIFFOpen(fileName, "r");
    TIFFGetField(tif, TIFFTAG_COMPRESSION, &compressed);

    // The other method won't work on specific slices of 3D images for now
    // so start slice must also be 0
    if(numWorkers > 1 || compressed > 1){
        #pragma omp parallel for
        for(w = 0; w < numWorkers; w++){
    
            uint8_t outCounter = 0;
            TIFF* tif = TIFFOpen(fileName, "r");
            while(!tif){
                tif = TIFFOpen(fileName, "r");
                if(outCounter == 3){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                    }
                    continue;
                }
                outCounter++;
            }
    
            void* buffer = malloc(x*stripSize*bytes);
    
    
            uint8_t counter = 0;
            while(!TIFFSetDirectory(tif, startSlice) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,0,counter+1);
                counter++;
                if(counter == 3){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                    }
                }
            }
            for (int64_t i = (w*batchSize); i < (w+1)*batchSize; i++)
            {
                if(i*stripSize >= y || err) break;
                //loading the data into a buffer
                int64_t cBytes = TIFFReadEncodedStrip(tif, i, buffer, stripSize*x*bytes);
                if(cBytes < 0){
                    #pragma omp critical
                    {
                        errBak = 1;
                        err = 1;
                        sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,i);
                    }
                    break;
                }
                if(!flipXY){
                    memcpy(tiff+((i*stripSize*x)*bytes),buffer,cBytes);
                    continue;
                }
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint8_t*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint16_t*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((float*)tiff)[((j*y)+(k+(i*stripSize)))] = ((float*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((double*)tiff)[((j*y)+(k+(i*stripSize)))] = ((double*)buffer)[j+(k*x)];
                            }
                        }
                                break;
                }
            }
            free(buffer);
            TIFFClose(tif);
        }
    }
    else{
        void* tiffC = NULL;
        FILE *fp = fopen(fileName, "rb");
        if(!fp){ 
			printf("File \"%s\" cannot be opened from Disk\n",fileName);
			err = 1;
			return err;
		}

        if(!tif){ 
			printf("File \"%s\" cannot be opened\n",fileName);
			err = 1;
			return err;
		}
        
		uint64_t offset = 0;
        uint64_t* offsets = NULL;
        TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
        if(!offsets){ 
			printf("Could not get offsets from the tiff file\n");
       		err = 1;
			return err;
		}
		offset = offsets[0];
        uint64_t zSize = x*y*bytes;
    
        fseek(fp, offset, SEEK_SET);


        TIFFClose(tif);
        
        if(!flipXY){
            fread(tiff, 1, zSize, fp);
        }
        else{
            uint64_t size = x*y*z*(bits/8);
            tiffC = malloc(size);
            fread(tiffC, 1, zSize, fp);
        }
        fclose(fp);
        if(flipXY){   
            for(uint64_t k = 0; k < z; k++){
                for(uint64_t j = 0; j < x; j++){
                    for(uint64_t i = 0; i < y; i++){
                        switch(bits){
                            case 8:
                                ((uint8_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint8_t*)tiffC)[j+(i*x)+(k*x*y)];
                                break;
                            case 16:
                                ((uint16_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint16_t*)tiffC)[j+(i*x)+(k*x*y)];
                                break;
                            case 32:
                                ((float*)tiff)[i+(j*y)+(k*x*y)] = ((float*)tiffC)[j+(i*x)+(k*x*y)];
                                break;
                            case 64:
                                ((double*)tiff)[i+(j*y)+(k*x*y)] = ((double*)tiffC)[j+(i*x)+(k*x*y)];
                                break;
                        }
                    }
                }
            }
            free(tiffC);
        }
    }

    if(err) {
        if(errBak) return readTiffParallel2DBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else printf(errString);
    }
	return err;
}


// Reading images saved by ImageJ
uint8_t readTiffParallelImageJ(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    uint8_t err = 0;
    FILE *fp = fopen(fileName, "rb");
    if(!fp){ 
		printf("File \"%s\" cannot be opened from Disk\n",fileName);
		err = 1;
		return err;
	}
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){ 
		printf("File \"%s\" cannot be opened\n",fileName);
		err = 1;
		return err;
	}
    uint64_t offset = 0;
    uint64_t* offsets = NULL;
    TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
    if(offsets) offset = offsets[0];

    TIFFClose(tif);

    fseek(fp, offset, SEEK_SET);

    uint64_t bytes = bits/8;
    //#pragma omp parallel for
    /*
    for(uint64_t i = 0; i < z; i++){
    uint64_t cOffset = x*y*bytes*i;
    //pread(fd,tiff+cOffset,x*y*bytes,offset+cOffset);
    read(fd,tiff+cOffset,x*y*bytes);
    }*/
    uint64_t chunk = 0;
    uint64_t tBytes = x*y*z*bytes;
    uint64_t bytesRead;
    uint64_t rBytes = tBytes;

    // Can probably read more than INT_MAX now that we use fread
    if(tBytes < INT_MAX) bytesRead = fread(tiff,1,tBytes,fp);
    else{
        while(chunk < tBytes){
            rBytes = tBytes-chunk;
            if(rBytes > INT_MAX) bytesRead = fread(tiff+chunk,1,INT_MAX,fp);
            else bytesRead = fread(tiff+chunk,1,rBytes,fp);
            chunk += bytesRead;
        }
    }
    fclose(fp);
    // Swap endianess for types greater than 8 bits
    // TODO: May need to change later because we may not always need to swap
    if(bits > 8){
        #pragma omp parallel for
        for(uint64_t i = 0; i < x*y*z; i++){
            switch(bits){
                case 16:
                    //((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] & 0xff) >> 8) | (((uint16_t*)tiff)[i] << 8);
                    //((uint16_t*)tiff)[i] = bswap_16(((uint16_t*)tiff)[i]);
                    ((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] << 8) & 0xff00) | ((((uint16_t*)tiff)[i] >> 8) & 0x00ff);
                    break;
                case 32:
                    //((num & 0xff000000) >> 24) | ((num & 0x00ff0000) >> 8) | ((num & 0x0000ff00) << 8) | (num << 24)
                    //((float*)tiff)[i] = bswap_32(((float*)tiff)[i]);
                    ((uint32_t*)tiff)[i] = ((((uint32_t*)tiff)[i] << 24) & 0xff000000 ) |
                        ((((uint32_t*)tiff)[i] <<  8) & 0x00ff0000 ) |
                        ((((uint32_t*)tiff)[i] >>  8) & 0x0000ff00 ) |
                        ((((uint32_t*)tiff)[i] >> 24) & 0x000000ff );
                    break;
                case 64:
                    //((double*)tiff)[i] = bswap_64(((double*)tiff)[i]);
                    ((uint64_t*)tiff)[i] = ( (((uint64_t*)tiff)[i] << 56) & 0xff00000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 40) & 0x00ff000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 24) & 0x0000ff0000000000UL ) |
                        ( (((uint64_t*)tiff)[i] <<  8) & 0x000000ff00000000UL ) |
                        ( (((uint64_t*)tiff)[i] >>  8) & 0x00000000ff000000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 24) & 0x0000000000ff0000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 40) & 0x000000000000ff00UL ) |
                        ( (((uint64_t*)tiff)[i] >> 56) & 0x00000000000000ffUL );
                    break;
            }

        }
    }
    // Find a way to do this in-place without making a copy
    if(flipXY){
        uint64_t size = x*y*z*(bits/8);
        void* tiffC = malloc(size);
        memcpy(tiffC,tiff,size);
        #pragma omp parallel for
        for(uint64_t k = 0; k < z; k++){
            for(uint64_t j = 0; j < x; j++){
                for(uint64_t i = 0; i < y; i++){
                    switch(bits){
                        case 8:
                            ((uint8_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint8_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 16:
                            ((uint16_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint16_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 32:
                            ((float*)tiff)[i+(j*y)+(k*x*y)] = ((float*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 64:
                            ((double*)tiff)[i+(j*y)+(k*x*y)] = ((double*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                    }
                }
            }
        }
        free(tiffC);
    }
	return err;
}


// tiff pointer guaranteed to be NULL or the correct size array for the tiff file
void* readTiffParallelWrapperHelper(const char* fileName, void* tiff, uint8_t flipXY, const std::vector<uint64_t> &zRange = {})
{
	TIFFSetWarningHandler(DummyHandler);
	TIFF* tif = TIFFOpen(fileName, "r");
	if(!tif) return NULL;

	uint64_t x = 1,y = 1,z = 1,bits = 1, startSlice = 0;
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
    z = getImageSizeZ(fileName);

	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
	uint64_t stripSize = 1;
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &stripSize);
	TIFFClose(tif);

	// Check if image is an imagej image with imagej metadata
	// Get the correct
	uint8_t imageJIm = 0;
	if(isImageJIm(fileName)){
		imageJIm = 1;
		uint64_t tempZ = imageJImGetZ(fileName);
		if(tempZ) z = tempZ;
	}

    if(zRange.size()){
        if(zRange.size() == 2){
            startSlice = zRange[0];
            z = zRange[1];
        }
        else{
            startSlice = zRange[0];
            z = zRange[0]+1;
        }
    }


	if(imageJIm){
		if(bits == 8){
			if(!tiff) tiff = (uint8_t*)malloc(x*y*z*sizeof(uint8_t));
			readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize,flipXY);
			return (void*)tiff;
		}
		else if(bits == 16){
			if(!tiff) tiff = (uint16_t*)malloc(x*y*z*sizeof(uint16_t));
			readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 32){
			if(!tiff) tiff = (float*)malloc(x*y*z*sizeof(float));
			readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 64){
			if(!tiff) tiff = (double*)malloc(x*y*z*sizeof(double));
			readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else{
			return NULL;
		}
	}
	else if(z <= 1){
		if(bits == 8){
			if(!tiff) tiff = (uint8_t*)malloc(x*y*z*sizeof(uint8_t));
			readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize,flipXY);
			return (void*)tiff;
		}
		else if(bits == 16){
			if(!tiff) tiff = (uint16_t*)malloc(x*y*z*sizeof(uint16_t));
			readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 32){
			if(!tiff) tiff = (float*)malloc(x*y*z*sizeof(float));
			readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 64){
			if(!tiff) tiff = (double*)malloc(x*y*z*sizeof(double));
			readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else{
			return NULL;
		}
	}
	else{
		if(bits == 8){
			if(!tiff) tiff = (uint8_t*)malloc(x*y*z*sizeof(uint8_t));
			readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 16){
			if(!tiff) tiff = (uint16_t*)malloc(x*y*z*sizeof(uint16_t));
			readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 32){
			if(!tiff) tiff = (float*)malloc(x*y*z*sizeof(float));
			readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else if(bits == 64){
			if(!tiff) tiff = (double*)malloc(x*y*z*sizeof(double));
			readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
			return (void*)tiff;
		}
		else{
			return NULL;
		}
	}

	// Should never get here but return NULL if we do
	return NULL;
}

void* readTiffParallelWrapper(const char* fileName)
{
	return readTiffParallelWrapperHelper(fileName,NULL,1);
}

void* readTiffParallelWrapperNoXYFlip(const char* fileName, const std::vector<uint64_t> &zRange)
{
	return readTiffParallelWrapperHelper(fileName,NULL,0,zRange);
}

// tTiff doesn't matter as tiff is set in the function
void readTiffParallelWrapperSet(const char* fileName, void* tiff){
	void* tTiff = readTiffParallelWrapperHelper(fileName,tiff,0);
}
