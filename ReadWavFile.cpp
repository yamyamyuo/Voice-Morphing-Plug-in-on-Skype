/*ReadWavFile.cpp
 * Read Sound data from given file end with ".wav"
 *     Sound should have only one sound strack
 */ 
 
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "SoundCompute.h"
#include "ReadWavFile.h"

static int ulaw2linear [] = 
      { -32124, -31100, -30076, -29052, -28028, -27004, -25980, -24956,
        -23932, -22908, -21884, -20860, -19836, -18812, -17788, -16764,
        -15996, -15484, -14972, -14460, -13948, -13436, -12924, -12412,
        -11900, -11388, -10876, -10364,  -9852,  -9340,  -8828,  -8316,
         -7932,  -7676,  -7420,  -7164,  -6908,  -6652,  -6396,  -6140,
         -5884,  -5628,  -5372,  -5116,  -4860,  -4604,  -4348,  -4092,
         -3900,  -3772,  -3644,  -3516,  -3388,  -3260,  -3132,  -3004,
         -2876,  -2748,  -2620,  -2492,  -2364,  -2236,  -2108,  -1980,
         -1884,  -1820,  -1756,  -1692,  -1628,  -1564,  -1500,  -1436,
         -1372,  -1308,  -1244,  -1180,  -1116,  -1052,   -988,   -924,
          -876,   -844,   -812,   -780,   -748,   -716,   -684,   -652,
          -620,   -588,   -556,   -524,   -492,   -460,   -428,   -396,
          -372,   -356,   -340,   -324,   -308,   -292,   -276,   -260,
          -244,   -228,   -212,   -196,   -180,   -164,   -148,   -132,
          -120,   -112,   -104,    -96,    -88,    -80,    -72,    -64,
           -56,    -48,    -40,    -32,    -24,    -16,     -8,      0,
         32124,  31100,  30076,  29052,  28028,  27004,  25980,  24956,
         23932,  22908,  21884,  20860,  19836,  18812,  17788,  16764,
         15996,  15484,  14972,  14460,  13948,  13436,  12924,  12412,
         11900,  11388,  10876,  10364,   9852,   9340,   8828,   8316,
          7932,   7676,   7420,   7164,   6908,   6652,   6396,   6140,
          5884,   5628,   5372,   5116,   4860,   4604,   4348,   4092,
          3900,   3772,   3644,   3516,   3388,   3260,   3132,   3004,
          2876,   2748,   2620,   2492,   2364,   2236,   2108,   1980,
          1884,   1820,   1756,   1692,   1628,   1564,   1500,   1436,
          1372,   1308,   1244,   1180,   1116,   1052,    988,    924,
           876,    844,    812,    780,    748,    716,    684,    652,
           620,    588,    556,    524,    492,    460,    428,    396,
           372,    356,    340,    324,    308,    292,    276,    260,
           244,    228,    212,    196,    180,    164,    148,    132,
           120,    112,    104,     96,     88,     80,     72,     64,
            56,     48,     40,     32,     24,     16,      8,      0
       };

static short alaw2linear[] = 
{
   -5504,  -5248,  -6016,  -5760,  -4480,  -4224,  -4992,  -4736,
   -7552,  -7296,  -8064,  -7808,  -6528,  -6272,  -7040,  -6784,
   -2752,  -2624,  -3008,  -2880,  -2240,  -2112,  -2496,  -2368,
   -3776,  -3648,  -4032,  -3904,  -3264,  -3136,  -3520,  -3392,
  -22016, -20992, -24064, -23040, -17920, -16896, -19968, -18944,
  -30208, -29184, -32256, -31232, -26112, -25088, -28160, -27136,
  -11008, -10496, -12032, -11520,  -8960,  -8448,  -9984,  -9472,
  -15104, -14592, -16128, -15616, -13056, -12544, -14080, -13568,
    -344,   -328,   -376,   -360,   -280,   -264,   -312,   -296,
    -472,   -456,   -504,   -488,   -408,   -392,   -440,   -424,
     -88,    -72,   -120,   -104,    -24,     -8,    -56,    -40,
    -216,   -200,   -248,   -232,   -152,   -136,   -184,   -168,
   -1376,  -1312,  -1504,  -1440,  -1120,  -1056,  -1248,  -1184,
   -1888,  -1824,  -2016,  -1952,  -1632,  -1568,  -1760,  -1696,
    -688,   -656,   -752,   -720,   -560,   -528,   -624,   -592,
    -944,   -912,  -1008,   -976,   -816,   -784,   -880,   -848,
    5504,   5248,   6016,   5760,   4480,   4224,   4992,   4736,
    7552,   7296,   8064,   7808,   6528,   6272,   7040,   6784,
    2752,   2624,   3008,   2880,   2240,   2112,   2496,   2368,
    3776,   3648,   4032,   3904,   3264,   3136,   3520,   3392,
   22016,  20992,  24064,  23040,  17920,  16896,  19968,  18944,
   30208,  29184,  32256,  31232,  26112,  25088,  28160,  27136,
   11008,  10496,  12032,  11520,   8960,   8448,   9984,   9472,
   15104,  14592,  16128,  15616,  13056,  12544,  14080,  13568,
     344,    328,    376,    360,    280,    264,    312,    296,
     472,    456,    504,    488,    408,    392,    440,    424,
      88,     72,    120,    104,     24,      8,     56,     40,
     216,    200,    248,    232,    152,    136,    184,    168,
    1376,   1312,   1504,   1440,   1120,   1056,   1248,   1184,
    1888,   1824,   2016,   1952,   1632,   1568,   1760,   1696,
     688,    656,    752,    720,    560,    528,    624,    592,
     944,    912,   1008,    976,    816,    784,    880,    848
};

unsigned int bingetu1 (FILE *f, bool *TF) {
	 int externalValue = getc (f);   // either -1 (EOF) or in the range 0..255
	  if (externalValue < 0) {
          std::cout<<"A byte."<<std::endl;
	      *TF = false;
		  return -2;
	   }
	   *TF = true;
	  return (unsigned int) externalValue;   
}	
   
int binget2LE(FILE *f, bool *TF){
   signed short s;
   if(fread(& s, sizeof(signed short), 1, f) != 1){
	   std::cout<<"A signed short integer error"<<std::endl;
	   *TF = false;
	   return -2;
	}
	*TF = true;
    return (int) s;	 
}

long binget3LE(FILE *f, bool *TF) {
	  unsigned char bytes [3];
	  if (fread (bytes, sizeof (unsigned char), 3, f) != 3){
	        std::cout<<"three bytes." << std::endl;
			*TF = false;
			return -2;
		}
	  uint32_t externalValue = ((uint32_t) bytes [2] << 16) | ((uint32_t) bytes [1] << 8) | (uint32_t) bytes [0];
	  if ((bytes [2] & 128) != 0)  			   // is the 24-bit sign bit on?
		  externalValue |= 0xFF000000;  	   // extend negative sign to 32 bits
		  *TF = true;
	  return (long) (int32_t) externalValue;   // first convert signedness, then perhaps extend sign to 64 bits!     
}

int binget4LE(FILE *f, bool *TF){
   long l;
   if (fread (& l, sizeof (long), 1, f) != 1){
       std::cout<<"A signed long integer error"<<std::endl;
	   *TF = false;
	   return -2;
	}
	*TF = true;
	return l;
}

double bingetr4LE (FILE *f, bool *TF) {
      float x;
	  if (fread (& x, sizeof (float), 1, f) != 1){
	       std::cout<<"A single-precision floating-point number."<<std::endl;
		   *TF = false;
		   return -2;
       }
      *TF = true;	   
	   return x;   
}

#ifndef WAVE_FORMAT_PCM
	#define WAVE_FORMAT_PCM  0x0001
#endif
#define WAVE_FORMAT_IEEE_FLOAT  0x0003
#define WAVE_FORMAT_ALAW  0x0006
#define WAVE_FORMAT_MULAW  0x0007
#define WAVE_FORMAT_DVI_ADPCM  0x0011

/*Audio encoding*/
#define Melder_LINEAR_8_UNSIGNED  2
#define Melder_LINEAR_16_LITTLE_ENDIAN  4
#define Melder_LINEAR_24_LITTLE_ENDIAN  6
#define Melder_LINEAR_32_LITTLE_ENDIAN  8
#define Melder_MULAW  9
#define Melder_ALAW  10
#define Melder_IEEE_FLOAT_32_LITTLE_ENDIAN  14


int checkWavFile(FILE *f, int *numOfChannels, int *encoding, double *sampleRate, long *startOfData, long *numOfSamples)
/**
    return value description:
	 0:   normal return
	-1:   format error
	-2:   data error
**/
{
    char data[8], ChunkID[4];
	bool formatChunkPresent = false, dataChunkPresent = false;
	int numberOfBitsPerSamplePoint = -1;
	long dataChunkSize = -1;
	bool TF = false;
	memset(data,0, sizeof(char)*8);
	memset(ChunkID, 0,sizeof(char)*4);
	/*check the header part of the audio file*/
	if(fread(data, 1, 4, f) < 4){
	    std::cout<<"File too small: no RIFF statement.\n";
	    return -1;
	}
	if(strncmp(data,"RIFF", 4)){
	    std::cout<<"Not a WAV file (RIFF statement expected).\n";
	    return -1;
	}
	if(fread(data, 1, 4, f) < 4){
	    std::cout<<"File too small: no size of RIFF chunk.\n";
	    return -1;
	}
	if(fread(data, 1, 4, f) < 4){
	    std::cout<<"File too small: no file type info (expected WAVE statement).\n";
	    return -1;
	}
	if(strncmp(data,"WAVE", 4)){
	    std::cout<<"Not a WAVE file (wrong file type info).\n";
	    return -1;
	}
	
	/*search for Format Chunck and Data Chunck*/
	while(fread(ChunkID, 1, 4, f) == 4)
	{
	    long chunkSize = binget4LE(f, &TF);    //long chunkSize = bingeti4LE (f);
		if(!TF) return -1;
		
		if(! strncmp(ChunkID, "fmt ", 4)){
		 // found a Format Chunk.
		 int winEncoding = binget2LE(f, &TF);   //winEncoding = bingeti2LE (f);
		 if(!TF) return -1;
		 
		 formatChunkPresent = true;
		 *numOfChannels = binget2LE(f, &TF);  //*numberOfChannels = bingeti2LE (f);
		 if(!TF) return -1;
		 
		 if(*numOfChannels < 1){
		     std::cout<<"Too few sound channels"<<std::endl;
			 std::cout<<"ReadWav.cpp Line 188."<<std::endl;
			 return -1;
		 }
		 *sampleRate = (double) binget4LE(f, & TF);
		 if(!TF) return -1;
		 if(*sampleRate <= 0.0){
		     std::cout<<"Wrong sampling frequency"<<std::endl;
			 std::cout<<"ReadWav.cpp Line 188."<<std::endl;
			 return -1;
		 }
		 (void) binget4LE (f, &TF);   // avgBytesPerSec
	     (void) binget2LE (f, &TF);   // blockAlign
		 
		 numberOfBitsPerSamplePoint = binget2LE(f, &TF);
		 if(!TF) return -1;	
		 
		 if(numberOfBitsPerSamplePoint == 0)
		    numberOfBitsPerSamplePoint = 16;     //default value
		 else if (numberOfBitsPerSamplePoint < 4){
		     std::cout<<"Too few bits per sample"<<std::endl;
			 std::cout<<"ReadWav.cpp Line 188.263"<<std::endl;
			 return -1;
		 }
		 else if (numberOfBitsPerSamplePoint > 32){
		     std::cout<<"Too few bits per sample"<<std::endl;
			 std::cout<<"ReadWav.cpp Line 188.268"<<std::endl;
			 return -1;		    
		 }
		 
		 switch(winEncoding){    // get encoding
		       case WAVE_FORMAT_PCM :
			       *encoding =  numberOfBitsPerSamplePoint > 24 ? Melder_LINEAR_32_LITTLE_ENDIAN :
						        numberOfBitsPerSamplePoint > 16 ? Melder_LINEAR_24_LITTLE_ENDIAN :
						        numberOfBitsPerSamplePoint > 8 ? Melder_LINEAR_16_LITTLE_ENDIAN :
						        Melder_LINEAR_8_UNSIGNED;
				   break;
			   case WAVE_FORMAT_IEEE_FLOAT :
			       *encoding = Melder_IEEE_FLOAT_32_LITTLE_ENDIAN;
				   break;
			   case WAVE_FORMAT_ALAW :
			       *encoding = Melder_ALAW;
				   break;
		       case WAVE_FORMAT_MULAW :
			        *encoding = Melder_MULAW ;
					break;
			   case WAVE_FORMAT_DVI_ADPCM :
					std::cout<<"Cannot read lossy compressed audio files (this one is DVI ADPCM).\n"
						"Please use uncompressed audio files. If you must open this file,\n"
						"please use an audio converter program to convert it first to normal (PCM) WAV format\n"
						"(Praat may have difficulty analysing the poor recording, though)."<<std::endl;
					std::cout<<"ReadWav.cpp Line 188.293"<<std::endl;
					return -2;
			   default:
			        std::cout<<"Unsupported Windows audio encoding %d."<< winEncoding <<std::endl;
                    std::cout<<"ReadWav.cpp Line 188.294"<<std::endl;
					return -2;							   
		 }    
		  
		 if(chunkSize & 1) chunkSize ++;
		 for(long i = 17; i <= chunkSize; i ++)
			if(fread(data, 1, 1, f) < 1){
			    std::cout<<"File too small: expected" << chunkSize << "bytes in fmt chunk, but found " << i << "."<<std::endl;   
			    std::cout<<"ReadWav.cpp Line 188.303"<<std::endl;
				return -1;
			} 
	  }else if(!strncmp(ChunkID, "data", 4))
	   {
	      // Found a Data Chunk.
		  dataChunkPresent = true;
		  dataChunkSize = chunkSize;
		  *startOfData = ftell(f);
		  if(chunkSize & 1) chunkSize ++;
		  if(chunkSize < 0) {// incorrect data chunk (sometimes -44); assume that the data run till the end of the file
		      fseek(f, 0, SEEK_END);
			  long endOfData = ftell(f);
			  dataChunkSize = chunkSize = endOfData - *startOfData;
			  fseek (f, *startOfData, SEEK_SET);
		  }
		  for(long i = 1;i <= chunkSize; ++ i){
			   if (fread (data, 1, 1, f) < 1){
			      	std::cout<<"File too small: expected " << chunkSize << " bytes in fmt chunck, but found " << i << "."<<std::endl;   
			        std::cout<<"ReadWav.cpp Line 188. 321"<<std::endl;
					return -1;
			   }
		  }
	  }else
		  break ;
   //    } else 
	  // {       //ignore the chunks
		 //  if (chunkSize & 1) chunkSize ++;
		 //  for (long i = 1; i <= chunkSize; i ++)
			//   if (fread (data, 1, 1, f) < 1){
			//      	std::cout<<"File too small: expected " << chunkSize << " bytes in fmt chunck, but found " << i << "."<<std::endl;   
			//        std::cout<<"ReadWav.cpp Line 188. 322"<<std::endl;
			//		return -1;
			//  }
		 //}  
	 }
		
    if(! formatChunkPresent) {
	   std::cout<< "Found no Format Chunck." <<std::endl;
	   std::cout<<"ReadWav.cpp Line 188."<<std::endl;
	   return -2;
	}
	if(! dataChunkPresent) {
	   std::cout<< "Found no Data Chunck." <<std::endl;
	   std::cout<<"ReadWav.cpp Line 188."<<std::endl;
	   return -2;
	}
	if(numberOfBitsPerSamplePoint != -1 && dataChunkSize != -1)
	    *numOfSamples = dataChunkSize / *numOfChannels / (numberOfBitsPerSamplePoint / 8);//(numberOfBitsPerSamplePoint + 7)
    else{ 
		 std::cout<<"ReadWav.cpp Line 188."<<std::endl;
	     return -2;
	 }
	return 0;
}

int readAudioToFloat(FILE  *f, int numberOfChannels, int encoding, double **buffer, long numberOfSamples)
/**
    return value description:
	 0:   normal return
	-1:   error
**/
{
	bool TF = false;
     switch(encoding){
	     case Melder_LINEAR_8_UNSIGNED:
			    for (long isamp = 1; isamp <= numberOfSamples; isamp ++)
					for (long ichan = 1; ichan <= numberOfChannels; ichan ++) {
						buffer[ichan][isamp] = bingetu1(f, &TF) * (1.0 / 128) - 1.0;
					    if(!TF) return -1;
					}
			    break;
		 case Melder_LINEAR_16_LITTLE_ENDIAN: {
				const int numberOfBytesPerSamplePerChannel = 2;
				if (numberOfChannels > (int) sizeof (double) / numberOfBytesPerSamplePerChannel) {
					for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
						for (long ichan = 1; ichan <= numberOfChannels; ichan ++) {
							buffer[ichan][isamp] = binget2LE (f, &TF) * (1.0 / 32768);
					        if(!TF) return -1;
							}
					}
				} else {      // optimize
					long numberOfBytes = numberOfChannels * numberOfSamples * numberOfBytesPerSamplePerChannel;
					unsigned char *bytes = (unsigned char *)&buffer[numberOfChannels][numberOfSamples] + sizeof(double) - numberOfBytes;
					if (fread (bytes, 1, numberOfBytes, f) < numberOfBytes)   /*** throw MelderError ();  ***/  // read 16-bit data into last quarter of buffer
					{
						std::cout<<"Error, read 16-bit data into last quarter of buffer"<<std::endl;
						std::cout<<"ReadWavFile.cpp: Line: 332."<<std::endl;
						return -1;
					}
					if (numberOfChannels == 1) {
						for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {  /**isamp = 1; isamp <= numberOfSamples;**/
							unsigned char byte1 = *bytes ++, byte2 = *bytes ++;
							long value = (int16_t)(((uint16_t) byte2 << 8) | (uint16_t) byte1);
							buffer[1][isamp] = value * (1.0 / 32768);
						}
					} else {
						for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {  /**isamp = 1; isamp <= numberOfSamples;**/
							for (long ichan = 1; ichan <= numberOfChannels; ichan ++) {  /**ichan = 1; ichan <= numberOfChannels;**/
								unsigned char byte1 = * bytes ++, byte2 = * bytes ++;
								long value = (int16_t) (((uint16_t) byte2 << 8) | (uint16_t) byte1);
								buffer[ichan][isamp] = value * (1.0 / 32768);
							}
						}
					}
				}
			}
			break; 
			
		 case Melder_LINEAR_24_LITTLE_ENDIAN: {
				const int numberOfBytesPerSamplePerChannel = 3;
				if (numberOfChannels > (int) sizeof (double) / numberOfBytesPerSamplePerChannel) {
					for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
						for (long ichan = 1; ichan <= numberOfChannels; ichan ++) {
							buffer[ichan][isamp] = binget3LE(f, &TF) * (1.0 / 8388608);
							if(!TF) return -1;
						}
					}
				} else { // optimize
					long numberOfBytes = numberOfChannels * numberOfSamples * numberOfBytesPerSamplePerChannel;
					unsigned char *bytes = (unsigned char *) & buffer [numberOfChannels - 1] [numberOfSamples] + sizeof (double) - numberOfBytes;
					if (fread (bytes, 1, numberOfBytes, f) < numberOfBytes) /*** throw MelderError (); ***/  // read 24-bit data into last three-eights of buffer
					{
						std::cout<<"Error, read 24-bit data into last quarter of buffer"<<std::endl;
						std::cout<<"ReadWavFile.cpp: Line: 332."<<std::endl;
						return -1;
					}
					if (numberOfChannels == 1) {
						for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
							unsigned char byte1 = * bytes ++, byte2 = * bytes ++, byte3 = * bytes ++;
							uint32_t unsignedValue = ((uint32_t) byte3 << 16) | ((uint32_t) byte2 << 8) | (uint32_t) byte1;
							if ((byte3 & 128) != 0) unsignedValue |= 0xFF000000;
							buffer[1][isamp] = (int32_t) unsignedValue * (1.0 / 8388608);
						}
					} else {
						for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
							for (long ichan = 1; ichan <= numberOfChannels; ichan ++) {
								unsigned char byte1 = * bytes ++, byte2 = * bytes ++, byte3 = * bytes ++;
								uint32_t unsignedValue = ((uint32_t) byte3 << 16) | ((uint32_t) byte2 << 8) | (uint32_t) byte1;
								if ((byte3 & 128) != 0) unsignedValue |= 0xFF000000;
								buffer[ichan][isamp] = (int32_t) unsignedValue * (1.0 / 8388608);
							}
						}
					}
				}
			}
			break;
		 case Melder_LINEAR_32_LITTLE_ENDIAN:
			    for (long isamp = 1; isamp <= numberOfSamples; isamp ++){ 
				     for (long ichan = 1; ichan <= numberOfChannels; ichan ++) 
							buffer[ichan][isamp] = binget4LE(f, &TF) * (1.0 / 32768 / 65536);
					        if(!TF) return -1;				
				}
			break;	
		 case Melder_MULAW:
				for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
				 	for (long ichan = 1; ichan <= numberOfChannels; ichan ++) 
						buffer[ichan][isamp] = ulaw2linear[bingetu1 (f, &TF)] * (1.0 / 32768);
					    if(!TF) return -1;	
				}
		    break;
		 case Melder_ALAW:
				for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
					for (long ichan = 1; ichan <= numberOfChannels; ichan ++) 
						buffer[ichan][isamp] = alaw2linear[bingetu1 (f, &TF)] * (1.0 / 32768);
					    if(!TF) return -1;	
				 }
		    break;
		 case Melder_IEEE_FLOAT_32_LITTLE_ENDIAN:
				for (long isamp = 1; isamp <= numberOfSamples; isamp ++) {
					for (long ichan = 1; ichan <= numberOfChannels; ichan ++)
						buffer[ichan][isamp] = bingetr4LE (f, &TF);
					    if(!TF) return -1;						
				}
			break;
	 }
	return 0;
}

Sound Sound_readFromSoundFile(std::string path)
  /* Read Data from Audio file */
{  
	if(path.length() == 0 || path == " ")
	    return NULL;

	FILE *fp = NULL;
    int numOfChannels = 0, encoding = 0;
	double sampleRate = 0;
	long startOfData = 0, numOfSamples = 0;
	int checkVal = 0, readval = 0;
	Sound me =- NULL;
	if((fp = fopen(path.c_str(), "rb")) == NULL)  return NULL;
	
	checkVal = checkWavFile(fp, &numOfChannels, &encoding, &sampleRate, &startOfData, &numOfSamples);

    if(checkVal != 0)  
		return NULL;
    me = Sound_createSimple( numOfSamples / sampleRate, sampleRate, numOfChannels);

    if(! me) 
		return NULL;
	fseek(fp, startOfData, SEEK_SET);

	readval = readAudioToFloat (fp, numOfChannels, encoding, my z, numOfSamples);

	if(readval == -1) 
		return NULL;

	fclose(fp);
	return me;
}

int Sound_writeTxt(Sound sound, FILE *fp){
    if(sound == NULL || fp == NULL)
		return 0;
	
    fprintf(fp, "File type = \"ooTextFile\"\n");
    fprintf(fp, "Object class = \"Sound 2\"\n\n");
	fprintf(fp, "xmin = %f\n", sound->xmin);
	fprintf(fp, "xmax = %18.17f\n", sound->xmax);
	fprintf(fp, "nx = %d \n", sound->nx);
	fprintf(fp, "dx = %17.16f \n", sound->dx);
	fprintf(fp, "x1 = %17.16f \n", sound->x1);
	fprintf(fp, "ymin = 1 \n");
	fprintf(fp, "ymax = 1 \n");
	fprintf(fp, "ny = 1 \n");
	fprintf(fp, "dx = 1 \n");
	fprintf(fp, "y1 = 1 \n");
	fprintf(fp, "z [] []: \n");
	fprintf(fp, "    z [1]:\n");
	
	for(long i = 1; i <= sound->nx; ++ i)
	     fprintf(fp, "%17.16f \n", sound->z[1][i]);//       z [1] [%d]: =      i,
		 
	fclose(fp);
	return 1;
}

