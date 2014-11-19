//ReadWavFile.h
#ifndef _READWAVFILE_H_
#define _READWAVFILE_H_

#include <iostream>
#include <string>
#include "stdint.h"
#include "Structure.h"

typedef int8_t int8;
typedef uint8_t uint8;
typedef int16_t int16;
typedef uint16_t uint16;
typedef int32_t int32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t uint64;

unsigned int bingetu1 (FILE *f, bool *TF);

int binget2LE(FILE *f, bool *TF);

long binget3LE(FILE *f, bool *TF);

int binget4LE(FILE *f, bool *TF);

double bingetr4LE (FILE *f, bool *TF);

int checkWavFile(FILE *f, int *numOfChannels, int *encoding, double *sampleRate, 
                  long *startOfData, long *numOfSamples);

int readAudioToFloat(FILE  *f, int numberOfChannels, int encoding, double **buffer, long numberOfSamples);

Sound Sound_readFromSoundFile(std::string path);

int Sound_writeTxt(Sound sound, FILE *fp);

#endif
