#ifndef _GET_DATA_TO_SOUND_h_
#define _GET_DATA_TO_SOUND_h_

#include "Structure.h"
#include <iostream>
#include <string>
#include "Structure.h"
#include "Pitch.h"
#include "RealTier.h"
#include "SoundCompute.h"
#include "ReadWavFile.h"
#include "Sound_to_Pitch.h"
#include "Pitch_to_PointProcess.h"
#include <fstream>
#include <iostream>
#include "NUM.h"


void Get_Data(double *srcdata, double *buff, int *Flag, int datalength);
void Convert_array_to_double(int *int_data, double *double_data, int datalength);
Sound Data_to_Sound(double *data,double samplingFrequency,int numOfChannels, double duration);
void copydata(double *srcdata, double *buff, int start, int end);
void copypitch(const Pitch me, Pitch him, int start, int end);
long pitch_xToIndex(Pitch me, double x);
PointProcess Sound_and_Pitch_to_Pulses (Sound me, Pitch him, double formantRatio,double new_pitch, 
	                                    double pitchRangeFactor, double durationFactor);
Sound Sound_and_Pitch_and_Pulses_changeGender_old (Sound me, Pitch him, PointProcess pulses, double formantRatio,double new_pitch, 
	                                    double pitchRangeFactor, double durationFactor);
PointProcess Sound_and_pitch_to_pulses_changespeaker(Sound me, Pitch him, 
						   double formantMultiplier , double pitchMultiplier,
						   double pitchRangeMultiplier, double durationMultiplier);
Sound Sound_and_pitch_and_pulses_changespeaker(Sound me, Pitch him, PointProcess pulses,
						   double formantMultiplier , double pitchMultiplier,
						   double pitchRangeMultiplier, double durationMultiplier);

double find_proper_endtime(PointProcess pulses, double time);
int find_pulses_vextor_index(PointProcess pulses, double time, int low, int high);
Pitch Pitch_create2(Pitch pitch0, int start, double tmin, double tmax, long nt, double dt, double t1, double ceiling, int maxCandidates);
void copycheckdata(double *from, double *to, int start, int copylength);
#endif
