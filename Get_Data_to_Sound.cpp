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
using namespace std;

#define FREQUENCY 22050.0
#define DURATION 0.1
#define DURATIONSTEP 0.04
#define FRAMESIZE FREQUENCY*0.1
#define FRAMESTEP FREQUENCY*0.04
#define MAX_T  0.02000000001   /* Maximum interval between two voice pulses (otherwise voiceless). */

int Pitch_to_Framestep(Pitch srcpitch, Sound me, Pitch thee, int framestep,int pretag)
{
	int i = srcpitch->nx, tag = 0;
	if(srcpitch->frame[srcpitch->nx].candidate[1].frequency==0)//unvoice
	{
		while(srcpitch->frame[srcpitch->nx].candidate[1].frequency==0)
		{
			tag++;
			i--;
		}
	}
	else if(srcpitch->frame[srcpitch->nx].candidate[1].frequency>0)//voice
	{
		while(srcpitch->frame[srcpitch->nx].candidate[1].frequency>0)
		{
			tag++;
			i--;
		}
	}
	if(srcpitch->frame[srcpitch->nx].candidate[1].frequency>0&&tag<3)//·ÅÈëÏÂÒ»²½´¦Àí
	{

	}
	else//Ö±½ÓÕâÒ»²½¾Í´¦ÀíÁË
	{

	}

	return tag;
}

void Get_Data(double *srcdata, double *buff, int *Flag, int datalength)
{
	int i,j;
   
	for(i = *Flag,j = 1; i < *Flag + FRAMESIZE; i++,j++)
	{
		buff[j] = srcdata[i%datalength]; 
	}
	int k = (int)FRAMESTEP;
	*Flag += k; 
}

void copydata(double *srcdata, double *buff, int start, int end)
{
	int i,j;
	for(i = start, j = 1; i <= end; i++, j++)
	{
		buff[j] = srcdata[i];
	}
}



void copypitch(const Pitch me, Pitch him, int start, int end)
{
	int i,j,k;
	for(i = start,j=1; i <= end; i++, j++)
	{
   	    him->frame[j].nCandidates = me->frame[i].nCandidates;
		him->frame[j].intensity = me->frame[i].intensity;

		for(k = 1 ; k <= me->frame[i].nCandidates; k++)
		{
			him->frame[j].candidate[k].frequency = me->frame[i].candidate[k].frequency;
			him->frame[j].candidate[k].strength = me->frame[i].candidate[k].strength;
		}

	}
}


void Convert_array_to_double(int *int_data, double *double_data, int datalength)
{
	int i;
	for(i = 0; i < datalength; i++)
	{
		double_data[i] = int_data[i]/FREQUENCY;//ÐèÒªÐÞ¸Ä
	}
}


Sound Data_to_Sound(double *data,double samplingFrequency,int numOfChannels, double duration)
{
	Sound sound = NULL;
	sound = Sound_createSimple(duration, samplingFrequency, numOfChannels);
	
	if(! sound) 
		return NULL;
	return sound;
}

long pitch_xToIndex(Pitch me, double x)
{
	return (long)floor((x-my x1)/my dx);
}

PointProcess Sound_and_Pitch_to_Pulses (Sound me, Pitch him, double formantRatio,double new_pitch, 
	                                    double pitchRangeFactor, double durationFactor){
     double samplingFrequency_old = 1 / my dx;
	 
	 if(my ny > 1){
		std::cout<<"Change Gender works only on mono sounds."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 if (my xmin != his xmin || my xmax != his xmax){
		std::cout<<"The Pitch and the Sound object must have the same starting times and finishing times."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 if(new_pitch < 0) {
		std::cout<<"The new pitch median must not be negative."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 Sound sound = (Sound)malloc(sizeof(structSound));
	 memcpy(sound, me, sizeof(structSound));

	 //Î´Ïû³ýÖ±Á÷·ÖÁ¿

	 if (formantRatio != 1) {
		// Shift all frequencies (inclusive pitch!)
		Sound_overrideSamplingFrequency (sound, samplingFrequency_old * formantRatio);
	}

	 Pitch pitch = Pitch_scaleTime_old(him, 1 / formantRatio); //¸ø³öÖØ²ÉÑùµÄÊý¾ÝÇø¼ä

	 PointProcess pulses = Sound_Pitch_to_PointProcess_cc (sound, pitch);


	 return pulses;
}

Sound Sound_and_Pitch_and_Pulses_changeGender_old (Sound me, Pitch him, PointProcess pulses, double formantRatio,double new_pitch, 
	                                    double pitchRangeFactor, double durationFactor){
     double samplingFrequency_old = 1 / my dx;
	 
	 if(my ny > 1){
		std::cout<<"Change Gender works only on mono sounds."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 if (my xmin != his xmin || my xmax != his xmax){
		std::cout<<"The Pitch and the Sound object must have the same starting times and finishing times."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 if(new_pitch < 0) {
		std::cout<<"The new pitch median must not be negative."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	 }

	 //Î´Ïû³ýÖ±Á÷·ÖÁ¿

	 Sound sound = (Sound)malloc(sizeof(structSound));
	 memcpy(sound, me, sizeof(structSound));
	 if (formantRatio != 1) {
		// Shift all frequencies (inclusive pitch!)
		Sound_overrideSamplingFrequency (sound, samplingFrequency_old * formantRatio);
	}

	 Pitch pitch = Pitch_scaleTime_old(him, 1 / formantRatio); //¸ø³öÖØ²ÉÑùµÄÊý¾ÝÇø¼ä

	 //PointProcess pulses = Sound_Pitch_to_PointProcess_cc (sound, pitch);

     PitchTier pitchTier = Pitch_to_PitchTier (pitch);

	 double median = Pitch_getQuantile (pitch, 0, 0, 0.5, kPitch_unit_HERTZ);

	 if (median != 0 && median != NUMundefined) {
		// Incorporate pitch shift from overriding the sampling frequency
		if (new_pitch == 0) 
			new_pitch = median / formantRatio;
		double factor = new_pitch / median;//Í³Ò»³ËÉÏÒ»¸öÊýÖµ
		PitchTier_multiplyFrequencies (pitchTier, sound->xmin, sound->xmax, factor);//ÆµÂÊÍ³Ò»³ËÉÏÁËfactor
		PitchTier_modifyRange_old (pitchTier, sound->xmin, sound->xmax, pitchRangeFactor, new_pitch);
		//factor=pitchRangeFactorºÍfmid=new_pitchÍ¨¹ýËã·¨f = fmid + (f - fmid) * factor;
	}else{
        std::cout<<"There were no voiced segments found."<<std::endl;
		std::cout<<"SoundCompute.cpp Line: 502"<<std::endl;
		return NULL;
	}

	 DurationTier duration = DurationTier_create (my xmin, my xmax);
	 RealTier_addPoint(duration, (my xmin + my xmax) / 2, formantRatio * durationFactor);//È¡ÖÐµã,Ïû³ýÀÛ»ýÎó²î£¬Ê¹É¾³ý²åÈë²»ÔÚÍ·²¿ÓëÎ²²¿

	 Sound thee = Sound_Point_Pitch_Duration_to_Sound (sound, pulses, pitchTier, duration,
		                   1.25 / Pitch_getMinimum(pitch, 0.0, 0.0, kPitch_unit_HERTZ, false),formantRatio);

	 // Resample to the original sampling frequencyÖØ²ÉÑù
	if (formantRatio != 1) 
		Sound_resample(thee, samplingFrequency_old, 10);

	return thee;
}

PointProcess Sound_and_pitch_to_pulses_changespeaker(Sound me, Pitch him, 
						   double formantMultiplier , double pitchMultiplier,
						   double pitchRangeMultiplier, double durationMultiplier)
{
	double samplingFrequency_old = 1 / my dx;
	/*if(my ny > 1){
		std::cout<<"Change Gender works only on mono sounds."<<std::endl;
		return NULL;
	}*/
	
	/*if(pitchMultiplier < 0){
		std::cout<<"The new pitch median must not be negative."<<std::endl;
		return NULL;	   
	}*/
	if(my xmin != his xmin || my xmax != his xmax){
	    std::cout<<"The Pitch and the Sound object must have the same start and end times."<<std::endl;
		std::cout<<"SoundCompute.cpp 418"<<std::endl;
		return NULL;
	}
	
	Sound sound = (Sound) malloc(sizeof(structSound));
	
	memcpy(sound, me, sizeof(structSound));
	Sound_substructMean(sound);

	if(formantMultiplier != 1)   //shift all frequencies (inclusive pitch!)
		Sound_overrideSamplingFrequency(sound, samplingFrequency_old * formantMultiplier);
   
	Pitch pitch = (Pitch) malloc(sizeof(structPitch));
	memcpy(pitch, him, sizeof(structPitch));
	Pitch_scaleDuration(pitch, 1 / formantMultiplier);
	Pitch_scalePitch(pitch, formantMultiplier);
 
    PointProcess pulses = Sound_Pitch_to_PointProcess_cc(sound, pitch);
	
	return pulses;
}


Sound Sound_and_pitch_and_pulses_changespeaker(Sound me, Pitch him, PointProcess pulses,
						   double formantMultiplier , double pitchMultiplier,
						   double pitchRangeMultiplier, double durationMultiplier)
{
	double samplingFrequency_old = 1 / my dx;
	/*if(my ny > 1){
		std::cout<<"Change Gender works only on mono sounds."<<std::endl;
		return NULL;
	}*/
	
	/*if(pitchMultiplier < 0){
		std::cout<<"The new pitch median must not be negative."<<std::endl;
		return NULL;	   
	}*/
	if(my xmin != his xmin || my xmax != his xmax){
	    std::cout<<"The Pitch and the Sound object must have the same start and end times."<<std::endl;
		std::cout<<"SoundCompute.cpp 418"<<std::endl;
		return NULL;
	}
	
	Sound sound = (Sound) malloc(sizeof(structSound));
	
	memcpy(sound, me, sizeof(structSound));
	Sound_substructMean(sound);

	if(formantMultiplier != 1)   //shift all frequencies (inclusive pitch!)
		Sound_overrideSamplingFrequency(sound, samplingFrequency_old * formantMultiplier);
   
	Pitch pitch = (Pitch) malloc(sizeof(structPitch));
	memcpy(pitch, him, sizeof(structPitch));
	Pitch_scaleDuration(pitch, 1 / formantMultiplier);
	Pitch_scalePitch(pitch, formantMultiplier);//ÆµÂÊ³ËÉÏformatMultiplier
 
    //PointProcess pulses = Sound_Pitch_to_PointProcess_cc(sound, pitch);
	PitchTier pitchTier = Pitch_to_PitchTier(pitch);

	double median = Pitch_getQuantile(pitch, 0, 0, 0.5, kPitch_unit_HERTZ);

	if(median != 0 && median != NUMundefined){
		PitchTier_multiplyFrequencies(pitchTier, sound->xmin, sound->xmax, pitchMultiplier/formantMultiplier);
		PitchTier_modifyExcursionRange(pitchTier, sound->xmin, sound->xmax, pitchRangeMultiplier, median);
	}else if (pitchMultiplier != 1){
		std::cout<<"Warning: Pitch has not been changed because the sound was entirely voiceless."<<std::endl;
	    std::cout<<"SoundCompute.cpp: Line 418"<<std::endl;

		return NULL;

	}

	DurationTier duration = DurationTier_create(my xmin, my xmax);
	RealTier_addPoint(duration, (my xmin + my xmax) / 2, formantMultiplier * durationMultiplier);

	Sound thee = Sound_Point_Pitch_Duration_to_Sound (sound, pulses, pitchTier, duration, MAX_T, formantMultiplier);

	if(formantMultiplier != 1)
	     thee = Sound_resample(thee, samplingFrequency_old, 10);
	return thee;
}

double find_proper_endtime(PointProcess pulses, double time)
{
	int i;
	for(i = pulses->nt; i > 0; i--)
	{
		if(pulses->t[i]<time)
			return pulses->t[i];
	}
	return 0;
}

int find_pulses_vextor_index(PointProcess pulses, double time, int low, int high)// low = 1, high = pulses->nt 
{

	int mid = (low + high)/2;
	if(pulses->t[mid] == time)//Ëã·¨ÓÐ´ýÑéÖ¤
	    return mid;
	else if(pulses->t[high] <time)
		return high;
	else if(pulses->t[mid] > time)
	    return find_pulses_vextor_index(pulses,time,low, mid - 1);
	else return find_pulses_vextor_index(pulses, time, mid + 1, high);
}

Pitch Pitch_create2(Pitch pitch0, int start, double tmin, double tmax, long nt, double dt, double t1, double ceiling, int maxCandidates){
   Pitch me = (Pitch)malloc(sizeof(structPitch));
   if(!me) return NULL;
   
   my xmin = tmin;
   my xmax = tmax;
   my nx = nt;
   my dx = dt;
   my x1 = t1;
   my ceiling = ceiling;
   my maxnCandidates = maxCandidates;
   my frame = (Pitch_Frame)malloc(sizeof(structPitch_Frame) * (nt + 1));
   
   /* Put one candidate in every frame (unvoiced, silent). */
   for(long i = 1,j = start; i <= nt; ++ i, j++)
	   Pitch_Frame_init(& my frame[i], pitch0->frame[j].nCandidates);
	   
	return me;
}

void copycheckdata(double *from, double *to, int start, int copylength)
{
	int i;
	for(i = 1; i <= copylength; i++)
	{
		to[start+i-1] = from[i];
	}

}
