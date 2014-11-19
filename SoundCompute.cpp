/*SoundCompute.cpp
 * implementation for the function of Pitch_Frame¡¢ Pitch¡¢ NUMfft defined in structure.h file 
 */
 
#include "math.h"
#include <iostream>
#include <string.h>
#include "stdlib.h"
#include "Structure.h"
#include "Sound_to_Pitch.h"
#include "NUM2.h"
#include "SystemData.h"
#include "SoundCompute.h"
#include "NUM.h"
#include "Pitch.h"
#include "RealTier.h"
#include "Pitch_to_PointProcess.h"
#include "Get_Data_to_Sound.h"

#define MAX_T  0.02000000001   /* Maximum interval between two voice pulses (otherwise voiceless). */

MorphingFactor Factor_create(double fm, double pm, double prm, double dm){
	MorphingFactor factor = (MorphingFactor)malloc(sizeof(structMorphingFactor));
	memset(factor, 0, sizeof(structMorphingFactor));

	factor->durationMultiplier = dm;
	factor->formantMultiplier = fm;
	factor->pitchMultiplier = pm;
	factor->pitchRangeMultiplier = prm;
	return factor;
}

LocationTuple Location_create(double xpos, double ypos){
    LocationTuple tuple = (LocationTuple)malloc(sizeof(structLocationTuple));
	memset(tuple, 0, sizeof(structLocationTuple));

	tuple->x = xpos;
	tuple->y = ypos;

	return tuple;
}


Sound Sound_create(double xmin, double xmax, long nx, double dx, double x1, long numOfChannels)
{ 
    if(xmin >= xmax || nx < 1 || dx < 0 || numOfChannels < 1)
	     return NULL;
		 
    Sound me = (Sound)malloc(sizeof(structSound));
	my xmin = xmin;
	my xmax = xmax;
	my nx = nx;
	my dx = dx;
	my x1 = x1;
	my ny = numOfChannels;
	my z = (double **)malloc(sizeof(double *) * (numOfChannels + 1));
	if(!(my z))
	    return NULL;
    for(long j = 0; j <= numOfChannels; ++ j){
	   my z[j] = (double *)malloc(sizeof(double) * (nx + 1));
	   memset(my z[j], 0, sizeof(double) * (nx + 1));
	  /* for(int i = 0; i <= nx; ++ i)
	       my z[j][i] = 0.0;*/
	}
	
	 return me;
}

Sound Sound_createSimple(double duration, double samplingFrequency, long numOfChannels){
    return Sound_create(0.0, duration, floor(duration * samplingFrequency + 0.5), 
					     1 / samplingFrequency, 0.5 / samplingFrequency, numOfChannels);
}

Sound Sound_upsample (Sound me)
{
   long nfft = 1;
   while(nfft < my nx + 2000)  nfft *= 2;
   Sound thee = Sound_create(my xmin, my xmax, my nx, my dx, my x1, my ny);
   double *data = (double *)malloc(sizeof(double) * (2 *nfft + 1));

   memset(data, 0.0, sizeof(double)*(2 *nfft + 1));
   for(long channel = 1; channel <= my ny; channel ++) {
	   memcpy(data + 1000, my z[channel], my nx * sizeof(double));
	   NUMrealft (data, nfft, 1);
	   long imin = (long) (nfft * 0.95);
	   for(long i = imin + 1; i <= nfft; i ++)
		   data [i] *= ((double) (nfft - i)/(nfft - imin));

	   data [2] = 0.0;
	   NUMrealft (data, 2 * nfft, -1);
	   double factor = 1.0 / nfft;
	   for (long i = 1; i <= thy nx; i ++) 
			thy z [channel] [i] = data [i + 2000] * factor;
   }
   return thee;
}

Sound Sound_resample(Sound me, double samplingFrequency, long precision)
{
	double upfactor = me->dx * samplingFrequency;
	if(fabs(upfactor - 2) < 1e-6)
		return Sound_upsample(me);
	if(fabs(upfactor - 1) < 1e-6){
		if(me == NULL) return NULL;
		Sound thee = Sound_create (my xmin, my xmax, my nx, my dx, my x1, my ny);
		memcpy(thee, me, sizeof(structSound));
		return thee;
	}

	long numberOfSamples = floor((my xmax - my xmin)*samplingFrequency + 0.5);
	if(numberOfSamples < 1){
		std::cout<<"The resampled Sound would have no samples."<<std::endl;
		std::cout<<"SoundCompute.cpp 74"<<std::endl;
		return NULL;
	}
	Sound filtered = NULL;
	if(upfactor < 1.0){                /* Need anti-aliasing filter? ¿¹»ìµþÂË²¨Æ÷´¦Àí*/
		long nfft = 1, antiTurnAround = 1000;
		while(nfft < my nx + antiTurnAround * 2) nfft *= 2;
		double *data = (double *)malloc(sizeof(double) * (nfft + 1));
		memset(data, 0.0, sizeof(double) * (nfft + 1));
		filtered = Sound_create (my xmin, my xmax, my nx, my dx, my x1, my ny);
		for (long channel = 1; channel <= my ny; channel ++) {
		    memset(data, 0.0, (sizeof(double) * (nfft + 1)));
		    memcpy(data + antiTurnAround, my z[channel], sizeof(double) * my nx);
		    NUMrealft (data, nfft, 1);   // go to the frequency domain
		    for (long i = floor (upfactor * nfft); i <= nfft; i ++) 
				data [i] = 0;                        /* Filter away high frequencies. */
		    data [2] = 0.0;
		    NUMrealft (data, nfft, -1);   // return to the time domain
			double factor = 1.0 / nfft;
			double *to = filtered -> z [channel];
			for (long i = 1; i <= my nx; i ++) 
				to [i] = data [i + antiTurnAround] * factor;
		}
		free(me);
		me = filtered;   // reference copy; remove at end
	}

	Sound thee = Sound_create (my xmin, my xmax, numberOfSamples, 1.0 / samplingFrequency,
			0.5 * (my xmin + my xmax - (numberOfSamples - 1) / samplingFrequency), my ny);
	for (long channel = 1; channel <= my ny; channel ++) {
		double *from = my z [channel];
		double *to = thy z [channel];
		if (precision <= 1) {
			for (long i = 1; i <= numberOfSamples; i ++) {
				double x = thy x1 + (i - 1) * thy dx;  /* Sampled_indexToX (thee, i); */
				double index = (x - my x1) / my dx + 1;  /* Sampled_xToIndex (me, x); */
				long leftSample = floor (index);
				double fraction = index - leftSample;
				to [i] = leftSample < 1 || leftSample >= my nx ? 0.0 :
					(1 - fraction) * from [leftSample] + fraction * from [leftSample + 1];
			}
		} else {
			for (long i = 1; i <= numberOfSamples; i ++) {
				double x = thy x1 + (i - 1) * thy dx;   /* Sampled_indexToX (thee, i); */
				double index = (x - my x1) / my dx + 1;   /* Sampled_xToIndex (me, x); */
				to [i] = NUM_interpolate_sinc (my z [channel], my nx, index, precision);
			}
		}
	}  
	return thee;
}

void Sound_getMaximumAndX (Sound me, double xmin, double xmax, long channel, int interpolation,
	                        double *return_maximum, double *return_xOfMaximum)
{
	long imin, imax, i, n = my nx;
	if(channel < 1 || channel > my ny){
	    std::cout<<"channel = "<< channel 
			<<". dosen't fit the conditions:channel >= 1 && channel <= my ny."<<std::endl; 
		std::cout<<"SoundCompute.cpp 142"<<std::endl;
	    return ;
	}
  /****	Melder_assert (channel >= 1 && channel <= my ny);  ****/
	double *y = my z[channel];
	double maximum, x;
	if (xmax <= xmin) { 
	    xmin = my xmin; 
		xmax = my xmax; 
	 }
	if (!Sampled_getWindowSamples (me, xmin, xmax, & imin, & imax)) {
		/*
		 * No samples between xmin and xmax.
		 * Try to return the greater of the values at these two points.
		 */
		double yleft = Sound_getValueAtX (me, xmin, channel,
			interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : 
							Vector_VALUE_INTERPOLATION_NEAREST);
		double yright = Sound_getValueAtX (me, xmax, channel,
			interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : 
							Vector_VALUE_INTERPOLATION_NEAREST);
		maximum = yleft > yright ? yleft : yright;
		x = yleft == yright ? (xmin + xmax) / 2 : yleft > yright ? xmin : xmax;
	} else {
		maximum = y[imin], x = imin;
		if (y[imax] > maximum) maximum = y[imax], x = imax;
		if (imin == 1) imin ++;
		if (imax == my nx) imax --;
		for (i = imin; i <= imax; i ++) {
			if (y[i] > y[i - 1] && y[i] >= y[i + 1]) {
				double i_real, localMaximum = NUMimproveMaximum (y, n, i, interpolation, & i_real);
				if (localMaximum > maximum) maximum = localMaximum, x = i_real;
			}
		}
		x = my x1 + (x - 1) * my dx;   /* Convert sample to x. */
		if (x < xmin) x = xmin; else if (x > xmax) x = xmax;
	}
	if (return_maximum) *return_maximum = maximum;
	if (return_xOfMaximum) *return_xOfMaximum = x;
}
	
void Sound_getMinimumAndX (Sound me, double xmin, double xmax, long channel, int interpolation,
	                       double *return_minimum, double *return_xOfMinimum)
/**void Vector_getMinimumAndX (Vector me, double xmin, double xmax, long channel, int interpolation,
	double *return_minimum, double *return_xOfMinimum)**/						   
{
	long imin, imax, n = my nx;
	if(channel < 1 || channel > my ny){
	    std::cout<<"channel = "<< channel 
			<<". dosen't fit the conditions:channel >= 1 && channel <= my ny."<<std::endl; 
		std::cout<<"SoundCompute.cpp 190"<<std::endl;
	    exit(0);
	}
	/****** Melder_assert (channel >= 1 && channel <= my ny); ******/
	double *y = my z[channel];
	double minimum, x;
	if (xmax <= xmin) {
     	xmin = my xmin; 
		xmax = my xmax; 
	}
	
	if (! Sampled_getWindowSamples (me, xmin, xmax, & imin, & imax)) {
		/*
		 * No samples between xmin and xmax. 
		 * Try to return the lesser of the values at these two points.
		 */
		double yleft = Sound_getValueAtX (me, xmin, channel,
			interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : 
							Vector_VALUE_INTERPOLATION_NEAREST);
		double yright = Sound_getValueAtX (me, xmax, channel,
			interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : 
							Vector_VALUE_INTERPOLATION_NEAREST);
		minimum = yleft < yright ? yleft : yright;
		x = yleft == yright ? (xmin + xmax) / 2 : yleft < yright ? xmin : xmax;
	} else {
		minimum = y[imin], x = imin;
		if (y[imax] < minimum) minimum = y[imax], x = imax;
		if (imin == 1) imin ++;
		if (imax == my nx) imax --;
		for (long i = imin; i <= imax; i ++) {
			if (y[i] < y[i - 1] && y[i] <= y[i + 1]) {
				double i_real, localMinimum = NUMimproveMinimum (y, n, i, interpolation, & i_real);
				if (localMinimum < minimum) minimum = localMinimum, x = i_real;
			}
		}
		x = my x1 + (x - 1) * my dx;        /* Convert sample to x. */
		if (x < xmin) x = xmin; else if (x > xmax) x = xmax;
	}
	if (return_minimum) *return_minimum = minimum;
	if (return_xOfMinimum) *return_xOfMinimum = x;
}

void Sound_getMaximumAndXAndChannel(Sound me, double xmin, double xmax, int interpolation,
	                       double *return_maximum, double *return_xOfMaximum, long *return_channelOfMaximum)
{
	double maximum, xOfMaximum;
	long channelOfMaximum = 1;
	Sound_getMaximumAndX (me, xmin, xmax, 1, interpolation, & maximum, & xOfMaximum);
	for (long channel = 2; channel <= my ny; channel ++) {
		double maximumOfChannel, xOfMaximumOfChannel; 
		Sound_getMaximumAndX (me, xmin, xmax, channel, interpolation, 
							 &maximumOfChannel, &xOfMaximumOfChannel);
		if (maximumOfChannel > maximum) {
			maximum = maximumOfChannel;
			xOfMaximum = xOfMaximumOfChannel;
			channelOfMaximum = channel;
		}
	}
	if (return_maximum) *return_maximum = maximum;
	if (return_xOfMaximum) *return_xOfMaximum = xOfMaximum;
	if (return_channelOfMaximum) *return_channelOfMaximum = channelOfMaximum;
}

void Sound_getMinimumAndXAndChannel (Sound me, double xmin, double xmax, int interpolation,
	                        double *return_minimum, double *return_xOfMinimum, long *return_channelOfMinimum)
/**void Vector_getMinimumAndXAndChannel (Vector me, double xmin, double xmax, int interpolation,
	double *return_minimum, double *return_xOfMinimum, long *return_channelOfMinimum)**/
{
	double minimum, xOfMinimum;
	long channelOfMinimum = 1;
	Sound_getMinimumAndX(me, xmin, xmax, 1, interpolation, & minimum, & xOfMinimum);
	for (long channel = 2; channel <= my ny; channel ++) {
		double minimumOfChannel, xOfMinimumOfChannel;
		Sound_getMinimumAndX(me, xmin, xmax, channel, interpolation, 
							&minimumOfChannel, &xOfMinimumOfChannel);
		if (minimumOfChannel < minimum) {
			minimum = minimumOfChannel;
			xOfMinimum = xOfMinimumOfChannel;
			channelOfMinimum = channel;
		}
	}
	if (return_minimum) *return_minimum = minimum;
	if (return_xOfMinimum) *return_xOfMinimum = xOfMinimum;
	if (return_channelOfMinimum) *return_channelOfMinimum = channelOfMinimum;
}

double Sound_getMaximum(Sound me, double xmin, double xmax, int interpolation) {   
	/**Vector_getMaximum (Vector me, double xmin, double xmax, int interpolation)**/
	double maximum;
	Sound_getMaximumAndXAndChannel(me, xmin, xmax, interpolation, & maximum, NULL, NULL);
	return maximum;
}

double Sound_getMinimum (Sound me, double xmin, double xmax, int interpolation) {  
   /**Vector_getMinimum (Vector me, double xmin, double xmax, int interpolation)**/
	double minimum;
	Sound_getMinimumAndXAndChannel (me, xmin, xmax, interpolation, & minimum, NULL, NULL);
	return minimum;
}

double Sound_getAbsoluteExtremum(Sound me, double xmin, double xmax, int interpolation){   /**double Vector_getAbsoluteExtremum (Vector me, double xmin, double xmax, int interpolation)**/
    double minimum = fabs (Sound_getMinimum (me, xmin, xmax, interpolation));      /**Vector_getMinimum**/
	double maximum = fabs (Sound_getMaximum (me, xmin, xmax, interpolation));      /**Vector_getMaximum**/
	return minimum > maximum ? minimum : maximum;
}

long Sampled_getWindowSamples (Sound me, double xmin, double xmax, long *ixmin, long *ixmax) {
	double rixmin = 1.0 + ceil ((xmin - my x1) / my dx);
	double rixmax = 1.0 + floor ((xmax - my x1) / my dx);
	*ixmin = rixmin < 1.0 ? 1 : (long) rixmin;
	*ixmax = rixmax > (double) my nx ? my nx : (long) rixmax;
	if (*ixmin > *ixmax) return 0;
	return *ixmax - *ixmin + 1;
}

// (Vector) Sound_getValueAtX () returns the average of all the interpolated channels.
double Sound_getValueAtX (Sound me, double x, long ilevel, int interpolation) {
	double leftEdge = my x1 - 0.5 * my dx, rightEdge = leftEdge + my nx * my dx;
	if (x <  leftEdge || x > rightEdge) return NUMundefined;
	if (ilevel > Vector_CHANNEL_AVERAGE) {
	   if(ilevel > my ny ){
	      std::cout<<"Errot, ilevel = " << ilevel <<" my n = "<<my ny
				   <<" doextn't fit the condition: ilevel <= my ny."<<std::endl;
		  std::cout<<" SoundCompute.cpp Line 316. "<<std::endl;
		  exit(0);
	   }
	/****	Melder_assert (ilevel <= my ny);   ****/
		return NUM_interpolate_sinc(my z[ilevel], my nx, ((x - my x1) / my dx + 1), //Sampled_xToIndex (me, x),
			interpolation == Vector_VALUE_INTERPOLATION_SINC70 ? NUM_VALUE_INTERPOLATE_SINC70 :
			interpolation == Vector_VALUE_INTERPOLATION_SINC700 ? NUM_VALUE_INTERPOLATE_SINC700 :
			interpolation);
	}
	double sum = 0.0;
	for (long channel = 1; channel <= my ny; channel ++) { 
		sum += NUM_interpolate_sinc(my z[channel], my nx, ((x - my x1) / my dx + 1),//Sampled_xToIndex (me, x),
			interpolation == Vector_VALUE_INTERPOLATION_SINC70 ? NUM_VALUE_INTERPOLATE_SINC70 :
			interpolation == Vector_VALUE_INTERPOLATION_SINC700 ? NUM_VALUE_INTERPOLATE_SINC700 :
			interpolation);
	}
	return sum / my ny;
}

VoicedUnvoiced GeVicedUnvoiced(Pitch pitch, PointProcess points, double periodsPerWindow, double minimumPitch)
{
	if(pitch == NULL && points == NULL)
		return NULL;

	double dt_window, halfdt_window;
	long k = 1;

	VoicedUnvoiced VUV = (VoicedUnvoiced)malloc(sizeof(structVoicedUnvoiced));

    VUV->nt = (long)ceil(double(pitch->nx / 2));
	VUV->vuv = (bool *)malloc(sizeof(bool) *(VUV->nt + 1));
	
	for(long i = 0; i < VUV->nt; ++ i)
		VUV->vuv[i] = false;

	dt_window = periodsPerWindow / minimumPitch;
	halfdt_window = dt_window / 2;

	if (dt_window < 0){
		std::cout<<"Analysis window too short."<<std::endl;
		std::cout<<"SoundCompute.cpp: Line 343"<<std::endl;
	       return NULL;		
	}

	for(long i = 1; i <= pitch->nx; i += 2){
		double start_dt = pitch->x1 + (i - 1) * pitch->dx - halfdt_window;
		double end_dt = pitch->x1 + (i - 1) * pitch->dx + halfdt_window;
		for(long j = 1; j <= points->nt; ){
			if(points->t[j] > end_dt){
				 k = j;  //not take overlaping of window into considerration 
			     break;
			}
			else if(points->t[j] >= start_dt && points->t[j] <= end_dt){
			    VUV->vuv[(i + 1)/2] = true;
			    ++ j;
			    continue;
			}
		    else{
				++ j;
				continue;
			}
		}
	}
	return VUV;
}

void Sound_substructMean(Sound me){
    for(long channel = 1; channel <= my ny; channel ++){
		double sum = 0.0;
		for(long i = 1 ; i < my nx; ++ i)
			sum += my z [channel] [i];
		for(long i = 1; i < my nx; ++ i)
			my z [channel] [i] -= sum / my nx;
	}
}

void Sound_overrideSamplingFrequency(Sound me, double rate) {
	my dx = 1 / rate;
	my x1 = my xmin + 0.5 * my dx;
	my xmax = my xmin + my nx * my dx;
}

Sound Sound_changespeaker(Sound me, double pitchMin, double pitchMax,  double formantMultiplier ,  // > 0
	                       double pitchMultiplier,  // > 0
						   double pitchRangeMultiplier,  // any number
						   double durationMultiplier  // > 0
						   )
{
	Pitch pitch = Sound_to_Pitch(me, 0.8/pitchMin, pitchMin, pitchMax);
	Sound thee = Sound_and_pitch_changespeaker(me, pitch, formantMultiplier, pitchMultiplier,
											   pitchRangeMultiplier, durationMultiplier);
	return thee;
}

Sound Sound_and_pitch_changespeaker(Sound me, Pitch him, 
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
	PitchTier pitchTier = Pitch_to_PitchTier(pitch);

	double median = Pitch_getQuantile(pitch, 0, 0, 0.5, kPitch_unit_HERTZ);

	if(median != 0 && median != NUMundefined){
		PitchTier_multiplyFrequencies(pitchTier, sound->xmin, sound->xmax, pitchMultiplier/formantMultiplier);
		PitchTier_modifyExcursionRange(pitchTier, sound->xmin, sound->xmax, pitchRangeMultiplier, median);
	}else if (pitchMultiplier != 1){
		std::cout<<"Warning: Pitch has not been changed because the sound was entirely voiceless."<<std::endl;
	    std::cout<<"SoundCompute.cpp: Line 418"<<std::endl;
		return NULL;//2013.3.14 Ìí¼Ó
	}

	DurationTier duration = DurationTier_create(my xmin, my xmax);
	RealTier_addPoint(duration, (my xmin + my xmax) / 2, formantMultiplier * durationMultiplier);

	Sound thee = Sound_Point_Pitch_Duration_to_Sound (sound, pulses, pitchTier, duration, MAX_T,formantMultiplier);

	if(formantMultiplier != 1)
	     thee = Sound_resample(thee, samplingFrequency_old, 10);
	return thee;
}

Sound Sound_changeGender_old (Sound me, double fmin, double fmax, double formantRatio,
	double new_pitch, double pitchRangeFactor, double durationFactor){
		Pitch pitch = Sound_to_Pitch(me, 0.8/fmin, fmin, fmax);
		Sound thee = Sound_and_Pitch_changeGender_old (me, pitch, formantRatio,
 		                 new_pitch, pitchRangeFactor, durationFactor);
		return thee;
}

Sound Sound_and_Pitch_changeGender_old (Sound me, Pitch him, double formantRatio,double new_pitch, 
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
	if (formantRatio != 1.0) 
		Sound_resample(thee, samplingFrequency_old, 10);

	return thee;
}

//Manipulation.cpp 331
Sound Sound_Point_Pitch_Duration_to_Sound (Sound me, PointProcess pulses,
										  PitchTier pitch, DurationTier duration, double maxT,double formantMultiplier)
{
	long ipointleft, ipointright;
	double deltat = 0, handledTime = my xmin;
	double startOfSourceNoise = 0, endOfSourceNoise = 0, startOfTargetNoise = 0, endOfTargetNoise = 0;
	double durationOfSourceNoise = 0, durationOfTargetNoise= 0 ;
	double startOfSourceVoice = 0, endOfSourceVoice = 0, startOfTargetVoice = 0, endOfTargetVoice = 0;
	double durationOfSourceVoice = 0, durationOfTargetVoice = 0;
	double startingPeriod = 0, finishingPeriod = 0, ttarget = 0, voicelessPeriod = 0;
	int mark = 1;//×÷ÎªÐÞÕûµÄ±ê¼Ç
//	double copytsource, copyperiod, copyisourcepulse, copyttarget;
//	double copytsource2, copyvoicelessperiod, copyttarget2;
//	int state = 0;

	if(duration->points->size == 0){
		std::cout<<"No duration points."<<std::endl;
		std::cout<<"SoundCompute.cpp 473"<<std::endl;
		return NULL;
	}

	//Create a Sound long enough to hold the longest possible duration-manipulated sound.
	Sound thee = Sound_create (my xmin, my xmin + 3 * (my xmax - my xmin), 3 * my nx, my dx, my x1,1);

	///Below, I'll abbreviate the voiced interval as "voice" and the voiceless interval as "noise".
    if (pitch && pitch -> points -> size)
		for (ipointleft = 1; ipointleft <= pulses->nt; ipointleft = ipointright + 1) {
	        //Find the beginning of the voice.
	        startOfSourceVoice = pulses -> t[ipointleft];    /* The first pulse of the voice. */
	        startingPeriod = 1.0 / RealTier_getValueAtTime (pitch, startOfSourceVoice);
	        startOfSourceVoice -= 0.5 * startingPeriod;      /* The first pulse is in the middle of a period. */

			//measure one noise.
			startOfSourceNoise = handledTime;
			endOfSourceNoise = startOfSourceVoice;
			durationOfSourceNoise = endOfSourceNoise - startOfSourceNoise;
			startOfTargetNoise = startOfSourceNoise + deltat;
			endOfTargetNoise = startOfTargetNoise + RealTier_getArea(duration, startOfSourceNoise, endOfSourceNoise);
			durationOfTargetNoise = endOfTargetNoise - startOfTargetNoise;

			//copy the noise
			voicelessPeriod = NUMrandomUniform (0.008, 0.012);
			ttarget = startOfTargetNoise + 0.5 * voicelessPeriod;

			while (ttarget < endOfTargetNoise) {
				
				double tsource;
				double tleft = startOfSourceNoise, tright = endOfSourceNoise;
				int i;
				for (i = 1; i <= 15; i ++) {
					double tsourcemid = 0.5 * (tleft + tright);
					double ttargetmid = startOfTargetNoise + RealTier_getArea (duration, startOfSourceNoise, tsourcemid);

 					if (ttargetmid < ttarget) 
						tleft = tsourcemid; 
				    else tright = tsourcemid;
				}
				tsource = 0.5 * (tleft + tright);
				/*
				if(ipointleft == 1 && mark == 1)
				{
					copyFall(me, tsource, tsource + voicelessPeriod, thee, ttarget);
					mark = 0;
				}*/
				copyBell (me, tsource, voicelessPeriod, voicelessPeriod, thee, ttarget);
				ttarget += voicelessPeriod;
				/*\
				if(ttarget >= endOfTargetNoise)
				{	
					copytsource2 = tsource;
					copyvoicelessperiod = voicelessPeriod;
					copyttarget2 = ttarget - voicelessPeriod;
					state = 1;
				}
				\*/
				voicelessPeriod = NUMrandomUniform (0.008, 0.012);
			}
			deltat += durationOfTargetNoise - durationOfSourceNoise;

			//find the end of the voice
			for (ipointright = ipointleft + 1; ipointright <= pulses -> nt; ipointright ++)

				if (pulses -> t [ipointright] - pulses -> t [ipointright - 1] > maxT)
					  break;

			ipointright --;
			endOfSourceVoice = pulses -> t [ipointright];   /* The last pulse of the voice. */
			finishingPeriod = 1.0 / RealTier_getValueAtTime (pitch, endOfSourceVoice);
			endOfSourceVoice += 0.5 * finishingPeriod;   /* The last pulse is in the middle of a period. */		
			 
			// Measure one voice.
			durationOfSourceVoice = endOfSourceVoice - startOfSourceVoice;
		
			// This will be copied to an interval with a different location and duration.
			startOfTargetVoice = startOfSourceVoice + deltat;
			endOfTargetVoice = startOfTargetVoice + 
							RealTier_getArea (duration, startOfSourceVoice, endOfSourceVoice);
			durationOfTargetVoice = endOfTargetVoice - startOfTargetVoice;
            
			ttarget = startOfTargetVoice + 0.5 * startingPeriod;
			while (ttarget < endOfTargetVoice) {
				double tsource, period;
				long isourcepulse;
				double tleft = startOfSourceVoice, tright = endOfSourceVoice;
				int i;
				for (i = 1; i <= 15; i ++) {
					double tsourcemid = 0.5 * (tleft + tright);
					double ttargetmid = startOfTargetVoice +	
								 RealTier_getArea (duration, startOfSourceVoice, tsourcemid);
					if (ttargetmid < ttarget) tleft = tsourcemid; else tright = tsourcemid;
				}
				tsource = 0.5 * (tleft + tright);
				period = 1.0 / RealTier_getValueAtTime (pitch, tsource);
				isourcepulse = PointProcess_getNearestIndex (pulses, tsource);
				copyBell2 (me, pulses, isourcepulse, period, period, thee, ttarget, maxT);
				/*
				if(ipointleft == 1 && mark == 1)
				{
					copyFall2(me, pulses, isourcepulse, period, period, thee, ttarget, maxT);
					mark = 0;
				}*/
				ttarget += period;
				/*\
				if(ttarget >= endOfTargetNoise && ipointright >= pulses->nt)
				{
					copytsource = tsource;
					copyperiod = period;
					copyisourcepulse = isourcepulse;
					copyttarget = ttarget - period;
					state = 2;
				}
				\*/
			}
			deltat += durationOfTargetVoice - durationOfSourceVoice;
			handledTime = endOfSourceVoice;
		}
		
		// Copy the remaining unvoiced part, if we are at the end.
		startOfSourceNoise = handledTime;
		endOfSourceNoise = my xmax;
		durationOfSourceNoise = endOfSourceNoise - startOfSourceNoise;
		startOfTargetNoise = startOfSourceNoise + deltat;
		endOfTargetNoise = startOfTargetNoise + 
						RealTier_getArea (duration, startOfSourceNoise, endOfSourceNoise);
		durationOfTargetNoise = endOfTargetNoise - startOfTargetNoise;
		voicelessPeriod = NUMrandomUniform (0.008, 0.012);
		ttarget = startOfTargetNoise + 0.5 * voicelessPeriod;
		while (ttarget < endOfTargetNoise) {
			double tsource;
			double tleft = startOfSourceNoise, tright = endOfSourceNoise;
			for (int i = 1; i <= 15; i ++) {
				double tsourcemid = 0.5 * (tleft + tright);
				double ttargetmid = startOfTargetNoise +
							   RealTier_getArea (duration, startOfSourceNoise, tsourcemid);
				if (ttargetmid < ttarget) 
					tleft = tsourcemid; 
				else tright = tsourcemid;
			}
			tsource = 0.5 * (tleft + tright);
			/*
			if(mark == 1)
			{
				copyFall(me, tsource, tsource + voicelessPeriod, thee, ttarget);
				mark = 0;
			}*/
			copyBell (me, tsource, voicelessPeriod, voicelessPeriod, thee, ttarget);
			ttarget += voicelessPeriod;
			/*\
			if(ttarget >= endOfTargetNoise)
			{		
				copytsource2 = tsource;
				copyvoicelessperiod = voicelessPeriod;
				copyttarget2 = ttarget - voicelessPeriod;
				state = 3;
			}
			\*/
			voicelessPeriod = NUMrandomUniform (0.008, 0.012);
		}


//		switch (state){
//			case 1:copyRise (me, copytsource2 - copyvoicelessperiod, copytsource2, thee, copyttarget2);break;
//			case 2:copyRise2 (me, pulses, copyisourcepulse, copyperiod, copyperiod, thee, copyttarget, maxT);break;
//			case 3:copyRise (me, copytsource2 - copyvoicelessperiod, copytsource2, thee, copyttarget2);break;}
		
		
		// Find the number of trailing zeroes and hack the sound's time domain.
		
		//thy xmax = my xmax;//thy xmin + RealTier_getArea (duration, my xmin, my xmax);
		
		//startOfTargetNoise;

		thy xmax = thy xmin + RealTier_getArea (duration, my xmin, my xmax);
		if (fabs (thy xmax - my xmax) < 1e-12) thy xmax = my xmax;  // Common situation.
		thy nx = (long)floor((thy xmax - thy x1)/thy dx) + 1;  //Sampled_xToLowIndex (thee, thy xmax);
		if (thy nx > 3 * my nx) thy nx = 3 * my nx;

		return thee;
}

//Manipulation.cpp 218
void copyRise (Sound me, double tmin, double tmax, Sound thee, double tmaxTarget) {
	long imin, imax, imaxTarget, distance, i;
	double dphase;
	imin = (long)ceil((tmin - my x1)/ my dx) + 1;   //Sampled_xToHighIndex (me, tmin);
	if (imin < 1) imin = 1;
	imax = (long)ceil((tmax - my x1) / my dx);   //Sampled_xToHighIndex (me, tmax) - 1;  
	      /* Not xToLowIndex: ensure separation of subsequent calls. */

	if (imax > my nx) imax = my nx;
	if (imax < imin) return;
	imaxTarget = (long)ceil((tmaxTarget - thy x1)/thy dx); //Sampled_xToHighIndex (thee, tmaxTarget) - 1;
	distance = imaxTarget - imax;
	dphase = NUMpi / (imax - imin + 1);
	for (i = imin; i <= imax; i ++) {
		long iTarget = i + distance;
		if (iTarget >= 1 && iTarget <= thy nx)
			thy z[1][iTarget] += my z[1][i] * 0.5 * (1 - cos (dphase * (i - imin + 0.5)));
	}
}

//Manipulation.cpp 236
void copyFall (Sound me, double tmin, double tmax, Sound thee, double tminTarget) {
	long imin, imax, iminTarget, distance, i;
	double dphase;
	imin = (long)ceil((tmin - my x1)/ my dx) + 1;  //Sampled_xToHighIndex (me, tmin);
	if (imin < 1) imin = 1;
	imax = (long)ceil((tmax - my x1) / my dx);   //Sampled_xToHighIndex (me, tmax) - 1;   
	  /* Not xToLowIndex: ensure separation of subsequent calls. */

	if (imax > my nx) imax = my nx;
	if (imax < imin) return;
	iminTarget = (long)ceil((tminTarget - thy x1)/thy dx) + 1;//Sampled_xToHighIndex (thee, tminTarget);
	distance = iminTarget - imin;
	dphase = NUMpi / (imax - imin + 1);
	for (i = imin; i <= imax; i ++) {
		long iTarget = i + distance;
		if (iTarget >= 1 && iTarget <= thy nx)
			thy z[1][iTarget] += my z[1][i] * 0.5 * (1 + cos(dphase * (i - imin + 0.5)));
	}
}

//Manipulation.cpp 254
void copyBell (Sound me, double tmid, double leftWidth, double rightWidth, Sound thee, double tmidTarget) {
	 copyRise (me, tmid - leftWidth, tmid, thee, tmidTarget);
	 copyFall (me, tmid, tmid + rightWidth, thee, tmidTarget);
}

//Manipulation.cpp 259
void copyBell2 (Sound me, PointProcess source, long isource, double leftWidth, double rightWidth,
	            Sound thee, double tmidTarget, double maxT)
{
	/*
	 * Replace 'leftWidth' and 'rightWidth' by the lengths of the intervals in the source (instead of target),
	 * if these are shorter.
	 */
	double tmid = source->t[isource];
	if (isource > 1 && tmid - source->t[isource - 1] <= maxT) {
		double sourceLeftWidth = tmid - source->t[isource - 1];
		if(sourceLeftWidth < leftWidth) leftWidth = sourceLeftWidth;
	}
	if (isource < source->nt && source->t[isource + 1] - tmid <= maxT) {
		double sourceRightWidth = source->t[isource + 1] - tmid;
		if (sourceRightWidth < rightWidth) rightWidth = sourceRightWidth;
	}
	copyBell(me, tmid, leftWidth, rightWidth, thee, tmidTarget);
}

//2013.3.16ÓÃÓÚ±ß½ç¼Ó´°´¦Àí
void copyFall2(Sound me, PointProcess source, long isource, double leftWidth, double rightWidth,
	            Sound thee, double tmidTarget, double maxT)
{
	double tmid = source->t[isource];
	if (isource > 1 && tmid - source->t[isource - 1] <= maxT) {
		double sourceLeftWidth = tmid - source->t[isource - 1];
		if(sourceLeftWidth < leftWidth) leftWidth = sourceLeftWidth;
	}
	if (isource < source->nt && source->t[isource + 1] - tmid <= maxT) {
		double sourceRightWidth = source->t[isource + 1] - tmid;
		if (sourceRightWidth < rightWidth) rightWidth = sourceRightWidth;
	}
	copyFall (me, tmid, tmid + rightWidth, thee, tmidTarget);
}

//2013.3.16ÓÃÓÚ±ß½ç¼Ó´°´¦Àí
void copyRise2 (Sound me, PointProcess source, long isource, double leftWidth, double rightWidth,
	            Sound thee, double tmidTarget, double maxT)
{
	double tmid = source->t[isource];
	if (isource > 1 && tmid - source->t[isource - 1] <= maxT) {
		double sourceLeftWidth = tmid - source->t[isource - 1];
		if(sourceLeftWidth < leftWidth) leftWidth = sourceLeftWidth;
	}
	if (isource < source->nt && source->t[isource + 1] - tmid <= maxT) {
		double sourceRightWidth = source->t[isource + 1] - tmid;
		if (sourceRightWidth < rightWidth) rightWidth = sourceRightWidth;
	}
	copyRise (me, tmid - leftWidth, tmid, thee, tmidTarget);
}
