
#include <iostream>
#include "Sound_to_Pitch.h"
#include "Structure.h"
#include "SoundCompute.h"
#include "Pitch.h"
#include "NUM.h"
#include "NUM2.h"

//timeStepStrategy  value
#define kTimeSoundAnalysisEditor_timeStepStrategy_AUTOMATIC 1
#define kTimeSoundAnalysisEditor_timeStepStrategy_FIXED 2
#define kTimeSoundAnalysisEditor_timeStepStrategy_VIEW_DEPENDENT 3 
 
//method  value
#define kTimeSoundAnalysisEditor_pitch_analysisMethod_AUTOCORRELATION  1
#define kTimeSoundAnalysisEditor_pitch_analysisMethod_CROSS_CORRELATION 2

//windows kinds
#define AC_HANNING  0
#define AC_GAUSS  1
#define FCC_NORMAL  2
#define FCC_ACCURATE  3

#define NUM_PEAK_INTERPOLATE_SINC70  3
#define NUM_PEAK_INTERPOLATE_SINC700  4


int Sampled_shortTermAnalysis (Sound me, double windowDuration, double timeStep, long *numberOfFrames, double *firstTime) {
	if (windowDuration <= 0.0)   return 0;  
	if (timeStep <= 0.0)  return 0;
	
	 double myDuration = my dx * my nx;
	if (windowDuration > myDuration){
		std::cout<<"Sound shorter than window length."<<std::endl; 
		std::cout<<"Sound_to_Pitch.cpp: Line 13"<<std::endl;
		return 0;
     }
	 
	*numberOfFrames = floor((myDuration - windowDuration) / timeStep) + 1; 
	if (*numberOfFrames < 1)  return -1;
	double ourMidTime = my x1 - 0.5 * my dx + 0.5 * myDuration; 
	double thyDuration = *numberOfFrames * timeStep;                  
	*firstTime = ourMidTime - 0.5 * thyDuration + 0.5 * timeStep;
	return 1;
}

Pitch Sound_to_Pitch_any (Sound me, double dt,     /*timeStepStradygy related*/
                         double minimumPitch,      /*Pitch settings realted*/
						 double periodsPerWindow,  /*kTimeSoundAnalysisEditor_pitch_analysisMethod  related*/
						 int maxnCandidates, 
						 int method,               /*method related*/
                         double silenceThreshold, double voicingThreshold, double octaveCost, double octaveJumpCost, 
						 double voicedUnvoicedCost, double ceiling)
{
	  NUMfft_Table fftTable = NUMfft_Table_create();
	  double duration, t1;
	  double dt_window;                       /* Window length in seconds. */
	  long nsamp_window, halfnsamp_window;   /* Number of samples per window. */
	  long nFrames, minimumLag, maximumLag;
	  long iframe, nsampFFT;
	  double interpolation_depth;
	  long nsamp_period, halfnsamp_period;   /* Number of samples in longest period. */
	  long brent_ixmax, brent_depth;
	  double brent_accuracy;                 /* Obsolete. */
	  double globalPeak;
	  
	   if (maxnCandidates < 2 || method < AC_HANNING && method > FCC_ACCURATE)
	   {
	       std::cout<<"Error: maxnCandidates: "<<maxnCandidates<<" method: "<<method<<"."<<std::endl;
		   std::cout<<"Sound_to_Pitch.cpp: Line 13. 69"<<std::endl;
		   return NULL;
	   }
	  
	   if (maxnCandidates < ceiling / minimumPitch) maxnCandidates = ceiling / minimumPitch;
 
	   if (dt <= 0.0) dt = periodsPerWindow / minimumPitch / 4.0;  /* e.g. 3 periods, 75 Hz: 10 milliseconds. */

		switch (method) {
			case AC_HANNING:
				brent_depth = NUM_PEAK_INTERPOLATE_SINC70;
				brent_accuracy = 1e-7;
				interpolation_depth = 0.5;
				break;
			case AC_GAUSS:
				periodsPerWindow *= 2;       /* Because Gaussian window is twice as long. */
				brent_depth = NUM_PEAK_INTERPOLATE_SINC700;
				brent_accuracy = 1e-11;
				interpolation_depth = 0.25;   /* Because Gaussian window is twice as long. */
				break;
			case FCC_NORMAL:
				brent_depth = NUM_PEAK_INTERPOLATE_SINC70;
				brent_accuracy = 1e-7;
				interpolation_depth = 1.0;
				break;
			case FCC_ACCURATE:
				brent_depth = NUM_PEAK_INTERPOLATE_SINC700;
				brent_accuracy = 1e-11;
				interpolation_depth = 1.0;
				break;
		}
		duration = my dx * my nx;
		if (minimumPitch < periodsPerWindow / duration) {
		     std::cout<<"To analyse this Sound, minimum pitch must not be less than "<< periodsPerWindow / duration<<" Hz."<<std::endl;
			 std::cout<<"Sound_to_Pitch.cpp: Line 31.103"<<std::endl;
			 return NULL;
		}
		
	   /*
		 * Determine the number of samples in the longest period.
		 * We need this to compute the local mean of the sound (looking one period in both directions),
		 * and to compute the local peak of the sound (looking half a period in both directions).
		 */
		nsamp_period = floor(1 / my dx / minimumPitch);
		halfnsamp_period = nsamp_period / 2 + 1;

		if (ceiling > 0.5 / my dx) ceiling = 0.5 / my dx;
		
	    // Determine window length in seconds and in samples.
		dt_window = periodsPerWindow / minimumPitch;
		nsamp_window = floor (dt_window / my dx);
		halfnsamp_window = nsamp_window / 2 - 1;
		if (halfnsamp_window < 2){
			std::cout<<"Analysis window too short."<<std::endl;
			std::cout<<"Sound_to_Pitch.cpp: Line 31.123"<<std::endl;
	        return NULL;		
		}
		nsamp_window = halfnsamp_window * 2;
		
	    // Determine the minimum and maximum lags.
		minimumLag = floor (1 / my dx / ceiling);
		if (minimumLag < 2) minimumLag = 2;
		maximumLag = floor (nsamp_window / periodsPerWindow) + 2;
		if (maximumLag > nsamp_window) maximumLag = nsamp_window;

		/*
		 * Determine the number of frames.
		 * Fit as many frames as possible symmetrically in the total duration.
		 * We do this even for the forward cross-correlation method,
		 * because that allows us to compare the two methods.
		 */  
	   if(!Sampled_shortTermAnalysis (me, method >= FCC_NORMAL ? 1 / minimumPitch + dt_window : dt_window, dt, & nFrames, & t1)){
           std::cout<<"The pitch analysis would give zero pitch frames."<<std::endl;   
           std::cout<<"Sound_to_Pitch.cpp: Line 31.142"<<std::endl;		   
		   return NULL;
	   }
	   	
	  // Create the resulting pitch contour. 
	    Pitch thee = Pitch_create (my xmin, my xmax, nFrames, dt, t1, ceiling, maxnCandidates);     
       
	   // Compute the global absolute peak for determination of silence threshold.
		globalPeak = 0.0;
		for (long channel = 1; channel <= my ny; channel ++) {
			double mean = 0.0;
			for (long i = 1; i <= my nx; i ++) {
				mean += my z [channel] [i];
			}
			mean /= my nx;
			for (long i = 1; i <= my nx; i ++) {
				double value = fabs (my z [channel] [i] - mean);
				if (value > globalPeak) globalPeak = value;
			}
		}
		if (globalPeak == 0.0)   return thee;
		
	   double **frame, *ac, *window, *windowR;	
	   
	   if (method >= FCC_NORMAL) {   /* For cross-correlation analysis. */			
		   // Create buffer for cross-correlation analysis.
		    frame = (double **)malloc(sizeof(double *) * (my ny + 1));
			for(long i = 1; i <= my ny; ++ i){
			   frame[i] = (double *)malloc(sizeof(double) * (nsamp_window + 1));
			   for(long j = 1; j <= nsamp_window; ++ j)
			      frame[i][j] = 0.0;
		    }   /****frame.reset (1, my ny, 1, nsamp_window);****/
				  
			brent_ixmax = nsamp_window * interpolation_depth;
		} else {   /* For autocorrelation analysis. */		   
		   /*
			* Compute the number of samples needed for doing FFT.
			* To avoid edge effects, we have to append zeroes to the window.
			* The maximum lag considered for maxima is maximumLag.
			* The maximum lag used in interpolation is nsamp_window * interpolation_depth.
			*/
			nsampFFT = 1; 
			while (nsampFFT < nsamp_window * (1 + interpolation_depth))  nsampFFT *= 2;
			
			// Create buffers for autocorrelation analysis.
		    frame = (double **)malloc(sizeof(double *) * (my ny + 1));
			for(long i = 1; i <= my ny; ++ i){
			   frame [i] = (double *)malloc(sizeof(double) * (nsampFFT + 1));
			   for(long j = 0; j <= nsampFFT; ++ j)
			      frame[i][j] = 0.0;
		    }  /****frame.reset (1, my ny, 1, nsampFFT);****/
			
			window = (double *)malloc(sizeof(double) * (nsamp_window + 1));
			for(long i = 0; i <= nsamp_window; ++ i)
			     window[i] = 0.0;
			/****window.reset (1, nsamp_window);****/		
			
			windowR = (double *)malloc(sizeof(double) * (nsampFFT + 1));
			ac = (double *)malloc(sizeof(double) * (nsampFFT + 1));
			for(long i = 0; i <= nsampFFT; ++ i)
			     windowR[i] = ac[i] = 0.0;
		     /****windowR.reset (1, nsampFFT); ac.reset (1, nsampFFT); ****/
			
			NUMfft_Table_init (fftTable, nsampFFT);
			
			/*
			* A Gaussian or Hanning window is applied against phase effects.
			* The Hanning window is 2 to 5 dB better for 3 periods/window.
			* The Gaussian window is 25 to 29 dB better for 6 periods/window.
			*/
			if (method == AC_GAUSS) { /* Gaussian window. */
				double imid = 0.5 * (nsamp_window + 1), edge = exp (-12.0);
				for (long i = 1; i <= nsamp_window; i ++)
					window[i] = (exp(-48.0*(i-imid)*(i-imid) /
						(nsamp_window + 1) / (nsamp_window + 1)) - edge) / (1 - edge);
			} else {  /* Hanning window*/
				for (long i = 1; i <= nsamp_window; i ++) 
					window [i] = 0.5 - 0.5 * cos (i * 2 * NUMpi / (nsamp_window + 1));
			}
			    
			// Compute the normalized autocorrelation of the window.
			for (long i = 1; i <= nsamp_window; i ++) windowR [i] = window [i];
			NUMfft_forward (fftTable, windowR);
			windowR [1] *= windowR [1];   // DC component
			for (long i = 2; i < nsampFFT; i += 2) {
				windowR [i] = windowR [i] * windowR [i] + windowR [i+1] * windowR [i+1];
				windowR [i + 1] = 0.0;   // power spectrum: square and zero
			}
			windowR [nsampFFT] *= windowR [nsampFFT];   // Nyquist frequency
			NUMfft_backward (fftTable, windowR);   // autocorrelation
			for (long i = 2; i <= nsamp_window; i ++) windowR [i] /= windowR [1];   // normalize
			windowR [1] = 1.0;   // normalize

			brent_ixmax = nsamp_window * interpolation_depth;
		}
		
	   double *r = (double *) malloc( sizeof(double) * (2 * (nsamp_window + 1) + 1) );
	   r += nsamp_window + 1;                                       //make "r" become a symetrical vectr 
	   long *imax = (long *) malloc( sizeof(long) * (maxnCandidates + 1));
	   double *localMean = (double *) malloc( sizeof(double) * (my ny + 1));
	   
	   for(iframe = 1; iframe <= nFrames; iframe ++){
	        Pitch_Frame pitchFrame = & thy frame[iframe];
			double t = thy x1 + (iframe - 1) *(thy dx), localPeak;
			long leftSample = (long) floor((t - my x1) / my dx) + 1;
			long rightSample = leftSample + 1;
			long startSample, endSample;
			
		   for(long channel = 1; channel <= my ny; ++ channel){   //Compute the local mean; look one longest period to both sides.
			    startSample = rightSample - nsamp_period;
				endSample = leftSample + nsamp_period;
				if ( startSample < 0 ) {
				    std::cout<<"StartSample < 1"<<std::endl;
					std::cout<<"Sound_to_Pitch.cpp: Line 31"<<std::endl;
					return NULL;
				}
				
				if (endSample > my nx){
				    std::cout<<"EndSample > my nx"<<std::endl;
					std::cout<<"Sound_to_Pitch.cpp: Line 31.262"<<std::endl;
					return NULL;
				}
				
				localMean[channel] = 0.0;
				for (long i = startSample; i <= endSample; i ++) {    
					localMean[channel] += my z[channel][i];
				}
				localMean[channel] /= 2 * nsamp_period;
		
				// Copy a window to a frame and subtract the local mean. We are going to kill the DC component before windowing.	 
				startSample = rightSample - halfnsamp_window;
				endSample = leftSample + halfnsamp_window;
				
				if ( startSample < 1 ) {
				    std::cout<<"StartSample < 1"<<std::endl;
					std::cout<<"Sound_to_Pitch.cpp: Line 31.281"<<std::endl;
					return NULL;
				}
				
				if (endSample > my nx){
				    std::cout<<"EndSample > my nx"<<std::endl;
					std::cout<<"Sound_to_Pitch.cpp: Line 31.287"<<std::endl;
					return NULL;
				}
			
	           if (method < FCC_NORMAL) {
					for (long j = 1, i = startSample; j <= nsamp_window; j ++)
						frame [channel] [j] = (my z [channel] [i ++] - localMean [channel]) * window [j];
					for (long j = nsamp_window + 1; j <= nsampFFT; j ++)
						frame [channel] [j] = 0.0;
				} else {
					for (long j = 1, i = startSample; j <= nsamp_window; j ++)
						frame [channel] [j] = my z [channel] [i ++] - localMean [channel];
				}
			}
          
		// Compute the local peak; look half a longest period to both sides.
            localPeak = 0.0;
			if ((startSample = halfnsamp_window + 1 - halfnsamp_period) < 1) startSample = 1;
			if ((endSample = halfnsamp_window + halfnsamp_period) > nsamp_window) endSample = nsamp_window;
			for (long channel = 1; channel <= my ny; channel ++) {
				for (long j = startSample; j <= endSample; j ++) {
					double value = fabs (frame [channel] [j]);
					if (value > localPeak) localPeak = value;
				}
			}
			pitchFrame->intensity = localPeak > globalPeak ? 1.0 : localPeak / globalPeak;  		
		
			// Compute the correlation into the array 'r'.		
		if (method >= FCC_NORMAL) {
			double startTime = t - 0.5 * (1.0 / minimumPitch + dt_window);
			long localSpan = maximumLag + nsamp_window, localMaximumLag, offset;
			if ((startSample = (long) floor ((startTime - my x1) / my dx)) + 1 < 1)
				 startSample = 1;
			if (localSpan > my nx + 1 - startSample) localSpan = my nx + 1 - startSample;
			localMaximumLag = localSpan - nsamp_window;
			offset = startSample - 1;
			double sumx2 = 0;                          /* Sum of squares. */
			for (long channel = 1; channel <= my ny; channel ++) {                         ///channel = 1; channel <= my ny
				double *amp = my z[channel] + offset;
				for (long i = 1; i <= nsamp_window; i ++) {                               ///i = 1; i <= nsamp_window
					double x = amp[i] - localMean[channel]; 
					sumx2 += x * x;
				}
			}
			double sumy2 = sumx2;                      /* At zero lag, these are still equal. */
			r[0] = 1.0;
			for (long i = 1; i <= localMaximumLag; i ++) {
				double product = 0.0;
				for (long channel = 1; channel <= my ny; channel ++) {                   ///channel = 1; channel <= my ny
					double *amp = my z[channel] + offset;
					double y0 = amp[i] - localMean[channel];
					double yZ = amp[i + nsamp_window] - localMean[channel];
					sumy2 += yZ * yZ - y0 * y0;
					for (long j = 1; j <= nsamp_window; j ++) {                          ///j = 1; j <= nsamp_window
						double x = amp[j] - localMean[channel];
						double y = amp[i + j] - localMean[channel];
						product += x * y;
					}
				}
				r[- i] = r[i] = product / sqrt (sumx2 * sumy2);
			}
		} else {			
			// The FFT of the autocorrelation is the power spectrum.		
	            for (long i = 1; i <= nsampFFT; i ++) 
					ac [i] = 0.0;
				for (long channel = 1; channel <= my ny; channel ++) {
					NUMfft_forward (fftTable, frame [channel]);   /* Complex spectrum. */
					ac [1] += frame [channel] [1] * frame [channel] [1];   /* DC component. */
					for (long i = 2; i < nsampFFT; i += 2) {
						ac [i] += frame [channel] [i] * frame [channel] [i] + frame [channel] [i+1] * frame [channel] [i+1]; /* Power spectrum. */
					}
					ac [nsampFFT] += frame [channel] [nsampFFT] * frame [channel] [nsampFFT];   /* Nyquist frequency. */
				}
				NUMfft_backward (fftTable, ac);   /* Autocorrelation. */

				/*
				 * Normalize the autocorrelation to the value with zero lag,
				 * and divide it by the normalized autocorrelation of the window.
				 */
				r [0] = 1.0;
				for (long i = 1; i <= brent_ixmax; i ++)
					r [- i] = r [i] = ac [i + 1] / (ac [1] * windowR [i + 1]);
		}
			
		// Create (too much) space for candidates
		Pitch_Frame_init (pitchFrame, maxnCandidates);

	    // Register the first candidate, which is always present: voicelessness.
		pitchFrame->nCandidates = 1;
		pitchFrame->candidate[1].frequency = 0.0;    /* Voiceless: always present. */
		pitchFrame->candidate[1].strength = 0.0;

		/*
		 * Shortcut: absolute silence is always voiceless.
		 * Go to next frame.
		 */
		if (localPeak == 0) continue;

		/*
		 * Find the strongest maxima of the correlation of this frame, 
		 * and register them as candidates.
		 */
		imax[1] = 0;
		for (long i = 2; i < maximumLag && i < brent_ixmax; i ++)
		    if (r[i] > 0.5 * voicingThreshold &&       /* Not too unvoiced? */
				r[i] > r[i-1] && r[i] >= r[i+1])       /* Maximum ? */
		{
			int place = 0;
		   // Use parabolic interpolation for first estimate of frequency,and sin(x)/x interpolation to compute the strength of this frequency.
			double dr = 0.5 * (r[i+1] - r[i-1]);
			double d2r = 2 * r[i] - r[i-1] - r[i+1];
			double frequencyOfMaximum = 1 / my dx / (i + dr / d2r);
			long offset = - brent_ixmax - 1;
			double strengthOfMaximum = /* method & 1 ? */
				NUM_interpolate_sinc (& r[offset], brent_ixmax - offset, 1 / my dx / frequencyOfMaximum - offset, 30)
				/* : r [i] + 0.5 * dr * dr / d2r */;
			   /* High values due to short windows are to be reflected around 1. */
			if (strengthOfMaximum > 1.0) strengthOfMaximum = 1.0 / strengthOfMaximum;

			// Find a place for this maximum.
			if (pitchFrame->nCandidates < thy maxnCandidates) { /* Is there still a free place? */
				place = ++ pitchFrame->nCandidates;
			} else {
			   /* Try the place of the weakest candidate so far. */
				double weakest = 2;
				for (int iweak = 2; iweak <= thy maxnCandidates; iweak ++) {   //iweak = 2; iweak <= thy maxnCandidates;
					/* High frequencies are to be favoured */
					/* if we want to analyze a perfectly periodic signal correctly. */
					double localStrength = pitchFrame->candidate[iweak].strength - octaveCost *
						NUMlog2 (minimumPitch / pitchFrame->candidate[iweak].frequency);
					if (localStrength < weakest) { 
					     weakest = localStrength; 
						 place = iweak; 
				      }
				}
				/* If this maximum is weaker than the weakest candidate so far, give it no place. */
				if (strengthOfMaximum - octaveCost * NUMlog2 (minimumPitch / frequencyOfMaximum) <= weakest)
					place = 0;
			}
			if (place) {              /* Have we found a place for this candidate? */
				pitchFrame->candidate[place].frequency = frequencyOfMaximum;
				pitchFrame->candidate[place].strength = strengthOfMaximum;
				imax [place] = i;
			}
		}
		
		// Second pass: for extra precision, maximize sin(x)/x interpolation ('sinc').
		for (long i = 2; i <= pitchFrame->nCandidates; i ++) { 
			if (method != AC_HANNING || pitchFrame->candidate[i].frequency > 0.0 / my dx) {
				double xmid, ymid;
				long offset = - brent_ixmax - 1;
				ymid = NUMimproveMaximum (& r[offset], brent_ixmax - offset, imax[i] - offset,
					pitchFrame->candidate[i].frequency > 0.3 / my dx ? NUM_PEAK_INTERPOLATE_SINC700 : brent_depth, & xmid);
				xmid += offset;
				pitchFrame->candidate[i].frequency = 1.0 / my dx / xmid;
				if (ymid > 1.0) ymid = 1.0 / ymid;
				pitchFrame->candidate[i].strength = ymid;
			}
		}
	}   /* Next frame. */
	
       Pitch_pathFinder (thee, silenceThreshold, voicingThreshold,octaveCost, octaveJumpCost,
			             voicedUnvoicedCost, ceiling, false);   
					   //false:  Melder_debug == 31 ? true : false   Melder_debug 31: Pitch analysis: formant pulling on
	return thee; 
}

Pitch Sound_to_Pitch_cc (Sound me,
	double dt, double minimumPitch, double periodsPerWindow, int maxnCandidates, int accurate,
	double silenceThreshold, double voicingThreshold,
	double octaveCost, double octaveJumpCost, double voicedUnvoicedCost, double ceiling)
{
	return Sound_to_Pitch_any (me, dt, minimumPitch, periodsPerWindow, maxnCandidates, 2 + accurate,
		silenceThreshold, voicingThreshold, octaveCost, octaveJumpCost, voicedUnvoicedCost, ceiling);
}

Pitch Sound_to_Pitch(Sound me, double timestep, double minimumPitch, double maximumPitch)
{
	return Sound_to_Pitch_any(me, timestep, minimumPitch, 3.0, 15, FALSE, 0.03, 0.45, 0.01, 0.35, 0.14, maximumPitch);
}

Pitch computePitch(Sound me, double fixedTimeStepStrategy, int method,  bool veryAccurate, double floor, double ceiling,
				   long maximumNumberOfCandidates,double silenceThreshold, double voicingThreshold, double octaveCost, 
				   double octaveJumpCost, double voicedUnvoicedCost)
{
    return Sound_to_Pitch_any(me, fixedTimeStepStrategy,
	                          floor,
							  method == kTimeSoundAnalysisEditor_pitch_analysisMethod_AUTOCORRELATION ? 3.0 : 1.0,
							  maximumNumberOfCandidates,
							  (method - 1) *2 + veryAccurate, ///veryAccurate + 2,    
						      silenceThreshold, voicingThreshold, octaveCost, octaveJumpCost, 
							  voicedUnvoicedCost, ceiling);
}

/*  Reference
static void computePitch_inside (TimeSoundAnalysisEditor me) {
	double margin = my pitch.veryAccurate ? 3.0 / my pitch.floor : 1.5 / my pitch.floor;
	forget (my pitch.data);
	try {
		autoSound sound = extractSound (me, my startWindow - margin, my endWindow + margin);
		double pitchTimeStep =
			my timeStepStrategy == kTimeSoundAnalysisEditor_timeStepStrategy_FIXED ? my fixedTimeStep :
			my timeStepStrategy == kTimeSoundAnalysisEditor_timeStepStrategy_VIEW_DEPENDENT ? (my endWindow - my startWindow) / my numberOfTimeStepsPerView :
			0.0;   // the default: determined by pitch floor
		    
		 my pitch.data = Sound_to_Pitch_any (sound.peek(), pitchTimeStep,
			                                 my pitch.floor,
											 my pitch.method == kTimeSoundAnalysisEditor_pitch_analysisMethod_AUTOCORRELATION ? 3.0 : 1.0,
											 my pitch.maximumNumberOfCandidates,
											(my pitch.method - 1) * 2 + my pitch.veryAccurate,
											 my pitch.silenceThreshold, my pitch.voicingThreshold,
											 my pitch.octaveCost, my pitch.octaveJumpCost, my pitch.voicedUnvoicedCost, my pitch.ceiling);
		my pitch.data -> xmin = my startWindow;
		my pitch.data -> xmax = my endWindow;
	} catch (MelderError) {
		Melder_clearError ();
	}
}*/
