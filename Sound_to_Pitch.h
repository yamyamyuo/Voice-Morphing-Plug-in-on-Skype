#ifndef _SOUND_TO_PITCH_H_
#define _SOUND_TO_PITCH_H_

#include "structure.h"

int Sampled_shortTermAnalysis (Sound me, double windowDuration, double timeStep, long *numberOfFrames, double *firstTime);

Pitch Sound_to_Pitch_any (Sound me,
	double dt,                 /* time step (seconds); 0.0 = automatic = periodsPerWindow / minimumPitch / 4 */
	double minimumPitch,       /* (Hz) */
	double periodsPerWindow,   /* ac3 for pitch analysis, 6 or 4.5 for HNR, 1 for FCC */
	int maxnCandidates,        /* maximum number of candidates per frame */
	int method,                /* 0 or 1 = AC, 2 or 3 = FCC, 0 or 2 = fast, 1 or 3 = accurate */

	double silenceThreshold,   /* relative to purely periodic; default 0.03 */
	double voicingThreshold,   /* relative to purely periodic; default 0.45 */
	double octaveCost,         /* favours higher pitches; default 0.01 */
	double octaveJumpCost,     /* default 0.35 */
	double voicedUnvoicedCost, /* default 0.14 */
	double maximumPitch);      /* (Hz) */
/*
	Function:
		acoustic periodicity analysis.    //ÉùÖÜÆÚ·ÖÎö
	Preconditions:
		minimumPitch > 0.0;
		maxnCandidates >= 2;
	Return value:
		the resulting pitch contour, or NULL in case of failure.    //ÓÉ´Ë²úÉúµÄ»ùÒôÂÖÀª
	Failures:
		Out of memory.
		Minimum frequency too low.
		Maximum frequency should not be greater than the Sound's Nyquist frequency.
	Description for method 0 or 1:
		There is a Hanning window (method == 0) or Gaussian window (method == 1)
		over the analysis window, in order to avoid phase effects.
		Zeroes are appended to the analysis window to avoid edge effects in the FFT.
		An FFT is done on the window, giving a complex spectrum.                    //´°ÉÏµÄFFT£¬¸ø³öÒ»¸´ºÏÆ×
		This complex spectrum is squared, thus giving the power spectrum.           //¸´ºÏÆ×Æ½·½µÃµ½¹¦ÂÊÆ×
		The power spectrum is FFTed back, thus giving the autocorrelation of         //¹¦ÂÊÆ×µÄ·´FFT£¬µÃµ½´°Ö¡µÄ×ÔÏà¹Ø
		the windowed frame. This autocorrelation is expressed relative to the power, 
		which is the autocorrelation for lag 0. The autocorrelation is divided by
		the normalized autocorrelation of the window, in order to bring
		all maxima of the autocorrelation of a periodic signal to the same height.
	General description:
		The maxima are found by sinc interpolation.
		The pitch values (frequencies) of the highest 'maxnCandidates' maxima
		are saved in 'result', together with their strengths (relative correlations).
		A Viterbi algorithm is used to find a smooth path through the candidates,
		using the last six arguments of this function.

		The 'maximumPitch' argument has no influence on the search for candidates.
		It is directly copied into the Pitch object as a hint for considering
		pitches above a certain value "voiceless".
*/

Pitch Sound_to_Pitch_cc (Sound me,
	                     double dt, 
						 double minimumPitch, 
						 double periodsPerWindow, 
						 int maxnCandidates, 
						 int accurate,  

						 double silenceThreshold, 
						 double voicingThreshold,
	                     double octaveCost, 
						 double octaveJumpCost, 
						 double voicedUnvoicedCost, 
						 double ceiling);
 
Pitch Sound_to_Pitch(Sound me, double timestep, double minimumPitch, double maximumPitch);

Pitch computePitch(Sound me, double fixedTimeStepStrategy = 0.0,  //auto
	               int method = 1,    /* = AUTOCORRELATION*/     /*kTimeSoundAnalysisEditor_pitch_analysisMethod*/
				   bool veryAccurate = false,				   
				   double floor = 75, double ceiling = 600, /*Pitch settings*/ 
				   long maximumNumberOfCandidates = 15,  /*default value,  >= 2*/
                   double silenceThreshold = 0.03, double voicingThreshold = 0.45, double octaveCost = 0.01, 
				   double octaveJumpCost = 0.35, double voicedUnvoicedCost = 0.14);


#endif
