/*Pitch.cpp
 *implementation for function of computing pitch
 */

#include <iostream>
#include "Pitch.h"
#include "SoundCompute.h"
#include "NUM.h"


#define doesUnitAllowNegativeValues(unit)  \
	( (unit) == kPitch_unit_HERTZ_LOGARITHMIC || (unit) == kPitch_unit_LOG_HERTZ ||  \
	  (unit) == kPitch_unit_SEMITONES_1 || (unit) == kPitch_unit_SEMITONES_100 ||  \
	  (unit) == kPitch_unit_SEMITONES_200 || (unit) == kPitch_unit_SEMITONES_440 )

void Pitch_Frame_init (Pitch_Frame me, int nCandidate)
{
    my nCandidates = nCandidate;
	my candidate = (Pitch_Candidate)malloc((nCandidate + 1) * sizeof(structPitch_Candidate));
	
	for(int i = 1; i <= nCandidate; ++ i){
	    my candidate[i].frequency = 0.0;
		my candidate[i].strength = 0.0;
	}
}

Pitch Pitch_create(double tmin, double tmax, long nt, double dt, double t1, double ceiling, int maxCandidates){
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
   for(long i = 0; i <= nt; ++ i)
       Pitch_Frame_init(& my frame[i], 1);
	   
	return me;
}

void Pitch_setCeiling (Pitch me, double ceiling) {
	my ceiling = ceiling;
}

int Pitch_getMaxnCandidates (Pitch me) {
	int result = 0;
	for (long i = 1; i <= my nx; i ++) {
		int nCandidates = my frame[i].nCandidates;
		if (nCandidates > result) result = nCandidates;
	}
	return result;
}

void Pitch_pathFinder (Pitch me, double silenceThreshold, double voicingThreshold, double octaveCost,
	                   double octaveJumpCost, double voicedUnvoicedCost, double ceiling, int pullFormants)
{
	 long maxnCandidates = Pitch_getMaxnCandidates (me);
	 long place;
     double maximum, value;  //volatile
	 double ceiling2 = pullFormants ? 2 * ceiling : ceiling;
	 double timeStepCorrection = 0.01 / my dx;
	 octaveJumpCost *= timeStepCorrection;
	 voicedUnvoicedCost *= timeStepCorrection;

	 my ceiling = ceiling;
	 
	 double **delta = (double **) malloc ( sizeof(double *) * (my nx + 1));
	 for(long i = 1; i <= my nx; ++ i)
	     delta[i] = (double *) malloc (sizeof(double) * maxnCandidates);	
	 long **psi = (long **) malloc ( sizeof(long *) * (my nx + 1));
	 for(long i = 1; i <= my nx; ++ i)
	      psi[i] = (long *) malloc(sizeof(long) * (maxnCandidates + 1));
     /// autoNUMmatrix <double> delta (1, my nx, 1, maxnCandidates);
	 /// autoNUMmatrix <long> psi (1, my nx, 1, maxnCandidates);

	 for (long iframe = 1; iframe <= my nx; iframe ++) {
		Pitch_Frame frame = & my frame [iframe];
		double unvoicedStrength = silenceThreshold <= 0 ? 0 :
			2 - frame->intensity / (silenceThreshold / (1 + voicingThreshold));
		unvoicedStrength = voicingThreshold + (unvoicedStrength > 0 ? unvoicedStrength : 0);
		for (long icand = 1; icand <= frame->nCandidates; icand ++) {
			Pitch_Candidate candidate = & frame->candidate [icand];
			int voiceless = candidate->frequency == 0 || candidate->frequency > ceiling2;
			delta [iframe] [icand] = voiceless ? unvoicedStrength :
				candidate->strength - octaveCost * NUMlog2 (ceiling / candidate->frequency);
		}
	 }

	/* Look for the most probable path through the maxima. 
	   There is a cost for the voiced/unvoiced transition, 
	   and a cost for a frequency jump. */
	for (long iframe = 2; iframe <= my nx; iframe ++) { 
		Pitch_Frame prevFrame = & my frame[iframe - 1], curFrame = & my frame[iframe];
		double *prevDelta = delta [iframe - 1], *curDelta = delta [iframe];
		long *curPsi = psi[iframe];
		for (long icand2 = 1; icand2 <= curFrame -> nCandidates; icand2 ++) {  /// icand2 = 1; icand2 <= curFrame -> nCandidates
			double f2 = curFrame -> candidate [icand2]. frequency;
			maximum = -1e30;
			place = 0;
			for (long icand1 = 1; icand1 <= prevFrame->nCandidates; icand1 ++) {  // icand1 = 1; icand1 <= prevFrame -> nCandidates
			    double f1 = prevFrame -> candidate [icand1]. frequency;
				double transitionCost;
				bool previousVoiceless = f1 <= 0 || f1 >= ceiling2;
				bool currentVoiceless = f2 <= 0 || f2 >= ceiling2;
				if (currentVoiceless) {
					if (previousVoiceless) 
						transitionCost = 0;   // both voiceless
					else 
						transitionCost = voicedUnvoicedCost;    // voiced-to-unvoiced transition
				} else {
					if (previousVoiceless) {
						transitionCost = voicedUnvoicedCost;    // unvoiced-to-voiced transition
/******** Melder_debug 30: pitch path finder: use octave jump cost across voiceless parts		
///                     if (Melder_debug == 30) {
						// Try to take into account a frequency jump across a voiceless stretch.
							long place1 = icand1;
							for (long jframe = iframe - 2; jframe >= 1; jframe --) {
								place1 = psi [jframe + 1] [place1];
								f1 = my frame [jframe]. candidate [place1]. frequency;
								if (f1 > 0 && f1 < ceiling) {
									transitionCost += octaveJumpCost * fabs (NUMlog2 (f1 / f2)) / (iframe - jframe);
									break;
								}
							}
						}          ***********/
					} else {
						transitionCost = octaveJumpCost * fabs (NUMlog2 (f1 / f2));   // both voiced
					}
				}
				value = prevDelta [icand1] - transitionCost + curDelta [icand2];
			//if (Melder_debug == 33) Melder_casual ("Frame %ld, current candidate %ld (delta %g), previous candidate %ld (delta %g), "
		    //	"transition cost %g, value %g, maximum %g", iframe, icand2, curDelta [icand2], icand1, prevDelta [icand1], transitionCost, value, maximum);
				if (value > maximum) {
					maximum = value;
					place = icand1;
				} else if (value == maximum) {
					/*** if (Melder_debug == 33) Melder_casual ("A tie in frame %ld, current candidate %ld, previous candidate %ld", iframe, icand2, icand1); ***/
				//	std::cout<<"A tie in frame "<<iframe<<" , current candidate "<< icand2<<", previous candidate "<< icand1 <<std::endl;
				//	std::cout<<"Pitch.cpp: Line: 58."<<std::endl;
				}
			}
			curDelta[icand2] = maximum;
			curPsi[icand2] = place;
		}
 	 } 
	 
	/* Find the end of the most probable path. */
	place = 1;
	maximum = delta [my nx] [place];
	for (long icand = 2; icand <= my frame [my nx]. nCandidates; icand ++) {
		if (delta [my nx] [icand] > maximum) {
			place = icand;
			maximum = delta [my nx] [place];
			}
	}

	 /* Backtracking: follow the path backwards. */
	 for (long iframe = my nx; iframe >= 1; iframe --) {
  /****** if (Melder_debug == 33) Melder_casual ("Frame %ld: swapping candidates 1 and %ld", iframe, place); ******/
		Pitch_Frame frame = & my frame[iframe];
		structPitch_Candidate help = frame->candidate[1];
		frame->candidate[1] = frame->candidate [place];
		frame->candidate[place] = help;
		place = psi[iframe][place]; 
	}

	 /* Pull formants: devoice frames with frequencies between ceiling and ceiling2. */
	 if (ceiling2 > ceiling) {
  /******  if (Melder_debug == 33) Melder_casual ("Pulling formants..."); *******/
		 for (long iframe = my nx; iframe >= 1; iframe --) {
			 Pitch_Frame frame = & my frame[iframe];
			 Pitch_Candidate winner = & frame->candidate[1];
			 double f = winner->frequency;
			 if (f > ceiling && f <= ceiling2) {
				 for (long icand = 2; icand <= frame->nCandidates; icand ++) {
				  	  Pitch_Candidate loser = & frame->candidate[icand];
					  if (loser->frequency == 0.0) {
						  structPitch_Candidate help = *winner;
						  *winner = *loser;
						  *loser = help;
					   	  break;
					  }
				  }
			 }
		 }
	 }
}

double Pitch_getValueAtTime (Pitch me, double x, int unit, int interpolate) {
	/// return Sampled_getValueAtX (me, time, Pitch_LEVEL_FREQUENCY, unit, interpolate);
	 long ilevel = Pitch_LEVEL_FREQUENCY;
	 if (x < my xmin || x > my xmax) return NUMundefined;
	 if (interpolate) {
	    double ireal = (x - my x1)/my dx + 1;
		long ileft = floor(ireal), inear, ifar;
		double phase = ireal - ileft;
		if (phase < 0.5) {
			inear = ileft, ifar = ileft + 1;
		} else {
			ifar = ileft, inear = ileft + 1;
			phase = 1.0 - phase;
		}
		if (inear < 1 || inear > my nx) return NUMundefined;   // x out of range?
		double fnear = v_getValueAtSample (me, inear, ilevel, unit);
		if (fnear == NUMundefined) return NUMundefined;   // function value not defined?
		if (ifar < 1 || ifar > my nx) return fnear;   // at edge? Extrapolate
		double ffar = v_getValueAtSample (me, ifar, ilevel, unit);
		if (ffar == NUMundefined) return fnear;   // neighbour undefined? Extrapolate
		return fnear + phase * (ffar - fnear);   // interpolate
	 }
	 
	 long isample = (long)((x - my x1) / my dx + 1.5);
	 return Sampled_getValueAtSample (me, isample, ilevel, unit);
}

/*     Reference:
double Sampled_getValueAtX (Sampled me, double x, long ilevel, int unit, int interpolate) {
	if (x < my xmin || x > my xmax) return NUMundefined;
	if (interpolate) {
		double ireal = Sampled_xToIndex (me, x);
		long ileft = floor (ireal), inear, ifar;
		double phase = ireal - ileft;
		if (phase < 0.5) {
			inear = ileft, ifar = ileft + 1;
		} else {
			ifar = ileft, inear = ileft + 1;
			phase = 1.0 - phase;
		}                                                     //±£Ö¤phase´óÐ¡ÔÚ[0.0, 0.5]Ö®¼ä
		if (inear < 1 || inear > my nx) return NUMundefined;   // x out of range?      x³¬³ö·¶Î§
		double fnear = my v_getValueAtSample (inear, ilevel, unit);
		if (fnear == NUMundefined) return NUMundefined;   // function value not defined?     //º¯ÊýÖµÎ´¶¨Òå
		if (ifar < 1 || ifar > my nx) return fnear;   // at edge? Extrapolate    //ÔÚ±ß½çÍâ²å
		double ffar = my v_getValueAtSample (ifar, ilevel, unit);
		if (ffar == NUMundefined) return fnear;   // neighbour undefined? Extrapolate     
		return fnear + phase * (ffar - fnear);   // interpolate                  //ÄÚ²åÖµ
	}
	return Sampled_getValueAtSample (me, Sampled_xToNearestIndex (me, x), ilevel, unit);
}*/

bool Pitch_isVoiced_i (Pitch me, long iframe) {
	return NUMdefined (Sampled_getValueAtSample (me, iframe, Pitch_LEVEL_FREQUENCY, kPitch_unit_HERTZ));
}

double Sampled_getValueAtSample (Pitch me, long isamp, long ilevel, int unit) {  /***(Sampled me, long isamp, long ilevel, int unit)***/
	if (isamp < 1 || isamp > my nx) return NUMundefined;
	return v_getValueAtSample (me, isamp, ilevel, unit);
}

//Sampled.cpp  110
double Sampled_getValueAtX (Pitch me, double x, long ilevel, int unit, int interpolate){
	if (x < my xmin || x > my xmax) 
		return NUMundefined;

	if (interpolate) {
		double ireal = (x - my x1) / my dx + 1;   //Sampled_xToIndex (me, x);
		long ileft = floor (ireal), inear, ifar;
		double phase = ireal - ileft;
		if (phase < 0.5) {
			inear = ileft, ifar = ileft + 1;
		} else {
			ifar = ileft, inear = ileft + 1;
			phase = 1.0 - phase;
		}
		if (inear < 1 || inear > my nx)
			return NUMundefined;   // x out of range?

		double fnear = v_getValueAtSample (me, inear, ilevel, unit);

		if (fnear == NUMundefined) 
			return NUMundefined;   // function value not defined?\

		if (ifar < 1 || ifar > my nx) 
			return fnear;   // at edge? Extrapolate

		double ffar = v_getValueAtSample (me, ifar, ilevel, unit);

		if (ffar == NUMundefined) 
			return fnear;   // neighbour undefined? Extrapolate

		return fnear + phase * (ffar - fnear);   // interpolate
	}
	long isamp = (long)floor((x - my x1)/my dx + 1.5);  //Sampled_xToNearestIndex (me, x)
	return Sampled_getValueAtSample (me, isamp, ilevel, unit);
}

//Sampled.cpp  133
long Sampled_countDefinedSamples (Pitch me, long ilevel, int unit) {
	long numberOfDefinedSamples = 0;
	for (long isamp = 1; isamp <= my nx; isamp ++) {
		double value = v_getValueAtSample (me, isamp, ilevel, unit);

		if (value == NUMundefined) 
			continue;
		numberOfDefinedSamples += 1;
	}
	return numberOfDefinedSamples;
}

double v_getValueAtSample (Pitch me, long iframe, long ilevel, int unit) {
	double f = my frame[iframe].candidate[1].frequency;
	if (f <= 0.0 || f >= my ceiling) return NUMundefined;   // frequency out of range (or NUMundefined)? Voiceless
	return v_convertStandardToSpecialUnit (ilevel == Pitch_LEVEL_FREQUENCY ? f : my frame[iframe].candidate[1].strength, ilevel, unit);}

double v_convertStandardToSpecialUnit (double value, long ilevel, int unit) {   /****structPitch :: v_convertStandardToSpecialUnit (double value, long ilevel, int unit)****/
	if (ilevel == Pitch_LEVEL_FREQUENCY) {
		return
			unit == kPitch_unit_HERTZ ? value :
			unit == kPitch_unit_HERTZ_LOGARITHMIC ? value <= 0.0 ? NUMundefined : log10 (value) :
			unit == kPitch_unit_MEL ? NUMhertzToMel (value) :
			unit == kPitch_unit_LOG_HERTZ ? value <= 0.0 ? NUMundefined : log10 (value) :
			unit == kPitch_unit_SEMITONES_1 ? value <= 0.0 ? NUMundefined : 12.0 * log (value / 1.0) / NUMln2 :
			unit == kPitch_unit_SEMITONES_100 ? value <= 0.0 ? NUMundefined : 12.0 * log (value / 100.0) / NUMln2 :
			unit == kPitch_unit_SEMITONES_200 ? value <= 0.0 ? NUMundefined : 12.0 * log (value / 200.0) / NUMln2 :
			unit == kPitch_unit_SEMITONES_440 ? value <= 0.0 ? NUMundefined : 12.0 * log (value / 440.0) / NUMln2 :
			unit == kPitch_unit_ERB ? NUMhertzToErb (value) :
			NUMundefined;
	} else {
		return
			unit == Pitch_STRENGTH_UNIT_AUTOCORRELATION ? value :
			unit == Pitch_STRENGTH_UNIT_NOISE_HARMONICS_RATIO ?
				value <= 1e-15 ? 1e15 : value > 1.0 - 1e-15 ? 1e-15 : (1.0 - value) / value :   /* Before losing precision. */
			unit == Pitch_STRENGTH_UNIT_HARMONICS_NOISE_DB ?
				value <= 1e-15 ? -150.0 : value > 1.0 - 1e-15 ? 150.0 : 10 * log10 (value / (1.0 - value)) :   /* Before losing precision. */
			NUMundefined;
	}
}

void Pitch_scaleDuration(Pitch me, double multiplier)
{
	if(multiplier != 1){    // keep xmin at the same value
		my dx *= multiplier;
		my x1 = my xmin + (my x1 - my xmin) * multiplier;
		my xmax = my xmin + (my xmax - my xmin) * multiplier;
	}
}

void Pitch_scalePitch(Pitch me, double multiplier)
{
	for(long i = 1; i < my nx; ++ i){
		double f = my frame[i].candidate[1].frequency;
		f *= multiplier;
		if(f < my ceiling)
			my frame[i].candidate[1].frequency = f;
	}
}

bool intersectRangeWithDomain(Pitch me, double *x1, double *x2)
{
	if(*x1 == *x2)  return false;
	if(*x1 < *x2){
	   if(*x1 < my xmin)  *x1 = my xmin;
	   if(*x2 > my xmax)  *x2 = my xmax;
	   if(*x2 <= *x1)  return false;
	}
	else 
	{
	   if(*x2 < my xmin)  *x1 = my xmin;
	   if(*x1 > my xmax)  *x2 = my xmax;
	   if(*x1 <= *x2)  return false;
	}
	return true;
}


long Sampled_getWindowSamples (Pitch me, double xmin, double xmax, long *ixmin, long *ixmax) {
	double rixmin = 1.0 + ceil ((xmin - my x1) / my dx);
	double rixmax = 1.0 + floor ((xmax - my x1) / my dx);
	*ixmin = rixmin < 1.0 ? 1 : (long) rixmin;
	*ixmax = rixmax > (double) my nx ? my nx : (long) rixmax;
	if (*ixmin > *ixmax) return 0;
	return *ixmax - *ixmin + 1;
}

double getQuantitle(Pitch me, double xmin, double xmax, double quantile, long ilevel, int unit)
{   //PraatÔ´ÂëÖÐ£º  Sampled.cpp 157  
	//double Sampled_getQuantile (Sampled me, double xmin, double xmax, double quantile, long ilevel, int unit) 

	double *values = (double *)malloc(sizeof(double) * (my nx + 1));
	if(xmin >= xmax){
		xmin = my xmin;
		xmax = my xmax;
	}

	if(!(intersectRangeWithDomain(me, & xmin, & xmax)))
		return NUMundefined;

	long imin, imax, numberOfDefinedSamples = 0;
	Sampled_getWindowSamples(me, xmin, xmax, & imin, & imax);

	for(long i = imin ; i < imax; ++ i){
	   double value = v_getValueAtSample(me, i, ilevel, unit);
	   if(NUMdefined(value))
		   values[++ numberOfDefinedSamples] = value;
	}

	double result = NUMundefined;
	if(numberOfDefinedSamples >= 1){
	    NUMsort_d(numberOfDefinedSamples, values);
		result = NUMquantile (numberOfDefinedSamples, values, quantile);
	}
	return result;
}

//Pitch.cpp 179
double Pitch_getQuantile (Pitch me, double tmin, double tmax, double quantile, int unit)
{
	double value = getQuantitle(me, tmin, tmax, quantile, Pitch_LEVEL_FREQUENCY, unit);
	//PraatÔ´ÂëÖÐ: Sampled_getQuantile
	if(value <= 0.0 && !doesUnitAllowNegativeValues(unit))
		value = NUMundefined;
	return value;
}

//Sound_extension.cpp 1501
Pitch Pitch_scaleTime_old (Pitch me, double scaleFactor){
	double dx = my dx, x1 = my x1, xmax = my xmax;
	if (scaleFactor != 1) {
		dx = my dx * scaleFactor;
		x1 = my xmin + 0.5 * dx;
		xmax = my xmin + my nx * dx;
	}
	
	Pitch thee = Pitch_create (my xmin, xmax, my nx, dx, x1, my ceiling, 2);
	for (long i = 1; i <= my nx; i++) {
		double f = my frame[i].candidate[1].frequency;
		thy frame[i].candidate[1].strength = my frame[i].candidate[1].strength;
		f /= scaleFactor;
		if (f < my ceiling)
			thy frame[i].candidate[1].frequency = f;
	}

	return thee;
}

//Pitch.cpp 227
double Pitch_getMinimum (Pitch me, double tmin, double tmax, int unit, int interpolate){
	double minimum;
	Pitch_getMinimumAndTime(me, tmin, tmax, unit, interpolate, & minimum, NULL);
	return minimum;
}

//Pitch.cpp 217
void Pitch_getMinimumAndTime (Pitch me, double tmin, double tmax, int unit, int interpolate,
	                          double *return_minimum, double *return_timeOfMinimum){

	Sampled_getMinimumAndX (me, tmin, tmax, Pitch_LEVEL_FREQUENCY, unit, interpolate, return_minimum, return_timeOfMinimum);
	if (!doesUnitAllowNegativeValues(unit) && return_minimum && *return_minimum <= 0.0)
	      *return_minimum = NUMundefined;     /* Not so unlikely. */
}

void Sampled_getMinimumAndX (Pitch me, double xmin, double xmax, long ilevel, int unit, int interpolate,
	double *return_minimum, double *return_xOfMinimum){
	long imin, imax, i;
	double minimum = 1e301, xOfMinimum = 0.0;
	if (xmin == NUMundefined || xmax == NUMundefined) {
		minimum = xOfMinimum = NUMundefined;
		goto end;
	}

	
    if (xmin >= xmax) {
		xmin = my xmin;
		xmax = my xmax;
	}//Ô´ÂëÖÐÎª£º Function_unidirectionalAutowindow (me, & xmin, & xmax);

	if (!IntersectRangeWithDomain(me, &xmin, &xmax) ) { //Ô´ÂëÖÐ£º Function_intersectRangeWithDomain (me, & xmin, & xmax)
		minimum = xOfMinimum = NUMundefined;   //requested range and logical domain do not intersect
		goto end;
	}
	if (! Sampled_getWindowSamples (me, xmin, xmax, & imin, & imax)) {
		/*
		 * No sample centres between xmin and xmax.
		 * Try to return the lesser of the values at these two points.
		 */
		double fleft =  Sampled_getValueAtX (me, xmin, ilevel, unit, interpolate);
		double fright = Sampled_getValueAtX (me, xmax, ilevel, unit, interpolate);
		if (NUMdefined (fleft) && fleft < minimum) minimum = fleft, xOfMinimum = xmin;
		if (NUMdefined (fright) && fright < minimum) minimum = fright, xOfMinimum = xmax;
	} else {
		for (i = imin; i <= imax; i ++) {
			double fmid = v_getValueAtSample (me, i, ilevel, unit);
			if (fmid == NUMundefined)
				continue;

			if (interpolate == FALSE) {
				if (fmid < minimum) 
					minimum = fmid, xOfMinimum = i;
			} else {
	     //Try an interpolation, possibly even taking into account a sample just outside the selection.
				double fleft = i <= 1 ? NUMundefined : v_getValueAtSample (me, i - 1, ilevel, unit);
				double fright = i >= my nx ? NUMundefined : v_getValueAtSample (me, i + 1, ilevel, unit);
				if (fleft == NUMundefined || fright == NUMundefined) {
					if (fmid < minimum)  
						 minimum = fmid, xOfMinimum = i;
				} else if (fmid < fleft && fmid <= fright) {
					double y [4], i_real, localMinimum;
					y [1] = fleft, y [2] = fmid, y [3] = fright;
					localMinimum = NUMimproveMinimum (y, 3, 2, NUM_PEAK_INTERPOLATE_PARABOLIC, & i_real);
					if (localMinimum < minimum)
						minimum = localMinimum, xOfMinimum = i_real + i - 2;
				}
			}
		}
		xOfMinimum = my x1 + (xOfMinimum - 1) * my dx;   /* From index plus phase to time. */
		/* Check boundary values. */
		if (interpolate) {
			double fleft = Sampled_getValueAtX (me, xmin, ilevel, unit, TRUE);
			double fright = Sampled_getValueAtX (me, xmax, ilevel, unit, TRUE);
			if (NUMdefined (fleft) && fleft < minimum) 
				minimum = fleft, xOfMinimum = xmin;
			if (NUMdefined (fright) && fright < minimum) 
				minimum = fright, xOfMinimum = xmax;
		}
		if (xOfMinimum < xmin) 
			xOfMinimum = xmin;
		if (xOfMinimum > xmax) 
			xOfMinimum = xmax;
	}
	if (minimum == 1e301)
		minimum = xOfMinimum = NUMundefined;

end:
	if (return_minimum)
		 *return_minimum = minimum;
	if (return_xOfMinimum) 
		 *return_xOfMinimum = xOfMinimum;
}

//Function.cpp 165
bool IntersectRangeWithDomain (Pitch me, double *x1, double *x2){
	//Ô´ÂëÖÐ  bool Function_intersectRangeWithDomain (Function me, double *x1, double *x2)  
	if (*x1 == *x2) 
		return false;

	if (*x1 < *x2) {
		if (*x1 < my xmin)
			*x1 = my xmin;   // intersect requested range with logical domain
		if (*x2 > my xmax)
			*x2 = my xmax;
		if (*x2 <= *x1) 
			return false;   // requested range and logical domain do not intersect
	} else {
		if (*x2 < my xmin) 
			*x2 = my xmin;   // intersect requested range with logical domain
		if (*x1 > my xmax) 
			*x1 = my xmax;
		if (*x1 <= *x2) 
			return false;   // requested range and logical domain do not intersect
	}

	return true;
}
