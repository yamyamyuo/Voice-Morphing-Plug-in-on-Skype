// Pitch_to_PointProcess.cpp
/*
 *implement the function for pith to pointProcess defined in Structure.cpp
 */

#include "Pitch_to_PointProcess.h"
#include <iostream>
#include "Structure.h"
#include "SoundCompute.h"
#include "Pitch.h"
#include "NUM.h"

static double findExtremum_3 (double *channel1_base, double *channel2_base, long d, long n, int includeMaxima, int includeMinima) {
	double *channel1 = channel1_base + d, *channel2 = channel2_base ? channel2_base + d: NULL;
	int includeAll = includeMaxima == includeMinima;
	long imin = 1, imax = 1, i, iextr;
	double minimum, maximum;
	if (n < 3) {
		if (n <= 0) return 0.0;   /* Outside. */
		else if (n == 1) return 1.0;
		else {    /* n == 2 */
			double x1 = channel2 ? 0.5 * (channel1[0] + channel2[0]) : channel1[0];
			double x2 = channel2 ? 0.5 * (channel1[1] + channel2[1]) : channel1[1];
			double xleft = includeAll ? fabs (x1) : includeMaxima ? x1 : - x1;
			double xright = includeAll ? fabs (x2) : includeMaxima ? x2 : - x2;
			if (xleft > xright) return 1.0;
			else if (xleft < xright) return 2.0;
			else return 1.5;
		}
	}
	minimum = maximum = channel2 ? 0.5 * (channel1[1] + channel2[1]) : channel1[1];
	for (i = 2; i <= n; i ++) {
		double value = channel2 ? 0.5 * (channel1[i] + channel2[i]) : channel1[i];
		if (value < minimum) { minimum = value; imin = i; }
		if (value > maximum) { maximum = value; imax = i; }
	}
	if (minimum == maximum) {
		return 0.5 * (n + 1.0);   /* All equal. */
	}
	iextr = includeAll ? ( fabs (minimum) > fabs (maximum) ? imin : imax ) : includeMaxima ? imax : imin;
	if (iextr == 1) return 1.0;
	if (iextr == n) return (double) n;
	/* Parabolic interpolation. */
	/* We do NOT need fabs here: we look for a genuine extremum. */
	double valueMid = channel2 ? 0.5 * (channel1[iextr] + channel2[iextr]) : channel1[iextr];
	double valueLeft = channel2 ? 0.5 * (channel1[iextr - 1] + channel2[iextr - 1]) : channel1[iextr - 1];
	double valueRight = channel2 ? 0.5 * (channel1[iextr + 1] + channel2[iextr + 1]) : channel1[iextr + 1];
	return iextr + 0.5 * (valueRight - valueLeft) / (2 * valueMid - valueLeft - valueRight);
}

static double Sound_findExtremum (Sound me, double tmin, double tmax, int includeMaxima, int includeMinima) {
	long imin = (long) floor ((tmin - my x1) / my dx) + 1; /** Sampled_xToLowIndex (me, tmin); **/
	long imax = (long) ceil ((tmin - my x1) / my dx) + 1; /** Sampled_xToHighIndex (me, tmax); **/
	if(! NUMdefined (tmin) || ! NUMdefined (tmax)) {
	    std::cout<<"tmin = "<<tmin<<", tmax = "<<tmax<<". \n"<<"Error, Pitch_to_PointProcess.cpp: Line 48"<<std::endl;
		exit(0);
	}
 /**** Melder_assert (NUMdefined (tmin));
	Melder_assert (NUMdefined (tmax)); ***/
	if (imin < 1) imin = 1;
	if (imax > my nx) imax = my nx;
	double iextremum = findExtremum_3 (my z[1], my ny > 1 ? my z[2] : NULL, imin - 1, imax - imin + 1, includeMaxima, includeMinima);
	if (iextremum)
		return my x1 + (imin - 1 + iextremum - 1) * my dx;
	else
		return (tmin + tmax) / 2;
}

static double Sound_findMaximumCorrelation (Sound me, double t1, double windowLength, double tmin2, 
											double tmax2, double *tout, double *peak) {
	double maximumCorrelation = -1.0, r1 = 0.0, r2 = 0.0, r3 = 0.0, r1_best, r3_best, ir;
	double halfWindowLength = 0.5 * windowLength;
	long ileft1 = (long) floor ((t1 - halfWindowLength - my x1) / my dx + 1.5);  
	long iright1 = (long) floor ((t1 + halfWindowLength - my x1) / my dx + 1.5);  
	long ileft2min = (long) floor ((tmin2 - halfWindowLength - my x1) / my dx) + 1; 
	long ileft2max = (long) ceil ((tmax2 - halfWindowLength - my x1) / my dx) + 1; 
	*peak = 0.0;        /* Default. */
	for (long ileft2 = ileft2min; ileft2 <= ileft2max; ileft2 ++) {
		double norm1 = 0.0, norm2 = 0.0, product = 0.0, localPeak = 0.0;
		if (my ny == 1) {
			for (long i1 = ileft1, i2 = ileft2; i1 <= iright1; i1 ++, i2 ++) {
				if (i1 < 1 || i1 > my nx || i2 < 1 || i2 > my nx) continue;
				double amp1 = my z[1][i1], amp2 = my z[1][i2];
				norm1 += amp1 * amp1;
				norm2 += amp2 * amp2;
				product += amp1 * amp2;
				if (fabs (amp2) > localPeak)
					localPeak = fabs (amp2);
			}
		} else {
			for (long i1 = ileft1, i2 = ileft2; i1 <= iright1; i1 ++, i2 ++) {
				if (i1 < 1 || i1 > my nx || i2 < 1 || i2 > my nx) continue;
				double amp1 = 0.5 * (my z[1][i1] + my z[2][i1]), amp2 = 0.5 * (my z[1][i2] + my z[2][i2]);
				norm1 += amp1 * amp1;
				norm2 += amp2 * amp2;
				product += amp1 * amp2;
				if (fabs (amp2) > localPeak)
					localPeak = fabs (amp2);
			}
		}
		r1 = r2;
		r2 = r3;
		r3 = product ? product / (sqrt (norm1 * norm2)) : 0.0;
		if (r2 > maximumCorrelation && r2 >= r1 && r2 >= r3) {
			r1_best = r1;
			maximumCorrelation = r2;
			r3_best = r3;
			ir = ileft2 - 1;
			*peak = localPeak;  
		}
	}
	
	// Improve the result by means of parabolic interpolation.
	if (maximumCorrelation > -1.0) {
		double d2r = 2 * maximumCorrelation - r1_best - r3_best;
		if (d2r != 0.0) {
			double dr = 0.5 * (r3_best - r1_best);
			maximumCorrelation += 0.5 * dr * dr / d2r;
			ir += dr / d2r;
		}
		*tout = t1 + (ir - ileft1) * my dx;
	}
	return maximumCorrelation;
}

long PointProcess_getLowIndex (PointProcess me, double t) {
	if (my nt == 0 || t < my t [1])
		return 0;
	if (t >= my t[my nt])   /* Special case that often occurs in practice. */
		return my nt;
	if(my nt == 1){
	    std::cout<<"my nt = 1.  Pitch_to_PointProcess.cpp Line 10"<<std::endl;
		std::cout<<"Pitch_to_PointProcess.cpp: Line 121."<<std::endl;
	    exit(0);
	}
	/**** Melder_assert (my nt != 1); ***/   /* May fail if t or my t [1] is NaN. */ 
	
	/* Start binary search. */
	long left = 1, right = my nt;
	while (left < right - 1) {
		long mid = (left + right) / 2;
		if (t >= my t[mid]) left = mid; 
		else right = mid;
	}
	/**** !!!???? Melder_assert (right == left + 1); ****/	
	return left;
}

void PointProcess_addPoint (PointProcess me, double t) {
	 if (t == NUMundefined){
		std::cout<<"Cannot add a point at an undefined time. Pitch_to_PointProcess.cpp: Line9."<<std::endl;
		exit(0);
	 }
	/**** Melder_throw ("Cannot add a point at an undefined time."); ****/
	 if (my nt >= my maxnt) {
	    // Create without change.
		 double *dum = (double *)malloc(sizeof(double) * (2 * my nt + 1));
	     /**** autoNUMvector <double> dum (1, 2 * my maxnt); ****/
		 memset(dum, 0, sizeof(double) * (2 * my nt + 1));
		 for(long i = 1; i <= my nt; ++ i)
			 dum[i] = my t[i];
		// memcpy(dum + sizeof(double), my t + sizeof(double), my nt * sizeof(double));
	     /**** NUMvector_copyElements (my t, dum.peek(), 1, my nt); ****/

  		//Change without error.
	  /* NUMvector_free (my t, 1);
		 my t = dum.transfer(); */
		 my t = dum;
		 my maxnt = 2 * my nt + 1;
	 }
	 if (my nt == 0 || t >= my t[my nt]) {   // special case that often occurs in practice
		 my t[++ my nt] = t;
	 } else {
		 long left = PointProcess_getLowIndex (me, t);
		 if (left == 0 || my t[left] != t) {
			 for (long i = my nt; i > left; i --) my t[i + 1] = my t[i];
		 	 my nt ++;
			 my t[left + 1] = t;
		 }
	 }
}
 
 void PointProcess_init(I, double tmin, double tmax, long initialMaxnt) {
	iam (PointProcess);
	my xmin = tmin;
	my xmax = tmax;
	if (initialMaxnt < 1) initialMaxnt = 1;
	my maxnt = initialMaxnt;
	my nt = 0;
	my t = (double *) malloc ( sizeof(double) * (my maxnt + 1));
	if(my t == NULL)
		return ;
	for(long i = 0; i <= my maxnt; ++ i)
		my t[i] = 0.0;
}
 
 PointProcess PointProcess_create(double tmin, double tmax, long initialMaxnt) {
    PointProcess me = (PointProcess)malloc(sizeof(structPointProcess));
	PointProcess_init(me, tmin, tmax, initialMaxnt);
	return me;
 } 
 
int Pitch_getVoicedIntervalAfter (Pitch me, double after, double *tleft, double *tright) {
	long ileft = (long)ceil((after - my x1) / my dx) + 1; /// Sampled_xToHighIndex (me, after);
	long iright = 0; 
	if (ileft > my nx) return 0;    /* Offright. */
	if (ileft < 1) ileft = 1;       /* Offleft. */

	/* Search for first voiced frame. */
	for (; ileft <= my nx; ileft ++)
		if (Pitch_isVoiced_i (me, ileft)) break;
	if (ileft > my nx) return 0;   /* Offright. */

	/* Search for last voiced frame. */
	for (iright = ileft; iright <= my nx; iright ++)
		if (! Pitch_isVoiced_i (me, iright)) break;
	iright --;

	*tleft = my x1 + (ileft - 1) * my dx - 0.5 * my dx;    /* The whole frame is considered voiced. */
	/**Sampled_indexToX (me, ileft) - 0.5 * my dx; ***/

	*tright =  my x1 + (iright - 1) * my dx + 0.5 * my dx;   /**** Sampled_indexToX (me, iright) + 0.5 * my dx; ****/
	if (*tleft >= my xmax - 0.5 * my dx) return 0;
	if (*tleft < my xmin) *tleft = my xmin;
	if (*tright > my xmax) *tright = my xmax;
	return 1;
}

 PointProcess Sound_Pitch_to_PointProcess_cc (Sound sound, Pitch pitch){
     PointProcess point = PointProcess_create (sound->xmin, sound->xmax, 10); 
	 double t = pitch->xmin;
	 double addedRight = -1e300;
	 double globalPeak = Sound_getAbsoluteExtremum(sound, sound->xmin, sound->xmax, 0), peak;
	 //¼ÆËãµÃµ½È«¾ÖµÄ¾ø¶ÔÖµ·åÖµ	 
	       /**globalPeak = Vector_getAbsoluteExtremum (sound, sound -> xmin, sound -> xmax, 0), peak;**/
	
	 // Cycle over all voiced intervals.
	/**** autoMelderProgress progress (L"Sound & Pitch: To PointProcess...");   ****/
	 for (;;) {
		 double tleft, tright;
		 if (! Pitch_getVoicedIntervalAfter (pitch, t, & tleft, & tright))
    		 break;
		// Go to the middle of the voice stretch. 
		 double tmiddle = (tleft + tright) / 2;
	/**** Melder_progress ((tmiddle - sound->xmin) / (sound->xmax - sound->xmin), L"Sound & Pitch to PointProcess"); ****/
		 double f0middle = Pitch_getValueAtTime (pitch, tmiddle, kPitch_unit_HERTZ, Pitch_LINEAR);     
		                                      /// #define Pitch_LINEAR 1   #define kPitch_unit_HERTZ 0
		
		// Our first point is near this middle.
		 if (f0middle == NUMundefined) {
		     std::cout<<"Sound_Pitch_to_PointProcess_cc: tleft: "<< tleft <<", tright:"<< tright 
					  <<", f0middle: "<< f0middle <<std::endl; 
			 std::cout<<"Pitch_to_PointProcess.cpp: Line 215"<<std::endl;
			 return NULL;
		  /****	 Melder_fatal ("Sound_Pitch_to_PointProcess_cc: tleft %ls, tright %ls, f0middle %ls",
			                Melder_double (tleft), Melder_double (tright), Melder_double (f0middle));  ****/
		 }
		 double tmax = Sound_findExtremum (sound, tmiddle - 0.5 / f0middle, tmiddle + 0.5 / f0middle, TRUE, TRUE);
		 //ÓïÒôÊý¾ÝÖÐµÄ¼«´óÖµ

		 if(!NUMdefined(tmax)){
		     std::cout<<"tmax is UnDefined!"<<std::endl;
			 std::cout<<tmax<<std::endl;
			 std::cout<<"Pitch_to_PointProcess.cpp: Line 215"<<std::endl;
			 return NULL;
		 }
	   /****	 Melder_assert (NUMdefined (tmax));    ****/
		 PointProcess_addPoint (point, tmax); 
  
		 double tsave = tmax;
		 for (;;) {
			 double f0 = Pitch_getValueAtTime(pitch, tmax, kPitch_unit_HERTZ, Pitch_LINEAR), correlation;
			 if (f0 == NUMundefined) break;
			 correlation = Sound_findMaximumCorrelation (sound, tmax, 1.0 / f0, tmax - 1.25 / f0, 
														 tmax - 0.8 / f0, &tmax, &peak);
			//¼ÆËãµÃµ½×î´óÏà¹ØÏµÊý£¬Çø¼äÄÚµÄ×î´óÖµ&tmaxÒÔ¼°¾Ö²¿·åÖµ&peak


			 if (correlation == -1) /*break*/ tmax -= 1.0 / f0;   /* This one period will drop out. */
			 //Ñ­»·ÖÕÖ¹Ìõ¼þ
			 if (tmax < tleft) {
				 if (correlation > 0.7 && peak > 0.023333 * globalPeak && tmax - addedRight > 0.8 / f0) {
			 		 PointProcess_addPoint (point, tmax);
				 }
				 break;
			 }
			 if (correlation > 0.3 && (peak == 0.0 || peak > 0.01 * globalPeak)) {
				 if (tmax - addedRight > 0.8 / f0) {  // do not fill in a short originally unvoiced interval twice
					 PointProcess_addPoint (point, tmax); 
				}
			 }
		 }
		 tmax = tsave;
		 for (;;) {
			 double f0 = Pitch_getValueAtTime (pitch, tmax, kPitch_unit_HERTZ, Pitch_LINEAR), correlation;
			 if (f0 == NUMundefined) break;
			 correlation = Sound_findMaximumCorrelation (sound, tmax, 1.0 / f0, tmax + 0.8 / f0, 
														 tmax + 1.25 / f0, &tmax, &peak);
			//¼ÆËãµÃµ½×î´óÏà¹ØÏµÊý£¬Çø¼äÄÚµÄ×î´óÖµÒÔ¼°¾Ö²¿·åÖµ

			 if (correlation == -1) /*break*/ tmax += 1.0 / f0;
			 //Ñ­»·ÖÕÖ¹Ìõ¼þ
			 if (tmax > tright) {
				 if (correlation > 0.7 && peak > 0.023333 * globalPeak) {
					 PointProcess_addPoint (point, tmax);
					 addedRight = tmax;
				 }
				 break;
			 }
			 if (correlation > 0.3 && (peak == 0.0 || peak > 0.01 * globalPeak)) {
				 PointProcess_addPoint (point, tmax); 
				 addedRight = tmax;
			 }
		 }
		 t = tright;
	 }
	 return point;
 }

 long PointProcess_getNearestIndex (PointProcess me, double t){
    if (my nt == 0)
		return 0;
	if (t <= my t [1])
		return 1;
	if (t >= my t [my nt])
		return my nt;
	/* Start binary search. */
	long left = 1, right = my nt;
	while (left < right - 1) {
		long mid = (left + right) / 2;
		if (t >= my t [mid]) left = mid; else right = mid;
	}

	if(right != left + 1){
	    std::cout<<"right != left + 1"<<std::endl;
		std::cout<<"Picth_to_Pointprocess.cpp 301."<<std::endl;
		return -1;
	}
	//Melder_assert (right == left + 1);
	return t - my t [left] < my t [right] - t ? left : right;
 }
