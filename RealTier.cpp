#include "iostream"
#include "Structure.h"
#include "RealTier.h"
#include "NUM.h"
#include "Pitch.h"
#include "math.h"

#include <stdio.h>
#include <assert.h>

template <typename T>

int compare(T x, T y){
   return x > y ? 1 : (x == y ? 0 : -1);
}

void PitchTier_modifyExcursionRange(PitchTier me, double tmin, double tmax, double multiplier, double fref_Hz) {
	if (fref_Hz <= 0 || multiplier == 1.0) 
		return;

	double fref_st = 12.0 * log (fref_Hz / 100.0) / NUMln2;
	for (long i = 1; i <= my points -> size; i++) {
		RealPoint point = (RealPoint)(my points->item[i]);   //ÀàÐÍ×ª»»
		double f = point -> value;
		if (point -> number < tmin || point -> number > tmax)
			continue;
		if (f > 0) {
			double f_st = fref_st + 12.0 * (log(f/fref_Hz)/log(2.0))*multiplier;  //Ô´ÂëÖÐÎª£ºlog2(f/fref_Hz)
			point -> value = 100.0 * exp (f_st * (NUMln2 / 12.0));
		}
	}
}

RealPoint RealPoint_create (double time, double value)
{
	RealPoint me = (RealPoint) malloc(sizeof(structRealPoint));
	my number = time;
	my value = value;
	return me;
}

void RealTier_init (RealTier me, double tmin, double tmax) {
	my xmin = tmin;
	my xmax = tmax;
	my points = SortedSetOfDouble_create ();
}

PitchTier PitchTier_create (double tmin, double tmax)
{
	PitchTier me = (PitchTier) malloc(sizeof(structPitchTier));
	RealTier_init (me, tmin,tmax);
	return me;
}

DurationTier DurationTier_create (double tmin, double tmax)
{
	DurationTier me = (DurationTier) malloc (sizeof(structDurationTier));
	RealTier_init (me, tmin,tmax);
	return me;
}

DurationTier PointProcess_upto_DurationTier (PointProcess me) {
   DurationTier thee = (DurationTier)malloc(sizeof(structDurationTier));
   for(long i = 1; i < my nt; ++ i)
	    RealTier_addPoint(thee, my t[i], 1.0);
   return thee;
}

void SortedSet_init(SortedSetOfDouble me, long initialCapacity)
{
	my _capacity = initialCapacity >= 1 ? initialCapacity : 1;
	my size = 0;
	my item = (SimpleDouble *)calloc(sizeof(structSimpleDouble *), (initialCapacity + 1));
   // my item --;   // base 1
}

SortedSetOfDouble SortedSetOfDouble_create (void)     //Collection.cpp 482
{
	SortedSetOfDouble me = (SortedSetOfDouble)calloc(sizeof(structSortedSetOfDouble), 1);
	SortedSet_init(me, 10);               //Collection.cpp 478
	return me;
}

long SortedSet_getposition(SortedSetOfDouble me, SimpleDouble data)    //Collection.cpp 398
{
	if(my size == 0)  return 1;          //empty set? then 'data' is going to be the first item 
	int pos = compare(data->number, my item[my size]->number);   // compare with the last item
	if(pos > 0) return my size + 1;     // insert at end
	if(pos == 0)  return 0;
	if(data->number < my item[1]->number)    // compare with the first item
		 return 1;        
	long left = 1, right = my size;

	while(left < right - 1){
	    long mid = (left + right) / 2;
		if(compare(data->number, my item[mid]->number) > 0)
			left = mid;
		else
			right = mid;
	}

	if(!compare(data->number, my item[left]->number) || !compare(data->number, my item[right]->number))
	     return 0;
	return right;
}
 // Collection.cpp 250
void SortedSet_insertItem(SortedSetOfDouble me, SimpleDouble data, long pos)
{ 
	if(my size >= my _capacity){
		SimpleDouble *dum = NULL;
		dum = (SimpleDouble *)calloc(sizeof(SimpleDouble), 2 * my _capacity*sizeof(SimpleDouble) + 1);
		memcpy(dum, my item, (my _capacity + 1)*sizeof(SimpleDouble));
		free(my item);
		my item = dum;	
      /*  
	      data lose when _capacity == 40 !!!???
	  my item = (SimpleDouble *)realloc(my item, 2 * my _capacity * sizeof(SimpleDouble) + 1); // some problem?? start position */

		my _capacity *= 2;
	}

	my size ++;
	for (long i = my size; i > pos; i --)
		my item [i] = my item[i - 1];
	my item[pos] = data;
}

//PitchTier.cpp 
void SortedSet_addItem(SortedSetOfDouble me, SimpleDouble data)   //Collection.cpp 267
{
    if(!data){
	  std::cout<<"Error, The data is NULL!"<<std::endl;
	  std::cout<<"PitchTier.cpp: Line 121"<<std::endl;
	  return ;
	}
	long index = SortedSet_getposition(me, data);
	
	if(index != 0)
		SortedSet_insertItem(me, data, index);                 //Collection.cpp 250
	else{
		//if(! my _dontOwnItems){
		//forget (data);   // could not insert; I am the owner, so I must dispose of the data
		std::cout<<"Cannot insert into the data."<<std::endl;
		std::cout<<"RealTier.cpp 121"<<std::endl;
		free(data);
		//}
	}
}

//Collection.cpp 478
void SortedSetOfDouble_init(SortedSetOfDouble me, long initialCapacity)
{   /*
	my _capacity = initialCapacity >= 1 ? initialCapacity : 1;
	my size = 0;
	my item = (SimpleDouble *) malloc (sizeof(structSimpleDouble *) * 2);*/
	SortedSet_init(me, initialCapacity);
}

void PitchTier_shiftFrequencies (PitchTier me, double tmin, double tmax, double shift, int unit)
{
	for(long i = 1; i < my points -> size; ++ i){
		RealPoint point = (RealPoint)my points->item[i];
		double frequency = point -> value;
		if(point->number < tmin || point->value > tmax)  continue;
		switch(unit){
			case kPitch_unit_HERTZ: {	
				frequency += shift;
				if (frequency <= 0.0){
					std::cout<<"The resulting frequency has to be greater than 0 Hz."<<std::endl;
				    std::cout<<"RealTier.cpp: Line 160"<<std::endl;
					return;
				}
			  } break;
		    case kPitch_unit_MEL: {
			      frequency = NUMhertzToMel (frequency) + shift;
			      if (frequency <= 0.0){
					std::cout<<"The resulting frequency has to be greater than 0 mel."<<std::endl;
				    std::cout<<"RealTier.cpp: Line 169"<<std::endl;
					return;				    
				  }
				  frequency = NUMmelToHertz(frequency);
				} break; 
			case kPitch_unit_LOG_HERTZ:
			        frequency = pow (10.0, log10 (frequency) + shift);
				    break; 
			case kPitch_unit_SEMITONES_1: 
			        frequency = NUMsemitonesToHertz(NUMhertzToSemitones(frequency) + shift);
				    break; 
			case kPitch_unit_ERB: {
					frequency = NUMhertzToErb(frequency) + shift;
					if (frequency <= 0.0){
					std::cout<<"The resulting frequency has to be greater than 0 mel."<<std::endl;
				    std::cout<<"RealTier.cpp: Line 184"<<std::endl;
					return;				    
				  }
				 frequency = NUMerbToHertz (frequency);
			}
		}
		point->value = frequency;
	}
}

void PitchTier_multiplyFrequencies (PitchTier me, double tmin, double tmax, double factor)
{
	if(factor <= 0.0){
		std::cout<<"Error, factor < 0.0"<<std::endl;
		std::cout<<"RealTier.cpp 194"<<std::endl;
		return;
	}
	//Ìí¼Óµ±factor==1Ê±²»ÐèÒªÐÞ¸ÄÖ±½Ó·µ»Ø
	if(factor == 1)
		return;

	for(long i = 1; i <= my points->size; ++ i){
		RealPoint point = (RealPoint)(my points->item[i]);
		if(point->number < tmin || point->number > tmax)
			continue;
		point->value *= factor;
	}
}

PitchTier Pitch_to_PitchTier (Pitch me)
{
	PitchTier thee = PitchTier_create(my xmin, my xmax);

	for(long i = 1; i <= my nx; ++ i){
		double frequency = my frame[i].candidate[1].frequency;
		//Count only voiced frames
		if(frequency > 0.0 && frequency < my ceiling){
		   double time = my x1 + (i - 1) * my dx;
		   RealTier_addPoint(thee, time, frequency);
		}
	}
	return thee;
}

void RealTier_addPoint(RealTier me, double t, double value)
{
	RealPoint point = RealPoint_create(t, value);
	SortedSet_addItem(my points, point);   
	///PraatÔ´ÂëÖÐµÄ Collection_addItem (my points, point.transfer());   Collection.cpp 267
}

//Pitch_PitchTier.cpp 126
Pitch Pitch_PuitchTier_to_Pitch(Pitch me, PitchTier tier){
	if(tier->points->size == 0){
	    std::cout<<"Error, tier size is equals to 0"<<std::endl;
		std::cout<<"RealTier.cpp 232"<<std::endl;
		return NULL;
	}

	Pitch thee = (Pitch)calloc(sizeof(structPitch), 1);
	if(!thee){
	    std::cout<<"Over flow"<<std::endl;
		std::cout<<"RealTier.cpp 232"<<std::endl;
		return NULL;	   
	}

	memcpy(thee, me, 1);

	for(long iframe = 1; iframe = my nx; ++ iframe){
	   Pitch_Frame frame = & thy frame[iframe];
	   Pitch_Candidate cand = & frame->candidate[1];
	   if(cand->frequency > 0.0 && cand->frequency <= my ceiling)
		   cand->frequency = RealTier_getValueAtTime(tier,my x1 + (iframe - 1) * my dx);//Sampled_indexToX (me, iframe)
	    cand -> strength = 0.9;
		frame -> nCandidates = 1;
	}

	return thee;
}

//RealTier.cpp 130
double RealTier_getValueAtTime(RealTier me, double t){
	long n = my points->size;
	if(n == 0)  return NUMundefined;
	RealPoint pointRight = (RealPoint)(me->points->item[1]);
	if(t <= pointRight->number)
		 return pointRight->value;   /* Constant extrapolation. */
	RealPoint pointLeft = (RealPoint)(me->points->item [n]);
	if (t >= pointLeft->number) 
		return pointLeft->value;   /* Constant extrapolation. */

	if(n < 1){
	   std::cout<<"RealTier points's size < 0"<<std::endl;
	   std::cout<<"RealTier.cpp  261"<<std::endl;
	   return 0;
	}

	long ileft = AnyTier_timeToLowIndex (me, t), iright = ileft + 1;
	if(ileft < 1 || iright > n){
	    std::cout<<"Error"<<std::endl;
	    std::cout<<"RealTier.cpp  261"<<std::endl;
	    return 0;
	}

	pointLeft = (RealPoint)(me->points->item[ileft]);
	pointRight = (RealPoint)(me->points->item[iright]);
    
	double tleft = pointLeft->number, fleft = pointLeft->value;
	double tright = pointRight->number, fright = pointRight->value;
	return t == tright ? fright   /* Be very accurate. */
		  : tleft == tright ? 0.5 * (fleft + fright)   /*Unusual, but possible; no preference.*/
		  : fleft + (t - tleft) * (fright - fleft) / (tright - tleft);   /*Linear interpolation.*/

}

//RealTier.cpp 171
double RealTier_getArea (RealTier me, double tmin, double tmax){
	long n = my points->size, imin, imax;      //long n = my f_getNumberOfPoints ()
	RealPoint *points = (RealPoint *) my points->item;
	if (n == 0) return NUMundefined;
	if (n == 1) return (tmax - tmin) * points[1]->value;
	imin = AnyTier_timeToLowIndex (me, tmin);
	if (imin == n) return (tmax - tmin) * points[n]->value;
	imax = AnyTier_timeToHighIndex (me, tmax);
	if (imax == 1) return (tmax - tmin) * points[1]->value;

	if(imin >= n || imax <= 1){
		std::cout<<"Error, imin = "<<imin<<", imax = "<<imax<<std::endl;
		std::cout<<"RealTier.cpp  265"<<std::endl;
		return -1;
	}
	/*
	 * Sum the areas between the points.
	 * This works even if imin is 0 (offleft) and/or imax is n + 1 (offright).
	 */
	double area = 0.0;
	for (long i = imin; i < imax; i ++) {
		double tleft, fleft, tright, fright;
		if (i == imin) tleft = tmin, fleft = RealTier_getValueAtTime (me, tmin);
		else tleft = points[i]->number, fleft = points[i]->value;
		if (i + 1 == imax) tright = tmax, fright = RealTier_getValueAtTime (me, tmax);
		else tright = points[i + 1]->number, fright = points[i + 1]->value;
		area += 0.5 * (fleft + fright) * (tright - tleft);
	}
	return area;  
}

//AnyTier.cpp 61
long AnyTier_timeToLowIndex (void * tc, double time){
	RealTier me = (RealTier)tc;
	if(my points->size == 0)  return 0;  // undefined
	long ileft = 1, iright = my points->size;
	AnyPoint *points = (AnyPoint *)(my points->item);

	double tleft = points[ileft]->number;
	if (time < tleft) return 0;   // offleft
	double tright = points[iright]->number;
	if (time >= tright) return iright;

	if(time < tleft || time > tright || tright < tleft){
		std::cout<<"get time Error"<<std::endl;
		std::cout<<"RealTier.cpp  392"<<std::endl;	
		return -1;
	}

	while (iright > ileft + 1) {
		long imid = (ileft + iright) / 2;
		double tmid = points[imid]->number;
		if (time < tmid) {
			iright = imid;
			tright = tmid;
		} else {
			ileft = imid;
			tleft = tmid;
		}
	}

	if(iright != ileft + 1 || ileft < 1 || iright > my points->size || 
	   time < points[ileft]->number || time > points[iright]->number){
	     std::cout<<"get time result error!"<<std::endl;
		 std::cout<<"RealTier.cpp  297"<<std::endl;
		 return -1;
	}

	return ileft;
}

//AnyTier.cpp 91
long AnyTier_timeToHighIndex (void *tc, double time){
	RealTier me = (RealTier)tc;
	if (my points->size == 0) return 0;   // undefined; is this right?
	long ileft = 1, iright = my points->size;
	AnyPoint *points = (AnyPoint *) my points->item;

	double tleft = points [ileft]->number;
	if (time <= tleft) return 1;
	double tright = points[iright]->number;
	if (time > tright) return iright + 1;   // offright

	if(time < tleft || time > tright || tright < tleft){
		std::cout<<"get time Error"<<std::endl;
		std::cout<<"RealTier.cpp  337"<<std::endl;
		return -1;
	}

	while (iright > ileft + 1) {
		long imid = (ileft + iright) / 2;
		double tmid = points [imid] -> number;
		if (time <= tmid) {
			iright = imid;
			tright = tmid;
		} else {
			ileft = imid;
			tleft = tmid;
		}
	}

	if(iright != ileft + 1 || ileft < 1 || iright > my points->size || 
	   time < points[ileft]->number || time > points[iright]->number){
	     std::cout<<"get time result error!"<<std::endl;
		 std::cout<<"RealTier.cpp  337"<<std::endl;
		 return -1;
	}

	return iright;
}

//RealTier.cpp 160 
double RealTier_getMinimumValue (RealTier me) {
	double result = NUMundefined;
	long n = my points -> size;
	for (long i = 1; i <= n; i ++) {
		RealPoint point = (RealPoint) my points->item [i];
		if (result == NUMundefined || point->value < result)
			result = point->value;
	}
	return result;
}

//Sound_extension.cpp 1489
void PitchTier_modifyRange_old (PitchTier me, double tmin, double tmax, double factor, double fmid){
	for (long i = 1; i <= my points -> size; i ++) {
		RealPoint point = (RealPoint) my points->item [i];
		double f = point->value;
		if (point->number < tmin || point->number > tmax) {
			continue;
		}
		f = fmid + (f - fmid) * factor;
		point->value = f < 0 ? 0 : f;
	}
}
