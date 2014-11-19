#include <iostream>
#include <string>
#include "Structure.h"
#include "Pitch.h"
#include "SoundCompute.h"
#include "ReadWavFile.h"
#include "Sound_to_Pitch.h"
#include "Pitch_to_PointProcess.h"
#include <fstream>
#include "Get_Data_to_Sound.h"
#include "math.h"
using namespace std;

#define FREQUENCY 22050.0
#define DURATION 0.1
#define DURATIONSTEP 0.04
#define FRAMESIZE FREQUENCY*0.1
#define FRAMESTEP FREQUENCY*0.04

#define formatRatio 1.0
#define pitchmedian 225.0
#define pitchRangeFactor 0.0
#define durationFactor 1.0
#define new_pitch 1.8

void test()
{
	double srcdata[120000];//*
	int count=120000;//*
	FILE* fd = NULL;

	memset(srcdata, 0, sizeof(double)*count);

	fd = fopen("D:/mf.txt","r");
	if(fd == NULL){
		std::cout<<"the file can't be opened"<<std::endl;
		return ;
	}
	int r;
	for(r=0; r < count; r++)
	{
		fscanf(fd,"%lf",&srcdata[r]);
	}
	fclose(fd);

	Sound sound0 = NULL, sound2 = NULL, changeSound = NULL;
	Sound checkSound = NULL;
	Pitch pitch0 = NULL, pitch2 = NULL;
	PointProcess pulses0 = NULL, pulses2 = NULL, pulses01 = NULL;

	std::string  pathfrom, sspathto, cspathto;
	FILE *csfp = NULL, *ssfp = NULL;

	int flag,i;
	int temp;
	double lasttime = 0, nexttime = 1; //¸Õ¿ªÊ¼³õÊ¼»¯Îª0  
	double pulselasttime = 0, pulsenexttime = 0;
	double xmin = 0, xmax=0;
	int presoundindex = 1,soundindex = 1, prepitchindex = 1, pitchindex = 1;//Êý×éµÄµÚÒ»¸öarray[0]²»Ê¹ÓÃ
	int vectorstart = 1, vectorend = 1;
	int checkindex = 1;
	int pointer = 1;
	int what = 0;

	for(i = 1,flag = 1; flag <90000; i++)//¸Ã½áÊøÌõ¼þÖ»ÊÇ³õ²½²âÊÔ
	{
		sound0 = Sound_createSimple(DURATION, FREQUENCY, 1);
		Get_Data(srcdata, sound0->z[1], &flag, count);

		if(i == 1)//Õâ²¿·ÖÖ»ÊÇÎªÁË²âÊÔ±äÉùµÄÊý¾Ý
		{
			checkSound = Sound_create(0,6,180000,sound0->dx,sound0->x1,1);//*
			for(int t = 0; t <= 180000; ++t)//*
			{
				checkSound->z[1][t] = 0.0;
			}
		}
		
		pitch0 = computePitch(sound0);

		
		if(nexttime > lasttime&&lasttime >= 0)// && soundindex > presoundindex    nexttime > DURATIONSTEP&&
		{
			printf("½øÈë±äÉù%d\n",i);

		
			//changeSound = Sound_and_Pitch_and_Pulses_changeGender_old (sound0, pitch0, pulses2, formatRatio ,pitchmedian, pitchRangeFactor,																				durationFactor);			
			//changeSound = Sound_and_pitch_and_pulses_changespeaker(sound0, pitch0, pulses2, formatRatio ,new_pitch, pitchRangeFactor,																					durationFactor);

			changeSound = Sound_and_Pitch_changeGender_old(sound0, pitch0, formatRatio ,pitchmedian, pitchRangeFactor,																				durationFactor);
			//changeSound = Sound_and_pitch_changespeaker(sound0, pitch0, formatRatio ,new_pitch, pitchRangeFactor,																					durationFactor);
			
			
			if(changeSound!=NULL)//there were no voiced segment found
			{
				what++;
				nexttime =changeSound->xmax;//ÕâÀïÔÙ¿¼ÂÇÒ»ÏÂ changeSound->nx/FREQUENCY;
				printf("±äÉùÕý³£\n");

				temp = floor(lasttime * FREQUENCY +0.5);
				if(temp == 0)
					temp = 1;
				for(int a = checkindex,b=0; a <= checkindex + changeSound->nx-temp; a++,b++)
				{
					checkSound->z[1][a] = changeSound->z[1][temp+b];
				}
				checkindex += changeSound->nx - temp + 1; //¼ì²é
				cout<<"flag = "<<flag<<endl;
				cout<<"checkindex = "<<checkindex<<endl;

				lasttime = nexttime - DURATIONSTEP;//ÔÙ¿¼ÂÇÒ»ÏÂ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
			}
			else
			{
				printf("changeSound = NULL,copyÔ­Êý¾Ý\n");
				temp = floor(lasttime * FREQUENCY +0.5);
				if(temp == 0)
					temp = 1;
				for(int a = checkindex,b=0; a <= checkindex + sound0->nx-temp; a++,b++)
				{
					checkSound->z[1][a] = sound0->z[1][temp+b];
				}
				checkindex += sound0->nx-temp+1; 
				lasttime = 0;
				flag += sound0->nx+1-FRAMESTEP;//Èç¹û¸Ä±äFRAMESTEP´óÐ¡£¬ÕâÀïÒ²Òª±ä£¬Òª²»ÒªÔÙ¼õÈ¥Ò»¸ötemp£¿£¿£¿£¿

			}


			free(changeSound);

		}
		else
		{
			printf("nexttime < lasttime£¬copyÔ­Êý¾Ý\n");
			//¸ü½øÒ»²½µÄ»°¾ÍÒªcopy´Ólasttimeµ½sound0->xmax,È»ºó±äÕâÒ»²¿·ÖµÄÉùÒô£¬¶ø²»ÊÇÏÂÃæÕâ¶Î
			temp = floor(lasttime * FREQUENCY +0.5);
			if(temp == 0)
				temp = 1;
			for(int a = checkindex,b=0; a <= checkindex + sound0->nx-temp; a++,b++)
			{
				checkSound->z[1][a] = sound0->z[1][temp+b];
			}
			checkindex += sound0->nx-temp+1; 
			lasttime = 0;
			flag += sound0->nx+1-FRAMESTEP;//Èç¹û¸Ä±äFRAMESTEP´óÐ¡£¬ÕâÀïÒ²Òª±ä    £¿£¿£¿£¿£¿-temp
		}


		free(sound0);	
		free(pitch0);	

	}
	std::cout<<"Please Enter the path to .xxx file: "<<std::endl;
	std::cin>>cspathto;

	csfp = fopen(cspathto.c_str(), "w+");
	Sound_writeTxt(checkSound, csfp);

	free(checkSound);
	cout<<"test end"<<endl;


}


double* test3(double *source, double frequency, double formatratio,double pitch_median,
								double pitchrangefactor, double durationfactor,double newpitch,bool button_info)
{
	Sound sound0 = NULL, changeSound = NULL;
	Pitch pitch0 = NULL;
	double changedata[882];
	int a,b;

	sound0 = Sound_createSimple(DURATION, frequency, 1);
	for(a = 1, b = 0; b < 2205; a++, b++)
	{
		sound0->z[1][a] = source[b]; 
	}
	pitch0 = computePitch(sound0);
	if(button_info == true)
		changeSound = Sound_and_Pitch_changeGender_old(sound0, pitch0, formatratio ,pitch_median, pitchrangefactor,																				durationfactor);
	else
		changeSound = Sound_and_pitch_changespeaker(sound0, pitch0, formatratio ,newpitch, pitchrangefactor,durationfactor);								
	if(changeSound!=NULL)
	{
		for(a = 0, b = 1324; a<882 && b<= changeSound->nx; a++,b++)
		{ 
			changedata[a] = changeSound->z[1][b];
		}
	}
	else
	{
		for(a = 0, b = 1324; a<882 && b<= sound0->nx; a++,b++)
		{
			changedata[a] = sound0->z[1][b];
		}
	}
	return changedata;
}


void test4()
{
	printf("hello world!\n");
}

extern "C" _declspec(dllexport) int add(int a, int b)
{
	printf("a+b³ÌÐò");
	return a+b;
}

int Add(int plus1, int plus2)
{
	int add_result = plus1 + plus2;
	return add_result;
}

int test2()
{
    Sound sound3 = NULL, changeSound3 = NULL;
	Pitch pitch3 = NULL;

	PointProcess pointprocess3 = NULL;
	VoicedUnvoiced VUV = NULL;

	std::string  pathfrom, sspathto, cspathto;
	std::cout<<"Please Enter the path of .wav file: "<<std::endl;
	std::cin>>pathfrom;
	sound3 = Sound_readFromSoundFile(pathfrom);
	
	if (!sound3) {
		std::cout<<"Error, No Sound read from file.!"<<std::endl;
		return 0;
	}
	pitch3 = computePitch(sound3);
	if (!pitch3) { 
		std::cout<<"Error, No Pitch computed from sound.!"<<std::endl;
		return 0;     
	}

   //changeSound3 = Sound_and_Pitch_changeGender_old(sound3, pitch3, 0.8, 0.0, 1.0, 1.0);
   changeSound3 = Sound_and_pitch_changespeaker(sound3, pitch3, 1.5, 1.5, 1, 1);
   
   FILE *csfp = NULL, *ssfp = NULL;
   do{
	   std::cout<<"Please Enter the path to .xxx file for storing the source soud file: "<<std::endl;
	   std::cin>>sspathto;
	   ssfp = fopen(sspathto.c_str(), "w+");
	   std::cout<<"Please Enter the path to .xxx file: "<<std::endl;
	   std::cin>>cspathto;
	   csfp = fopen(cspathto.c_str(), "w+");
	   if(csfp == NULL || ssfp == NULL)
           std::cout<<"open file "<<sspathto.c_str()<<" or "<<cspathto.c_str()<<" error. the file don't exist\nre-input file path:"<<std::endl;
   }while(csfp == NULL || ssfp == NULL);
	    


	Sound_writeTxt(sound3, ssfp);
    Sound_writeTxt(changeSound3, csfp);

}

void test3(double source[], double changedata[], double frequency, double formatratio,double pitch_median,
								double pitchrangefactor, double newpitch,bool button_info)
{
	
	Sound sound0 = NULL, changeSound = NULL;
	Pitch pitch0 = NULL;
	double duration;
	int flag;


	flag = (int)changedata[0];
	memcpy(changedata, source, flag*sizeof(double));
	
	if(flag <= 2)
		return ;

	duration = flag/frequency;
	sound0 = Sound_createSimple(duration, frequency, 1);
	int flag2 = sound0->nx - 1;

	memcpy(sound0->z[1] + 1, source, flag*sizeof(double));
	
	/*
	for(a = 1, b = 0; a < flag2 && b < flag; a++, b++)
	{
		sound0->z[1][a] = source[b];
	}*/

	
	pitch0 = computePitch(sound0);
	if(pitch0 == NULL)
		return;

	if(button_info)
		changeSound = Sound_and_Pitch_changeGender_old(sound0, pitch0, formatratio ,pitch_median, pitchrangefactor,1.0);
	else
		changeSound = Sound_and_pitch_changespeaker(sound0, pitch0, formatratio ,newpitch, pitchrangefactor,1.0);								
	if(changeSound != NULL)
		memcpy(changedata, changeSound->z[1] + 1, flag*sizeof(double));
	
	free(sound0);
	free(pitch0);
	free(changeSound);
	
}


//ËÍºó40ÃëµÄtestº¯Êý
void test40(double source[], double changedata[], double frequency, double formatratio,double pitch_median,
								double pitchrangefactor, double newpitch,bool button_info)
{

	Sound sound0 = NULL, changeSound = NULL;
	Pitch pitch0 = NULL;
	double duration;
	int flag;
	int start = floor(0.06*frequency) + 1;

	flag = (int)changedata[0];
	
	if(flag <= 2)
	{
		memcpy(changedata, source + start , floor(0.04*frequency)*sizeof(double));//copy the last 40ms
		return ;
	}
	duration = flag/frequency;
	sound0 = Sound_createSimple(duration, frequency, 1);

	if(sound0->nx > flag)
		memcpy(sound0->z[1] + 1, source, flag*sizeof(double));
	else
		memcpy(sound0->z[1] + 1, source, sound0->nx*sizeof(double));

	pitch0 = computePitch(sound0);
	if(pitch0 == NULL)
	{
		memcpy(changedata, source + start , floor(0.04*frequency)*sizeof(double));//copy the last 40ms
		return;
	}
	if(button_info)
		changeSound = Sound_and_Pitch_changeGender_old(sound0, pitch0, formatratio ,pitch_median, pitchrangefactor,1.0);
	else
		changeSound = Sound_and_pitch_changespeaker(sound0, pitch0, formatratio ,newpitch, pitchrangefactor,1.0);	

	if(changeSound != NULL)
		memcpy(changedata, changeSound->z[1] + start,  floor(0.04*frequency)*sizeof(double));
	else
	{
		memcpy(changedata, source + start ,  floor(0.04*frequency)*sizeof(double));//copy the last 40ms
	}
	
	 
	free(sound0);
	free(pitch0);
	free(changeSound);
	
}

Sound test_multi(double source[], double changedata[], double frequency, double formatratio,double pitch_median,
								double pitchrangefactor, double newpitch,bool button_info)
{
	
	Sound sound0 = NULL, changeSound = NULL,checkSound = NULL;
	Pitch pitch0 = NULL;
	double duration;
	int flag;
	int a;
	double max;

	flag = (int)changedata[0];
	max = source[0];

	for(a = 0; a < flag; a++)
	{
		if(fabs(source[a]) >= fabs(max))
			max = source[a];
	}
	for(a = 0; a < flag; a++)
	{
		source[a] /= fabs(max);
	}

	memcpy(changedata, source, flag*sizeof(double));
	
	if(flag <= 2)
		return NULL;

	duration = flag/frequency;
	sound0 = Sound_createSimple(duration, frequency, 1);
	int flag2 = sound0->nx - 1;

	memcpy(sound0->z[1] + 1, source, flag*sizeof(double));

	pitch0 = computePitch(sound0);
	if(pitch0 == NULL)
		return NULL;

	if(button_info)
		changeSound = Sound_and_Pitch_changeGender_old(sound0, pitch0, formatratio ,pitch_median, pitchrangefactor,1.0);
	else
		changeSound = Sound_and_pitch_changespeaker(sound0, pitch0, formatratio ,newpitch, pitchrangefactor,1.0);								
	if(changeSound != NULL)
	{
		
		checkSound = Sound_create(0,6,180000,sound0->dx,sound0->x1,1);//*
		for(int t = 0; t <= 180000; ++t)//*
		{
			checkSound->z[1][t] = 0.0;
		}
		memcpy(changedata, changeSound->z[1] + 1, flag*sizeof(double));
		for(a = 0; a < flag; a++)
		{
			changedata[a] *= fabs(max);
		}
		memcpy(checkSound->z[1], changedata, flag*sizeof(double));
	}
	else
	{	
		for(a = 0; a < flag; a++)
		{
			changedata[a] *= fabs(max);
		}
		return NULL;
	}



	free(sound0);
	free(pitch0);
	free(changeSound);
	return checkSound;
}

void testoftest3()
{
	double srcdata[64000];//*
	int count=64000;//*
	FILE* fd = NULL;

	memset(srcdata, 0, sizeof(double)*count);

	fd = fopen("D:/mf2.txt","r");
	if(fd == NULL){
		std::cout<<"the file can't be opened"<<std::endl;
		return ;
	}
	int r;
	for(r=0; r < count; r++)
	{
		fscanf(fd,"%lf",&srcdata[r]);
	}
	fclose(fd);

	double changedata[63000];
	changedata[0] = 63000;

	Sound sound = NULL;

	sound = test_multi(srcdata, changedata, FREQUENCY, formatRatio, pitchmedian, pitchRangeFactor, new_pitch, false);
	
	std::string  pathfrom, sspathto, cspathto;
	FILE *csfp = NULL, *ssfp = NULL;

	std::cout<<"Please Enter the path to .xxx file: "<<std::endl;
	std::cin>>cspathto;

	csfp = fopen(cspathto.c_str(), "w+");
	Sound_writeTxt(sound, csfp);
	
	

}


int main(int argc, char *argv[])
{
	//test2();
	//test();
   testoftest3();
   return 0;
}

