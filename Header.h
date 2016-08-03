#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <iomanip> 
#include <math.h>
#include <list>
#include <sstream>
#include "fft.h"
#include <functional>
#include <set>

#define TRUE        1
#define FALSE       0

#define HEADER_SIZE  44
#define SAMPLE_RATE  44100
#define LP_CUTOFF    2500
#define HP_CUTOFF    500
#define MAX_VALUE    32768
#define OFFSET       4
#define NEIGHBOURS   3

using namespace std;

class wavHdr{
public:
    char                RIFF[4];        // RIFF Header     
    unsigned long       ChunkSize;      // RIFF Chunk Size  
    char                WAVE[4];        // WAVE Header      
    char                fmt[4];         // FMT header       
    unsigned long       Subchunk1Size;  // Size of the fmt chunk                                
    unsigned short      AudioFormat;    // Audio format
    unsigned short      NumOfChan;      // Number of channels                
    unsigned long       SamplesPerSec;  // Sampling Frequency in Hz                             
    unsigned long       bytesPerSec;    // Bytes per second 
    unsigned short      blockAlign;     // 2=16-bit mono, 4=16-bit stereo 
    unsigned short      bitsPerSample;  // Number of bits per sample      
    char                Subchunk2ID[4]; // "data"  string   
    unsigned long       Subchunk2Size;  // Sampled data length    
};

class WaveBother
{
public:
    WaveBother();
   ~WaveBother();

	void   addTrack(string name);
	void   findTrack(string name);
	multimap<double, string> database;	
	int loadDB();
	string dbPrefix;
private:
	
	string suspPrefix;
	string dbFilename;
	


    int length;

	short searchMax(short* data); 

	int findMaxAmp(int i, int k, vector<double> frequencies);

	set<double> WaveBother::fingerPrint(string name);

	void   loadWavs();
    void   normalize(short* signal);   
};
