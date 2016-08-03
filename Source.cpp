#include "Header.h"

WaveBother::WaveBother()
{
	dbPrefix   = "C:\\Media\\";
	suspPrefix = "D:\\Media\\";
	dbFilename = "C:\\Media\\fingerprintBase.spct";
}

WaveBother::~WaveBother()
{
}

int WaveBother::loadDB()
{
	string name;

	ifstream source(dbFilename);

	double fPrint;
	
	cout << "Loading database...";

	while(!source.fail())
	{
		source >> fPrint;
		getline(source, name);
		database.insert(pair<double, string>(fPrint, name));
	}

	if(database.size() == 0)
	{
        cout << endl << "Database is empty." <<endl;
		
		return 0;
	}
	else
		cout << "Done." << endl;
	
	source.close();

	return 1;
}

/*
void WaveBother::normalize(short* signal)
{
	short max = searchMax(signal);
    float ratio = MAX_VALUE / abs(max);
	
    for(int i = 0; i < length; i++)
        signal[i] = signal[i] * ratio;
}


int WaveBother::findTempo(short* data, float optRatio, float optTreshold)
{
	short* filteredData = new short[length]; 

	for(int i = 0; i < length; i++)
		filteredData[i] = data[i];

	normalize(filteredData);

	int max       = searchMax(filteredData);
	int threshold = abs(0.90 * max);

	vector<int>             peaks;
	vector<double>         tempos;
	map<double, int>         vars;
	map<double, int>::iterator it;

    for(int i = 0; i < length;)
    {
		if(abs(filteredData[i]) > threshold)
		{
			peaks.push_back(i);
			i += 10000;
		}

		i++;
	}

	double sum = 0;

	for(int i = 1; i < peaks.size(); i++)
	{
		tempos.push_back((60 * 44100)/(peaks[i] - peaks[i - 1]));

		while(tempos[i - 1] < 90) 
			tempos[i - 1] = ceil(tempos[i - 1] * 2);

		while(tempos[i - 1] > 180) 
			tempos[i - 1] = ceil(tempos[i - 1] / 2);
				
		sum+= tempos[i-1];

		vars.insert(pair<double,int>(tempos[i - 1], 0));

		it = vars.find(tempos[i - 1]);
		it->second++;
	}

	if(vars.size() == 0)
		return 0;

	int maxtempo = it->second;
	int tempo = 0;

	for(auto i = vars.cbegin(); i != vars.cend(); ++i)
	{
		if(i->second > maxtempo)
		{
			maxtempo = i->second;
			tempo    = i->first;
		}
	}

	float ratio = (float)maxtempo/tempos.size();
	
	if(ratio >= optRatio)
		return tempo;
	else
		return 0;
}
*/

void WaveBother::addTrack(string name)
{ 
	ofstream spectre(dbPrefix + "fingerprintBase.spct", ios_base::app); 

	set<double> fingerprints = fingerPrint(suspPrefix + name);
	set<double>::iterator it;

	for(it = fingerprints.begin(); it != fingerprints.end(); ++it)
	{
		spectre << setprecision(16) << *it << " " << name << endl;
	}
	
	spectre.close();

	cout << name + " added to base." << endl; 
}

void WaveBother::findTrack(string name)
{
	set<double> fingerprints = fingerPrint(dbPrefix + name);
	
	multimap<double, string>::iterator it;
	
	map<string, int> occasions;
	map<string, int>::iterator iter;

	clock_t startTime = clock();
	
	for(it = database.begin(); it != database.end(); ++it)
	{
		if(fingerprints.find(it->first) != fingerprints.end())
		{
			occasions[it->second]++;
		}
	}

	cout << name << " " << double(clock() - startTime) / (double)CLOCKS_PER_SEC<< " seconds." << endl;

	int max = 0;
	int k = 0;
	string j;

	for(iter = occasions.begin(); iter != occasions.end(); ++iter)
	{
		 //cout << iter->first << " " << iter->second << endl;

		if(iter->second > max)
		{
			max = iter->second;
			j = iter->first;
		}
	}
	
	cout << j << " " << occasions[j] << endl << endl;
}

set<double> WaveBother::fingerPrint(string name)
{
	int numOfSamples;
	
	vector<double>         frequencies;
	vector<pair<int, int>> freqTime;
	set<double>            diffs;

	ifstream source(name,  ios_base::binary);

	wavHdr header;
	
	source.read((char*)&header, HEADER_SIZE);

	numOfSamples = (header.Subchunk2Size * 8)/header.bitsPerSample;
	this->length = numOfSamples;
	
    short* mono = new short[numOfSamples];
	short* outputBuff = new short[numOfSamples];
	short* pout       = outputBuff;
	source.read((char*)mono, header.Subchunk2Size);

	// Downsampling and Hamming window application
	//
	if(header.SamplesPerSec = 44100)
	{
		for(int i = 0; i < numOfSamples; i += 2)
		{
			*pout++ = (mono[i] + mono[i + 1]) / 2;
			outputBuff[i/2] *= 0.54 - 0.46 * cos((M_PI * i) / 1023);
		}
	}

	numOfSamples /= 2;
	//////////////////////////////////////////////

	int*    tempPos = new int[10];
	double  tempSum = 0;

	complex *pSignal = new complex[2048];

	// Starting FFT Algorithm
	//
	for(int i = 0; i < numOfSamples / 2048; i++)
	{
		//Computing Fourier coeffs
		//
		for(int k = 0; k < 2048; k++)
			pSignal[k] = outputBuff[i * 2048 + k];

		CFFT::Forward(pSignal, 2048);

		for(int k = 0; k < 1024; k++)
			frequencies.push_back(10 * log10(sqrt(pSignal[k].re() * pSignal[k].re() + pSignal[k].im() * pSignal[k].im())));
		//////////////////////////////////////////////////

		// Dividing spectre into bands
		//
		tempPos[0] = findMaxAmp(100, 105, frequencies);
		tempPos[1] = findMaxAmp(105, 110, frequencies);
		tempPos[2] = findMaxAmp(110, 117, frequencies);
		tempPos[3] = findMaxAmp(117, 130, frequencies);
		tempPos[4] = findMaxAmp(130, 150, frequencies);
		tempPos[5] = findMaxAmp(150, 180, frequencies);
		tempPos[6] = findMaxAmp(180, 220, frequencies);
		tempPos[7] = findMaxAmp(220, 300, frequencies);
		tempPos[8] = findMaxAmp(300, 400, frequencies);
		tempPos[9] = findMaxAmp(400, 500, frequencies);
		//////////////////////////////////////////////////

		// Average amplitude of bands
		//
		for(int j = 0; j < 10; j++)
			tempSum += frequencies[tempPos[j]];

		tempSum /= 10;
		/////////////////////////////////////////////////


		// Filling time-frequency vector
		//
		for(int j = 0; j < 10; j++)
		{
			if((frequencies[tempPos[j]] > tempSum))
			{
				freqTime.push_back(pair<int,int>(i, tempPos[j]));
			}
		}
		////////////////////////////////////////////////

		tempSum  = 0;
		frequencies.clear();
	}
	////////////////////////////////////////////////////
	free(pSignal);
	stringstream ss;

	// Points union logic
	//
	for(int i = 0; i < freqTime.size() - OFFSET * NEIGHBOURS; i++)
	{
		ss << freqTime[i].second;

		for(int l = 0; l < NEIGHBOURS; l++)
		{
				ss << freqTime[i + OFFSET * l].first - freqTime[i].first
				   << freqTime[i + OFFSET * l].second;
		}

		diffs.insert(stod(ss.str()));

		ss.str("");
	}
	///////////////////////////////////////////////////

	free(mono);
	free(outputBuff);

	source.close();

	return diffs;
}

short WaveBother::searchMax(short* data)
{
	short max = data[0];

    for(int i = 0; i < length; i++)
	{
        if (abs(data[i]) > abs(max))
            max = data[i];
	}

    return max;
}

int WaveBother::findMaxAmp(int i, int k, vector<double> frequencies)
{
	double max = frequencies[i];
	int    pos = i;

	while(i < k)
	{
		if(frequencies[i] >= max)
		{
			max = frequencies[i];
			pos = i;
		}

		i++;
	}

	return pos;
}

int main(int argc, char *argv[])
{
	WaveBother bd;
	string t;

	/*
	for(int i = 1; i <= 900; i++)
	{
		t = "part (" + to_string(i);
		//t += i;
		t += ").wav";
		bd.addTrack(t);
	}
	*/
	//bd.database.get_allocator().allocate(30000000);
	if(bd.loadDB() == 0)
	{
		getchar();
		
		return 0;
	}

	string selection;
	string name;

	while(TRUE)
	{
		cout << "[1] For add track to the base" << endl;
		cout << "[2] For find track in the base" << endl;
		cout << ":";

		cin >> selection;

		if((selection == "1") || (selection == "2"))
		{
			cout << "Put .WAV file (44100 or 22050 Mono PCM) into C:\\Media\\ directory." << endl;
			cout << "Enter the name: ";
			cin >> name;

			if(selection == "2")
				bd.findTrack(name);
			else
				bd.addTrack(name);
		}
	}

	getchar();
    return 0;
}