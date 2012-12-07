#ifndef CWavH
#define CWavH

class CWav
{
private:
	char* 	myData;
public:
	int 	myChunkSize;
	int		mySubChunk1Size;
	short 	myFormat;
	short 	myChannels;
	int   	mySampleRate;
	int   	myByteRate;
	short 	myBlockAlign;
	short 	myBitsPerSample;
	int		myDataSize;

	short *my_signal;
	int mySignalSize;

public:

	void readInput(char *filename);

};

#endif
