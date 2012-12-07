#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#include "CWav.h"

using namespace std;

void CWav::readInput(char *filename)
{
	
	ifstream inFile( filename, ios::in | ios::binary);
	
	//printf("Reading wav file...\n"); // for debugging only
	
	inFile.seekg(4, ios::beg);
	inFile.read( (char*) &myChunkSize, 4 ); // read the ChunkSize
	
	inFile.seekg(16, ios::beg);
	inFile.read( (char*) &mySubChunk1Size, 4 ); // read the SubChunk1Size
	
	inFile.seekg(20, ios::beg);
	inFile.read( (char*) &myFormat, sizeof(short) ); // read the file format.  This should be 1 for PCM
	
	//inFile.seekg(22, ios::beg);
	inFile.read( (char*) &myChannels, sizeof(short) ); // read the # of channels (1 or 2)
	
	//inFile.seekg(24, ios::beg);
	inFile.read( (char*) &mySampleRate, sizeof(int) ); // read the samplerate
	
	//inFile.seekg(28, ios::beg);
	inFile.read( (char*) &myByteRate, sizeof(int) ); // read the byterate
	
	//inFile.seekg(32, ios::beg);
	inFile.read( (char*) &myBlockAlign, sizeof(short) ); // read the blockalign
	
	//inFile.seekg(34, ios::beg);
	inFile.read( (char*) &myBitsPerSample, sizeof(short) ); // read the bitspersample
	
	inFile.seekg(40, ios::beg);
	inFile.read( (char*) &myDataSize, sizeof(int) ); // read the size of the data
	
	
	// read the data chunk
	myData = new char[myDataSize];
	inFile.seekg(44, ios::beg);
	inFile.read(myData, myDataSize);
	
	inFile.close(); // close the input file
	
	my_signal = NULL;
	
	if ( myBitsPerSample == 8 )
	{
		mySignalSize = myDataSize;
		my_signal = new short[mySignalSize];
		for ( int i = 0; i < myDataSize; i++ )
			my_signal[i] = (short)( (unsigned char) myData[i] );
		
	}
	else if ( myBitsPerSample == 16 ){
		mySignalSize = myDataSize / 2;
		my_signal = new short[mySignalSize];
		short val;
		for ( int i = 0; i < myDataSize; i+=2 )
		{
			val = (short)( (unsigned char) myData[i] );
			val += (short)( (unsigned char) myData[i+1] ) * 256;
			my_signal[i/2] = val;
		}
	}
}
