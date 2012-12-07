/******************************************************************************
 *
 *     Program:       testtone
 *     
 *     Description:   This program generates a two-second 440 Hz test tone and
 *                    writes it to a .wav file.  The sound file has 16-bit
 *                    samples at a sample rate of 44.1 kHz, and is monophonic. 
 *
 *     Author:        Leonard Manzara
 *
 *     Date:          November 21, 2009
 *
 ******************************************************************************/


/*  HEADER FILES  ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "CWav.h"
using namespace std;

#define DEBUG_MODE

/*  CONSTANTS  ***************************************************************/
#define PI					3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY			440.0

/*  Test tone duration in seconds  */
#define DURATION			2.0				

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE			44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE		16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC			1
#define STEREOPHONIC		2


/*  FUNCTION PROTOTYPES  *****************************************************/
void createTestTone(double frequency, double duration,
					int numberOfChannels, int bitsPerSample, int sampleRate, char *filename);
void writeWaveFileHeader(int channels, int numberSamples, int bitsPerSample,
						 double sampleRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

void readInput(char *filename);

void manipulate(CWav *original, int manipulateChoice);
void manipulate(CWav *original, CWav *impulse, char *outputFile);

void convolve(float x[], int N, float h[], int M, float y[], int P);
void convolve(CWav *original, CWav *impulse, float y[], int P);

/******************************************************************************
 *
 *	function:	main
 *
 *	purpose:	Creates the test tone and writes it to the
 *               specified .wav file.
 *			
 *   arguments:	argv[1]:  the filename for the output .wav file
 *                       
 ******************************************************************************/

int main (int argc, char *argv[])
{
	char *outputFilename = "sine.wav";
	
	/*  Create the sine wave test tone, using the specified
	 frequency, duration, and number of channels, writing to
	 a .wav file with the specified output filename  */
	
	//createTestTone(FREQUENCY, DURATION, MONOPHONIC, BITS_PER_SAMPLE, SAMPLE_RATE, outputFilename);
		
	CWav *inputSignal = new CWav();
	inputSignal->readInput("testCase1.wav");
	
	//manipulate(inputSignal, 2);
	
	CWav *impulse = new CWav();
	impulse->readInput("smallCave.wav");
	//impulse->readInput("Pool.wav");
	
	cout << "Input Signal: " << inputSignal->mySignalSize << ", Impulse Size: " << impulse->mySignalSize << endl;
	
	manipulate(inputSignal, impulse, "output_smallCave.wav");
	
	/*  End of program  */
	return 0;
}



/******************************************************************************
 *
 *       function:       createTestTone
 *
 *       purpose:        Calculates and writes out a sine test tone to file
 *
 *		arguments:		frequency:  frequency of the test tone in Hz
 *						duration:  length of the test tone in seconds
 *                       numberOfChannels:  number of audio channels
 *						filename:  name of the file to create
 *                       
 *       internal
 *       functions:      writeWaveFileHeader, fwriteShortLSB
 *
 *       library
 *       functions:      ceil, pow, fopen, fprintf, sin, rint, fclose
 *
 ******************************************************************************/

void createTestTone(double frequency, double duration,
					int numberOfChannels, int bitsPerSample, int sampleRate, char *filename)
{
	int i;
	
	/*  Calculate the number of sound samples to create,
	 rounding upwards if necessary  */
	int numberOfSamples = (int)ceil(duration * sampleRate);
	
	/*  Calculate the maximum value of a sample  */
	int maximumValue = (int)pow(2.0, (double)bitsPerSample - 1) - 1;
	
	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
		return;
	}
	

	/*  Write the WAVE file header  */
	writeWaveFileHeader(numberOfChannels, numberOfSamples, bitsPerSample,
						sampleRate, outputFileStream);
	
	/*  Create the sine tone and write it to file  */
	/*  Since the frequency is fixed, the angular frequency
	 and increment can be precalculated  */
	double angularFrequency = 2.0 * PI * frequency;
	double increment = angularFrequency / sampleRate;
	
	int *testTone = new int[numberOfSamples];
	
	for (i = 0; i < numberOfSamples; i++) {
		/*  Calculate the sine wave in the range -1.0 to + 1.0  */
		double value = sin(i * increment);
		
		/*  Convert the value to a 16-bit integer, with the
		 range -maximumValue to + maximumValue.  The calculated
		 value is rounded to the nearest integer  */
		short int sampleValue = rint(value * maximumValue);
		
		testTone[i] = (int)sampleValue;
		
		/*  Write out the sample as a 16-bit (short) integer
		 in little-endian format  */
		fwriteShortLSB(sampleValue, outputFileStream);
		
		/*  If stereo output, duplicate the sample in the right channel  */
		if (numberOfChannels == STEREOPHONIC)
			fwriteShortLSB(sampleValue, outputFileStream);
	}
	
	
	
	/*  Close the output file stream  */
	fclose(outputFileStream);
}


void manipulate(CWav *original, int manipulateChoice)
{
	float *impulse_response = new float[1000];
	int impulse_size;

	for ( int j = 0; j < 1000; j++ )
		impulse_response[j] = 0;
	

	switch (manipulateChoice) {
		case 0:
			/*  Create an "identity" impulse response.  The output should be
			 the same as the input when convolved with this  */
			impulse_response[0] = 1;
			impulse_size = 1;			
			break;
		case 1: // turn down the volume
			impulse_response[0] = 0.1;
			impulse_response[1] = 0.2;
			impulse_response[2] = 0.1;
			impulse_size = 3;	
			break;
		default:
			for ( int j = 0; j < 1000; j++ )
				impulse_response[j] = 1.0/50;
			impulse_size = 1000;	
			break;
	}


	int output_size = original->mySignalSize + impulse_size - 1;
	
	float *output_signal = new float[output_size];
	
	float *signalInFloat = new float[original->mySignalSize];
	for ( int i = 0; i < original->mySignalSize; i++ )
		signalInFloat[i] = original->my_signal[i];
		
	convolve(signalInFloat, original->mySignalSize, impulse_response, impulse_size, output_signal, output_size);

	
	int i;
	
	/*  Calculate the number of sound samples to create,
	 rounding upwards if necessary  */
	int numberOfSamples = output_size;
	
	
	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen("artificial.wav", "wb");
	
	/*  Write the WAVE file header  */
	writeWaveFileHeader(original->myChannels, numberOfSamples, original->myBitsPerSample,
						original->mySampleRate, outputFileStream);
	
	float maxValInResult = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if ( output_signal[i] > maxValInResult )
			maxValInResult = output_signal[i];

	float maxValInInput = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if ( original->my_signal[i] > maxValInInput )
			maxValInInput = original->my_signal[i];
			
	for (i = 0; i < numberOfSamples; i++ )
		fwriteShortLSB((short)(output_signal[i] / maxValInResult * maxValInInput), outputFileStream);
	
	
	/*  Close the output file stream  */
	fclose(outputFileStream);
	

}

void manipulate(CWav *original, CWav *impulse, char *outputFile)
{
	int output_size = original->mySignalSize + impulse->mySignalSize - 1;
	
	float *output_signal = new float[output_size];
		
	convolve(original, impulse, output_signal, output_size);
		
	int i;
	
	/*  Calculate the number of sound samples to create,
	 rounding upwards if necessary  */
	int numberOfSamples = output_size;
	
	
	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen(outputFile, "wb");
	
	/*  Write the WAVE file header  */
	writeWaveFileHeader(original->myChannels, numberOfSamples, original->myBitsPerSample,
						original->mySampleRate, outputFileStream);
	
	
	float maxValInResult = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if ( output_signal[i] > maxValInResult )
			maxValInResult = output_signal[i];

	float maxValInInput = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if ( original->my_signal[i] > maxValInInput )
			maxValInInput = original->my_signal[i];
			
	for (i = 0; i < numberOfSamples; i++ )
		fwriteShortLSB((short)(output_signal[i] / maxValInResult * maxValInInput), outputFileStream);
	
	
	/*  Close the output file stream  */
	fclose(outputFileStream);
	
}

/******************************************************************************
 *
 *       function:       writeWaveFileHeader
 *
 *       purpose:        Writes the header in WAVE format to the output file.
 *
 *		arguments:		channels:  the number of sound output channels
 *						numberSamples:  the number of sound samples
 *                       outputRate:  the sample rate
 *						outputFile:  the output file stream to write to
 *                       
 *       internal
 *       functions:      fwriteIntLSB, fwriteShortLSB
 *
 *       library
 *       functions:      ceil, fputs
 *
 ******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples, int bitsPerSample,
						 double sampleRate, FILE *outputFile)
{
	/*  Calculate the total number of bytes for the data chunk  */
	int dataChunkSize = channels * numberSamples * (bitsPerSample / 8);
	
	/*  Calculate the total number of bytes for the form size  */
	int formSize = 36 + dataChunkSize;
	
	/*  Calculate the total number of bytes per frame  */
	short int frameSize = channels * (bitsPerSample / 8);
	
	/*  Calculate the byte rate  */
	int bytesPerSecond = (int)ceil(sampleRate * frameSize);
	
	/*  Write header to file  */
	/*  Form container identifier  */
	fputs("RIFF", outputFile);
	
	/*  Form size  */
	fwriteIntLSB(formSize, outputFile);
	
	/*  Form container type  */
	fputs("WAVE", outputFile);
	
	/*  Format chunk identifier (Note: space after 't' needed)  */
	fputs("fmt ", outputFile);
	
	/*  Format chunk size (fixed at 16 bytes)  */
	fwriteIntLSB(16, outputFile);
	
	/*  Compression code:  1 = PCM  */
	fwriteShortLSB(1, outputFile);
	
	/*  Number of channels  */
	fwriteShortLSB((short)channels, outputFile);
	
	/*  Output Sample Rate  */
	fwriteIntLSB((int)sampleRate, outputFile);
	
	/*  Bytes per second  */
	fwriteIntLSB(bytesPerSecond, outputFile);
	
	/*  Block alignment (frame size)  */
	fwriteShortLSB(frameSize, outputFile);
	
	/*  Bits per sample  */
	fwriteShortLSB(bitsPerSample, outputFile);
	
	/*  Sound Data chunk identifier  */
	fputs("data", outputFile);
	
	/*  Chunk size  */
	fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
 *
 *       function:       fwriteIntLSB
 *
 *       purpose:        Writes a 4-byte integer to the file stream, starting
 *                       with the least significant byte (i.e. writes the int
 *                       in little-endian form).  This routine will work on both
 *                       big-endian and little-endian architectures.
 *
 *       internal
 *       functions:      none
 *
 *       library
 *       functions:      fwrite
 *
 ******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];
	
    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
 *
 *       function:       fwriteShortLSB
 *
 *       purpose:        Writes a 2-byte integer to the file stream, starting
 *                       with the least significant byte (i.e. writes the int
 *                       in little-endian form).  This routine will work on both
 *                       big-endian and little-endian architectures.
 *
 *       internal
 *       functions:      none
 *
 *       library
 *       functions:      fwrite
 *
 ******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];
	
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}


void convolve(float x[], int N, float h[], int M, float y[], int P) {
	int n, m;
	/*	Make sure the output buffer is the right size: P = N + M - 1	*/ 
	if (P != (N + M - 1)) {
		printf("Output signal vector is the wrong size\n"); printf("It is %-d, but should be %-d\n", P, (N + M - 1)); printf("Aborting convolution\n"); return;
	}
	
	/* Clear the output buffer y[] to all zero values */ 
	for (n = 0; n < P; n++)
		y[n] = 0.0;
	/* Do the convolution */ /* Outer loop: process each input value x[n] in turn */ 
		/* Inner loop: process x[n] with each sample of h[] */ 
	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++)
			y[n+m] += (x[n] * h[m]);
	}
}

void convolve(CWav *original, CWav *impulse, float y[], int P) {
	float *h = new float[impulse->mySignalSize];
	for ( int i = 0; i < impulse->mySignalSize; i++ )
		h[i] = (float)impulse->my_signal[i];
		
	float *x = new float[original->mySignalSize];
	for ( int i = 0; i < original->mySignalSize; i++ )
		x[i] = (float)original->my_signal[i];
	
	convolve(x, original->mySignalSize, h, impulse->mySignalSize, y, P);
}

void readInput(char *filename)
{
	CWav *waveObject = new CWav();
	waveObject->readInput(filename);
	
	cout << "Showing the first 100 samples in " << filename << endl;
	
	for ( int i = 0; i < waveObject->mySignalSize; i++ )
		cout <<  i << " : " << waveObject->my_signal[i] << endl;
	
}
