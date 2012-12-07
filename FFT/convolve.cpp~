#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "CWav.h"

using namespace std;

#define DEBUG_MODE

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

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr


void outputToFile(double data[], CWav* original , int numberOfSamples, char* outputFile);
void four1(double data[], int nn, int isign);
void writeWaveFileHeader(int channels, int numberSamples, int bitsPerSample, double sampleRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
void getSignal(CWav *input, double x[]);
size_t fwriteShortLSB(short int data, FILE *stream);
void outputToFile(double data[], CWav* original , int numberOfSamples, char* outputFile);


//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

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

void getSignal(CWav *input, double x[]) {
	for ( int i = 0; i < input->mySignalSize; i++ ) 
		x[i] = ((double)input->my_signal[i])/32678.0;
}

void outputToFile(double output_signal[], CWav* original , int numberOfSamples, char* outputFile) {

	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen(outputFile, "wb");
	
	/*  Write the WAVE file header  */
	writeWaveFileHeader(original->myChannels, numberOfSamples, original->myBitsPerSample,
						original->mySampleRate, outputFileStream);
	
	int i;
	float maxValInResult = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if ( output_signal[i] > maxValInResult )
			maxValInResult = output_signal[i];

	float maxValInInput = -1;
	for (i = 0; i < numberOfSamples; i++ )
		if (original->my_signal[i] > maxValInInput )
			maxValInInput = original->my_signal[i];
			
	for (i = 0; i < numberOfSamples; i++ )
		fwriteShortLSB((short)(output_signal[i] / maxValInResult * maxValInInput), outputFileStream);
	
	
	/*  Close the output file stream  */
	fclose(outputFileStream);
}
//convolve inputfile IRfile outputfile as inputs to main
int main(int argc, char* argv[])
{
	if (argc != 4) {
		cout << "USAGE: ./convolve inputfile IRfile outputfile" << endl;
		return 0;
	}
		cout << argv[0] << endl;
		cout << "Input File:" << argv[1] << endl;
		cout << "IRfile:" << argv[2] << endl;
		cout << "outputfile File:" <<argv[3] << endl;


	char *outputFilename = argv[3];
	
	/*  Create the sine wave test tone, using the specified
	 frequency, duration, and number of channels, writing to
	 a .wav file with the specified output filename  */
	
	//createTestTone(FREQUENCY, DURATION, MONOPHONIC, BITS_PER_SAMPLE, SAMPLE_RATE, outputFilename);
		
	CWav *inputSignal = new CWav();
	inputSignal->readInput(argv[1]);
	
	//manipulate(inputSignal, 2);
	
	CWav *impulse = new CWav();
	impulse->readInput(argv[2]);
		
	cout << "Input Signal: " << inputSignal->mySignalSize << ", Impulse Size: " << impulse->mySignalSize << endl;
	
	double h[impulse->mySignalSize];
	double x[inputSignal->mySignalSize];

	getSignal(impulse, h);	
	getSignal(inputSignal, x);	
	

	int sizeH = impulse->mySignalSize;
	int sizeX = inputSignal->mySignalSize;

	cout << "SIZES(H,X)" << endl;
	cout << sizeH << endl;
	cout << sizeX << endl;
	              
	int maxSize;
	if(sizeH >= sizeX) {
		maxSize = sizeH;
	}
	else {
		maxSize = sizeX;
	}
	cout << "maxSize: " << maxSize << endl;

	int power = 0;	
	int pow2 = 0;
	
	pow2 = (int) log2(maxSize) + 1;
	pow2 = pow(2,pow2);

	cout << "POW :" << pow2 << endl;

	int i = 0;

	//set hComplex with 0's
	double hComplex[2 * pow2];
	for(i = 0; i < 2 * pow2; i++) {
		hComplex[i] = 0.0;			
	}	

	//set xComplex with 0's

	double *xComplex = new double[2 * pow2];
 	
	for(i = 0; i < 2 * pow2; i++) {
		xComplex[i] = 0.0;
	}		
	
	//padding the complex number with 0 and the real number with original value for h
	for(i = 0; i < sizeH; i++) {
		hComplex[2*i] = h[i];
	}

	//padding the complex number with 0 and the real number with original value for x
	for(i = 0; i < sizeX; i++) {
		xComplex[2*i] = x[i];
	}


	four1(hComplex, pow2, 1);
	four1(xComplex, pow2, 1);

	double *yComplex = new double[2 * pow2];		
	for(i = 0; i < pow2 ; i++) {
		yComplex[i*2] = xComplex[i] * hComplex[i] - xComplex[i+1] * hComplex[i+1];
		yComplex[i*2+1] = xComplex[i+1] * hComplex[i] + xComplex[i] * hComplex[i+1];
	}

	four1(yComplex-1, pow2, -1);
	outputToFile(yComplex, inputSignal, pow2, outputFilename);
}
