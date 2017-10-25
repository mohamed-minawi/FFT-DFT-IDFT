#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <complex>
#include <valarray>
#include <vector>
#include <ctime>

using namespace std;

#define DFTpoints 128
#define SamplingRate 5000
#define pi 3.14159265

typedef complex<double> Complex;
typedef valarray<Complex> ComplexArray;

struct polarc
{
	float magnitude, phase;
};
struct rectangle {
	float real, imaj;
};

void DFT(vector<float> arrayx, rectangle arr[], polarc arrx[]);
void FFT(ComplexArray & x, polarc arrx[]);
void FFTrecursive(ComplexArray &);

void writetofile(polarc fftd[], polarc dftd[], float, float);
void signaldatasampling(vector<Complex> &, vector<float> &);
void signaldatasamplingwindow(vector<Complex> &, vector<float> &);

int main()
{
	// Time Variables
	clock_t dftbegin, dftend, fftbegin, fftend;

	// FFT DataStructures
	vector<Complex> dataFFT(DFTpoints);
	polarc polarFFT[DFTpoints];
		
	// DFT DataStructures
	vector<float> dataDFT(DFTpoints);
	rectangle rectDFT[DFTpoints]; 
	polarc polarDFT[DFTpoints];

	// Signal Sampling
	signaldatasampling(dataFFT, dataDFT);

	//FFT DataStructure
	ComplexArray x(&dataFFT[0], DFTpoints);

	// Fourier Transforms
	dftbegin = clock();
	DFT(dataDFT, rectDFT, polarDFT);
	dftend = clock();

	fftbegin = clock();
	FFT(x, polarFFT);
	fftend = clock();

	// Writing Data to Files
	writetofile(polarFFT, polarDFT, difftime(dftbegin, dftend)/CLOCKS_PER_SEC, difftime(fftbegin, fftend)/ CLOCKS_PER_SEC);

	system("pause");
	return 0;
}

void DFT(vector<float> data, rectangle arr[], polarc arrx[])
{
	for (int i = 0; i < DFTpoints; i++)
	{
		arr[i].imaj = 0;
		arr[i].real = 0;
		for (int j = 0; j < DFTpoints; j++)
		{
			arr[i].real += data[j] * cos(2 * pi * i * j / DFTpoints);
			arr[i].imaj -= data[j] * sin(2 * pi * i * j / DFTpoints);
		}

	}

	for (int i = 0; i < DFTpoints; i++)
	{
		arrx[i].magnitude = pow((pow(arr[i].imaj, 2) + pow(arr[i].real, 2)), 0.5);
		arrx[i].phase = atan(arr[i].imaj / arr[i].real) * 180.0 / pi;

	}
}
void FFT(ComplexArray & x, polarc arrx[])
{
	FFTrecursive(x);

	for (int i = 0; i < DFTpoints; i++)
	{
		arrx[i].magnitude = pow((pow(x[i].imag(), 2) + pow(x[i].real(), 2)), 0.5);
		arrx[i].phase = atan(x[i].imag() / x[i].real()) * 180.0 / pi;

	}
}
void FFTrecursive(ComplexArray & x)
{
	int size = x.size();
	if (size <= 1) return;

	ComplexArray evenarr = x[slice(0, size / 2, 2)];
	ComplexArray  oddarr = x[slice(1, size / 2, 2)];

	FFTrecursive(evenarr);
	FFTrecursive(oddarr);

	for (int i = 0; i < size / 2; ++i)
	{
		Complex W = polar(1.0, -2 * pi * i / size) * oddarr[i];
		x[i] = evenarr[i] + W;
		x[i + size / 2] = evenarr[i] - W;
	}
}

void writetofile(polarc fftd[], polarc dftd[], float DFTime, float FFTime)
{
	ofstream outputDFT(to_string(DFTpoints) + "DFTMagnitude.txt");
	ofstream outputFFT(to_string(DFTpoints) + "FFTMagnitude.txt");

	outputDFT << "The time for DFT is: " << DFTime << endl;
	outputFFT << "The time for FFT is: " << FFTime << endl;
	for (int i = 0; i < DFTpoints; ++i)
		outputDFT << dftd[i].magnitude << endl, outputFFT << fftd[i].magnitude << endl;

	outputDFT.close();
	outputFFT.close();
}
void signaldatasampling(vector<Complex> & fft, vector<float> & dft)
{
	ofstream output(to_string(SamplingRate) + "SamplingRate" + to_string(DFTpoints) + "DFTpoints.txt");
	for (int i = 0; i < DFTpoints; i++)
		fft[i] = dft[i] = sin(2000 * pi*i / SamplingRate) + 3 * sin(4000 * pi*i / SamplingRate), output << dft[i] << endl;
	output.close();
}
void signaldatasamplingwindow(vector<Complex> & fft, vector<float> & dft)
{
	ofstream output(to_string(SamplingRate) + "SamplingRate" + to_string(DFTpoints) + "DFTpointsWindow.txt");
	for (int i = 0; i < DFTpoints; i++)
		fft[i] = dft[i] = (sin(2000 * pi*i / SamplingRate) + 3 * sin(4000 * pi*i / SamplingRate))*(0.54 - 0.46*cos(2 * pi*i / DFTpoints)), output << dft[i] << endl;
	output.close();
}
