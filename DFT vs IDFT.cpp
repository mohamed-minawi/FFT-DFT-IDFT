#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#define DFTpoints 8000
#define pi 3.14159265

using namespace std;
struct polar
{
	float magnitude, phase;
};

struct rectangle {
	float real, imaj;
};

void calculateresult(float arrayx[]);
void rectangulardft(float arrayx[], rectangle arr[]);
void polardft(rectangle arr[], polar arrx[]);
void inversedft(float[], rectangle[]);
void readaudiofile(float[]);
void newaudiofile(float[]);

int main()
{
	float wave[DFTpoints], resultinversedft[DFTpoints];
	polar polararray[DFTpoints];
	rectangle rectanglearray[DFTpoints];
	if (DFTpoints == 256 || DFTpoints == 128)
	{
		cout << "normal" << endl;
		calculateresult(wave);
		rectangulardft(wave, rectanglearray);
		polardft(rectanglearray, polararray);
		inversedft(resultinversedft, rectanglearray);
	}
	else if (DFTpoints == 8000)
	{
		cout << "audio" << endl;
		readaudiofile(wave);
		rectangulardft(wave, rectanglearray);
		polardft(rectanglearray, polararray);
		inversedft(resultinversedft, rectanglearray);
		newaudiofile(resultinversedft);
	}
	
	system("pause");
	return 0;
}

void calculateresult(float arrayx[])
{
	ofstream output("SamplesDFTfunction"+to_string(DFTpoints)+".txt");
	for (int i = 0; i < DFTpoints; i++)
		arrayx[i] = (sin((i * pi* 2.0) / 40.0) + 2 * sin((i *2.0*pi) / 16.0))*exp(-1 * pow(((i - 128) / 64.0), 2.0)),output << arrayx[i] << endl;
	
	output.close();
}
void rectangulardft(float arrayx[], rectangle arr[])
{
	ofstream routput("Real" + to_string(DFTpoints) + ".txt");
	ofstream ioutput("Imaj" + to_string(DFTpoints) + ".txt");

	for (int i = 0; i < DFTpoints; i++)
	{
		arr[i].imaj = 0;
		arr[i].real = 0;
		for (int j = 0; j < DFTpoints; j++)
		{
			arr[i].real += arrayx[j] * cos(2 * pi * i * j / DFTpoints);
			arr[i].imaj -= arrayx[j] * sin(2 * pi * i * j / DFTpoints);
		}
		
		routput << arr[i].real << endl;
		ioutput << arr[i].imaj << endl;
	} 
	routput.close();
	ioutput.close();
}
void polardft(rectangle arr[], polar arrx[])
{
	ofstream moutput("Magnitude" + to_string(DFTpoints) + ".txt");
	ofstream poutput("Phase" + to_string(DFTpoints) + ".txt");

	for (int i = 0; i < DFTpoints; i++)
	{
		arrx[i].magnitude = pow((pow(arr[i].imaj,2) + pow(arr[i].real, 2)), 0.5);
		arrx[i].phase = atan(arr[i].imaj/arr[i].real) * 180.0/pi;
		moutput << arrx[i].magnitude << endl;
		poutput << arrx[i].phase << endl;
	}
	moutput.close();
	poutput.close();
}
void inversedft(float idft[], rectangle dft[])
{
	rectangle* temp = new rectangle[DFTpoints];
	ofstream ioutput("InverseDFT" + to_string(DFTpoints) + ".txt");

	for (int i = 0; i < DFTpoints; i++)
	{
		temp[i].imaj = dft[i].imaj / DFTpoints;
		temp[i].real = dft[i].real / DFTpoints;
		idft[i] = 0;
	}

	for (int i = 0; i < DFTpoints; i++)
	{
		for (int j = 0; j < DFTpoints; j++)
		{
			idft[i] += temp[j].real * cos(2 * pi * i * j / DFTpoints);
			idft[i] -= temp[j].imaj * sin(2 * pi * i * j / DFTpoints);
		}
		
		ioutput << idft[i] << endl;
	}
	ioutput.close();
	delete[] temp;
}
void readaudiofile(float arrayx[])
{
	string x = "audiodata.txt";
	ifstream inputdata(x);
	float temp;
	
	for (int i = 0; i < 8000; i++)
	{
		inputdata >> temp;
		arrayx[i] = temp;
	}
	inputdata.close();

}
void newaudiofile(float arrayx[])
{
	string x = "newaudiodata.bin";
	ofstream newaudio(x, ios::out | ios::trunc | ios::binary);

	for (int i = 0; i < DFTpoints; i++)
		newaudio << char((arrayx[i]*128));
	
	newaudio.close();

}
