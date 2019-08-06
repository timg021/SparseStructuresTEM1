// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fftw3.h>
#include <XArray3D.h>

int main()
{
    std::cout << "Hello World!\n";
	xar::XArray3D<double> aaa(4, 4, 4, 1.0);


	fftw_plan_dft_r2c_3d(2, 2, 2, &aaa[0][0][0], &aaa[0][0][0], FFTW_ESTIMATE);

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
