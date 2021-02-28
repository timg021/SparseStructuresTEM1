#pragma once
#include <string>
#include <vector>

void SaveXYZfile(std::string fileout, float ctblength, float ctblengthz, int natom, int* Znum, float* x, float* y, float* z, float* occ, float* wobble);
int AddIce(float iceThick, float ctblength, int natom, int** pZnum, float** px, float** py, float** pz, float** pocc, float** pwobble, float wobbleaver, unsigned long* piseed);

