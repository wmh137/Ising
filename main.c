#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
int main(int argc, char** argv) {
	int NumProcs;
	int MyID;
	char direction[2] = "+-";
	//setting begin
	FILE* ini;
	ini=fopen("mTH.ini", "rb");
	if (!ini) return 0;
	unsigned int issize;
	fread(&issize, sizeof(unsigned int), 1, ini);
	const unsigned int size = issize * issize;
	unsigned int step;
	fread(&step, sizeof(unsigned int), 1, ini);
	printf("%d\n", step);
	unsigned int Tlen;
	fread(&Tlen, sizeof(unsigned int), 1, ini);
	unsigned int Hlen;
	fread(&Hlen, sizeof(unsigned int), 1, ini);
	float T;
	float* TL = (float*)malloc(Tlen * sizeof(float));
	fread(TL, sizeof(float), Tlen, ini);
	float H;
	float* HL = (float*)malloc(Hlen * sizeof(float));
	fread(HL, sizeof(float), Hlen, ini);
	//setting end
	unsigned int i, j, ind;
	unsigned int* up = (unsigned int*)malloc(size * sizeof(unsigned int));
	unsigned int* down = (unsigned int*)malloc(size * sizeof(unsigned int));
	unsigned int* left = (unsigned int*)malloc(size * sizeof(unsigned int));
	unsigned int* right = (unsigned int*)malloc(size * sizeof(unsigned int));
	for (i = 0; i < issize; i++) {
		for (j = 0; j < issize; j++) {
			ind = i + j * issize;
			up[ind] = i + issize * ((j + issize - 1) % issize);
			down[ind] = i + issize * ((j + 1) % issize);
			left[ind] = (i + issize - 1) % issize + j * issize;
			right[ind] = (i + 1) % issize + j * issize;
		}
	}
	clock_t beg_t, end_t;
	float t;
	//parallel begin//
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	srand(1000);
	rand();
	char* ising = (char*)malloc(size * sizeof(char));
	unsigned int Ti, Hi;
	FILE* moutput, * output;
	float tp;
	float dE;
	int* count = (int*)malloc(Hlen * sizeof(int));
	if (!count) return 0;
	int* counti = (int*)malloc(Hlen * step * sizeof(int));
	if (!counti) return 0;
	for ( Ti = MyID; Ti < 2 * Tlen; Ti += NumProcs) {
		T = TL[Ti / 2];
		printf("Process%d T=%.5f%c\n", MyID, T, direction[Ti % 2]);
		fflush(stdout);
		char fname[64], fnamei[64], mfname[64];
		sprintf(fname, "%d//cdat//T-%.5f%c.dat", issize, T, direction[Ti % 2]);
		sprintf(fnamei, "%d//cdat//T-%.5f%c_i.dat", issize, T, direction[Ti % 2]);
		sprintf(mfname, "%d//Mat//MT-%.5f%c.dat", issize, T, direction[Ti % 2]);
		*strchr(fname, '.') = '_';
		*strchr(fnamei, '.') = '_';
		*strchr(mfname, '.') = '_';
		for (i = 0; i < size; i++) ising[i] = rand() % 2 * 2 - 1;
		memset(count, 0, Hlen * sizeof(int));
		memset(counti, 0, Hlen * step * sizeof(int));
		//calculate begin//
		beg_t = clock();
		for (Hi = 0; Hi < Hlen; Hi++) {
			if (Ti % 2) H = HL[Hi];
			else H = HL[Hlen - Hi - 1];
			printf("Process%d T=%.5f H=%.5f%c\n", MyID, T, H, direction[Ti % 2]);
			fflush(stdout);
			for (i = 0; i < step; i++) {
				for (ind = 0; ind < size; ind++) {
					dE = (ising[up[ind]] + ising[down[ind]] + ising[right[ind]] + ising[left[ind]] + H) * ising[ind] * 2;
					if (dE < 0)ising[ind] = -ising[ind];
					else {
						tp = exp(-dE / T);
						if ((float)rand() / RAND_MAX < tp) ising[ind] = -ising[ind];
					}
					
				}
			}
			for (ind = 0; ind < size; ind++) {
				if (Ti % 2)count[Hi] += ising[ind];
				else count[Hlen - Hi - 1] += ising[ind];
			}
		}
		end_t = clock();
		t = (float)(end_t - beg_t) / CLOCKS_PER_SEC;
		printf("Process%d T=%.5f calculation complete cost time=%fs\n", MyID, T, t);
		fflush(stdout);
		//calculation end//
		//output begin//
		output = fopen(fname, "wb");
		fwrite(&Hlen, sizeof(unsigned int), 1, output);
		fwrite(HL, sizeof(float), Hlen, output);
		fwrite(count, sizeof(int), Hlen, output);
		fclose(output);
		output = fopen(fnamei, "wb");
		fwrite(counti, sizeof(unsigned int), Hlen * step, output);
		fclose(output);
		moutput = fopen(mfname, "wb");
		fwrite(ising, sizeof(char), size, moutput);
		fclose(moutput);
		printf("Process%d T=%.5f output complete\n", MyID, T);
		fflush(stdout);
		//output end//
	}
	free(ising);
	MPI_Finalize();
	//parallel end//
	return 0;
}