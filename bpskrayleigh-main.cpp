
// Study : Point - to - Point communication over an AWGN channel with Rayleigh fading

// Source modulation scheme : BPSK or 2-PAM

// Reads: Degree Distribution of LDPC codes from .dat file

// Code by:
// Aimal Rehman
// COMSATS IIT Lahore

//___________________________________________________________________________________//


#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include "stdlib.h"
#include "Random_Fns.h"
#include "LDPCfnsBPSK.h"
#include "filehandlingfns.h"
#include <limits.h>
#include <stdio.h>


using std::string;

// Nomrally DIstributed Random Numbers
float drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

float random_normal()
 /* normal distribution, centered on 0, std dev 0.707 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

//GENERATE BINARY DATA OF LENGTH N AND RETURN THE POINTER
void gendata(char *data, int N)
{
	int ii;

	for (ii=0; ii<N; ii++)
	{
		data[ii] = char (ran2(&idum)>=0.5)?1:0;			//equally likely data
		//data[ii] = 0;
	}

}

// BPSK MODULATION
void modulation(float *mod, int n, float c1, char *data)
{
	int ii;
	for (ii=0; ii<n; ii++){
		mod[ii] = c1*(1 - 2*data[ii]); // c1 is the avg. source power constraint !
	}
}

void addnoise(float *y, int n, double c1, float *data1, float var)
{
	int i;
	float *z;
	float tempvar, tempvar2;

	z = new float [n];


	tempvar = 0.0;
	for (i=0; i<n; i++)
	{
		z[i] = randn(&idum);
		tempvar += z[i]*z[i];

	}
	tempvar = tempvar / n;


	for (i=0; i<n; i++)
	{
		z[i] = z[i] * sqrt(var/tempvar);		//sqrt(var / tempvar) makes the variance of z exactly equal to var
		tempvar2 += z[i]*z[i];
		y[i] = c1*(data1[i]);
		y[i] += z[i];

	}

	delete(z);
}

void genGains(float *Rgains, int n){

// generate channel gains : Rayleigh coefficients

int i;

  float rands[n];
  float randi[n];

  for (i=0; i<n; i++){
  rands[i] = 0.707*random_normal();
  randi[i] = 0.707*random_normal(); }

 // Generate two arrays for Gaussian dist. numbers

  for (i=0; i<n; i++){
	Rgains[i] = sqrt(rands[i]*rands[i] + randi[i]*randi[i]);
  }

}

void genZeros(char *X, int n){
	// function to generate all zero codeword !

	int ii;
	for (ii=0; ii<n; ii++){
		X[ii] = 0;
	}
}

void genhalfZeros(char *X, int n){
	// function to generate half- zero codeword !

	int ii;
	for (ii=0; ii<n/2; ii++){
		X[ii]=0;
	}
	for (ii=n/2; ii<n; ii++){
		X[ii] = 1;
	}
}



int main(void){

	int totalsims=40; // Total simulations

	//BER computation
	string filename;
	FILE *fid;
	fid = fopen("ber.txt","w"); // here Bit Error Rate against every simulation is stored


	// ----Initialization ---- //
	char *X, *synd;
	float *Xmod;
	float *Yd;
	float *Gains;
	float *BER, BERavg;
	float erronous, nn;

	//degree distribution parameters for LDPC code
	double *lamda, *rho;
	int *degL, *degR;
	int ndegL, ndegR;
	double R;

	int n = 1000000; // number of bits in a codeword or block length

	int k, m;
	int ii, jj;
	int errors;
	int minerrors;
	int itergap;


	//----------------Reading the degree distributions-------------//

	R = read_degdist("Urbanke91.dat",   &lamda, &rho, &degL, &degR, &ndegL, &ndegR);

	k = (int) (n*R); // length of information bits in a block

	m = n - k;		// length of the parity bits part in a block

	//---------------------Channel parameters--------------------//

	float pathgain;

	float dB = 0.1; // input power for simulation

	float Power = pow(10, (dB/10)); // power in linear scale

	//Allocating memories

	X = new char [n];
	synd = new char [m];
	Xmod = new float[n];
	Yd = new float [n];	// noisy signal at receiver


	//-------------- Parity check matrix creation----------------//

	struct Hmatrix H;

	normalizeprob(lamda, ndegL);

	H = createHmatrix(nb, kb, lamda, degL, ndegL, rho, degR, ndegR, ELIM, 100, 0);

	removelength4cycles(H);

	init_random(&idum);

	// --------- Memory allocation ------------ //
	BER = new float[totalsims];
	Gains = new float[totalsims];

	genGains(Gains, totalsims);  // gains generated


	//------Decoding variables------//
	float **Lcb, **Lbc;  // SPA edges' messages
	float *Lch;	  			 // channel LLRs
	float *LLR;	  			 // bit nodes' LLRs

	//Allocting memories for the decoding variables
	Lcb = mallocLcb(H);
	Lbc = mallocLbc(H);

	LLR = new float [H.n]; // LLRs for variable nodes

	Lch = new float [H.n]; // channel LLRs

	itergap = 50;  // inner iterations of sum product algorithm !!!

	int sim;

	errors=0;
	int var = 1; // variance of the noise in channel
	//pathgain = 1; // Set this, if only AWGN channel is to be studied
	float loss = 1;

	//-------------------  Decoding ------------------------------ //
	for (sim=1; sim<=totalsims; sim++)
	{

		fprintf(fid,"--------------------Simulation Number : %d ---------------\n",sim);
		printf("\n ------ Simulation Number %d ------- \n", sim);

		//--------Generating random codeword-------------------------------//
		gendata(X, n);

		//--------------syndrome generation-------------------------------//
		calcsyndrome(synd, X, H);

		modulation(Xmod, n, sqrt(Power), X);

		//pathgain= Gains[sim-1]; // set this, if Rayleigh fading gain is needed

		addnoise(Yd, nb, loss*pathgain, Xmod, var);

		//--------- initializing check to node messages-----------------------//

		initLcb(Lcb, H);
		int burst;

		//----------------Decoding the LDPC code------------------------//

		printf("\n ---- Decoding LDPC Codes ---- \n");

		minerrors=n;

		calcLch(Lch, Yd, loss*pathgain*sqrt(Power), var, n);

		errors = counterrors(X, Lch, n);
		printf("\nErrors before decoding: %d\n", errors);

		//-------------- Outer decoding iterations---------------------//

		for (int loop=0; loop<20; loop++)
		{
		printf("\n\n\n ### OUTER ITERATION NUMBER: %d   ####\n\n", loop+1);

		fprintf(fid,"\n\n\n ### OUTER ITERATION NUMBER: %d   ####\n\n", loop+1);

		SPA(Lcb, Lbc, Lch, H, synd, itergap); //SPA: full-fledged SPA ; for approx, use SPAapprox

		Lcb2LLR(LLR, H, Lcb, Lch, INT);

		errors = counterrors(X, LLR, n);

		if (errors < minerrors)
		{
			minerrors = errors;
			burst = 0;
			}
		else
			burst++;

		if (burst>=3 && loop>=5)
		{
		printf(" \n Breaking outerloop ! \n");
		break;
	}

		printf("\nIter # %d: Errors: %d",(loop+1), errors); // in a single decoding iteration set, errors = errors
		fprintf(fid,"\nIter # %d: Errors= %d", (loop+1), errors);

		if (errors==0){
			break;
			}

		} // Outer decoding loop ends here !

		erronous = minerrors; // minimum errors decoded made possible

		nn=n; // just to convert to floating point for division purpose

		BER[sim-1] = erronous/nnb; // calculation of BER in this simulation

		printf("\n BER in %d th simulations = %lf\n", sim, BER[sim-1]);

	} // Simulation ends here !

		BERavg=0;

		for (int k=0; k<totalsims; k++){
			BERavg += BER[k];
		}

		float tsims=totalsims; // again, just convert to float

		printf("\n\nOverall Bit Error Rate: %lf", BERavg/(tsims));
		printf("\n Power (dB) : %lf", dB);

		fprintf(fid,"\n\n Bit Error Rate : %lf", BERavg/(tsims));
		fprintf(fid, "\n Power (dB) : %lf", dB);
		fclose(fid);

		return 0;

}
