#include <fstream>
#include <string>
#include "stdlib.h"
using std::string;
using std::ifstream;


void printy(float *y, int n, char filename[])
{
	int ii;
	FILE *fid;

	fid = fopen(filename, "w");
		for (ii = 0; ii < n; ii++)
		fprintf(fid, "%3f ", y[ii]);
	fclose(fid);

}



//prints the matrix H to the file given by filename in the "alist" format
/*void printmatrix(struct Hmatrix H, char filename[])
{
	FILE *fid;
	int ii,jj;
	char buffer[7];

	fid = fopen(filename,"w");
	
	itoa(H.n,buffer,10);
	fprintf(fid,"%s ",buffer);
	itoa(H.m,buffer,10);
	fprintf(fid,"%s\n",buffer);

	//printing maximum degrees
	itoa(H.maxdegL, buffer, 10);
	fprintf(fid, "%s ", buffer);
	itoa(H.maxdegR, buffer, 10);
	fprintf(fid, "%s\n", buffer);

	//---- printing the degrees ----

	for (ii = 0; ii < H.n; ii++)
	{
		itoa(H.degL[ii], buffer, 10);
		if (ii < (H.n - 1))
			fprintf(fid, "%s ", buffer);
		else
			fprintf(fid, "%s\n", buffer);
	}

	
	itoa(H.degR[0], buffer, 10);
	fprintf(fid, "%s", buffer);

	for (ii = 1; ii < H.m; ii++)
	{
		itoa(H.degR[ii], buffer, 10);
		fprintf(fid, " %s", buffer);
	}
	
	fprintf(fid, "\n");

	//------------- printing the edges --------------------
	//from the cw perspective
	for (ii = 0; ii < H.n; ii++)
	{
		
		
		for (jj = 0; jj < H.degL[ii]; jj++)
		{
			itoa(H.edgesL[ii][jj] + 1, buffer, 10);
			if (jj == 0)
				fprintf(fid, "%s", buffer);
			else
				fprintf(fid, " %s", buffer);
		}

		for (jj = H.degL[ii]; jj < H.maxdegL; jj++)
		{
			fprintf(fid, " 0"); 
		}
		fprintf(fid, "\n");

	}

	//from the bit node perspective

	for (ii = 0; ii < H.m; ii++)
	{
		for (jj = 0; jj < H.degR[ii]; jj++)
		{
			itoa(H.edgesR[ii][jj] + 1, buffer, 10);
			if (jj == 0)
				fprintf(fid, "%s",buffer);
			else
				fprintf(fid, " %s", buffer);
		}


		for (jj = H.degR[ii]; jj < H.maxdegR; jj++)
		{
			fprintf(fid, " 0"); 
		}
		fprintf(fid, "\n");

	}

	fclose(fid);
}
*/
//parses the buffer until the a ' ' or a '\n' character and returns the integer
int seekint(char *buffer, long *countb)
{
	char string[7];
	int counts = 0;

	strcpy(string,"\n");
	
	//skipping whitespace if any
	while (buffer[*countb] == ' ' || buffer[*countb] == '\n' || buffer[*countb] == 13)
		(*countb)++;
	



	while (buffer[*countb] != ' ' && buffer[*countb] != '\n' && buffer[*countb] != 9 )
	{
		string[counts] = buffer[*countb];
		(*countb)++;
		counts++;
	}
	(*countb)++;
	return (atoi(string));
}

//READS PARITY CHECK MATRIX IN A-LIST FORMAT FROM FILENAME AND RETURNS THE MATRIX
/*struct Hmatrix readfile(char filename[])
{
	struct Hmatrix G;
	FILE * fid;
	char *buffer;
	long size;
	long pointer = 0;
	int ii, jj;
	

	fid = fopen(filename,"rb");
	if (fid == NULL)
	{
		printf("File %s not found",filename);
		exit(0);
	}

	//obtain file size
	fseek(fid, 0, SEEK_END);
	size = ftell(fid);
	rewind(fid);

	//allocate memory to contain the while file
	
	buffer = new char [size];
	if (buffer == NULL)
	{
		printf("Not enough memory to allocate buffer\n");
		exit(1);
	}

	fread(buffer,1,size,fid);
	fclose(fid);
	//free(buffer);
	
	//parsing the buffer and putting into the G structure
	
	//-------------- reading the size, n and k
	G.n = seekint(buffer, &pointer);
	G.m = seekint(buffer, &pointer);

	G.maxdegL = seekint(buffer, &pointer);
	G.maxdegR = seekint(buffer, &pointer);

	//allocating memory for the degree vectors 
	G.degL = new int [G.n];
	G.degR = new int [G.m];

	//reading the degree of the nodes
	for (ii = 0; ii < G.n; ii++)
	{
		G.degL[ii] = seekint(buffer, &pointer);
	}

	for (ii = 0; ii < G.m; ii++)
	{
		G.degR[ii] = seekint(buffer, &pointer);
	}

	
	//allocating memory for the edge connections matrices
	G.edgesL = new int * [G.n];
	G.edgesR = new int * [G.m];

	for (ii=0; ii < G.n; ii++)
		G.edgesL[ii] = new int [G.degL[ii]];

	for (ii=0; ii < G.m; ii++)
		G.edgesR[ii] = new int [G.degR[ii]];

	// Reading edges from the left nodes perspective
	for (ii = 0; ii < G.n; ii++)
	{
		for (jj = 0; jj < G.degL[ii]; jj++)
		{
			G.edgesL[ii][jj] = seekint(buffer, &pointer) - 1;
		}

		while (jj < G.maxdegL)
		{
			seekint(buffer, &pointer);
			jj++;
		}
	}


	// Reading edges from the information nodes perspective
	for (ii = 0; ii < G.m; ii++)
	{
		for (jj = 0; jj < G.degR[ii]; jj++)
		{
			G.edgesR[ii][jj] = seekint(buffer, &pointer) - 1;
		}

		while (jj < G.maxdegR)
		{
			seekint(buffer, &pointer);
			jj++;
		}
	}


	delete(buffer);


	return G;

}
*/
double read_degdistLDPC(string filename,
 double **lamda1, 
 double **lamda2, 
 double **rho, 
 int **degL1, 
 int **degL2,
 int **degR, 
 int *ndegL1,
 int *ndegL2,
 int *ndegR)
{
	
	char * fn=(char *)filename.c_str();
	int dg;
	double tmp, R;
	int maxdegL, maxdegR, nzdegL1, nzdegL2, nzdegR;
	int count=0;
	string par;
	ifstream in(fn);
	string rd;
	in>>maxdegR;	
	in>>maxdegL;
	in>>nzdegR;
	in>>nzdegL1;
	in>>nzdegL2;
	*ndegL1 = nzdegL1;
	*ndegL2 = nzdegL2;
    *ndegR = nzdegR;
	
	lamda1[0] = new double[nzdegL1];
    lamda2[0] = new double[nzdegL2];
	degL1[0] = new int [nzdegL1];
	degL2[0] = new int [nzdegL2];
	rho[0] = new double[nzdegR];
	degR[0] = new int [nzdegR];

	while(1)
	{
		in>>dg>>tmp;
		if(dg>0)
		{
			lamda1[0][count]=tmp;
			degL1[0][count] = dg;
			count++;
		}
		else
			break;
	}
	count = 0;
	
	while(1)
	{
		in>>dg>>tmp;
		if(dg>0)
		{
			lamda2[0][count]=tmp;
			degL2[0][count] = dg;
			count++;
		}
		else
			break;
	}
	count = 0;
	
    while(1)
	{
		in>>dg>>tmp;
		if(dg>0)
		{
			rho[0][count]=tmp;
			degR[0][count] = dg;
			count++;
		}
		else
			break;
	}
	in>>R;
	in.close();
	return R;
	
}

//open code profile
double read_degdist(string filename, double **lamda, double **rho, int **degL, int **degR, int *ndegL, int *ndegR)
{
	
	char * fn=(char *)filename.c_str();
	int dg;
	double tmp, R;
	int maxdegL, maxdegR, nzdegL, nzdegR;
	int count=0;
	string par;
	ifstream in(fn);
	string rd;
	in>>maxdegR;	
	in>>maxdegL;
	in>>nzdegR;
	in>>nzdegL;
	*ndegL = nzdegL;
	*ndegR = nzdegR;
	
	lamda[0] = new double[nzdegL];
	degL[0] = new int [nzdegL];
	rho[0] = new double[nzdegR];
	degR[0] = new int [nzdegR];

	while(1)
	{
		in>>dg>>tmp;
		if(dg>0)
		{
			lamda[0][count]=tmp;
			degL[0][count] = dg;
			count++;
		}
		else
			break;
	}
	count = 0;
	while(1)
	{
		in>>dg>>tmp;
		if(dg>0)
		{
			rho[0][count]=tmp;
			degR[0][count] = dg;
			count++;
		}
		else
			break;
	}
	in>>R;
	in.close();
	return R;
	
}
