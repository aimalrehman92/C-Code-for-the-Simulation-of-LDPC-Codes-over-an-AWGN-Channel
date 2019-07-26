#define INF 10.0
#define PRUNE 333
#define ELIM 444
#define SYST 111
#define IRA 222
#define EXT 555
#define INT 666

#include <fstream>

// This code generates the structrure of Hmatrix  (Wifi)

//STRUCTURE FOR THE PARITY CHECK MATRIX
struct Hmatrix
{
	int n, m;		//n = number of left nodes, m = no of right nodes

	int maxdegL, maxdegR;
	
	int **edgesL;		//connections from the left nodes perspective
	int **edgesR;		//connections from the right nodes perspective
	
	int *degL;			//vector of the degree of information nodes
	int *degR;		    //vector of the degree of codeword nodes

};


void normalizeprob
(
double prob[],		//the probability mass function to be normalized
int ldeg			//length of the vector
)
{
	int ii;
	double sum = 0.0;

	for (ii = 0; ii < ldeg; ii++)
	{
		sum+= prob[ii];
	}

	for (ii = 0; ii < ldeg; ii++)
	{
		prob[ii] = prob[ii] / sum;
	}
}


//CONVERTS THE DISTRIBUTION FROM THE EDGE PERSPECTIVE TO THE NODE PERSPECTIVE
void edge2node
(
double probout[],		//Output probability distribution from the node perspective
double probin[],		//Input probability distribution from the edge perspective
int deg[],				//The degrees			
int ldeg				//length of the vectors
)
{
	int ii;
	for (ii = 0; ii < ldeg; ii++)
	{
		probout[ii] = probin[ii]/deg[ii];
	}
	normalizeprob(probout, ldeg);
}



// PARTITION THE COLUMNS ACCORDING TO THE SPECIFIED PROPORTIONS.  It
//   may not be possible to do this exactly.  Returns a pointer to an
//  array of integers containing the numbers of columns corresponding 
//   to the entries in the distribution passed. 

int *column_partition
( double prob[],	//containing the proportion
  int ldeg,		//length of degree distribution
  int n			// Total number of columns to partition 
)

{
  int *part;
  int cur, used;
  int i, j;
  double *bias;

  part = new int [ldeg];
  bias = new double [ldeg];

  used = 0;
  for (i = 0; i<ldeg; i++)
  { cur = floor(prob[i]*n);
    part[i] = cur; 
    bias[i] = bias[i] + prob[i]*n - cur; 
    used += cur; 
  }

  if (used>n) 
  { printf("\nerror in column partition\n");
  }
  
  while (used<n)
  { cur = 0;
    for (j = 1; j<ldeg; j++) 
    { if (bias[j]>bias[cur])
      { cur = j;
      }
    }
    part[cur] += 1;
    used += 1;
    bias[cur] = bias[cur] -1.0;
  }

  delete (bias);
  return part;
}

//CALCULATES TOTAL NUMBER OF EDGES
int calcnoofedges
(int *noofnodes,			//nuumber of nodes with degree from the vector deg
 int deg[],
 int length					//length of the vectors
 )
{
	int ii;
	int totaledges=0;

	for (ii = 0; ii < length; ii++)
	{
		totaledges += noofnodes[ii]*deg[ii];
	}
	return totaledges;
}

//GENERATES RANDOM INTGER PERMUTATIONS
//generates N distinct random integers between 0 and max-1 and returns the result; if N = max, returns a radnom permutation
int *genperm(int max, int N)
{
	int *reqarray;
	int *intarray;
	int ii, lowbound;
	int randindx;
	int temp;

	intarray = new int [max];
	reqarray = new int[N];

	for (ii = 0; ii < max; ii++)
		intarray[ii] = ii;

	for (lowbound = 0; lowbound < N; lowbound++)
	{
		randindx = lowbound + floor((double)(max - 0.01 - lowbound)*ran2(&idum));		//generate a random number between lowbound and max-1
		
		//exchange tempnumbers[randindx] and tempnumbers[lowbound]
		temp = intarray[lowbound];
		intarray[lowbound] = intarray[randindx];
		intarray[randindx] = temp;
	}

	for (ii = 0; ii < N; ii++)
		reqarray[ii] = intarray[ii];

	delete (intarray);
	return reqarray;

}

//SORT IN ASCENDING ORDER
void sortascending(int *numbers, int length)
{
	int ii, jj;
	int min, indexmin;
	int temp;

	for (ii = 0; ii < length; ii++)
	{
		min = numbers[ii];
		indexmin = ii;
		for (jj = ii+1; jj < length; jj++)
		{
			if (numbers[jj] < min)
			{
				min = numbers[jj];
				indexmin = jj;
			}

		}

		temp = numbers[ii];
		numbers[ii] = min;
		numbers[indexmin] = temp;
	}

}


//CREATES PARITY CHECK MATRIX H AND RETURNS IT
struct Hmatrix createHmatrix
(
int n,				
int k,
double lamda[],			//degree distribution (from the edge perspective) of the left nodes
int degL[],			//indices of the nonzero elements of L(x) 
int ndegL,			//number of non-zero entries in L(x)
double rho[],			//  |Same as above
int degR[],			//  |for the right nodes 
int ndegR,			//  |
int flagcycle2,	    // if flag == PRUNE, removes the length 2 cycles by pruning the edges, else if flag == ELIM, removes them by rearranging edges	
int maxpass1,		//maximum number of passes for the search of length-4 cycle elimination (within the existing edge distribution)
int maxpass2		//maximum number of passes after the search within the existing edge distribution has been exhausted -- adds new edges to random nodes
)
{
	struct Hmatrix H;
	int *noofnodesL, *noofnodesR;
	int noofedgesL, noofedgesR;
	int ii, jj, temp;
	int edgecount, counte;
	int *randomedges, *edgeleftnodenumber;
	int **edgerightnodenumber;
	int N;
	int *edgecountR, rightnode;
	int isbefore, kk;
	int *perm;
	double *L, *R;

	L = new double [ndegL];
	R = new double [ndegR];



	edge2node(L,lamda,degL,ndegL);		//converting from edge perspective to node perspective 
	edge2node(R,rho,degR,ndegR);


	H.n = n;
	H.m = n-k;

	noofnodesL = column_partition(L, ndegL, H.n);		//finding the number of nodes with deg degL
	noofnodesR = column_partition(R, ndegR, H.m);		//-------------------------------------degR
	
	noofedgesL = calcnoofedges(noofnodesL, degL, ndegL);		//number of edges on the left
	noofedgesR = calcnoofedges(noofnodesR, degR, ndegR);		//number of edges on the right

	H.degL = new int[H.n];
	H.degR = new int[H.m];
	
	perm = genperm(H.n, H.n);
	edgecount = 0;
	//-----------------assigning the degrees to each nodes
	//left nodes
	for (ii = 0; ii < ndegL; ii++)
	{
		for (jj = 0; jj < noofnodesL[ii]; jj++)
		{
			H.degL[perm[edgecount]] = degL[ii];	//degrees assigned randomly
		//	H.degL[edgecount] = degL[ii];		//degrees in ascending order
			edgecount++;
		}
	}
	delete(perm);
	perm = genperm(H.m, H.m);

	edgecount = 0;
	//right nodes
	for (ii = 0; ii < ndegR; ii++)
	{
		for (jj = 0; jj < noofnodesR[ii]; jj++)
		{
			H.degR[perm[edgecount]] = degR[ii];
			edgecount++;
		}
	}
	//----------------------------------------------------
	
	
	//------------ making sure the number of edges on the left and right are equal --------------
	if (noofedgesL < noofedgesR)
	{
		//remove a total of N = (nofedgesR - noofedgesL) random edges from the right
		N = noofedgesR - noofedgesL;
		randomedges = genperm(noofedgesR, N);
		sortascending(randomedges, N);
		edgecount = 0;
		counte = 0;

		for (ii = 0; ii < H.m; ii++)
		{
			temp = H.degR[ii];
			for (jj = 0; jj < temp; jj++)
			{
				if (edgecount == randomedges[counte])
				{
					H.degR[ii]--;
					counte++;
				
				}
				edgecount++;
			}
		}

		edgecount = 0;
		for (ii = 0; ii < H.m; ii++)
		{
			edgecount+= H.degR[ii];
		}
		
		if (edgecount != noofedgesL)
			printf ("Number of edges on left not equal to on right\n");
		noofedgesR = noofedgesL;

	}
	
	else if (noofedgesR < noofedgesL)
	{
		//remove a total of N = (nofedgesL - noofedgesR) random edges from the left
		N = noofedgesL - noofedgesR;
		randomedges = genperm(noofedgesL, N);
		sortascending(randomedges, N);
		edgecount = 0;
		counte = 0;

		for (ii = 0; ii < H.n; ii++)
		{
			temp = H.degL[ii];
			for (jj = 0; jj < temp; jj++)
			{
				if (edgecount == randomedges[counte])
				{
					H.degL[ii]--;
					counte++;
				
				}
				edgecount++;
			}
		}

		edgecount = 0;
		for (ii = 0; ii < H.n; ii++)
		{
			edgecount+= H.degL[ii];
		}

		if (edgecount != noofedgesR)
			printf ("Number of edges on left not equal to on right\n");
		noofedgesL =  noofedgesR;
	}

	delete(randomedges);

	
	
	
	//----------------- Creating the actual graph ---------------------------
	
	randomedges = new int[noofedgesL];
	edgerightnodenumber = new int *[noofedgesL];
	edgeleftnodenumber = new int [noofedgesL];
	edgecountR = new int[H.m];
	//allocating memory for the edge connection vectors
	H.edgesL = new int* [H.n];
	H.edgesR = new int* [H.m];
	int *intarray;
	int randedge, randedgeindx;
	int iscycle, rn1, rn2, ll, mm, pass;
	int *edgetemp;

	intarray = new int [noofedgesL];

	for (ii = 0; ii < noofedgesL; ii++)
		intarray[ii] = ii;

	for (ii=0; ii<H.n; ii++)
	{
		H.edgesL[ii] = new int [H.degL[ii]];
	}
	
	for (ii=0; ii<H.m; ii++)
	{
		H.edgesR[ii] = new int [H.degR[ii]];
	}

	for (ii = 0; ii < noofedgesL; ii++)
	{
		edgerightnodenumber[ii] = new int [2];
	}
	
	//generating a random permutation of the edges
	randomedges = genperm(noofedgesL, noofedgesL);

	//associating each edge with its corresponding right node number
	edgecount = 0;
	for (ii = 0; ii < H.m; ii++)
	{
		for (jj = 0; jj < H.degR[ii]; jj++)
		{
			edgerightnodenumber[edgecount][0] = ii;
			edgerightnodenumber[edgecount][1] = jj;
			edgecount++;
		}
	
	}

	//associating each edge with its corresponding left node number
	edgecount = 0;
	for (ii = 0; ii < H.n; ii++)
	{
		for (jj = 0; jj < H.degL[ii]; jj++)
		{
			edgeleftnodenumber[edgecount] = ii;
			edgecount++;
		}
	
	}

	//initializing the right edge count to an all zero vector
	for (ii = 0; ii < H.m; ii++)
		edgecountR[ii] = 0;
	
	edgecount = 0;
	for (ii = 0; ii < H.n; ii++)
	{
		for (jj = 0; jj < H.degL[ii]; jj++)
		{
			
			for (pass = 0; pass < (maxpass1 + maxpass2); pass++)
			{
				iscycle = 0;
				if (pass < maxpass1)		//if within first round of pass
				{
					randedgeindx = edgecount + floor((double)(noofedgesL - 0.01 - edgecount)*ran2(&idum));		//generate a random number between edgecount and noofedgesL
					randedge = intarray[randedgeindx];
					rn1 = edgerightnodenumber[randedge][0];
				}

				else		//if after first round can't find any edge which eliminates the cycle
				{
					//pick a random node from the right
					rn1 = floor(double(H.m - 0.01)*ran2(&idum));
				}

				//check to see if length 2 cycles
				for (kk = 0; kk < jj; kk++)
				{
					if (H.edgesL[ii][kk] == edgerightnodenumber[randedge][0])
					{
						iscycle = 1;
						break;
					}
						
				}

				//check to see if length 4 cycle
				if (!iscycle)
				{
					for (kk = 0; kk < jj; kk++)
					{
						if (iscycle)
							break;
						rn2 = H.edgesL[ii][kk];

						for (ll = 0; ll < edgecountR[rn1]; ll++)
						{
							if (iscycle)
							break;
						
							for (mm = 0; mm < edgecountR[rn2]; mm++)
							{
								if (H.edgesR[rn1][ll] == H.edgesR[rn2][mm])	//length 4 cycle detected
								{
									iscycle = 1;
									break;
								}
							}
						}
					}

				}//end of if(!iscycle)

				if (!iscycle)
					break;
			}//end of pass loop

			
			//exchange
			if (pass < maxpass1 || pass == (maxpass1 + maxpass2))
			{
				if (pass == (maxpass1 + maxpass2))
				{
					rn1 = edgerightnodenumber[randedge][0];
				}

				temp = intarray[edgecount];
				intarray[edgecount] = intarray[randedgeindx];
				intarray[randedgeindx] = temp;

				H.edgesR[rn1][edgecountR[rn1]] = ii;
				edgecountR[rn1]++;
				edgecount++;
			}

			else
			{
				edgetemp = new int [edgecountR[rn1]];
				
				for (ll = 0; ll < edgecountR[rn1]; ll++)
				{
					edgetemp[ll] = H.edgesR[rn1][ll];
				}
				delete(H.edgesR[rn1]);
				H.degR[rn1]++;
				H.edgesR[rn1] = new int [H.degR[rn1]];

				for (ll = 0; ll < edgecountR[rn1]; ll++)
				{
					H.edgesR[rn1][ll] = edgetemp[ll];
				}
				H.edgesR[rn1][edgecountR[rn1]] = ii;
				edgecountR[rn1]++;
				delete(edgetemp);
				
			}
			
			H.edgesL[ii][jj] = rn1;
	

		}//end of jj loop
	sortascending(H.edgesL[ii], H.degL[ii]);
	}//end of ii loop

	for (ii = 0; ii < H.m; ii++)
	{
		if (edgecountR[ii] != H.degR[ii])
			H.degR[ii] = edgecountR[ii];
	}
	
	
	
	
	/*
	//Connecting the edges on the two sides and recording the corresponding node numbers
	edgecount = 0;		//edge count
	for (ii = 0; ii < H.n; ii++)
	{
		for (jj = 0; jj < H.degL[ii]; jj++)
		{
			rightnode = edgerightnodenumber[randomedges[edgecount]][0];
			
			H.edgesL[ii][jj] = rightnode;
			H.edgesR[rightnode][edgecountR[rightnode]] = ii;
			edgecountR[rightnode]++;
			edgecount++;
		}
		sortascending(H.edgesL[ii], H.degL[ii]);
	}

*/
	
	
	

	
	//----------------------------   trying to eliminate remaining cycles of length 2

	int newrightnode, indxnewleftnode, indxrightnode, newleftnode;
	int removededgecount = 0;

	for (ii = 0; ii < H.n; ii++)
	{
		
		for (jj = 0; jj < H.degL[ii]; jj++)
		{
			rightnode = H.edgesL[ii][jj];
			for (kk = 0; kk < jj; kk++)
			{
				
				if (H.edgesL[ii][kk] == rightnode)		//cycle present
				{
					if (flagcycle2 == ELIM)
					{
						for (pass = 0; pass < 10; pass++)	//for a maximum of 10 passes
						{
							iscycle = 0;
							randedge = floor(ran2(&idum)*(noofedgesL-0.001));
							newrightnode = edgerightnodenumber[randedge][0];

							//checking to verify that the new edge does not itself form a cycle

							//check if the new right node is already connected to the left node ii or not
							for (ll = 0; ll < H.degL[ii]; ll++)
							{
								if (H.edgesL[ii][ll] == newrightnode)
								{
									iscycle = 1;
									break;
								}
							}

							if (!iscycle)
							{
								//check to see if the left node of newrightnode at randedge is connected to rightnode or not
								newleftnode = H.edgesR[newrightnode][edgerightnodenumber[randedge][1]];

								for (ll = 0; ll < H.degL[newleftnode]; ll++)
								{
									if (H.edgesL[newleftnode][ll] == newrightnode)
									{
										indxnewleftnode = ll;
									}
								
									if (H.edgesL[newleftnode][ll] == rightnode)
									{
										iscycle = 1;
										break;
									}
						
								}

							}

							if (!iscycle)
							{
								for (ll = 0; ll < H.degR[rightnode]; ll++)
								{
									if (H.edgesR[rightnode][ll] == ii)
									{
										indxrightnode = ll;
										break;
									}
								}
								break;
							}
						}//end of for pass loop

					}//end of if flagcycle2==ELIM

					else		//if not prune make iscycle = 1 so that the next if jumps directly to else
					{
						iscycle = 1;
					}

					
					if (!iscycle)		//if cycle was removed exchange the edges
					{
						H.edgesL[ii][jj] = newrightnode;
						H.edgesR[newrightnode][edgerightnodenumber[randedge][1]] = ii;
						H.edgesL[newleftnode][indxnewleftnode] = rightnode;
						H.edgesR[rightnode][indxrightnode] = newleftnode;

						sortascending(H.edgesL[newleftnode], H.degL[newleftnode]);
						sortascending(H.edgesR[rightnode], H.degR[rightnode]);
						sortascending(H.edgesR[newrightnode], H.degR[newrightnode]);

					}

					else		//otherwise remove the edge
					{
						for (ll = jj; ll < H.degL[ii]-1; ll++)
						{
							H.edgesL[ii][ll] = H.edgesL[ii][ll+1];
						}
						H.degL[ii]--;

						
						for (ll = 0; ll < H.degR[rightnode]; ll++)
						{
							if (H.edgesR[rightnode][ll] == ii)
							{
								indxrightnode = ll+1;
								break;
							}
						}
						
						for (ll = indxrightnode; ll < H.degR[rightnode]-1; ll++)
						{
							H.edgesR[rightnode][ll] = H.edgesR[rightnode][ll+1];
						}
						H.degR[rightnode]--;
						removededgecount++;
						jj--;

					}

					


				}//end of if (H.edgesL[ii][kk] == rightnode)
			
			
			}//end of kk loop
		}//end of jj loop
	
		sortascending(H.edgesL[ii], H.degL[ii]);

	}

	
	
		//finding the maximum degrees
	H.maxdegL = H.degL[0];
	H.maxdegR = H.degR[0];

	for (ii = 0; ii < H.n; ii++)
	{
		if (H.degL[ii] > H.maxdegL)
			H.maxdegL = H.degL[ii];
	}

	for (ii = 0; ii < H.m; ii++)
	{
		if (H.degR[ii] > H.maxdegR)
			H.maxdegR = H.degR[ii];
	}	
		
		
	//deleting allocated memories
	delete (randomedges);
	delete(edgerightnodenumber);
	delete(edgeleftnodenumber);
	delete(edgecountR);
	delete(L);
	delete(R);

	return H;

}

int countlength4cycles(struct Hmatrix H)
{
	int ii, jj, kk, ll,mm;
	int rn1, rn2, ln1;
	int noofcycles = 0;
	


	int *edgecounter;
	edgecounter = new int [H.m];



	
	noofcycles = 0;
	for (ii = 0; ii < H.m; ii++)
	{
		edgecounter[ii] = 0;
	}
	
	for (ln1 = 0; ln1 < H.n; ln1++)
	{
	
		for (jj = 0; jj < H.degL[ln1]; jj++)
		{
			edgecounter[H.edgesL[ln1][jj]]++;
		}


		for (jj = 0; jj < H.degL[ln1]; jj++)
		{
			rn1 = H.edgesL[ln1][jj];
			for (kk = jj+1; kk < H.degL[ln1]; kk++)
			{
				rn2 = H.edgesL[ln1][kk];
			
				//if rn1 and rn2 have any other common left node (beside node ln1) then there's a length 4 cycle
				for (ll = edgecounter[rn1]; ll < H.degR[rn1]; ll++)
				{
					
					for (mm = edgecounter[rn2]; mm < H.degR[rn2]; mm++)
					{
							
						if (H.edgesR[rn1][ll] == H.edgesR[rn2][mm] && H.edgesL[rn1][ll]!=ln1 && H.edgesL[rn2][mm]!=ln1)		//length 4 cycle detected
						{
							noofcycles++;
							
						}
					
					}//end of mm loop
				}//end of ll loop

			}//end of kk loop
			
				
		}//end of jj loop
		

	}//end of ln1 loop

return noofcycles;


}


void deletematrix(struct Hmatrix H)
{
	int ii;

	for (ii = 0; ii < H.n; ii++)
	{
		delete(H.edgesL[ii]);
	}

	for (ii = 0; ii < H.m; ii++)
	{
		delete(H.edgesR[ii]);
	}
	delete(H.edgesL);
	delete(H.edgesR);
	delete(H.degL);
	delete(H.degR);
}


struct Hmatrix assignmatrix(struct Hmatrix Hin)
{
	struct Hmatrix Hout;
	int ii, jj;

	Hout.n = Hin.n;
	Hout.m = Hin.m;
	Hout.maxdegL = Hin.maxdegL;
	Hout.maxdegR = Hin.maxdegR;

	Hout.degL = new int [Hout.n];
	Hout.degR = new int [Hout.m];
	Hout.edgesL = new int * [Hout.n];
	Hout.edgesR = new int * [Hout.m];



	for (ii = 0; ii < Hout.n; ii++)
	{
		Hout.degL[ii] = Hin.degL[ii];
		Hout.edgesL[ii] = new int [Hout.degL[ii]];
	}

	for (ii = 0; ii < Hout.m; ii++)
	{
		Hout.degR[ii] = Hin.degR[ii];
		Hout.edgesR[ii] = new int [Hout.degR[ii]];
	}

	

	for (ii = 0; ii < Hout.n; ii++)
	{
		for (jj = 0; jj < Hout.degL[ii]; jj++)
		{
			Hout.edgesL[ii][jj] = Hin.edgesL[ii][jj];
		}
	}
	
	for (ii = 0; ii < Hout.m; ii++)
	{
		for (jj = 0; jj < Hout.degR[ii]; jj++)
		{
			Hout.edgesR[ii][jj] = Hin.edgesR[ii][jj];
		}
	}

	return Hout;
}


struct Hmatrix createcyclefreeH(int n, int k, double *lamda, int *degL, int ndegL, 
								double *rho, int *degR, int ndegR)
{
	struct Hmatrix Htemp;
	struct Hmatrix H;
	int noofcycles, minnoofcycles, ii;

	Htemp = createHmatrix(n,k,lamda,degL,ndegL,rho,degR, ndegR, ELIM, 1, 0);	//create a random parity check matrix
	
	H = assignmatrix(Htemp);
	noofcycles = countlength4cycles(Htemp);
	minnoofcycles = noofcycles;
	
	for (ii = 0; ii < 0; ii++)
	{
		if (minnoofcycles == 0)
			break;

		printf("%d\n",ii);
		deletematrix(Htemp);
		Htemp = createHmatrix(n,k,lamda,degL,ndegL,rho,degR, ndegR, ELIM, 100, 0);	//create a random parity check matrix
		noofcycles = countlength4cycles(Htemp);

		if (noofcycles < minnoofcycles)
		{
			deletematrix(H);
			minnoofcycles = noofcycles;
			H = assignmatrix(Htemp);
		}

	}
	return H;
}


void removelength4cycles(struct Hmatrix H)
{
	int ii, jj, kk, ll,mm;
	int rn1, rn2, ln1;
	int noofcycles = 0, cumnoofcycles = 0, pass;
	int noofedges=0, edgecount, randedge, newrightnode, newleftnode, indxnewleftnode, iscycle2, pass2;
	int **edgerightnodenumber;

	int *edgecounter;
	edgecounter = new int [H.m];


	for (ii = 0; ii < H.n; ii++)
	{
		for (jj = 0; jj < H.degL[ii]; jj++)
			noofedges++;
	}

	edgerightnodenumber = new int*[noofedges];

	for (ii = 0; ii < noofedges; ii++)
	{
		edgerightnodenumber[ii] = new int [2];
	}
	

	//assigining the node numbers to the edges
	edgecount = 0;
	for (ii = 0; ii < H.m; ii++)
	{
		for (jj = 0; jj < H.degR[ii]; jj++)
		{
			edgerightnodenumber[edgecount][0] = ii;
			edgerightnodenumber[edgecount][1] = jj;
			edgecount++;
		}
	
	}


	for (pass = 0; pass < 10; pass++)
	{
		noofcycles = 0;
		for (ii = 0; ii < H.m; ii++)
		{
			edgecounter[ii] = 0;
		}
		
		for (ln1 = 0; ln1 < H.n; ln1++)
		{
			//printf("%d\n",ln1);
			
			for (jj = 0; jj < H.degL[ln1]; jj++)
			{
				edgecounter[H.edgesL[ln1][jj]]++;
			}


			for (jj = 0; jj < H.degL[ln1]; jj++)
			{
				rn1 = H.edgesL[ln1][jj];
nextkk:
				for (kk = jj+1; kk < H.degL[ln1]; kk++)
				{
					rn2 = H.edgesL[ln1][kk];
				
					//if rn1 and rn2 have any other common left node (beside node ln1) then there's a length 4 cycle
					//for (ll = edgecounter[rn1]; ll < H.degR[rn1]; ll++)
					for (ll = 0; ll < H.degR[rn1]; ll++)
					{
						if (H.edgesR[rn1][ll] > ln1)
						{	
							//for (mm = edgecounter[rn2]; mm < H.degR[rn2]; mm++)
							for (mm = 0; mm < H.degR[rn2]; mm++)
							{
								if (H.edgesR[rn2][mm] > ln1)
							
								{
							
								if (H.edgesR[rn1][ll] == H.edgesR[rn2][mm] && H.edgesL[rn1][ll]!=ln1 && H.edgesL[rn2][mm]!=ln1)		//length 4 cycle detected
								{
									noofcycles++;
									cumnoofcycles++;
																
									for (pass2 = 0; pass2 < 10; pass2++)
									{
										iscycle2 = 0;
										randedge = floor(ran2(&idum)*(noofedges-0.001));
										newrightnode = edgerightnodenumber[randedge][0];

										for (ii = 0; ii < H.degR[newrightnode]; ii++)
										{
											if (H.edgesL[newrightnode][ii] == ln1)
											{
												iscycle2 = 1;
												break;
											}
										}
									
										if (!iscycle2)
										{
											//check to see if the newleftnode is connected to rn2 or not
											newleftnode = H.edgesR[newrightnode][edgerightnodenumber[randedge][1]];

											for (ii = 0; ii < H.degL[newleftnode]; ii++)
											{
												if (H.edgesL[newleftnode][ii] == newrightnode)
												{
													indxnewleftnode = ii;
												}
									
												if (H.edgesL[newleftnode][ii] == rn2)
												{
													iscycle2 = 1;
													break;
												}
											}
						
										}
										if (!iscycle2)
											break;

									}//end of pass 2 loop



									//exchanging the edge from ln1 -> rn1 with the edge randedge if not a 2 length cycle	
									if (!iscycle2)
									{
										H.edgesL[ln1][kk] = newrightnode;
										H.edgesR[newrightnode][edgerightnodenumber[randedge][1]] = ln1;
									
										for (ii = 0; ii < H.degR[rn1]; ii++)
										{
											if (H.edgesR[rn2][ii] == ln1)
											{
												H.edgesR[rn2][ii] = newleftnode;
												break;
											}
										}
									
									
										H.edgesL[newleftnode][indxnewleftnode] = rn2;
								

										sortascending(H.edgesR[newrightnode], H.degR[newrightnode]);
										sortascending(H.edgesL[newleftnode], H.degL[newleftnode]);
										sortascending(H.edgesR[rn2], H.degR[rn2]);
										goto nextkk;
									}
								
								}//end of if cycle detected
								}//end of  if (H.edgesR[rn2][mm] > ln1)
							
							
							}
					
						}//end of if H.edgesR[rn1][ll] > ln1 
					}

				}//end of kk loop
			
				
			}//end of jj loop
		sortascending(H.edgesL[ln1], H.degL[ln1]);

		ln1 = ln1;

		}//end of ln1 loop

		
		if (noofcycles == 0)
			break;
	
	}//end of pass loop	


}


//---------------------------------- Phi -------------------------------------------
//Implements the phi function in LDPC SPA decoding
//-----------------------------------------------------------------------------------
float phi(float x)
{
	float res;
	res = log((exp(x)+1)/(exp(x)-1));

	if (x < 6e-20)
		return INF;

	else if (x > INF)
		return 0.0;

	else
		return res;
}

//--------------------------------- Sign ----------------------------------------------
//returns the sign of in. 1 if positive (or zero), -1 if negative
//--------------------------------------------------------------------------------------
char signof(float in)
{
	char out;
	if (in >= 0)
		out = +1;
	else
		out = -1;

	return out;
}

float sspalookup(float a)
{
	float b;
	if (a>=0 && a <= 0.5)
    b= 1.1000;
	else if ((a > 0.5) && (a <= 1.0))
    b= 0.5000;
	else if ((a > 1.0) && (a <= 1.8))
    b=0.1250;
	else if ((a > 1.8) && (a <= 2.5))
    b=0.0250;
	else if ((a > 2.5) && (a <= 3.0))
    b= 0.0100;
	else if (a > 3.0)
    b= 0.0025;   
    return b;
}

float sspalookupi(float a)
{
	float b;
	if ((a>=0) && (a <=0.006))
    b=3.25;
	else if ((a>0.006) && (a<=0.018))
    b=2.50;
	else if ((a>0.018)&& (a<=0.075))
    b=2.00;
	else if ((a>0.075) && (a<=0.400))
    b=1.00;
	else if ((a>0.400) && (a<=1.00))
    b=0.50;
	else if (a > 1.000)
    b=0.25;
    return b;
}

float philookup(float a)
{
	float b;
	
if  (a <= 0.15) 
    b= 3.00; 

else if ((a > 0.15) && (a <= 0.20))
    b= 2.40;
    
else if ((a > 0.20) && (a <= 0.40))
    b= 2.00;
    
else if ((a > 0.40) && (a <= 0.70))
    b= 1.30;
    
else if ((a > 0.70) && (a <= 1.2))
    b= 0.80; 
    
else if ((a > 1.20) && (a <= 1.60))
    b= 0.50;

else if ((a > 1.60) && (a <= 2.50))
    b= 0.25;
    
else if ((a > 2.50))
    b= 0.1;
    
    return b;
    
}


//LDPC SPA DECODING IN THE LOG DOMAIN
//-------------------------------------------------------------------------------------------------
void SPA(
		 float **LcbA,	//Check to bit message ; to be passed from outside the function
		 float **LbcA,	//Bit to check message; to be passed from outside the function
		 float *LchA,	//Channel LLR : log P(0)/P(1);
		 struct Hmatrix H1,  
		 char *syndromeA, 
		 int noofiter)
{
	//float **Lbc;		//messages from bit to check nodes
	int ii, jj, kk;
	float LLRtempA;
	int nodeindxRA, nodeindxLA;
	int *edgeindxRA, *edgeindxLA;
	char signtemp;


	
	//-----------for bit to check messages
	//Lbc = new float* [H.m];
	//if (Lbc == NULL)
	//	printf("Out of memory for Lbc in SPA \n");
	
	//for (ii=0; ii < H.m; ii++)
	//{
	//	Lbc[ii] = new float [H.degR[ii]];
	//	if (Lbc[ii] == NULL)
	//		printf("Out of memory for Lbc[%d] in SPA\n",ii);
	//}

	

	//for the degree indices
	edgeindxRA = new int [H1.m];
	edgeindxLA = new int [H1.n];

	if (edgeindxRA == NULL || edgeindxLA == NULL)
		printf("Out of memory in SPA\n");

	//--------------------------------------------------------------------------

	
	// ------------------ Main Iterations ---------------------------------------
	for (kk = 0; kk < noofiter;  kk ++)
	{
		//printf("Iter # = %d\n", kk);
		//------------passing message from the bit nodes to the check nodes
		//printf("Passing message from bit to check nodes\n");
		for (ii = 0; ii < H1.m; ii++)
		{
			edgeindxRA[ii] = 0;
		}
		
		
		for (ii = 0; ii < H1.n; ii++)
		{
			LLRtempA = LchA[ii];
			for (jj=0; jj < H1.degL[ii]; jj++)
			{
				LLRtempA += LcbA[ii][jj];
			}

			for (jj=0; jj < H1.degL[ii]; jj++)
			{
				 nodeindxRA = H1.edgesL[ii][jj];
				 LbcA[nodeindxRA][edgeindxRA[nodeindxRA]] = LLRtempA - LcbA[ii][jj];
				 edgeindxRA[nodeindxRA]++;
				 
			}

		}
		
		
		//printf("Passing message from check to bit nodes\n");
		//----------- passing messages from the check nodes to the bit nodes
		for (ii = 0; ii < H1.n; ii++)
		{
			edgeindxLA[ii] = 0;
		}
		
				
		for (ii = 0; ii < H1.m; ii++)
		{
			signtemp = 1 - 2*syndromeA[ii];		//if syndrome[ii] = 1, invert the sign otherwise don't
			LLRtempA  = 0.0;
	
	
			for (jj=0; jj < H1.degR[ii]; jj++)
			{
				LLRtempA += phi(fabs(LbcA[ii][jj]));
				signtemp *= signof(LbcA[ii][jj]);
			}
			
			for (jj=0; jj < H1.degR[ii]; jj++)
			{
				 nodeindxLA = H1.edgesR[ii][jj];
				 
				 LcbA[nodeindxLA][edgeindxLA[nodeindxLA]] = (signtemp / signof(LbcA[ii][jj]))*(phi(LLRtempA - phi(fabs(LbcA[ii][jj]))));
				 edgeindxLA[nodeindxLA]++;
			}
		}

	}
			
	//deleting memory
	delete(edgeindxRA);
	delete(edgeindxLA);
	
	
}


void SPAapprox(
		 float **LcbB,	//Check to bit message ; to be passed from outside the function
		 float **LbcB,	//Bit to check message; to be passed from outside the function
		 float *LchB,	//Channel LLR : log P(0)/P(1);
		 struct Hmatrix H2,  
		 char *syndromeB, 
		 int noofiter)
{
	//float **Lbc;		//messages from bit to check nodes
	int ii, jj, kk;
	float LLRtempB;
	int nodeindxRB, nodeindxLB;
	int *edgeindxRB, *edgeindxLB;
	char signtemp;


	
	//-----------for bit to check messages
	//Lbc = new float* [H.m];
	//if (Lbc == NULL)
	//	printf("Out of memory for Lbc in SPA \n");
	
	//for (ii=0; ii < H.m; ii++)
	//{
	//	Lbc[ii] = new float [H.degR[ii]];
	//	if (Lbc[ii] == NULL)
	//		printf("Out of memory for Lbc[%d] in SPA\n",ii);
	//}

	

	//for the degree indices
	edgeindxRB = new int [H2.m];
	edgeindxLB  = new int [H2.n];

	if (edgeindxRB == NULL || edgeindxLB == NULL)
		printf("Out of memory in SPA\n");

	//--------------------------------------------------------------------------

	
	// ------------------ Main Iterations ---------------------------------------
	for (kk = 0; kk < noofiter;  kk ++)
	{
		//printf("Iter # = %d\n", kk);
		//------------passing message from the bit nodes to the check nodes
		//printf("Passing message from bit to check nodes\n");
		for (ii = 0; ii < H2.m; ii++)
		{
			edgeindxRB[ii] = 0;
		}
		
		
		for (ii = 0; ii < H2.n; ii++)
		{
			LLRtempB = LchB[ii];
			for (jj=0; jj < H2.degL[ii]; jj++)
			{
				LLRtempB += LcbB[ii][jj];
			}

			for (jj=0; jj < H2.degL[ii]; jj++)
			{
				 nodeindxRB = H2.edgesL[ii][jj];
				 LbcB[nodeindxRB][edgeindxRB[nodeindxRB]] = LLRtempB - LcbB[ii][jj];
				 edgeindxRB[nodeindxRB]++;
				 
			}

		}
		
		
		//printf("Passing message from check to bit nodes\n");
		//----------- passing messages from the check nodes to the bit nodes
		for (ii = 0; ii < H2.n; ii++)
		{
			edgeindxLB[ii] = 0;
		}
		
				
		for (ii = 0; ii < H2.m; ii++)
		{
			signtemp = 1 - 2*syndromeB[ii];		//if syndrome[ii] = 1, invert the sign otherwise don't
			LLRtempB  = 0.0;
	
	
			for (jj=0; jj < H2.degR[ii]; jj++)
			{
				LbcB[ii][jj] = LbcB[ii][jj]/2;
				LLRtempB += sspalookup(fabs(LbcB[ii][jj]));
				signtemp *= signof(LbcB[ii][jj]);
			}
			
			for (jj=0; jj < H2.degR[ii]; jj++)
			{
				 nodeindxLB = H2.edgesR[ii][jj];
				 
				LcbB[nodeindxLB][edgeindxLB[nodeindxLB]] = (signtemp / signof(LbcB[ii][jj]))*(phi(LLRtempB - phi(fabs(LbcB[ii][jj]))));
				edgeindxLB[nodeindxLB]++;
			}
		}

	}
			
	//deleting memory
	delete(edgeindxRB);
	delete(edgeindxLB);
	
}

//CALCULATES THE SYNDROME AND RETURNS THE POINTER
void calcsyndrome(char *output, char *data, struct Hmatrix H)
{
	int ii,jj;
	int temp;


	for (ii = 0; ii<H.m; ii++)
	{
		temp=0;
		for (jj=0; jj<H.degR[ii]; jj++)
		{
			temp+=data[H.edgesR[ii][jj]];
		}
		output[ii] = temp%2;
		
	}

}



//INITIALIZES THE CHECK TO BIT MESSAGES TO ZERO
void initLcb(float **Lcb, struct Hmatrix H)
{
	int ii, jj;

	//Initializing check to bit messages to zero

	for (ii=0; ii < H.n; ii++)
	{
		for (jj = 0; jj < H.degL[ii]; jj++)
			Lcb[ii][jj] = 0.0;
	}

}


//ALLOCATES MEMORY FOR THE CHECK TO BIT MESSAGES
float **mallocLcb(struct Hmatrix H)
{
	int ii;
	float **Lcb;

	Lcb = new float* [H.n];
	for (ii=0; ii < H.n; ii++)
	{
		Lcb[ii] = new float [H.degL[ii]];
	}
	return Lcb;
}


float **mallocLbc(struct Hmatrix H)
{
	
	int ii;
	float **Lbc;
	
	Lbc = new float* [H.m];
	if (Lbc == NULL)
		printf("Out of memory for Lbc in mallocLbc \n");
	
	for (ii=0; ii < H.m; ii++)
	{
		Lbc[ii] = new float [H.degR[ii]];
		if (Lbc[ii] == NULL)
			printf("Out of memory for Lbc[%d] in mallocLbc\n",ii);
	}
	return Lbc;
	
}


float **mallocbitlog(int length){
	
	int ii;
	float **bitlog;
	
	bitlog = new float* [length];
	
	for (ii=0; ii<length; ii++){
		
		bitlog[ii] = new float [2]; // 2 because we have 2 bits/symbol, its 4-PAm signaling.
	
	}
		
	return bitlog;

}


void calcLch(float *Lch, float *y, float c1, float nvar, int length)
{
	// calculates the channel LLRs for AWGN channel
	int ii;
	for (ii=0; ii<length; ii++){
	Lch[ii] = 2*c1*y[ii]/nvar;
	}
	
}


//CALCULATES THE FINAL LLR 
//from the channel LLRs Lch and the check to bit messages Lch
void Lcb2LLR(float *LLR, struct Hmatrix H, float **Lcb, float *Lch, int flag)
{
	int ii, jj;

	
	for (ii=0; ii< H.n; ii++)
	{
		if (flag == EXT)	//extrinsic information
			LLR[ii] = 0.0;
		else
			LLR[ii] = Lch[ii];

		for (jj=0; jj < H.degL[ii]; jj++)
		{
			LLR[ii] += Lcb[ii][jj];
			
		}
	}
}

//COUNTS NUMBER OF ERRORS IN THE DECODED BIT STREAM
int counterrors(char *data, float *LLR, int length)
{
	int errors = 0;
	int ii;
	char decodeddata;

	for (ii=0; ii < length; ii++)
	{
		decodeddata = LLR[ii] < 0.0;
		if (data[ii] != decodeddata){
			errors++;
	
	}
}
	return errors;

}

//Hard Thresholds THE LLR vector

void hardthreshold(char *decodeddata, float *LLR, int length)
{
	int ii;
	
	for (ii = 0; ii < length; ii++)
	{
		if (LLR[ii] < 0 )
			decodeddata[ii] = 1;
		else
			decodeddata[ii] = 0;
	}
	
}
