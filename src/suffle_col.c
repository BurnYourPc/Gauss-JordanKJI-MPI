#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

int find_max(float *mat,int col,int k,int n,int piv[]){    //function that finds the max of a column
	int i;
	int loc=k;
	float max=abs(mat[piv[k]+col*(n+1)]);
	for(i=k+1; i<n; i++){
		if (abs(mat[piv[i]+col*(n+1)])>max){
			max=abs(mat[piv[i]+col*(n+1)]);
			loc=i;
		}
	}
	return loc;
}

int give_sign(int x){
if (x > 0) return 1;
if (x < 0) return 0;
return 1;
}

int main(int argc, char* argv[]) {

	if(argc<2 || argc>3){
		printf("la8os  orismata\n");
		return 0;
	}

	int i,k,j,l,start,end,max_loc,temp,num_comm,rank, size;
	int n=atoi(argv[1]);
	
	double t1, t2;
	float *mat,num;
	int N=n+1,counter=0,numcols=0;
	int piv[n];
	char cn;
	for(i=0; i<n; i++){
		piv[i]=i;   // Create the pivot array
	}
	float message[n+1];

	div_t root,mycol,block;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Compute the number of collumns that belong to each proccess
	for(i=0;i<N;i+=size){
		if(rank+counter*size<N){
			numcols++;
			counter++;
		}
	}

	int vec_rank[size];
	int vec_size[size];
	MPI_Comm vec_comm[size];    //Create the communicators' array
	vec_comm[0]=MPI_COMM_WORLD;
	vec_rank[0]=rank;
	vec_comm[1]=MPI_COMM_WORLD;
	vec_rank[1]=rank;

	for(i=1; i<size-1; i++){
		if(rank<=i && rank!=0){
			num_comm=0;
		}else{
			num_comm=1;
		}
		MPI_Comm_split(vec_comm[i], num_comm, vec_rank[i], (vec_comm+i+1));   // Create new communicator
		MPI_Comm_rank(vec_comm[i+1], (vec_rank+i+1));
	}

	for(i=1; i<size-1; i++){
		if(rank<=i && rank!=0){
			MPI_Comm_free((vec_comm+i+1));     // Delete useless communicators
		}
	}
	
	counter=0;
	if(rank==0){
		if(argc==3){
			float *matrix = (float *)malloc(((n+1)*(n+1)) * sizeof(float));    // Bind memory for the matrix
			mat = (float *)malloc((N*numcols) * sizeof(float));  // Bind memory for the collumns
			FILE *myfile;
			float myvariable;
			myfile=fopen(argv[2], "r");

			// Read the matrix from a txt
			for(i = 0; i < n; i++){
				for (j = 0 ; j < N; j++){
					fscanf(myfile,"%f",&myvariable);
					matrix[i+j*N]=myvariable;
					if(j%size==0){
						mycol=div(j,size);
						mat[i+mycol.quot*N]=matrix[i+j*N];
					}
				}
			}
			fclose(myfile);
			for(j=0; j<N; j++){
				if(j%size != 0){
					MPI_Send(matrix + j*N, N, MPI_FLOAT, j%size, j, MPI_COMM_WORLD);
				}
			}
			free(matrix);
		}else if(argc==2){
			float *matrix = (float *)malloc(((n+1)*(n+1)) * sizeof(float));    // Bind memory for the matrix
			mat = (float *)malloc((N*numcols) * sizeof(float));     // Bind memory for the collumns
			for(j=0; j<N; j++){
				for(i=0; i<n; i++){
					matrix[i+j*N]=(float)drand48()*10.0;   // Create a random matrix
					if(j%size==0){
						mycol=div(j,size);
						mat[i+mycol.quot*N]=matrix[i+j*N];
					}
				}
				if(j%size != 0){
					MPI_Send(matrix + j*N, N, MPI_FLOAT, j%size, j, MPI_COMM_WORLD);  // Send collumns to other proccesses
				}
			}
			free(matrix);
		}
	}else{
		mat = (float *)malloc((N*numcols) * sizeof(float));  // Bind memory for the collumns
		for(i=0;i<N; i+=size){
			if(rank+counter*size<N){
				MPI_Recv(&message, N, MPI_FLOAT, 0, rank+counter*size, MPI_COMM_WORLD,&status);   // Receive collumns
				for(k=0; k<n; k++){
					mat[k+counter*N]=message[k];
				}
				counter++;
			}
		}
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	// Start the computations
	for(k=0; k<n; k++){
		root=div(k,size);
		block=div(n,size);
		if(k%size==rank){     // Check if the collumn belongs to proccess
			mycol=div(k,size);
			max_loc=find_max(mat,mycol.quot,k,n,piv);    // Find max location
			temp=piv[max_loc];
			piv[max_loc]=piv[k];
			piv[k]=temp;
			mat[n+mycol.quot*N]=(float)max_loc*1.0;
			if(mycol.quot==numcols-1){    // Check if the k-collumn is the last collumn of the proccess
				MPI_Bcast(mat + (mycol.quot)*N, N, MPI_FLOAT, vec_rank[rank], vec_comm[rank]);  // Send k-collumn to other proccesses
				continue;
			}else{
				MPI_Bcast(mat + (mycol.quot)*N, N, MPI_FLOAT, rank, MPI_COMM_WORLD);    // Send k-collumn to other proccesses

				// Execute computations
				for(j=mycol.quot+1;j<numcols;j++){
					mat[piv[k]+j*N]=mat[piv[k]+j*N]/mat[piv[k]+mycol.quot*N];
					for(i=0;i<n;i++){
						if(i!=piv[k]){
							mat[i+j*N]=mat[i+j*N]-(mat[i+mycol.quot*N])*(mat[piv[k]+j*N]);
						}	
					}
				}
			}

		}else if(k%size<rank || root.quot<numcols-1){
			if(k%size==0 || root.quot!=block.quot-1 ){   // Check if k-collumn is among last p collumns (p is the number of proccesses)
				MPI_Bcast(&message, N, MPI_FLOAT, k%size, MPI_COMM_WORLD);  // Receive k-collumn
			}else{
				MPI_Bcast(&message, N, MPI_FLOAT, 1, vec_comm[k%size]);   // Receive k-collumn
			}
			
			max_loc=(int)message[n];
			temp=piv[max_loc];
			piv[max_loc]=piv[k];
			piv[k]=temp;
			root=div(k,size);

			//Execute computations
			for(j=root.quot+give_sign(k%size-rank);j<numcols;j++){
				mat[piv[k]+j*N]=mat[piv[k]+j*N]/message[piv[k]];
				for(i=0;i<n;i++){
					if(i!=piv[k]){
						mat[i+j*N]=mat[i+j*N]-(message[i])*(mat[piv[k]+j*N]);
					}	
				}
			}
		}else{
			free(mat);
			if(rank>1){
				MPI_Comm_free((vec_comm+rank));   //Delete useless communicator
			}
			MPI_Finalize();  // Kill useless procceesses
			return 0;
		}
	}

	// Print the results
	if(n%size==rank){
		t2 = MPI_Wtime();
		printf( "Elapsed time is %f\n", t2 - t1 );      //Print executional time
		printf("Print the solution? [Y/N]\n");
		scanf("%c", &cn);
		if(cn=='Y'){
			for (l=0;l<n;l++){
				printf("%lf\n",mat[piv[l]+(numcols-1)*N]);
			}
		}
	}
	free(mat);
	for(i=2; i<size; i++){
		MPI_Comm_free((vec_comm+i));     //Delete the last communicator
	}
	MPI_Finalize();
	return 0;
}
