#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

int find_max(float *mat,int col,int k,int n,int piv[]){  //function that finds the max of a column
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


int main(int argc, char* argv[]) {

	if(argc<2 || argc>3){
		printf("la8os  orismata\n");
		return 0;
	}
	
	int i,k,j,l,start,end,max_loc,temp,numcols,mycol,num_comm,rank, size;
	int n=atoi(argv[1]);
	double t1, t2;
	float *mat,num;
	int N=n+1;
	char cn;
	int counter=-1;
	int piv[n];
	for(i=0; i<n; i++){
		piv[i]=i;    // Create the pivot array
	}
	float message[n+1];
	div_t root;

	MPI_Status status;
	
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank==0){
		printf("Print the solution? [Y/N]\n");
		scanf("%c", &cn);
		MPI_Send(&cn, 1, MPI_CHAR, size-1, 101, MPI_COMM_WORLD);
	}else if(rank==size-1){
		MPI_Recv(&cn, N, MPI_CHAR, 0, 101, MPI_COMM_WORLD, &status);
	}

	MPI_Comm vec_comm[size];  //Create the communicators' array
	int vec_rank[size];
	int vec_size[size];
	vec_comm[0]=MPI_COMM_WORLD;
	vec_rank[0]=rank;

	for(i=0; i<size-1; i++){
		if(rank<=i){
			num_comm=0;
		}else{
			num_comm=1;
		}
		if(rank>=i){
			MPI_Comm_split(vec_comm[i], num_comm, vec_rank[i], (vec_comm+i+1));  // Create new communicator
			MPI_Comm_rank(vec_comm[i+1], (vec_rank+i+1));
		}
		if(rank==i){
			MPI_Comm_free((vec_comm+i+1));   // Delete useless communicators
		}
	}
	
	//Define start and end column for each proccess
	start=(rank)*(n/size);
	end=start+(n/size)-1;

	// Bind memory for the collumns
	if(rank==size-1){
		end++;
		numcols=end-start+1;
		mat = (float *)malloc((N*numcols-1) * sizeof(float));
	}else{
		numcols=end-start+1;
		mat = (float *)malloc((N*numcols) * sizeof(float));
	}

	
	if(rank==0){
		float *matrix = (float *)malloc(((n+1)*(n+1)) * sizeof(float));
		mat = (float *)malloc((N*numcols) * sizeof(float));
		
		// Read the matrix from a txt
		if(argc==3){
		
			FILE *myfile;
			float myvariable;
			myfile=fopen(argv[2], "r");
	
			for(i = 0; i < n; i++){
				for (j = 0 ; j < N; j++){
					fscanf(myfile,"%f",&myvariable);
					matrix[i+j*N]=myvariable;
					if(j<=end){
						mat[i+j*N]=myvariable;
					}
				}
			}
			fclose(myfile);

		// Create a random matrix
		}else if(argc==2){
			srand48(clock());
			for(j=0; j<N; j++){
				for(i=0; i<n; i++){
					matrix[i+j*N]=(float)drand48()*10.0;
					if(j<=end){
						mat[i+j*N]=matrix[i+j*N];
					}
				}
			}
		}

		// Send collumns to other proccesses
		for(j=end+1; j<N; j++){
			root=div(j,(n/size));
			if(root.quot>size-1){
				MPI_Send(matrix + j*N, N, MPI_FLOAT, root.quot-1, j, MPI_COMM_WORLD);
			}else{
				MPI_Send(matrix + j*N, N, MPI_FLOAT, root.quot, j, MPI_COMM_WORLD);
			}
		}
		free(matrix);
	// Receive collumns
	}else{
		mat = (float *)malloc((N*numcols) * sizeof(float));
		for(j=start;j<=end; j++){
			MPI_Recv(&message, N, MPI_FLOAT, 0, j, MPI_COMM_WORLD,&status);
			for(k=0; k<n; k++){
				mat[k+(j-start)*N]=message[k];
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	// Start the computations
	for(k=0; k<n; k++){
		if(k>=start && k<=end){  //  Check if the collumn belongs to proccess
			mycol=k-start;
			max_loc=find_max(mat,mycol,k,n,piv); // Find max location
			temp=piv[max_loc];
			piv[max_loc]=piv[k];
			piv[k]=temp;
			mat[n+mycol*N]=(float)max_loc*1.0;
	
			if(rank != size-1){ // Checks if the proccess isn't the last
				MPI_Bcast(mat + mycol*N, N, MPI_FLOAT, 0, vec_comm[rank]);  // Send the k-collumn
			}
			if(k==end){  //Check if the proccess must end computations
				continue;
			}else{ // Execute computations
				for(j=mycol+1;j<numcols;j++){
					mat[piv[k]+j*N]=mat[piv[k]+j*N]/mat[piv[k]+mycol*N];
					for(i=0;i<n;i++){
						if(i!=piv[k]){
							mat[i+j*N]=mat[i+j*N]-(mat[i+mycol*N])*(mat[piv[k]+j*N]);
						}	
					}
				}
			}
		}else if(k<start){  // Check if the proccess' zone is after k-collumn
			root=div(k,(n/size));
			MPI_Bcast(&message, N, MPI_FLOAT, 0, vec_comm[root.quot]);  // Receive k-collumn

			max_loc=(int)message[n];
			temp=piv[max_loc];
			piv[max_loc]=piv[k];
			piv[k]=temp;
			for(j=0;j<numcols;j++){  //Execute computations
				mat[piv[k]+j*N]=mat[piv[k]+j*N]/message[piv[k]];
				for(i=0;i<n;i++){
					if(i!=piv[k]){
						mat[i+j*N]=mat[i+j*N]-(message[i])*(mat[piv[k]+j*N]);
					}	
				}
			}
		}else{
			free(mat);
			if(rank !=0){
				MPI_Comm_free((vec_comm+rank));  //Delete useless communicator
			}
			MPI_Finalize();    // Kill useless procceesses
			return 0;

		}
	}

	// Print the results
	t2 = MPI_Wtime();
	if(rank==size-1){
		if(cn=='Y'){
			for (l=0;l<n;l++){
				printf("%lf\n",mat[piv[l]+(numcols-1)*N]);
			}
		}
		printf( "Elapsed time is %f\n", t2 - t1 );  //Print executional time	
	}
	free(mat);
	MPI_Comm_free((vec_comm+rank));  //Delete the last communicator
	MPI_Finalize();
	return 0;
}
