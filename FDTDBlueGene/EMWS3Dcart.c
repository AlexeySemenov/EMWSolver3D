#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <mpix.h>
#include <omp.h>

#include "Phys.h"
#include "FieldFloat.h"

int step;
float dx, dy, dz, dt;
float SigE, SigM;

float** InterfaceEx;
float** InterfaceEy;
float** InterfaceEz;

float** InterfaceHx;
float** InterfaceHy;
float** InterfaceHz;

float* bufx;
float* bufy;

float* rbufx;
float* rbufy;

int sliceSize, addSize;
MPI_Request requestE1;
MPI_Request requestE2;

MPI_Request requestH1;
MPI_Request requestH2;
//======================Additional functions===============================
float*** Malloc3darray(int x, int y, int z)
{
	float*** ar;
	int i, j, k;
	ar = (float***) malloc(x*sizeof(float**));
	for (i = 0; i < x; i++)
		ar[i] = (float**) malloc(y*sizeof(float*));
	for (i = 0; i < x; i++)
		for(j = 0; j < y; j++)
			ar[i][j] = (float*) malloc(z*sizeof(float));

	for(i = 0; i < x; i++)
		for(j = 0; j < y; j++)
			for(k = 0; k < z; k++)
				ar[i][j][k] = 0;
	return ar;

}

float** Malloc2darray(int x, int y)
{
	float** ar;
	int i, j, k;
	ar = (float**) malloc(x*sizeof(float*));
	for (i = 0; i < x; i++)
		ar[i] = (float*) malloc(y*sizeof(float));


	for(i = 0; i < x; i++)
		for(j = 0; j < y; j++)
			ar[i][j] = 0;
	return ar;

}

//======================Update FDTD===============================
void UpdateMainDomain()
{

	int i, j, k;
	#pragma omp parallel for
	for(i = 1; i < gridX-1; i++)
		for(j = 1; j < gridY-1; j++)
			for(k = 1; k < sliceSize-1; k++)
			{
				Ex[i][j][k] = Ex[i][j][k] + (dt/(Eps[i][j][k]*EPS_Z))*
					((Hz[i][j][k] - Hz[i][j-1][k])/dy - (Hy[i][j][k] - Hy[i][j][k-1])/dz);
				Ey[i][j][k] = Ey[i][j][k] + (dt/(Eps[i][j][k]*EPS_Z))*
					((Hx[i][j][k] - Hx[i][j][k-1])/dz - (Hz[i][j][k] - Hz[i-1][j][k])/dx);
				Ez[i][j][k] = Ez[i][j][k] + (dt/(Eps[i][j][k]*EPS_Z))*
					((Hy[i][j][k] - Hy[i-1][j][k])/dx - (Hx[i][j][k] - Hx[i][j-1][k])/dy);
			}
	#pragma omp parallel for
	for(i = 1; i < gridX-1; i++)
		for(j = 1; j < gridY-1; j++)
			for(k = 1; k < sliceSize-1; k++)
			{
				Hx[i][j][k] = Hx[i][j][k] + (dt/(Mu[i][j][k]*MU_Z))*
					((Ey[i][j][k+1] - Ey[i][j][k])/dz - (Ez[i][j+1][k] - Ez[i][j][k])/dy);
				Hy[i][j][k] = Hy[i][j][k] + (dt/(Mu[i][j][k]*MU_Z))*
					((Ez[i+1][j][k] - Ez[i][j][k])/dx - (Ex[i][j][k+1] - Ex[i][j][k])/dz);
				Hz[i][j][k] = Hz[i][j][k] + (dt/(Mu[i][j][k]*MU_Z))*
					((Ex[i][j+1][k] - Ex[i][j][k])/dy - (Ey[i+1][j][k] - Ey[i][j][k])/dx);
			}
}

void UpdateInterfaceZ(int _rank, int _num)
{
	int i, j;
	int es;

	es = sliceSize-1;
	if(_rank != 0)
	{
		#pragma omp parallel for
		for(i = 1; i < gridX-1; i++)
			for(j = 1; j < gridY-1; j++)
			{
				Ex[i][j][0] = Ex[i][j][0] + (dt/(Eps[i][j][0]*EPS_Z))*
						((Hz[i][j][0] - Hz[i][j-1][0])/dy - (Hy[i][j][0] - InterfaceHy[i][j])/dz);
				Ey[i][j][0] = Ey[i][j][0] + (dt/(Eps[i][j][0]*EPS_Z))*
						((Hx[i][j][0] - InterfaceHx[i][j])/dz - (Hz[i][j][0] - Hz[i-1][j][0])/dx);
				Ez[i][j][0] = Ez[i][j][0] + (dt/(Eps[i][j][0]*EPS_Z))*
					((Hy[i][j][0] - Hy[i-1][j][0])/dx - (Hx[i][j][0] - Hx[i][j-1][0])/dy);

				Hx[i][j][0] = Hx[i][j][0] + (dt/(Mu[i][j][0]*MU_Z))*
					((Ey[i][j][1] - Ey[i][j][0])/dz - (Ez[i][j+1][0] - Ez[i][j][0])/dy);
				Hy[i][j][0] = Hy[i][j][0] + (dt/(Mu[i][j][0]*MU_Z))*
					((Ez[i+1][j][0] - Ez[i][j][0])/dx - (Ex[i][j][1] - Ex[i][j][0])/dz);
				Hz[i][j][0] = Hz[i][j][0] + (dt/(Mu[i][j][0]*MU_Z))*
					((Ex[i][j+1][0] - Ex[i][j][0])/dy - (Ey[i+1][j][0] - Ey[i][j][0])/dx);
			}
	}

	if(_rank != _num-1)
	{
		#pragma omp parallel for
		for(i = 1; i < gridX-1; i++)
			for(j = 1; j < gridY-1; j++)
			{
				Ex[i][j][es] = Ex[i][j][es] + (dt/(Eps[i][j][es]*EPS_Z))*
					((Hz[i][j][es] - Hz[i][j-1][es])/dy - (Hy[i][j][es] - Hy[i][j][es-1])/dz);
				Ey[i][j][es] = Ey[i][j][es] + (dt/(Eps[i][j][es]*EPS_Z))*
					((Hx[i][j][es] - Hx[i][j][es-1])/dz - (Hz[i][j][es] - Hz[i-1][j][es])/dx);
				Ez[i][j][es] = Ez[i][j][es] + (dt/(Eps[i][j][es]*EPS_Z))*
					((Hy[i][j][es] - Hy[i-1][j][es])/dx - (Hx[i][j][es] - Hx[i][j-1][es])/dy);


				Hx[i][j][es] = Hx[i][j][es] + (dt/(Mu[i][j][es]*MU_Z))*
							((InterfaceEy[i][j] - Ey[i][j][es])/dz - (Ez[i][j+1][es] - Ez[i][j][es])/dy);
				Hy[i][j][es] = Hy[i][j][es] + (dt/(Mu[i][j][es]*MU_Z))*
							((Ez[i+1][j][es] - Ez[i][j][es])/dx - (InterfaceEx[i][j] - Ex[i][j][es])/dz);
				Hz[i][j][es] = Hz[i][j][es] + (dt/(Mu[i][j][es]*MU_Z))*
							((Ex[i][j+1][es] - Ex[i][j][es])/dy - (Ey[i+1][j][es] - Ey[i][j][es])/dx);
			}
	}
}

//======================Send===============================
void SendE(int _rank, int _num, MPI_Comm _comm)
{
	int i, j;

	MPI_Request request1;
	MPI_Request request2;

	#pragma omp parallel for
	for(i = 0; i < gridX; i++)
		for(j = 0; j < gridY; j++)
		{
			bufx[i*gridY + j] = InterfaceEx[i][j];
			bufy[i*gridY + j] = InterfaceEy[i][j];
		}

	//printf("(%04d): BSendEx.\n", _rank);
	MPI_Isend(bufx, gridX*gridY, MPI_FLOAT, _rank+1, 0, _comm, &request1);
	//printf("(%04d): SendEx.\n", _rank);
	MPI_Isend(bufy, gridX*gridY, MPI_FLOAT, _rank+1, 1, _comm, &request2);
	//printf("(%04d): SendEy.\n", _rank);
}
void SendH(int _rank, int _num, MPI_Comm _comm)
{
	int i, j;

	MPI_Request request1;
	MPI_Request request2;
	
#pragma omp parallel for
for(i = 0; i < gridX; i++)
		for(j = 0; j < gridY; j++)
		{
			bufx[i*gridY + j] = InterfaceHx[i][j];
			bufy[i*gridY + j] = InterfaceHy[i][j];
		}
	//printf("(%04d): BSendHx.\n", _rank);
	MPI_Isend(bufx, gridX*gridY, MPI_FLOAT, _rank-1, 2, _comm, &request1);
	//printf("(%04d): SendHx.\n", _rank);
	MPI_Isend(bufy, gridX*gridY, MPI_FLOAT, _rank-1, 3, _comm, &request2);
	//printf("(%04d): SendHy.\n", _rank);
}

void RecvE(int _rank, int _num, MPI_Comm _comm)
{
	int i, j;

	MPI_Irecv(rbufx, gridX*gridY, MPI_FLOAT, _rank-1, 0, _comm, &requestE1);
	MPI_Irecv(rbufy, gridX*gridY, MPI_FLOAT, _rank-1, 1, _comm, &requestE2);

}

void CopyE(void)
{
	int i, j;
	#pragma omp parallel for
	for(i = 0; i < gridX; i++)
		for(j = 0; j < gridY; j++)
		{
			InterfaceEx[i][j] = rbufx[i*gridY + j];
			InterfaceEy[i][j] = rbufy[i*gridY + j];
		}
}

void RecvH(int _rank, int _num, MPI_Comm _comm)
{
	int i, j;
	
	MPI_Irecv(rbufx, gridX*gridY, MPI_FLOAT, _rank+1, 2, _comm, &requestH1);
	MPI_Irecv(rbufy, gridX*gridY, MPI_FLOAT, _rank+1, 3, _comm, &requestH2);


	
}

void CopyH(void)
{
	int i, j;
	#pragma omp parallel for
	for(i = 0; i < gridX; i++)
		for(j = 0; j < gridY; j++)
		{
			InterfaceHx[i][j] = rbufx[i*gridY + j];
			InterfaceHy[i][j] = rbufy[i*gridY + j];
		}
}
//======================Main===============================
int main(int argc, char** argv)
{
	int rank, numprocs;
	int i, j, k;
	int t = 0;
	float startTime;
	int res; 

    	MPI_Comm comm;// = MPI_COMM_WORLD;

    	MPI_Init(&argc, &argv);
    	
    	res = MPIX_Cart_comm_create(&comm);
    	if (res != MPI_SUCCESS)
    	{
        	fprintf(stderr, "MPIX_Cart_comm_create() failed "
                	"(return code: %d)\n", res);
        	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    	}

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm,&numprocs);

	MPI_Status status;

	/*if (rank == 0)
    	{
        	const char *val, *name;

        	name = "OMP_NUM_THREADS";
       	 printf("%s: %s\n", name, (val = getenv(name)) ? val: "");

        	name = "BG_MAPPING";
        	printf("%s: %s\n", name, (val = getenv(name)) ? val: "");
    	}*/

	/*int rankInDim;
	MPI_Comm_rank(comm, &rankInDim);

	int rank_source, rank_desta, rank_destb, rank_destc, rank_destd;
	MPI_Cart_shift(comm, 0,1,&rank_source, &rank_desta);
	MPI_Cart_shift(comm, 1,1,&rank_source, &rank_destb);
	MPI_Cart_shift(comm, 2,1,&rank_source, &rank_destc);
	MPI_Cart_shift(comm, 3,1,&rank_source, &rank_destd);
	fprintf(stderr,"I am known in the world as %d my adjacents are -> %d - %d - %d - %d \n", rankInDim, rank_desta, rank_destb, rank_destc, rank_destd);
	*/

	gridX = 512;
	gridY = 512;
	gridZ = 1024;

	t = 0;

	step = 20;

	dx = LAMBDA/step;
	dy = LAMBDA/step;
	dz = LAMBDA/step;

	dt = dx/(2.0*CC);

	SigE = 0;
	SigM = 0;

	sliceSize = gridZ/numprocs;
	addSize = gridZ - sliceSize*numprocs;
	if(rank == numprocs-1)
		sliceSize += addSize;

	printf("node:%d - %d \n.", rank, sliceSize);
	Ex = Malloc3darray(gridX, gridY, sliceSize);
	Ey = Malloc3darray(gridX, gridY, sliceSize);
	Ez = Malloc3darray(gridX, gridY, sliceSize);

	Hx = Malloc3darray(gridX, gridY, sliceSize);
	Hy = Malloc3darray(gridX, gridY, sliceSize);
	Hz = Malloc3darray(gridX, gridY, sliceSize);

	Eps = Malloc3darray(gridX, gridY, sliceSize);
	Mu = Malloc3darray(gridX, gridY, sliceSize);

	InterfaceEx = Malloc2darray(gridX, gridY);
	InterfaceEy = Malloc2darray(gridX, gridY);
	InterfaceEz = Malloc2darray(gridX, gridY);

	InterfaceHx = Malloc2darray(gridX, gridY);
	InterfaceHy = Malloc2darray(gridX, gridY);
	InterfaceHz = Malloc2darray(gridX, gridY);

	bufx = (float*) malloc(gridX*gridY*sizeof(float));
	bufy = (float*) malloc(gridX*gridY*sizeof(float));

	rbufx = (float*) malloc(gridX*gridY*sizeof(float));
	rbufy = (float*) malloc(gridX*gridY*sizeof(float));

	
	for(i = 0; i < gridX; i++)
		for(j = 0; j < gridY; j++)
			for(k = 0; k < sliceSize; k ++)
			{
				Eps[i][j][k] = 1;
				Mu[i][j][k] = 1;
			}
	printf("(%04d): Initialized!\n", rank);
	MPI_Barrier(comm);
	if(rank == 0) startTime = MPI_Wtime();
	while (t < 10)
	{
//==================Data irecv==========================
		if(rank == 0)
		{
			RecvH(rank, numprocs, comm);
		}
		else if(rank != numprocs-1)
		{
			RecvE(rank, numprocs, comm);
			RecvH(rank, numprocs, comm);
		}
		else
		{
			RecvE(rank, numprocs, comm);
		}
//==================Data isend==========================
		if(rank == 0)
		{
			SendE(rank, numprocs, comm);
		}
		else if(rank != numprocs-1)
		{
			SendE(rank, numprocs, comm);
			SendH(rank, numprocs, comm);
		}
		else
		{
			SendH(rank, numprocs, comm);
		}
//=================Update main domain		
		if(rank == 0)
			printf("(%04d, %d): Solving..\n", rank, t);
		UpdateMainDomain();
		MPI_Barrier(comm);
		if(rank == 0)
			printf("(%04d, %d): Main domain updated.\n", rank, t);
//==================Wait for data==============
		if(rank == 0)
		{
			//printf("(%04d, %d): Waiting H...\n", rank, t);
			MPI_Wait(&requestH1, &status);
			MPI_Wait(&requestH2, &status);
			//printf("(%04d, %d): Done H.\n", rank, t);

		}
		else if(rank != numprocs-1)
		{
			//printf("(%04d, %d): Waiting E H...\n", rank, t);
			MPI_Wait(&requestE1, &status);
			MPI_Wait(&requestE2, &status);
			MPI_Wait(&requestH1, &status);
			MPI_Wait(&requestH2, &status);
			//printf("(%04d, %d): Done E H.\n", rank, t);
		}
		else
		{
			//printf("(%04d, %d): Waiting E...\n", rank, t);
			MPI_Wait(&requestE1, &status);
			MPI_Wait(&requestE2, &status);
			//printf("(%04d, %d): Done E.\n", rank, t);
		}
//============Copy Data=======================================
		if(rank == 0)
		{
			CopyH();
		}
		else if(rank != numprocs-1)
		{
			CopyE();
			CopyH();
		}
		else
		{
			CopyE();
		}
		//printf("(%04d, %d): Data received.\n", rank, t);
		UpdateInterfaceZ(rank, numprocs);
		if(rank == 0)
		printf("(%04d, %d): done.\n", rank, t);
		MPI_Barrier(comm);
		t++;
	}
    
	

#pragma omp parallel
    if(rank ==0)
		printf("Hail to the king:%lf\n", MPI_Wtime()-startTime);

    MPI_Finalize();
    return 0;
}
