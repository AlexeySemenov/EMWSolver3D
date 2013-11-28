#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#include "Phys.h"
#include "Field.h"


int step;
double dx, dy, dz, dt;
double SigE, SigM;

double** InterfaceEx;
double** InterfaceEy;
double** InterfaceEz;

double** InterfaceHx;
double** InterfaceHy;
double** InterfaceHz;

double* bufx;
double* bufy;

double* rbufx;
double* rbufy;

int sliceSize, addSize;
MPI_Request requestE1;
MPI_Request requestE2;

MPI_Request requestH1;
MPI_Request requestH2;
//======================Additional functions===============================
double*** Malloc3darray(int x, int y, int z)
{
	double*** ar;
	int i, j, k;
	ar = (double***) malloc(x*sizeof(double**));
	for (i = 0; i < x; i++)
		ar[i] = (double**) malloc(y*sizeof(double*));
	for (i = 0; i < x; i++)
		for(j = 0; j < y; j++)
			ar[i][j] = (double*) malloc(z*sizeof(double));

	for(i = 0; i < x; i++)
		for(j = 0; j < y; j++)
			for(k = 0; k < z; k++)
				ar[i][j][k] = 0;
	return ar;

}

double** Malloc2darray(int x, int y)
{
	double** ar;
	int i, j, k;
	ar = (double**) malloc(x*sizeof(double*));
	for (i = 0; i < x; i++)
		ar[i] = (double*) malloc(y*sizeof(double));


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
	MPI_Isend(bufx, gridX*gridY, MPI_DOUBLE, _rank+1, 0, _comm, &request1);
	//printf("(%04d): SendEx.\n", _rank);
	MPI_Isend(bufy, gridX*gridY, MPI_DOUBLE, _rank+1, 1, _comm, &request2);
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
	MPI_Isend(bufx, gridX*gridY, MPI_DOUBLE, _rank-1, 2, _comm, &request1);
	//printf("(%04d): SendHx.\n", _rank);
	MPI_Isend(bufy, gridX*gridY, MPI_DOUBLE, _rank-1, 3, _comm, &request2);
	//printf("(%04d): SendHy.\n", _rank);
}

void RecvE(int _rank, int _num, MPI_Comm _comm)
{
	int i, j;

	MPI_Irecv(rbufx, gridX*gridY, MPI_DOUBLE, _rank-1, 0, _comm, &requestE1);
	MPI_Irecv(rbufy, gridX*gridY, MPI_DOUBLE, _rank-1, 1, _comm, &requestE2);

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
	
	MPI_Irecv(rbufx, gridX*gridY, MPI_DOUBLE, _rank+1, 2, _comm, &requestH1);
	MPI_Irecv(rbufy, gridX*gridY, MPI_DOUBLE, _rank+1, 3, _comm, &requestH2);


	
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
	double startTime; 

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	MPI_Status status;

	gridX = 512;
	gridY = 512;
	gridZ = 4096;

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

	bufx = (double*) malloc(gridX*gridY*sizeof(double));
	bufy = (double*) malloc(gridX*gridY*sizeof(double));

	rbufx = (double*) malloc(gridX*gridY*sizeof(double));
	rbufy = (double*) malloc(gridX*gridY*sizeof(double));

	
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
		printf("Calc time:%lf\n Nproc = %d \n", MPI_Wtime()-startTime, numprocs);

    MPI_Finalize();
    return 0;
}
