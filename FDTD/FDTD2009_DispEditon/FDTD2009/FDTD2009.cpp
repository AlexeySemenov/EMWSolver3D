#pragma once



#include "stdhead.h"
#include "Field.h"
#include "FieldTE.h"
#include "FieldTM.h"
#include "Solver.h"
#include "Solver3d.h"

using namespace std;

int main(int argc, char* argv[])
{
	clock_t wrkTime;
	cout << "<<FDTD2009 Solver>>" << endl;
	#pragma omp parallel
	  {
		printf("Numder of threads - %d. ", omp_get_num_threads());
		printf("My number is %d\n", omp_get_thread_num());
		
	  }
	system("PAUSE");
	CSolver sol("TM_UPML", 240, 640, 0);

	sol.InitSimpleFDTD_TM_UPML();
	sol.TimeStepLoop(5000);
	//CSolver3d sol(1000,1000,1);
	//sol.Init3DFDTD_UPML();
	//sol.TimeStepLoop(0);

	//CFieldTM EM(100, 100);
	//double z1 = 50;

	//double dx,dy,dz;
	//dx = dy = dz = 0.000000000005;

	//z1 = z1*dz;

	//double eps1 = 1;
	//double eps2 = 5;

	//double H0 = 1;

	//double teta = M_PI/4;

	//double k1 = OMEGA*sqrt(eps1);
	//double k2 = OMEGA*sqrt(eps2);

	//double beta1x;
	//double beta2x;
	//double beta1y;
	//double beta2y;

	//beta1x = k1*cos(teta);
	//beta1y = k1*sin(teta);
	//beta2y = beta1y;

	//beta2x = sqrt(k2*k2 - beta2y*beta2y);

	//double G;
	//double tmp1, tmp2;
	//tmp1=beta1x/(OMEGA*eps1);
	//tmp2=beta2x/(OMEGA*eps2);
	//G = (tmp1 - tmp2)/(tmp1  + tmp2);

	//double T;
	//T = 1 + G;



	//for(int i = 0; i < EM.gridX; i++)
	//	for(int j = 0; j < EM.gridY; j++)
	//	{
	//		//if(i>500) EM.Ez[i][j] = cos(kx*j*dx)*Wplus*cos(kz1*i*dz);
	//		//else EM.Ez[i][j] = cos(kx*i*dx)*Wminus*cos(-kz0*j*dz)+cos(kx*i*dx)*A*cos(kz0*j*dz);
	//		if (i<50)EM.Ez[i][j] =H0*(1+G*cos(2*beta1x*i*dx))*cos(-beta1x*i*dx - beta1y*j*dy);
	//		else EM.Ez[i][j] = H0*T*cos(-beta2x*i*dx - beta2y*j*dy);
	//	}
	//	EM.WriteToNetCDF(21745,0,0,0,0);





	wrkTime = clock();
	cout << "Solve Time: "<<wrkTime/CLOCKS_PER_SEC << endl;

	system("PAUSE");
	return 1;
}