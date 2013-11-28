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
	cout << "<<FDTD2011 Solver>>" << endl;
	#pragma omp parallel
	  {
		printf("Numder of threads - %d. ", omp_get_num_threads());
		printf("My number is %d\n", omp_get_thread_num());
		
	  }
	system("PAUSE");
	//CSolver sol("TE", 420, 1080, 0);
	//sol.InitSimpleFDTD_TE_UPML();

	
	//sol.TimeStepLoopTE(4000);

	/*CSolver3d* sol;
	double epsStp = 2; 
	for(int i = 1; i <= 20; i++)
	{
		sol = new CSolver3d(100,60,420);
		sol->taskID = i;
		sol->Init3DFDTD_UPML(0, epsStp + i);
		sol->TimeStepLoop(450);

		delete sol;
	}

	for(int i = 1; i <= 20; i++)
	{
		sol = new CSolver3d(100,60,420);
		sol->taskID = i+100;
		sol->Init3DFDTD_UPML(20, epsStp + i);
		sol->TimeStepLoop(450);

		delete sol;
	}

	for(int i = 1; i <= 20; i++)
	{
		sol = new CSolver3d(100,60,420);
		sol->taskID = i+200;
		sol->Init3DFDTD_UPML(30, epsStp + i);
		sol->TimeStepLoop(450);

		delete sol;
	}*/

	CSolver3d* sol;
	double blockStp = 2;
	//for(int j = 4; j <= 4; j++)
	//for(int i = 2; i <= 2; i++)
	{
		sol = new CSolver3d(90,50,980);
		sol->taskID = 8002;//i+ j*100+9000;
		sol->Init3DFDTD_UPML(60, 20);
		sol->TimeStepLoop(600);

		delete sol;
	}


	//CSolver sol("1D", 1000,1000,1000);
	//sol.Init1DSimpleFDTD();

	//sol.TimeStepLoop1D(1000);




	wrkTime = clock();
	cout << "Solve Time: "<<wrkTime/CLOCKS_PER_SEC << endl;

	system("PAUSE");
	return 1;
}