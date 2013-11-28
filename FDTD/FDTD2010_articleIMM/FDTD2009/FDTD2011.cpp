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
	CSolver sol("TM_UPML", 1040, 840, 0);
	sol.InitSimpleFDTD_TM_UPML();

	
	sol.TimeStepLoopTM(10000);

	//CSolver sol2("TM_UPML", 420, 1080, 0);
	//sol2.InitSimpleFDTD_TM_UPML();

	//sol2.TimeStepLoopTM(4000);




	wrkTime = clock();
	cout << "Solve Time: "<<wrkTime/CLOCKS_PER_SEC << endl;

	system("PAUSE");
	return 1;
}