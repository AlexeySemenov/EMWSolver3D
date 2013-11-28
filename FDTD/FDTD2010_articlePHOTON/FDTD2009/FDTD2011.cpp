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
	CSolver sol("TE", 1080, 420, 0);
	sol.InitSimpleFDTD_TE_UPML();	
	sol.srcwidth = 120;
	sol.TimeStepLoopTE(10000,"K:\\Hz\\120");
	sol.DeleteSimpleFDTD_TE_UPML();

	CSolver sol2("TE", 1080, 420, 0);
	sol2.InitSimpleFDTD_TE_UPML();	
	sol2.srcwidth = 180;
	sol2.TimeStepLoopTE(10000,"K:\\Hz\\180");
	sol2.DeleteSimpleFDTD_TE_UPML();

	CSolver sol3("TE", 1080, 420, 0);
	sol3.InitSimpleFDTD_TE_UPML();	
	sol3.srcwidth = 240;
	sol3.TimeStepLoopTE(10000,"K:\\Hz\\240");
	sol3.DeleteSimpleFDTD_TE_UPML();

	CSolver sol4("TE", 1080, 420, 0);
	sol4.InitSimpleFDTD_TE_UPML();	
	sol4.srcwidth = 300;
	sol4.TimeStepLoopTE(10000,"K:\\Hz\\300");
	sol4.DeleteSimpleFDTD_TE_UPML();

	CSolver sol5("TE", 1080, 420, 0);
	sol5.InitSimpleFDTD_TE_UPML();	
	sol5.srcwidth = 360;
	sol5.TimeStepLoopTE(10000,"K:\\Hz\\360");
	sol5.DeleteSimpleFDTD_TE_UPML();

	CSolver sol6("TE", 1080, 420, 0);
	sol6.InitSimpleFDTD_TE_UPML();	
	sol6.srcwidth = 420;
	sol6.TimeStepLoopTE(10000,"K:\\Hz\\420");
	sol6.DeleteSimpleFDTD_TE_UPML();

	CSolver sol7("TE", 1080, 420, 0);
	sol7.InitSimpleFDTD_TE_UPML();	
	sol7.srcwidth = 480;
	sol7.TimeStepLoopTE(10000,"K:\\Hz\\480");
	sol7.DeleteSimpleFDTD_TE_UPML();
	//CSolver sol2("TM_UPML", 420, 1080, 0);
	//sol2.InitSimpleFDTD_TM_UPML();

	//sol2.TimeStepLoopTM(4000);




	wrkTime = clock();
	cout << "Solve Time: "<<wrkTime/CLOCKS_PER_SEC << endl;

	system("PAUSE");
	return 1;
}