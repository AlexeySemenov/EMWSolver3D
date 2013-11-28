#include "stdhead.h"
#include "Solver3d.h"

using namespace std;

CSolver3d::CSolver3d(void)
{
	data = new CField(100, 100, 100);
}
CSolver3d::CSolver3d(int SizeX, int SizeY, int SizeZ)
{
	taskID = 0;
	data = new CField(SizeX, SizeY, SizeZ);
		Fz = new double**[data->gridX];
		Gz = new double**[data->gridX];

		Fx = new double**[data->gridX];
		Gx = new double**[data->gridX];

		Fy = new double**[data->gridX];
		Gy = new double**[data->gridX];

		Bx = new double**[data->gridX];
		By = new double**[data->gridX];
		Bz = new double**[data->gridX];

		Cb = new double**[data->gridX];
		Db = new double**[data->gridX];

		CaFz = new double**[data->gridX];
		CbFz = new double**[data->gridX];
		CaGz = new double**[data->gridX];
		CbGz = new double**[data->gridX];
		CaEz = new double**[data->gridX];
		CbEz = new double**[data->gridX];
		CcEz = new double**[data->gridX];

		CaFx = new double**[data->gridX];
		CbFx = new double**[data->gridX];
		CaGx = new double**[data->gridX];
		CbGx = new double**[data->gridX];
		CaEx = new double**[data->gridX];
		CbEx = new double**[data->gridX];
		CcEx = new double**[data->gridX];

		CaFy = new double**[data->gridX];
		CbFy = new double**[data->gridX];
		CaGy = new double**[data->gridX];
		CbGy = new double**[data->gridX];
		CaEy = new double**[data->gridX];
		CbEy = new double**[data->gridX];
		CcEy = new double**[data->gridX];

		DaBx = new double**[data->gridX];
		DbBx = new double**[data->gridX];
		DaBy = new double**[data->gridX];
		DbBy = new double**[data->gridX];
		DaBz = new double**[data->gridX];
		DbBz = new double**[data->gridX];
		DaHx = new double**[data->gridX];
		DbHx = new double**[data->gridX];
		DcHx = new double**[data->gridX];
		DaHy = new double**[data->gridX];
		DbHy = new double**[data->gridX];
		DcHy = new double**[data->gridX];
		DaHz = new double**[data->gridX];
		DbHz = new double**[data->gridX];
		DcHz = new double**[data->gridX];

		Eps = new double**[data->gridX];
		Mju = new double**[data->gridX];
		Sig = new double**[data->gridX];
		SigS = new double**[data->gridX];

		Ew = new double**[data->gridX];

		for(int i = 0; i < data->gridX; i++)
		{
			Fz[i] = new double*[data->gridY];
			Gz[i] = new double*[data->gridY];
			Fx[i] = new double*[data->gridY];
			Gx[i] = new double*[data->gridY];
			Fy[i] = new double*[data->gridY];
			Gy[i] = new double*[data->gridY];

			Bx[i] = new double*[data->gridY];
			By[i] = new double*[data->gridY];
			Bz[i] = new double*[data->gridY];

			Cb[i] = new double*[data->gridY];
			Db[i] = new double*[data->gridY];


			CaFz[i] = new double*[data->gridY];
			CbFz[i] = new double*[data->gridY];
			CaGz[i] = new double*[data->gridY];
			CbGz[i] = new double*[data->gridY];
			CaEz[i] = new double*[data->gridY];
			CbEz[i] = new double*[data->gridY];
			CcEz[i] = new double*[data->gridY];

			CaFx[i] = new double*[data->gridY];
			CbFx[i] = new double*[data->gridY];
			CaGx[i] = new double*[data->gridY];
			CbGx[i] = new double*[data->gridY];
			CaEx[i] = new double*[data->gridY];
			CbEx[i] = new double*[data->gridY];
			CcEx[i] = new double*[data->gridY];

			CaFy[i] = new double*[data->gridY];
			CbFy[i] = new double*[data->gridY];
			CaGy[i] = new double*[data->gridY];
			CbGy[i] = new double*[data->gridY];
			CaEy[i] = new double*[data->gridY];
			CbEy[i] = new double*[data->gridY];
			CcEy[i] = new double*[data->gridY];

			DaBx[i] = new double*[data->gridY];
			DbBx[i] = new double*[data->gridY];
			DaBy[i] = new double*[data->gridY];
			DbBy[i] = new double*[data->gridY];
			DaBz[i] = new double*[data->gridY];
			DbBz[i] = new double*[data->gridY];
			DaHx[i] = new double*[data->gridY];
			DbHx[i] = new double*[data->gridY];
			DcHx[i] = new double*[data->gridY];
			DaHy[i] = new double*[data->gridY];
			DbHy[i] = new double*[data->gridY];
			DcHy[i] = new double*[data->gridY];
			DaHz[i] = new double*[data->gridY];
			DbHz[i] = new double*[data->gridY];
			DcHz[i] = new double*[data->gridY];

			Eps[i] = new double*[data->gridY];
			Mju[i] = new double*[data->gridY];
			Sig[i] = new double*[data->gridY];
			SigS[i] = new double*[data->gridY];

			Ew[i] = new double*[data->gridY];
		}
		for(int i = 0; i < data->gridX; i++)
			for(int j = 0; j < data->gridY; j++)
			{
				Fz[i][j] = new double[data->gridZ];
				Gz[i][j] = new double[data->gridZ];
				Fx[i][j] = new double[data->gridZ];
				Gx[i][j] = new double[data->gridZ];
				Fy[i][j] = new double[data->gridZ];
				Gy[i][j] = new double[data->gridZ];

				Bx[i][j] = new double[data->gridZ];
				By[i][j] = new double[data->gridZ];
				Bz[i][j] = new double[data->gridZ];

				Cb[i][j] = new double[data->gridZ];
				Db[i][j] = new double[data->gridZ];


				CaFz[i][j] = new double[data->gridZ];
				CbFz[i][j] = new double[data->gridZ];
				CaGz[i][j] = new double[data->gridZ];
				CbGz[i][j] = new double[data->gridZ];
				CaEz[i][j] = new double[data->gridZ];
				CbEz[i][j] = new double[data->gridZ];
				CcEz[i][j] = new double[data->gridZ];

				CaFx[i][j] = new double[data->gridZ];
				CbFx[i][j] = new double[data->gridZ];
				CaGx[i][j] = new double[data->gridZ];
				CbGx[i][j] = new double[data->gridZ];
				CaEx[i][j] = new double[data->gridZ];
				CbEx[i][j] = new double[data->gridZ];
				CcEx[i][j] = new double[data->gridZ];

				CaFy[i][j] = new double[data->gridZ];
				CbFy[i][j] = new double[data->gridZ];
				CaGy[i][j] = new double[data->gridZ];
				CbGy[i][j] = new double[data->gridZ];
				CaEy[i][j] = new double[data->gridZ];
				CbEy[i][j] = new double[data->gridZ];
				CcEy[i][j] = new double[data->gridZ];

				DaBx[i][j] = new double[data->gridZ];
				DbBx[i][j] = new double[data->gridZ];
				DaBy[i][j] = new double[data->gridZ];
				DbBy[i][j] = new double[data->gridZ];
				DaBz[i][j] = new double[data->gridZ];
				DbBz[i][j] = new double[data->gridZ];
				DaHx[i][j] = new double[data->gridZ];
				DbHx[i][j] = new double[data->gridZ];
				DcHx[i][j] = new double[data->gridZ];
				DaHy[i][j] = new double[data->gridZ];
				DbHy[i][j] = new double[data->gridZ];
				DcHy[i][j] = new double[data->gridZ];
				DaHz[i][j] = new double[data->gridZ];
				DbHz[i][j] = new double[data->gridZ];
				DcHz[i][j] = new double[data->gridZ];

				Eps[i][j] = new double[data->gridZ];
				Mju[i][j] = new double[data->gridZ];
				Sig[i][j] = new double[data->gridZ];
				SigS[i][j] = new double[data->gridZ];

				Ew[i][j] = new double[data->gridZ];
			}

}





void CSolver3d::Init3DFDTD_UPML(int bOff, double bEps)
{
	//Init timestep parameters
	S = 0.99;
	S2D = S / sqrt(2.0);
	Nx = data->gridX;
	Ny = data->gridY;
	Nz = data->gridZ;
	//dx = 0.0001;
	//dy = 0.0001;
	double Nl = 40;
	double Nlz = 120;
	double leps = bEps;
	dz = LAMBDA/Nlz;
	dx = LAMBDA/Nl;
	dy = LAMBDA/Nl;
	dxy = dx;
	dt = S* dz / (sqrt(3.0)* CC);

	if (dt > 1.0/(2.0*FREQ))
	{cout << "Warning::Computed dt is greater than Nyquist: " << dt << endl; system("PAUSE");}


	emiTime = -1;

	int Nxh = Nx/2;
	int Nyh = Ny/2;
	int Nzh = Nz/2;
	res = Nl;
	 _a = 2*Nl;//0.02 / dz;//2*res;
	 _b = Nl;////0.005 / dz;
	 L0 = Nlz;//0.02 / dz;
	 int blockSize = Nlz/10;
	 L1 = L0;
	 L2 = 2.5*res;
	 _w = 2.5;
	 int blockOffsetZ = bOff;
	int _x1 = Nxh-_a/2;
	int _x2 = Nxh+_a/2;
	int _y1 = Nyh-_b/2;

	int _y2 = Nyh+_b/2;
	startZ = Nzh - L0/2;

	int aThick = _a/2 - blockSize/2;
	int bThick = _b/2 - blockSize/2;
	int cThick = 0;//5;
	int ofs = 0;
	double gamma0 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z - M_PI*M_PI/(_a*_a));
	double gamma1 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z*leps - M_PI*M_PI/(_a*_a));
	double k0 = OMEGA*OMEGA*EPS_Z*MU_Z;
	//double Nlz = Nl;
	double Na = _a;
	double Sz = CC * dt / dz;
	double g0A = (1.0 / Sz) * sin(Sz * M_PI * k0 / (Nlz * gamma0)) ;
	//g0A *= g0A;
	double g0B = (2.0 * Na / (gamma0 * Na * dx * Nlz)) * sin(M_PI / (2.0 * Na));
	g0B *= g0B;
	double g1A = sin(OMEGA*dt/2.0) / dt;
	double g1B = sin(M_PI * dx / (2.0 * _a * dx)) / dx;
	g1A *= g1A;
	g1B *= g1B;
	gamma0_num = gamma0 * (Nlz / M_PI) * asin(sqrt(g0A - g0B)); 
	gamma0_num = (2.0/dz) * asin(dz * sqrt(EPS_Z*MU_Z * g1A - g1B));
	gamma1_num = (2.0/dz) * asin(dz * sqrt(leps*EPS_Z*MU_Z * g1A - g1B));

	double LambdaK = 2*_a*dz;
	double A1 = 1;
	//double gamma1 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z*leps - M_PI*M_PI/(_a*_a));
	double regz = cos(gamma1 * L0 * dz);
	double imgz = sin(gamma1 * L0 * dz * (1 + (gamma1/gamma0)*(gamma1/gamma0))/2.0);
	double modF = cos(gamma1*L0*dx)*cos(gamma1*L0*dx) + (gamma1*0.5/gamma0 + gamma0*0.5/gamma1)*(gamma1*0.5/gamma0 + gamma0*0.5/gamma1);
	double ReF = (cos(gamma0 * L0 * dz) * regz + sin(gamma0 * L0 * dz) * imgz)/(regz*regz + imgz*imgz);
	double ImF = (-cos(gamma0 * L0 * dz) * imgz + sin(gamma0 * L0 * dz) * regz)/(regz*regz + imgz*imgz);
	cout << "LAMBDA K = " << LambdaK << endl << "Lambda = " << LAMBDA << endl;
	cout << "k0 = " << k0*_a*dx <<endl;
	cout << "dz = " << dz << endl;
	cout << "dt = " << dt << endl;
	cout << "gridX = " << Nx << endl;
	cout << "gridY = " << Ny << endl;
	cout << "gridZ = " << Nz << endl;
	cout << "gamma0 = " << gamma0 << endl;
	cout << "gamma1 = " << gamma1 << endl;
	cout << "gamma0 numerical = " << gamma0_num << endl;
	cout << "gamma1 numerical = " << gamma1_num << endl;
	cout << "T0 = " << 2*M_PI/OMEGA/dt << endl;
	cout << "ReF = " << ReF << endl;
	cout << "ImF = " << ImF << endl;
	system("pause");


	for(int i = 0; i < data->gridX; i++)
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
		{

			Eps[i][j][k] = 1.0;
			Mju[i][j][k] = 1.0;
			Sig[i][j][k] = 0.0;
			SigS[i][j][k] = 0.0;
			
			if ((i>=_x1)&&(i<=_x2))
				if ((j>=_y1)&&(j<_y2))
				{
					if((k>=startZ) && (k <=(startZ + L0)))
					{
						Eps[i][j][k] = 2;

						if ((i>=_x1+aThick)&&(i<=_x2-aThick))
							if ((j>=_y1+bThick)&&(j<_y2-bThick))
							{
								if((k>=startZ+blockOffsetZ) && (k <=(startZ + blockSize +blockOffsetZ)))
								{
									Eps[i][j][k] = leps;
									//Sig[i][j][k] = 0.567;

									/*if ((i>=_x1+3*aThick)&&(i<=_x2-3*aThick))
										if ((j>=_y1+3*bThick)&&(j<_y2-3*bThick))
											if((k>=startZ+3*cThick) && (k <=(startZ + L0-3*cThick)))
												Eps[i][j][k] = leps;*/
								}
							}
					}

					//if((i>startZ+ofs*8) && (i <=(startZ+(ofs)*8+8)))
					//{
					//	Eps[i][j][k] = 21;
					//}
	
				}


			

			Gz[i][j][k] = 0.0;
			Fz[i][j][k] = 0.0;
			Gx[i][j][k] = 0.0;
			Fx[i][j][k] = 0.0;
			Gy[i][j][k] = 0.0;
			Fy[i][j][k] = 0.0;
			Bx[i][j][k] = 0.0;
			By[i][j][k] = 0.0;
			Bz[i][j][k] = 0.0;
		}



	//Init UPML
		double mu_r_x_1, mu_r_x_2, eps_r_x_1,eps_r_x_2;
		
		double boundaryWidth, boundaryHeight;

		mu_r_x_1 = 1;
		mu_r_x_2 = 1;
		eps_r_x_1 = 1;
		eps_r_x_2 = 1;
		pmlUD = 0;
		pmlLR = 0;
		pmlTB = 0;

		boundaryWidth = (double)pmlLR*dx;
		boundaryHeight = (double)pmlUD*dy;
		

		ka_max = 1.0;
		exponent = 4;
		R_err = 1e-16;
		double eta_1 = sqrt(MU_Z*mu_r_x_1/EPS_Z/eps_r_x_1);
		double eta_2 = sqrt(MU_Z*mu_r_x_2/EPS_Z/eps_r_x_2);
		sigma_max_1= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryWidth);
		sigma_max_2= -log(R_err) * (exponent + 1.0) / (2.0 * IMP_Z * boundaryHeight);;
		boundaryFactor = sigma_max_1 / ( dx * (pow(boundaryWidth,exponent)) * (exponent + 1));

		tw = 26.53e-12;
		t0 = 4*tw;

	//Init FDTDSolver Consts
		
	for(int i = 0; i < data->gridX; i++)
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
		{
			Cb[i][j][k] = (dt / (Eps[i][j][k]* EPS_Z)) / (1.0 + (Sig[i][j][k] * dt)/(2 * Eps[i][j][k]* EPS_Z));

			Db[i][j][k] = (dt / (Mju[i][j][k]* MU_Z)) / (1.0 + (SigS[i][j][k] * dt)/(2 * Mju[i][j][k]* MU_Z));


			CaFz[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFz[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGz[i][j][k] = (2*EPS_Z * K_x(i) - Sig_x(i,false)*dt)/
				(2*EPS_Z * K_x(i) + Sig_x(i,false)*dt);
			CbGz[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_x(i) + Sig_x(i,false)*dt);

			CaEz[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,false)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);
			CbEz[i][j][k] = (2*EPS_Z * K_z(k) + Sig_z(k,true)*dt)/
				(2*EPS_Z* K_y(j) + Sig_y(j,false)*dt);
			CcEz[i][j][k] = (2*EPS_Z*K_z(k)-Sig_z(k,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,false)*dt);

			CaFx[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFx[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGx[i][j][k] = (2*EPS_Z * K_y(j) - Sig_y(j,false)*dt)/
				(2*EPS_Z * K_y(i) + Sig_y(j,false)*dt);
			CbGx[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_y(j) + Sig_y(j,false)*dt);

			CaEx[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,false)*dt)/
				(2*EPS_Z*K_z(k) + Sig_z(k,false)*dt);
			CbEx[i][j][k] = (2*EPS_Z * K_x(i) + Sig_x(i,true)*dt)/
				(2*EPS_Z* K_z(k) + Sig_z(k,false)*dt);
			CcEx[i][j][k] = (2*EPS_Z*K_x(i)-Sig_x(i,true)*dt)/
				(2*EPS_Z*K_z(k) + Sig_z(k,false)*dt);

			CaFy[i][j][k] = (2*EPS_Z * Eps[i][j][k] - Sig[i][j][k]*dt)/
				(2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);
			CbFy[i][j][k] = 2 * dt / (2*EPS_Z * Eps[i][j][k] + Sig[i][j][k]*dt);

			CaGy[i][j][k] = (2*EPS_Z * K_z(k) - Sig_z(k,false)*dt)/
				(2*EPS_Z * K_z(k) + Sig_z(k,false)*dt);
			CbGy[i][j][k] = 2*EPS_Z/(2*EPS_Z * K_z(k) + Sig_z(k,false)*dt);

			CaEy[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,false)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);
			CbEy[i][j][k] = (2*EPS_Z * K_y(j) + Sig_y(j,true)*dt)/
				(2*EPS_Z* K_x(i) + Sig_x(i,false)*dt);
			CcEy[i][j][k] = (2*EPS_Z*K_y(j)-Sig_y(j,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,false)*dt);

			DaBx[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			DbBx[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			DaHx[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,true)*dt)/
				(2*EPS_Z*K_z(k) + Sig_z(k,true)*dt);
			DbHx[i][j][k] = (2*EPS_Z*K_x(i) + Sig_x(i,false)*dt)/
				((2*EPS_Z*K_z(k) + Sig_z(k,true)*dt)*MU_Z*Mju[i][j][k]);
			DcHx[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,false)*dt)/
				((2*EPS_Z*K_z(k) + Sig_z(k,true)*dt)*MU_Z*Mju[i][j][k]);


			DaBy[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,true)*dt)/
				(2*EPS_Z*K_z(k) - Sig_z(k,true)*dt);
			DbBy[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_z(k) + Sig_z(k,true)*dt);
			DaHy[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);
			DbHy[i][j][k] = (2*EPS_Z*K_y(j) + Sig_y(j,false)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)*MU_Z*Mju[i][j][k]);
			DcHy[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,false)*dt)/
				((2*EPS_Z*K_x(i) + Sig_x(i,true)*dt)*MU_Z*Mju[i][j][k]);

			DaBz[i][j][k] = (2*EPS_Z*K_x(i) - Sig_x(i,true)*dt)/
				(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);
			DbBz[i][j][k] = 2*EPS_Z*dt/(2*EPS_Z*K_x(i) + Sig_x(i,true)*dt);
			DaHz[i][j][k] = (2*EPS_Z*K_y(j) - Sig_y(j,true)*dt)/
				(2*EPS_Z*K_y(j) + Sig_y(j,true)*dt);
			DbHz[i][j][k] = (2*EPS_Z*K_z(k) + Sig_z(k,false)*dt)/
				((2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)*MU_Z*Mju[i][j][k]);
			DcHz[i][j][k] = (2*EPS_Z*K_z(k) - Sig_z(k,false)*dt)/
				((2*EPS_Z*K_y(j) + Sig_y(j,true)*dt)*MU_Z*Mju[i][j][k]);
			

		}


	//system("PAUSE");
}

void CSolver3d::FDTD3DSolver_UPML(int timestep)
{
	double rbc, rbc1, tmp1, tmp2;
	//double kvec = 2*M_PI/LAMBDA;
	//TFSF part
	int Nxh = data->gridX/2;
	int Nyh = data->gridY/2;
	int Nzh = data->gridZ/2;
	//double gamma0 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z - M_PI*M_PI/(_a*_a));
	double expS = 1.0;
//	if(timestep <=120)
//		expS = 1.0/exp(double(120-timestep));
	int _x1 = Nxh-_a/2;
	int _x2 = Nxh+_a/2;
	int _y1 = Nyh-_b/2;
	int _y2 = Nyh+_b/2;
	int sourceOffset = 0;
	//tmp1 = sqrt(OMEGA*OMEGA/(CC*CC) - M_PI*M_PI/(_a*_a*dx*dx));

	//if (timestep == 0)
	//{
	//	for(int i = _x1; i <_x2; i++)
	//		for(int j = _y1; j <_y2; j++)
	//		{
				//data->Ey[i][j][sourceOffset] = expS*sin(M_PI*(i-_x1)/((double)_a))*cos(OMEGA*(timestep-1)*dt);
				//data->Hx[i][j][sourceOffset] = -expS*gamma0*sin(M_P`I*(i-_x1+0.5)/((double)_a))*cos(OMEGA*(timestep-0.5)*dt)/(OMEGA*MU_Z);
				//data->Hz[i][j][sourceOffset] = -expS*(M_PI/_a)*cos(M_PI*(i-_x1+0.5)/((double)_a))*sin(OMEGA*(timestep-0.5)*dt)/(OMEGA*MU_Z);

				//data->Ey[i][j][1] = data->Ey[i][j][0];
				//data->Hx[i][j][1] = data->Hx[i][j][0];
				//data->Hz[i][j][1] = data->Hz[i][j][0];
	//		}


	//}
	
	double Ey0 = M_PI/2.274;
	//double Hx0 = (-dt / (MU_Z * dz)) * Ey0 * sin(gamma0_num * dz / 2.0) / sin(OMEGA * dt / 2.0);
	//double Hz0 = (dt / (MU_Z * dx)) * Ey0 * sin(M_PI * dx / (2.0 * _a * dx)) / sin(OMEGA * dt / 2.0);
	//Update E - Fields
	for(int i = 1; i < (data->gridX - 1); i++)
		for(int j = 1; j < (data->gridY - 1); j++)
			#pragma omp parallel for 
			for(int k = 1; k < (data->gridZ - 1); k++)
		{
			double Fx_pr, Fy_pr, Fz_pr, Gx_pr, Gy_pr, Gz_pr; 

			Fx_pr = Fx[i][j][k];
			Gx_pr = Gx[i][j][k];
			Fx[i][j][k] =  CaFx[i][j][k]*Fx[i][j][k] + CbFx[i][j][k]*((data->Hz[i][j][k] - data->Hz[i][j-1][k])/dy - 
					(data->Hy[i][j][k] - data->Hy[i][j][k-1])/dz);
			Gx[i][j][k] = CaGx[i][j][k]*Gx[i][j][k] + CbGx[i][j][k]*(Fx[i][j][k]-Fx_pr);
			data->Ex[i][j][k] = CaEx[i][j][k]*data->Ex[i][j][k] + CbEx[i][j][k]*Gx[i][j][k] -
				CcEx[i][j][k]*Gx_pr;

			Fy_pr = Fy[i][j][k];
			Gy_pr = Gy[i][j][k];
			Fy[i][j][k] =  CaFy[i][j][k]*Fy[i][j][k] + CbFy[i][j][k]*((data->Hx[i][j][k] - data->Hx[i][j][k-1])/dz - 
					(data->Hz[i][j][k] - data->Hz[i-1][j][k])/dx);
			Gy[i][j][k] = CaGy[i][j][k]*Gy[i][j][k] + CbGy[i][j][k]*(Fy[i][j][k]-Fy_pr);
			data->Ey[i][j][k] = CaEy[i][j][k]*data->Ey[i][j][k] + CbEy[i][j][k]*Gy[i][j][k] -
				CcEy[i][j][k]*Gy_pr;
			
			Fz_pr = Fz[i][j][k];
			Gz_pr = Gz[i][j][k];
			Fz[i][j][k] =  CaFz[i][j][k]*Fz[i][j][k] + CbFz[i][j][k]*((data->Hy[i][j][k] - data->Hy[i-1][j][k])/dx - 
					(data->Hx[i][j][k] - data->Hx[i][j-1][k])/dy);
 			Gz[i][j][k] = CaGz[i][j][k]*Gz[i][j][k] + CbGz[i][j][k]*(Fz[i][j][k]-Fz_pr);
			data->Ez[i][j][k] = CaEz[i][j][k]*data->Ez[i][j][k] + CbEz[i][j][k]*Gz[i][j][k] -
				CcEz[i][j][k]*Gz_pr;



			}

			//TFSF Correction
			for(int i = 1; i < (data->gridX - 1); i++)
				for(int j = 1; j < (data->gridY - 1); j++)
				{
					#pragma omp parallel for
					for(int k = 1; k < (data->gridZ - 1); k++)
					{
						double Eyinc1 = Ey0 * cos(OMEGA * (timestep + 0.5) * dt - gamma0_num * k * dz) * sin(M_PI * (i - _x1) * dx / (_a * dx));
						double Eyinc0 = Ey0 * cos(OMEGA * (timestep - 0.5) * dt - gamma0_num * k * dz) * sin(M_PI * (i - _x1) * dx / (_a * dx));
						Ew[i][j][k] = Eyinc1;
						data->Ey[i][j][k] = data->Ey[i][j][k] - (1.0 - 1.0/Eps[i][j][k])*
							(Eyinc1 - Eyinc0) - (dt * Sig[i][j][k] / (Eps[i][j][k]*EPS_Z))*Eyinc1;
					}
				}

	//=======Waveguide
	//for(int i = _x1+1; i <_x2-1; i++)
	//	for(int j = _y1+1; j <_y2-1; j++)
	//	{
	//		data->Ey[i][j][sourceOffset] = expS*sin(M_PI*(i-_x1)/((double)_a))*sin(OMEGA*timestep*dt);
	//		//data->Ey[i][j][1] = data->Ey[i][j][0];
	//	}	
	//
	for(int i = _x1-1; i <= _x2+1; i++)
		for(int k = sourceOffset; k <data->gridZ; k++)
		{
			data->Ez[i][_y1][k] = 0;
			data->Ez[i][_y2][k] = 0;

			data->Ex[i][_y1][k] = 0;
			data->Ex[i][_y2][k] = 0;		

		}
		for(int j = _y1-1; j <=_y2+1; j++)
			for(int k = sourceOffset; k <data->gridZ; k++)
			{
				data->Ez[_x1][j][k] = 0;
				data->Ez[_x2][j][k] = 0;

				data->Ey[_x1][j][k] = 0;
				data->Ey[_x2][j][k] = 0;
			}

	/*for(int i = Nxh-10; i <= Nxh+10; i++)
		for(int j = Nyh-10; j <= Nyh+10; j++)
			data->Ez[i][j][90] = ((2*M_PI)/LAMBDA)*cos(OMEGA*timestep*dt);*/
	//Update H - Field
	for(int i = 1; i < (data->gridX - 1); i++)
		for(int j = 1; j < (data->gridY - 1); j++)
		{
			#pragma omp parallel for
			for(int k = 1; k < (data->gridZ - 1); k++)
			{
				double Bx_pr, By_pr, Bz_pr;
			Bx_pr = Bx[i][j][k];
			Bx[i][j][k] = DaBx[i][j][k]* Bx[i][j][k] + DbBx[i][j][k]*
				((data->Ez[i][j][k] - data->Ez[i][j+1][k])/dy + (data->Ey[i][j][k+1] - data->Ey[i][j][k])/dz);

			data->Hx[i][j][k] = DaHx[i][j][k]*data->Hx[i][j][k]+ DbHx[i][j][k]*Bx[i][j][k]-
				DcHx[i][j][k]*Bx_pr;

			By_pr = By[i][j][k];
			By[i][j][k] = DaBy[i][j][k]* By[i][j][k] + DbBy[i][j][k]*
				((data->Ez[i+1][j][k] - data->Ez[i][j][k])/dx+(data->Ex[i][j][k] - data->Ex[i][j][k+1])/dz);

			data->Hy[i][j][k] = DaHy[i][j][k]*data->Hy[i][j][k]+ DbHy[i][j][k]*By[i][j][k]-
				DcHy[i][j][k]*By_pr;

			Bz_pr = Bz[i][j][k];
			Bz[i][j][k] = DaBz[i][j][k]* Bz[i][j][k] + DbBz[i][j][k]*
				((data->Ey[i][j][k] - data->Ey[i+1][j][k])/dx+(data->Ex[i][j+1][k] - data->Ex[i][j][k])/dy);

			data->Hz[i][j][k] = DaHz[i][j][k]*data->Hz[i][j][k]+ DbHz[i][j][k]*Bz[i][j][k]-
				DcHz[i][j][k]*Bz_pr;


			}
		}
		//!!!!WAVEGUIDE H == ZERO
		/*for(int i = 1; i < (data->gridX - 1); i++)
			for(int j = 1; j < (data->gridY - 1); j++)
			{
				#pragma omp parallel for
				for(int k = 1; k < (data->gridZ - 1); k++)
				{
					double Hxinc1 = Hx0 * cos(OMEGA * (timestep + 1) * dt - gamma0_num * k * dz) * sin(M_PI * (i - _x1) * dx / (_a * dx));
					double Hxinc0 = Hx0 * cos(OMEGA * (timestep) * dt - gamma0_num * k * dz) * sin(M_PI * (i - _x1) * dx / (_a * dx));
					double Hzinc1 = Hz0 * sin(OMEGA * (timestep + 1) * dt - gamma0_num * k * dz) * cos(M_PI * (i - _x1) * dx / (_a * dx));
					double Hzinc0 = Hz0 * sin(OMEGA * (timestep) * dt - gamma0_num * k * dz) * cos(M_PI * (i - _x1) * dx / (_a * dx));
					data->Hx[i][j][k] = data->Hx[i][j][k] - (1.0 - 1.0/Mju[i][j][k])*
						(Hxinc1 - Hxinc0) - (dt * SigS[i][j][k] / (Mju[i][j][k] * MU_Z)) * Hxinc1;
					data->Hz[i][j][k] = data->Hz[i][j][k] - (1.0 - 1.0/Mju[i][j][k])*
						(Hzinc1 - Hzinc0) - (dt * SigS[i][j][k] / (Mju[i][j][k] * MU_Z)) * Hzinc1 ;
				}
			}*/
		//=====================Waveguide
		//for(int i = _x1; i <_x2; i++)
		//	for(int j = _y1; j <_y2; j++)
		//	{
		//		//data->Hx[i][j][sourceOffset] = -expS*gamma0*sin(M_PI*(i-_x1+0.5)/((double)_a))*cos(OMEGA*(timestep+0.5)*dt)/(OMEGA*MU_Z);
		//		//data->Hz[i][j][sourceOffset] = -expS*(M_PI/_a)*cos(M_PI*(i-_x1+0.5)/((double)_a))*sin(OMEGA*(timestep+0.5)*dt)/(OMEGA*MU_Z);
		//		//data->Hz[i][j][1] = data->Hz[i][j][0];
		//		//data->Hx[i][j][1] = data->Hx[i][j][0];
		//		//data->Hx[i][j][1] = data->Hx[i][j][0];
		//		//data->Hz[i][j][1] = data->Hz[i][j][0];
		//	}
		for(int i = _x1-1; i <= _x2+1; i++)
			for(int k = sourceOffset; k <data->gridZ; k++)
			{
				data->Hz[i][_y1-1][k] = 0;
				data->Hz[i][_y2][k] = 0;

				data->Hx[i][_y1-1][k] = 0;
				data->Hx[i][_y2][k] = 0;		

			}
		for(int j = _y1-1; j <=_y2+1; j++)
			for(int k = sourceOffset; k <data->gridZ; k++)
			{
				data->Hz[_x1-1][j][k] = 0;
				data->Hz[_x2][j][k] = 0;

				data->Hy[_x1-1][j][k] = 0;
				data->Hy[_x2][j][k] = 0;
			}
		
				
}

void CSolver3d::TimeStepLoop(int N)
{
	for(int t = 0; t < N; t++)
	{
		cout << "Solving Step: " << t+1 << " of " << N; 
		FDTD3DSolver_UPML(t);
		//HardSourceMIncidence(50,50, t);
		//Simple1DFDTDSolver(t);
		//FDTD1DIncidence(0, t);
		cout << " ...done." << endl;
		//data1D->InitNetCDF(t, dx);
	   // dataTE->InitNetCDF(t, dx);
		//if(t>=250000)
		//if(t%400 ==0) dataTM->WriteToNetCDF(t/100,20,20,20,20);

		//if(t%2 ==0)dataTM->WriteToNetCDF(t/2,20,20,20,20);
		//int i = 0;
		//Debugout(DaBz, i); i++;
		//Debugout(DbBz, i); i++;
		//Debugout(DaHz, i); i++;
		//Debugout(DbHz, i); i++;
		//Debugout(DcHz, i); i++;

		//Debugout(DaBx, i); i++;
		//Debugout(DbBx, i); i++;
		//Debugout(DaHx, i); i++;
		//Debugout(DbHx, i); i++;
		//Debugout(DcHx, i); i++;

		//Debugout(DaBy, i); i++;
		//Debugout(DbBy, i); i++;
		//Debugout(DaHy, i); i++;
		//Debugout(DbHy, i); i++;
		//Debugout(DcHy, i); i++;

		//data->WriteToBinAmp(t, 5, 5, 5, 5, 0, 0);
		if((t>=350))
		{
			//data->WriteToBinEx(t 5, 5, 5, 5, 0, 0);
			//if(t%5 == 0)
			data->WriteToBinEy(taskID, t, 0, 0, 0, 0, 0, 0, Ew);
			//data->WriteToBinEz(t, 5, 5, 5, 5, 0, 0);
			//if(t%5 == 0)
			//data->WriteToBinHx(t/5, 0, 0, 0, 0, 0, 0);
			//data->WriteToBinHy(t, 5, 5, 5, 5, 0, 0);
			//data->WriteToBinHz(t, 5, 5, 5, 5, 0, 0);
		}
		//Debugout();
	}

}

void CSolver3d::Debugout(double*** A, int index)
{
	string fileName;
	stringstream s;
	s << "K:\\Ez\\debug["<<index<<"].dat";
	fileName = s.str();
	ofstream file(fileName.c_str(), ios::out | ios::binary);
	double tmp;

	if(!file.is_open())
		cout << "ERROR::Cannot open file for write!";

	cout << "Writing to file... ";
	for(int i = 0; i < data->gridX; i++)
		for(int j = 0; j < data->gridY; j++)
			for(int k = 0; k < data->gridZ; k++)
			{
				tmp = A[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
	cout << data->gridX*data->gridY*data->gridZ*sizeof(double) << " BYTE written succesfully" << endl;
	file.close();
}

double CSolver3d::source(double t)
{
	const double E0 = 1.0;
	return E0*sin(OMEGA*dt*t);
}

double CSolver3d::PlaneWaveLookUpd(int i, int j)
{
	return (1.0/S2DS1D)*(CosSource_Angle*(i-Oia)+SinSource_Angle*(j-Oja))+Offset;
}

double CSolver3d::interp1(double* y, double xi)
{
	double xi2,temp ;
	xi2 = xi - (int)xi;
	temp = (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1];
	return (1-xi2)*y[int(xi)] + xi2*y[int(xi)+1]	;
}

double CSolver3d::Sig_x(double x,bool offset)
{
	double tmp;
	double x1;
	double x2;
	double dist = 0;
	//offset = !offset;

	if((x > pmlLR-1)&&(x < (data->gridX - pmlLR))) tmp = 0;
	else 
	{
		if(x <= pmlLR-1)
		{
			dist = pmlLR - x;
			if(offset)
			{
				x1 = (dist + 0.0) * dx;       // upper bounds for point i
				x2 = (dist - 1.0) * dx;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dx;       // upper bounds for point i
				x2 = (dist - 0.5) * dx;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)));   //   polynomial grading
		}
		if(x >= (data->gridX - pmlLR)) 
		{
			dist = x - (data->gridX - pmlLR);
			if(offset)
			{
				x1 = (dist + 1.0) * dx;       // upper bounds for point i
				x2 = (dist + 0.0) * dx;       // lower bounds for point i
			}
			else
			{
				x1 = (dist + 0.5) * dx;       // upper bounds for point i
				x2 = (dist - 0.5) * dx;       // lower bounds for point i
				if(dist == 0)
				{
					x1 = 0;
					x2 = 0;
				}
			}
			tmp = boundaryFactor * (pow(x1,(exponent+1)) - pow(x2,(exponent+1)) );   //   polynomial grading
		}
	}
		//gradientK = 1.0 + (kmax - 1.0)* (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );

	if(tmp<0) cout << "ERROR IN PML" <<endl;
	return tmp;
}
double CSolver3d::Sig_y(double y, bool offset)
{
	double tmp;
	double y1;
	double y2;
	double dist = 0;
	//offset = !offset;
	if((y > pmlUD-1)&&(y < (data->gridY - pmlUD)))
		tmp = 0;
	else 
	{
		if(y <= pmlUD-1)
		{
			dist = pmlUD - y;
			if(offset)
			{
				y1 = (dist + 0.0) * dy;       // upper bounds for point i
				y2 = (dist - 1.0) * dy;       // lower bounds for point i
			}
			else
			{
				y1 = (dist + 0.5) * dy;       // upper bounds for point i
				y2 = (dist - 0.5) * dy;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(y1,(exponent+1)) - pow(y2,(exponent+1)) );   //   polynomial grading
		}
		if(y >= (data->gridY - pmlUD)) 
		{
			dist = y - (data->gridY - pmlUD);
			if(offset)
			{
				y1 = (dist + 1.0) * dy;       // upper bounds for point i
				y2 = (dist + 0.0) * dy;       // lower bounds for point i
			}
			else
			{
				y1 = (dist + 0.5) * dy;       // upper bounds for point i
				y2 = (dist - 0.5) * dy;       // lower bounds for point i
				if(dist == 0)
				{
					y1 = 0;
					y2 = 0;
				}
			}
			tmp = boundaryFactor * (pow(y1,(exponent+1)) - pow(y2,(exponent+1)) );   //   polynomial grading
		}
	}
		//gradientK = 1.0 + (kmax - 1.0)* (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );
	if(tmp<0) cout << "ERROR IN PML" <<endl;
	return tmp;
}


double CSolver3d::Sig_z(double z, bool offset)
{
		double tmp;
	double z1;
	double z2;
	double dist = 0;

	if((z > pmlTB-1)&&(z < (data->gridZ - pmlTB))) 
		tmp = 0;
	else 
	{
		if(z <= pmlTB-1)
		{
			dist = pmlTB - z;
			if(offset)
			{
				z1 = (dist + 0.0) * dz;       // upper bounds for point i
				z2 = (dist - 1.0) * dz;       // lower bounds for point i
			}
			else
			{
				z1 = (dist + 0.5) * dz;       // upper bounds for point i
				z2 = (dist - 0.5) * dz;       // lower bounds for point i
			}
			tmp = boundaryFactor * (pow(z1,(exponent+1)) - pow(z2,(exponent+1)));   //   polynomial grading
		}
		if(z >= (data->gridZ - pmlTB)) 
		{
			dist = z - (data->gridZ - pmlTB);
			if(offset)
			{
				z1 = (dist + 1.0) * dz;       // upper bounds for point i
				z2 = (dist + 0.0) * dz;       // lower bounds for point i
			}
			else
			{
				z1 = (dist + 0.5) * dz;       // upper bounds for point i
				z2 = (dist - 0.5) * dz;       // lower bounds for point i
				if(dist == 0)
				{
					z1 = 0;
					z2 = 0;
				}
			}
			tmp = boundaryFactor * (pow(z1,(exponent+1)) - pow(z2,(exponent+1)) );   //   polynomial grading
		}
	}
		//gradientK = 1.0 + (kmax - 1.0)* (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );

	if(tmp<0) cout << "ERROR IN PML" <<endl;
	return tmp;
}

double CSolver3d::K_x(double x)
{
	//double tmp;
	//if( (x > dx*pmlLR)&&(x < dx*(data->gridX - pmlLR)))tmp = 1.0;
	//else if(x <= dx*pmlLR) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR - x)/((double) pmlLR*dx) ,exponent);
	//else if(x >= dx*(data->gridX - pmlLR)) tmp = 1.0 + (ka_max - 1.0)*pow( (dx*pmlLR + x - dx*data->gridX)/((double) pmlLR*dx) ,exponent);
	//if (tmp!=1)cout << "K_X = " << tmp <<endl;
	return 1;
}
double CSolver3d::K_y(double y)
{
	//double tmp;
	//if( (y > dy*pmlUD)&&(y < dy*(data->gridY - pmlUD)))tmp = 1.0;
	//else if(y <= dy*pmlUD) tmp = 1.0 + (ka_max - 1.0)*pow( (dy*pmlUD - y)/((double) pmlUD*dy) ,exponent);
	//else if(y >= dy*(data->gridY - pmlUD)) tmp = (1.0 + (ka_max - 1.0)*pow( (dy*pmlUD + y - dy*data->gridY)/((double) pmlUD*dy) ,exponent));
	//if (tmp!=1)cout << "K_Y = " << tmp <<endl;
	return 1;
}
double CSolver3d::K_z(double z)
{
	return 1;
}
CSolver3d::~CSolver3d(void)
{

		
		for(int i = 0; i < data->gridX; i++)
			for(int j = 0; j < data->gridY; j++)
			{
				delete[] Fz[i][j];
				delete[] Gz[i][j];
				delete[] Fx[i][j];
				delete[] Gx[i][j];
				delete[] Fy[i][j];
				delete[] Gy[i][j];

				delete[] Bx[i][j];
				delete[] By[i][j];
				delete[] Bz[i][j];

				delete[] Cb[i][j];
				delete[] Db[i][j];


				delete[] CaFz[i][j];
				delete[] CbFz[i][j];
				delete[] CaGz[i][j];
				delete[] CbGz[i][j];
				delete[] CaEz[i][j];
				delete[] CbEz[i][j];
				delete[] CcEz[i][j];

				delete[] CaFx[i][j];
				delete[] CbFx[i][j];
				delete[] CaGx[i][j];
				delete[] CbGx[i][j];
				delete[] CaEx[i][j];
				delete[] CbEx[i][j];
				delete[] CcEx[i][j];

				delete[] CaFy[i][j];
				delete[] CbFy[i][j];
				delete[] CaGy[i][j];
				delete[] CbGy[i][j];
				delete[] CaEy[i][j];
				delete[] CbEy[i][j];
				delete[] CcEy[i][j];

				delete[] DaBx[i][j];
				delete[] DbBx[i][j];
				delete[] DaBy[i][j];
				delete[] DbBy[i][j];
				delete[] DaBz[i][j];
				delete[] DbBz[i][j];
				delete[] DaHx[i][j];
				delete[] DbHx[i][j];
				delete[] DcHx[i][j];
				delete[] DaHy[i][j];
				delete[] DbHy[i][j];
				delete[] DcHy[i][j];
				delete[] DaHz[i][j];
				delete[] DbHz[i][j];
				delete[] DcHz[i][j];

				delete[] Eps[i][j];
				delete[] Mju[i][j];
				delete[] Sig[i][j];
				delete[] SigS[i][j];
				delete[] Ew[i][j];

			}
for(int i = 0; i < data->gridX; i++)
		{
			delete[] Fz[i];
			delete[] Gz[i];
			delete[] Fx[i];
			delete[] Gx[i];
			delete[] Fy[i];
			delete[] Gy[i];

			delete[] Bx[i];
			delete[] By[i];
			delete[] Bz[i];

			delete[] Cb[i];
			delete[] Db[i];


			delete[] CaFz[i];
			delete[] CbFz[i];
			delete[] CaGz[i];
			delete[] CbGz[i];
			delete[] CaEz[i];
			delete[] CbEz[i];
			delete[] CcEz[i];

			delete[] CaFx[i];
			delete[] CbFx[i];
			delete[] CaGx[i];
			delete[] CbGx[i];
			delete[] CaEx[i];
			delete[] CbEx[i];
			delete[] CcEx[i];

			delete[] CaFy[i];
			delete[] CbFy[i];
			delete[] CaGy[i];
			delete[] CbGy[i];
			delete[] CaEy[i];
			delete[] CbEy[i];
			delete[] CcEy[i];

			delete[] DaBx[i];
			delete[] DbBx[i];
			delete[] DaBy[i];
			delete[] DbBy[i];
			delete[] DaBz[i];
			delete[] DbBz[i];
			delete[] DaHx[i];
			delete[] DbHx[i];
			delete[] DcHx[i];
			delete[] DaHy[i];
			delete[] DbHy[i];
			delete[] DcHy[i];
			delete[] DaHz[i];
			delete[] DbHz[i];
			delete[] DcHz[i];

			delete[] Eps[i];
			delete[] Mju[i];
			delete[] Sig[i];
			delete[] SigS[i];
			delete[] Ew[i];
		}
		delete[] Fz;
		delete[] Gz;

		delete[] Fx;
		delete[] Gx;

		delete[] Fy;
		delete[] Gy;

		delete[] Bx;
		delete[] By;
		delete[] Bz;

		delete[] Cb;
		delete[] Db;

		delete[] CaFz;
		delete[] CbFz;
		delete[] CaGz;
		delete[] CbGz;
		delete[] CaEz;
		delete[] CbEz;
		delete[] CcEz;

		delete[] CaFx;
		delete[] CbFx;
		delete[] CaGx;
		delete[] CbGx;
		delete[] CaEx;
		delete[] CbEx;
		delete[] CcEx;

		delete[] CaFy;
		delete[] CbFy;
		delete[] CaGy;
		delete[] CbGy;
		delete[] CaEy;
		delete[] CbEy;
		delete[] CcEy;

		delete[] DaBx;
		delete[] DbBx;
		delete[] DaBy;
		delete[] DbBy;
		delete[] DaBz;
		delete[] DbBz;
		delete[] DaHx;
		delete[] DbHx;
		delete[] DcHx;
		delete[] DaHy;
		delete[] DbHy;
		delete[] DcHy;
		delete[] DaHz;
		delete[] DbHz;
		delete[] DcHz;

		delete[] Eps;
		delete[] Mju;
		delete[] Sig;
		delete[] SigS;
		delete[] Ew;
		delete data;// = new CField(SizeX, SizeY, SizeZ);
}
