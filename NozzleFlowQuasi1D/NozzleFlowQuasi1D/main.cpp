/*
Based on book: 
COMPUTATIONAL FLUID DYNAMICS: The Basics with Applications,
John D. Anderson, Jr.
This program intends to solve quasi-1d nozzle flow problem
described in chapter 7.
*/
#include<iostream>
#include<cstdio>
#include<algorithm>
#include<cstring>
#include<cmath>
using namespace std;
typedef double Float;
const int GRIDL = 31;
const Float C = 0.5;//Courant number
const Float gamma = 1.4;//specific heats for standard day condition
class NDGrid {//all values are nondimensional
public:
	Float dx, dt;
	//"nondimensional" values described in page 297
	Float rho[GRIDL];//density field
	Float V[GRIDL];//velocity field
	Float T[GRIDL];//temperature field
	Float A[GRIDL];//diameter of nozzle at grid points
	void print(void) {
		printf("dx=%.4lf dt=%.4lf\n", dx, dt);
		printf("     x       A       rho      V       T       M\n");
		for (int i = 0; i < GRIDL; i++) {
			printf("%8.1lf%8.3lf%8.3lf%8.3lf%8.3lf%8.3lf\n", dx*i, A[i], rho[i], V[i], T[i], V[i] / sqrt(T[i]));
		}
	}
	void init(void) {
		dx = 3.0 / (GRIDL - 1);
		for (int i = 0; i < GRIDL; i++) {
			Float x = dx*i;
			A[i] = (1 + 2.2*(x - 1.5)*(x - 1.5)) / 1.0;
			rho[i] = 1.0 - 0.3146*x;
			T[i] = 1.0 - 0.2314*x;
			V[i] = (0.1 + 1.09*x)*sqrt(T[i]);
		}
		calc_dt();
	}
	void calc_dt(void) {
		dt = 1e10;
		for (int i = 0; i < GRIDL; i++) {
			Float a = sqrt(T[i]);
			double now = C*dx / (a + V[i]);
			dt = min(now, dt);
		}
	}
	void continuity_dt_fwd(Float rho_dt[]) {
		//see equation 7.51 in page 299
		for (int i = 0; i + 1 < GRIDL; i++) {
			rho_dt[i] =
				-rho[i] * (V[i + 1] - V[i]) / dx
				- rho[i] * V[i] * (log(A[i + 1]) - log(A[i])) / dx
				- V[i] * (rho[i + 1] - rho[i]) / dx;
		}
		//linear extrapolation for last point
		rho_dt[GRIDL - 1] = 2 * rho_dt[GRIDL - 2] - rho_dt[GRIDL - 3];
	}
	void continuity_dt_rwd(Float rho_dt[]) {
		//see equation 7.57 in page 300
		for (int i = 1; i < GRIDL; i++) {
			rho_dt[i] =
				-rho[i] * (V[i] - V[i - 1]) / dx
				- rho[i] * V[i] * (log(A[i]) - log(A[i - 1])) / dx
				- V[i] * (rho[i] - rho[i - 1]) / dx;
		}
		//linear extrapolation for first point
		rho_dt[0] = 2 * rho_dt[1] - rho_dt[2];
	}
	void momentum_dt_fwd(Float V_dt[]) {
		//see equation 7.52 in page 300
		for (int i = 0; i + 1 < GRIDL; i++) {
			V_dt[i] =
				-V[i] * (V[i + 1] - V[i]) / dx
				- ((T[i + 1] - T[i]) / dx + T[i] / rho[i] * (rho[i + 1] - rho[i]) / dx) / gamma;
		}
		V_dt[GRIDL - 1] = 2 * V_dt[GRIDL - 2] - V_dt[GRIDL - 3];
	}
	void momentum_dt_rwd(Float V_dt[]) {
		//see equation 7.58 in page 300
		for (int i = 1; i < GRIDL; i++) {
			V_dt[i] =
				-V[i] * (V[i] - V[i - 1]) / dx
				- ((T[i] - T[i - 1]) / dx + T[i] / rho[i] * (rho[i] - rho[i - 1]) / dx) / gamma;
		}
		V_dt[0] = 2 * V_dt[1] - V_dt[2];
	}
	void energy_dt_fwd(Float T_dt[]) {
		//see equation 7.53 in page 300
		for (int i = 0; i + 1 < GRIDL; i++) {
			T_dt[i] =
				-V[i] * (T[i + 1] - T[i]) / dx
				- (gamma - 1)*T[i] * ((V[i + 1] - V[i]) / dx + V[i] * (log(A[i + 1]) - log(A[i])) / dx);
		}
		T_dt[GRIDL - 1] = 2 * T_dt[GRIDL - 2] - T_dt[GRIDL - 3];
	}
	void energy_dt_rwd(Float T_dt[]) {
		//see equation 7.59 in page 300
		for (int i = 1; i < GRIDL; i++) {
			T_dt[i] =
				-V[i] * (T[i] - T[i - 1]) / dx
				- (gamma - 1)*T[i] * ((V[i] - V[i - 1]) / dx + V[i] * (log(A[i]) - log(A[i - 1])) / dx);
		}
		T_dt[0] = 2 * T_dt[1] - T_dt[2];
	}
	void MacCormack_Step(void) {
		static Float rho_dt0[GRIDL], V_dt0[GRIDL], T_dt0[GRIDL];
		continuity_dt_fwd(rho_dt0);
		momentum_dt_fwd(V_dt0);
		energy_dt_fwd(T_dt0);
		NDGrid H = (*this);
		for (int i = 0; i < GRIDL; i++) {
			H.rho[i] = rho[i] + dt*rho_dt0[i];
			H.V[i] = V[i] + dt*V_dt0[i];
			H.T[i] = T[i] + dt*T_dt0[i];
		}
		H.rho[0] = 1;
		H.T[0] = 1;
		static Float rho_dt1[GRIDL], V_dt1[GRIDL], T_dt1[GRIDL];
		H.continuity_dt_rwd(rho_dt1);
		H.momentum_dt_rwd(V_dt1);
		H.energy_dt_rwd(T_dt1);
		//for (int i = 0; i < GRIDL; i++) { cout << V_dt0[i] << " "; }cout << endl;
		//cout << H.V[15] << " " << H.T[15] << endl;
		//cout << rho_dt1[15] << endl;
		for (int i = 0; i < GRIDL; i++) {
			rho[i] = rho[i] + dt*(rho_dt0[i] + rho_dt1[i]) / 2.0;
			V[i] = V[i] + dt*(V_dt0[i] + V_dt1[i]) / 2.0;
			T[i] = T[i] + dt*(T_dt0[i] + T_dt1[i]) / 2.0;
		}
		rho[0] = 1;
		T[0] = 1;
	}
};
NDGrid G;
int main() {
	G.init();
	G.print();
	for (int t = 0; t < 1400; t++) {
		//printf("iteration %d\n", t);
		G.MacCormack_Step();
		//G.print();
		//printf("\n\n");
	}
	G.print();
	system("pause");
	return 0;
}