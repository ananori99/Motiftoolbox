
/*
 * FitzHugh Nagumo oscillator according to Schwabedal and Pikovsky, PRE 2010
 *
 * Modified on Sept 8, 2014 to accomodate for the Ana/Aaron version of the equations used to study Post Inhibitory Rebound
 */


#include <math.h>
#include <stdio.h>


#define N_EQ1		2
#define NP1		6
#define N_EQ3		(3*N_EQ1)
#define NP3		(3*NP1)
#define N_EQ4		(4*N_EQ1)
#define NP4		(4*NP1)
#define NC		6


#define COUPLING_THRESHOLD	0. //
#define THRESHOLD_SLOPE		100. //

extern "C" {

double boltzmann(const double V, const double V_0, const double k)
{
	return 1./(1.+exp(-k*(V-V_0))); // this boltzman is unchanged, but will have to be applied more in modified eqns - Ana
}


void derivs_one(double* y, double* dxdt, const double* p)
{
	dxdt[0] = p[5]*(y[0]-y[0]*y[0]*y[0])-y[1]+p[0]*(boltzman(cell2,0.3,20) + boltzman(cell3, 0.3,20))*boltzman(self,-0.3,-20); // x' = m (x-x^3)-y+I // I did not remove m, but I is now multipled
	// by the boltz dependent on all 3 cells - Ana
	dxdt[1] = p[1]*(y[0] + 0.574 - 0.01*y[1]); // y' = e (Bfun(x, x_0, k_2) - y) // this e is the time scale diffrence 
	// parameter? new eqn : e*(x +a -b*y)
	//FIND OUT WHAT a,b REPRESENT IN MATLAB CODE AND SWITCH THAT IN
	// for now, I put in constants from my prior PIR summer work. Just to check the code will run...
}


void integrate_one_rk4(double* y, const double* params, double* output, const double dt, const unsigned N, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ1], y2[N_EQ1], k1[N_EQ1], k2[N_EQ1], k3[N_EQ1], k4[N_EQ1];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<N_EQ1; j++) output[j] = y[j];
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			derivs_one(y, k1, params);
			for(k=0; k<N_EQ1; k++) y1[k] = y[k]+k1[k]*dt2; 			
			derivs_one(y1, k2, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_one(y2, k3, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k3[k]*dt; 			
			derivs_one(y2, k4, params);
			for(k=0; k<N_EQ1; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<N_EQ1; j++) output[N_EQ1*i+j] = y[j];
	}
}


void derivs_three(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_factor[3];
	bm_factor[0] = boltzmann(y[0],       COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[1] = boltzmann(y[N_EQ1],   COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[2] = boltzmann(y[2*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);

	dxdt[0] +=       (p[4]-y[0])*            (kij[0]*bm_factor[1] + kij[1]*bm_factor[2]);
	dxdt[N_EQ1] +=   (p[4+NP1]-y[N_EQ1])*    (kij[2]*bm_factor[0] + kij[3]*bm_factor[2]);
	dxdt[2*N_EQ1] += (p[4+2*NP1]-y[2*N_EQ1])*(kij[4]*bm_factor[0] + kij[5]*bm_factor[1]);
}


void derivs_four(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_factor[4];
	bm_factor[0] = boltzmann(y[0],       COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[1] = boltzmann(y[N_EQ1],   COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[2] = boltzmann(y[2*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[3] = boltzmann(y[3*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);
	derivs_one(y+3*N_EQ1, dxdt+3*N_EQ1, p+3*NP1);

	dxdt[0] +=       (p[4]-y[0])*            (kij[0]*bm_factor[1] + kij[1]* bm_factor[2] + kij[2]* bm_factor[3]);
	dxdt[N_EQ1] +=   (p[4+NP1]-y[N_EQ1])*    (kij[3]*bm_factor[0] + kij[4]* bm_factor[2] + kij[5]* bm_factor[3]);
	dxdt[2*N_EQ1] += (p[4+2*NP1]-y[2*N_EQ1])*(kij[6]*bm_factor[0] + kij[7]* bm_factor[1] + kij[8]* bm_factor[3]);
	dxdt[3*N_EQ1] += (p[4+3*NP1]-y[3*N_EQ1])*(kij[9]*bm_factor[0] + kij[10]*bm_factor[1] + kij[11]*bm_factor[2]);
}


void integrate_three_rk4(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ3], y2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3];

	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<3; j++)
		output[j] = y[j*N_EQ1];

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_three(y, k1, params, kij);

			for(k=0; k<N_EQ3; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_three(y1, k2, params, kij);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k2[k]*dt2;

			derivs_three(y2, k3, params, kij);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_three(y2, k4, params, kij);

			for(k=0; k<N_EQ3; k++)
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<3; j++)
			output[3*i+j] = y[j*N_EQ1]; 					
	}
}


void integrate_four_rk4(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ4], y2[N_EQ4], k1[N_EQ4], k2[N_EQ4], k3[N_EQ4], k4[N_EQ4];

	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<4; j++)
		output[j] = y[j*N_EQ1]; 	

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_four(y, k1, params, kij);

			for(k=0; k<N_EQ4; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_four(y1, k2, params, kij);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k2[k]*dt2;

			derivs_four(y2, k3, params, kij);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_four(y2, k4, params, kij);

			for(k=0; k<N_EQ4; k++)
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<4; j++)
			output[4*i+j] = y[j*N_EQ1]; 					
	}
}



void derivs_n(const unsigned num_osci, double* y, double* dxdt, const double* p, const double g_inh)
{
	unsigned i, j;
	double bm_factor[num_osci], bmf_sum;
	for(i=0; i<num_osci; i++)
	{
		bm_factor[i] = boltzmann(y[i*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);
		derivs_one(y+i*N_EQ1, dxdt+i*N_EQ1, p);
	}
	for(i=0; i<num_osci; i++)
	{
		bmf_sum = 0.;
		for(j=0; j<i; j++)
		{
			bmf_sum += bm_factor[j];
		}
		for(j=i+1; j<num_osci; j++)
		{
			bmf_sum += bm_factor[j];
		}
		dxdt[i*N_EQ1] += g_inh*(p[4]-y[i*N_EQ1])*bmf_sum;
	}
}



void step_n_em(const unsigned num_osci, double* y, const double* params, const double g_inh, const double dt, const unsigned stride, double* noise)
{
	unsigned i, j, k;
	unsigned n_eq = num_osci*N_EQ1;
	double sqdt=sqrt(dt), right_hand_side[n_eq];

	for(i=0; i<stride; i++)
	{
		derivs_n(num_osci, y, right_hand_side, params, g_inh);

		for(j=0; j<n_eq; j++)
		{
			y[j] += right_hand_side[j]*dt; 			
		}
		for(j=0; j<num_osci; j++)
		{
			y[N_EQ1*j] += noise[i*num_osci+j]*sqdt; 			
		}
	}
}

}
