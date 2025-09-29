#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"pcg_random.h"

typedef struct conf{
		int spin;
		int left;
		int right;
}Ising;

void initialize(int N, Ising *S);
void mcstep(double T, int N, Ising *S);
double magnetization(int N, Ising *S);
double energy(int N, Ising *S);

int main(int argc, char **argv)
{
		int N;
		int i;
		int j;
		int time;
		Ising *S;
		double T;
		double mag=0;
		double tot_Erg=0;
		double capacity=0.0;
		double kai=0.0; 

		init_rnd(gus());

		N=atoi(argv[1]);

		S=(Ising *)malloc(sizeof(Ising) * N);

		initialize(N, S);

		for(T=3.01; T>0.01; T-=0.01)
		{
				mag=0.0;
				tot_Erg=0.0;
				kai=0.0;
				capacity=0.0;

				for(time=0; time<200000; time++)
						mcstep(T, N, S);

				for(time=0; time<200000; time++)
				{
						mcstep(T, N, S);
						mag+=magnetization(N, S);
						tot_Erg+=energy(N, S);
						kai+=pow(magnetization(N,S), 2);
						capacity+=pow(energy(N,S), 2);
				}

				mag=mag/200000;
				tot_Erg=tot_Erg/200000;
				kai=kai/200000;
				capacity=capacity/200000;

				kai=(kai-mag*mag)/T;
				capacity=(capacity - tot_Erg*tot_Erg)/(T*T);

				kai=kai*N;
				capacity=capacity*N;

				printf("%lf\t%lf\t%lf\t%lf\t%lf\n", T, mag, tot_Erg, kai, capacity);
		}

		free(S);
		return 0;
}

void initialize(int N, Ising *S)
{
		int i;

		for(i=0; i<N; i++)
		{
				if(drnd()<0.5)
						S[i].spin=1;
				else 
						S[i].spin=-1;
		}

		// boundary codnition 
		for(i=0; i<N; i++)
		{
				if(i==0)
				{
						S[i].left=N-1;
						S[i].right=i+1;
				}
				else if(i==N-1)
				{
						S[i].left=i-1;
						S[i].right=0;
				}
				else
				{
						S[i].left=i-1;
						S[i].right=i+1;
				}
		}

}

void mcstep(double T, int N, Ising *S)
{
		int i,j;
		double E1, E2;
		double E_flip;

		for(j=0; j<N; j++)
		{
				i=(int)(drnd()*N);
				E1=-S[i].spin*(S[S[i].left].spin + S[S[i].right].spin);
				E2=-(-S[i].spin)*(S[S[i].left].spin + S[S[i].right].spin);

				E_flip=E2-E1;

				if(E_flip<0 || drnd() < exp(-E_flip/T))
						S[i].spin=-S[i].spin;
		}
}

double magnetization(int N, Ising *S)
{
		int i;
		double mag;
		mag=0;

		for(i=0; i<N; i++)
				mag+=S[i].spin;

		mag=mag/N;

		return fabs(mag);
}

double energy(int N, Ising *S)
{
		int i;
		double E_tot;
		E_tot=0.0;

		for(i=0; i<N; i++)
				E_tot+=-S[i].spin*(S[S[i].left].spin + S[S[i].right].spin);

		E_tot=E_tot/N;

		return E_tot;
}

