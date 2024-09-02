#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void conjugate_gradient(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N]);
void initialisation(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N],double beta_sq,double T_top,double T_bottom,double T_left,double T_right);
double dot_product(int n,double A[n],double B[n]);
void print_final_temp(int N,int n,double T[N],double T_top,double T_bottom,double T_left,double T_right);
void M_inverse(int N,int n,double M_inverse_r[N],double A_p[N],double A_n[N],double A_w[N],double A_e[N],double A_s[N],double r[N]);


void main()
{
int N_x,N_y;
double T_top,T_bottom,T_left,T_right,L,H,dx,dy,beta,beta_sq;

//---------------reading inputs----------
FILE *fp1;
fp1 = fopen("input.txt","r");
FILE *fp2;
fp2 = fopen("T.csv","w");

fscanf(fp1,"Number of grids in x-direction = %d\n",&N_x);
fscanf(fp1,"Number of grids in y-direction = %d\n",&N_y);
fscanf(fp1,"Temperature at top = %lf\n",&T_top);
fscanf(fp1,"Temperature at bottom = %lf\n",&T_bottom);
fscanf(fp1,"Temperature at left = %lf\n",&T_left);
fscanf(fp1,"Temperature at right = %lf\n",&T_right);
fscanf(fp1,"Width of the domain = %lf\n",&L);
fscanf(fp1,"Height of the domain = %lf\n",&H);

//printf("%lf\n",L);

//---------------Declaring 1-d array size for lexicographic arrangement ---------------
int N = N_x * N_y; //---------------All grid points (128 x 128)--------------
int N_sq = N * N;	//----------size of matrix A (128 x 128 X 28 x 128)
double T[N+1],b[N+1],A_p[N+1],A_e[N+1],A_w[N+1],A_n[N+1],A_s[N+1]; //------declaring 5 diagonals of A matrix (Ap , Ae , Aw , An , As)------

dx = L / N_x;
dy = H / N_y;


beta = dx / dy;
beta_sq = beta * beta;
//printf("%d",N_sq);

initialisation(N+1,N_x,T,b,A_p,A_e,A_w,A_n,A_s,beta_sq,T_top,T_bottom,T_left,T_right);

conjugate_gradient(N+1,N_x,T,b,A_p,A_e,A_w,A_n,A_s);

print_final_temp(N+1,N_x,T,T_top,T_bottom,T_left,T_right);








}
void initialisation(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N],double beta_sq,double T_top,double T_bottom,double T_left,double T_right)
{FILE *fp1;
fp1 = fopen("initial.csv","w");
    //printf("%lf",beta_sq);
    
  for(int i = 1; i< N ; i++)
  	{
  	T[i] = 0;
  	A_p[i] = 2*(1 + beta_sq);
  	A_n[i] = -beta_sq;
  	A_s[i] = -beta_sq;
  	A_e[i] = -1;
  	A_w[i] = -1;
  	b[i] = 0;
  	
  	
  	//--------------Left boundary--------------
  	if(i % n == 1)
  		{
  		A_p[i] += 1;
  		A_w[i] = 0;
  		b[i] = 2*T_left;
  		}
  	//--------------Right boundary--------------
  	if(i % n == 0)
  		{
  		A_p[i] += 1;
  		A_e[i] = 0;
  		b[i] = 2*T_right;
  		}
  		
  	//--------------Bottom boundary--------------
  	if(i >= N-n)
  		{
  		A_p[i] += beta_sq;
  		A_s[i] = 0;
  		b[i] = 2*T_bottom; 
  		}
  		//--------------Top boundary--------------
  	if(i <= n)
  		{
  		A_p[i] += beta_sq;
  		A_n[i] = 0;
  		b[i] = 2*T_top; 
  		}
  		
  	}
 fprintf(fp1,"A_p\t A_e \t A_w \t A_n \t A_s \t T \t b\n");
  for(int i = 1; i< N ; i++)fprintf(fp1,"%.3lf\t %.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf\n",A_p[i],A_e[i],A_w[i],A_n[i],A_s[i],T[i],b[i]);
 
}

void conjugate_gradient(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N])
{double Ax[N],r[N],r_new[N],d[N],Ad[N],M_inverse_r[N],M[N],M_inverse_r_new[N],alpha,beta,residue;
int iteration = 1;
FILE *fp2;
fp2 = fopen("ILU_P_CG_Iter_vs_residue.dat","w");

//----------------------Jacobi Preconditioner-----------------------



    //--------------calculating AT or Ax--------------
    for(int i =1;i<N;i++)
    {
        Ax[i] = A_n[i]*T[i-n] + A_w[i]*T[i-1] + A_p[i]*T[i] + A_e[i]*T[i+1] + A_s[i]*T[i+n]; 
        r[i] = b[i] - Ax[i];
        //d[i] = r[i]/M[i]; //----------------M^-1--------------
    }
M_inverse(N,n,d,A_p,A_n,A_w,A_e,A_s,r);
fprintf(fp2,"r\t b \t Ax \n");
for(int i =1;i<N;i++)fprintf(fp2,"%.3lf\t %.3lf \t %.3lf \n",r[i],b[i],Ax[i]);    

M_inverse(N,n,M_inverse_r,A_p,A_n,A_w,A_e,A_s,r);

do{
//--------------calculating Ad[i]--------------
for(int i =1;i<N;i++)
    {
        Ad[i] = A_n[i]*d[i-n] + A_w[i]*d[i-1] + A_p[i]*d[i] + A_e[i]*d[i+1] + A_s[i]*d[i+n];
         //---------------M^-1--------------
    }

//--------------calculating Alpha--------------
alpha = dot_product(N,r,M_inverse_r) / dot_product(N,d,Ad);
    //printf("\nAlpha = %lf",alpha);
    
//--------------calculating new temp--------------
for(int i =1;i<N;i++)T[i] += alpha*d[i];

if(iteration % 50 == 0){
for(int i =1;i<N;i++)
{
Ax[i] = A_n[i]*T[i-n] + A_w[i]*T[i-1] + A_p[i]*T[i] + A_e[i]*T[i+1] + A_s[i]*T[i+n]; 
r_new[i] = b[i] - Ax[i];
//M_inverse_r_new[i] = r_new[i] / M[i];
}
M_inverse(N,n,M_inverse_r_new,A_p,A_n,A_w,A_e,A_s,r_new);
}
//--------------calculating new residue--------------
else for(int i =1;i<N;i++)
{r_new[i] = r[i] - alpha*Ad[i];
//M_inverse_r_new[i] = r_new[i] / M[i]; //--------------Jacobi precondition--------------new M^-1--------------
}
M_inverse(N,n,M_inverse_r,A_p,A_n,A_w,A_e,A_s,r_new);
//--------------calculating beta--------------
beta = dot_product(N,r_new,M_inverse_r_new) / dot_product(N,r,M_inverse_r);

residue = 0;
//--------------calculating new direction--------------
for(int i =1;i<N;i++)
	{
	M_inverse_r[i] = M_inverse_r_new[i];
	d[i] = M_inverse_r_new[i] + beta*d[i];
	residue += pow(r_new[i],2);
	r[i] = r_new[i];
	}
residue = sqrt(residue);
printf("\nIteration = %d ,Residue = %lf ",iteration,residue); 
fprintf(fp2,"%d \t%lf \n",iteration,residue); 
iteration++;
if(iteration > 1000)break;
}while (residue > 1e-6);



}

double dot_product(int n,double A[n],double B[n])
{double res = 0;
for(int i=1;i<n;i++) res += A[i]*B[i];
return res;
}

void exact_temp(int x,int y,double T_exact[x][y])
{
for(int i = 1;i<x;i++)
	{
	for(int j =1; j < y;j++)
		{
		}
	}
}
void print_final_temp(int N,int n,double T[N],double T_top,double T_bottom,double T_left,double T_right)
{FILE *fp2;
fp2 = fopen("T.csv","w");
for(int i = 1; i <= n+1 ; i++)fprintf(fp2,"%.4lf \t",T_top);
fprintf(fp2,"\n%.4lf \t",T_left);
for(int i = 1; i < N-n ; i++)
	{if (i%n != 0)fprintf(fp2,"%.4lf \t",0.25*(T[i]+T[i+1]+T[i+n]+T[i+n+1]));
	else fprintf(fp2,"%.4lf\n%.4lf \t",T_right,T_left);
	}
for(int i = N-n+1; i <= N ; i++)fprintf(fp2,"%.4lf \t",T_bottom);

}
void M_inverse(int N,int n,double M_inverse_r[N],double A_p[N],double A_n[N],double A_w[N],double A_e[N],double A_s[N],double r[N])
{
double L_n[N],L_w[N],L_p[N],U_e[N],U_s[N],R[N];
FILE *fp1;
fp1 = fopen("ILU.csv","w");

//printf("------------\nN = %d\n----------------",N);
for(int i=1; i < N;i++)
	{
	L_n[i] = A_n[i];
	L_w[i] = A_w[i];
	if(i == 1) L_p[i] = A_p[i];
	else if( i <= n) L_p[i] = A_p[i] - (L_w[i] * U_e[i-1]); 
		else L_p[i] = A_p[i] - (L_w[i] * U_e[i-1]) - (L_n[i] * U_s[i-n]);
	U_e[i] = A_e[i] / L_p[i];
	U_s[i] = A_s[i] / L_p[i];
	}



for(int i=1; i < N;i++){
if(i == 1) R[i] = r[i] / L_p[i] ;
	else if( i <= n) R[i] = (r[i] - (L_w[i] * R[i-1])) / L_p[i]; 
		else R[i] = (r[i] - (L_w[i] * R[i-1]) - (L_n[i] * R[i-n]))/L_p[i];
}

for(int i= N-1; i >= 1;i--)
{
if(i == N-1) M_inverse_r[i] = R[i];
	else if( i >= N-n-1) M_inverse_r[i] = R[i] - (U_e[i] * M_inverse_r[i+1]); 
		else M_inverse_r[i] = R[i] - (U_e[i] * M_inverse_r[i+1]) - (U_s[i] * M_inverse_r[i+n]);


}
fprintf(fp1,"L_p\t L_n \t L_w \t U_e \t U_s \t R \t delta\n");
  for(int i = 1; i< N ; i++)fprintf(fp1,"%.3lf\t %.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf\n",L_p[i],L_n[i],L_w[i],U_e[i],U_s[i],R[i],M_inverse_r[i]);

}
