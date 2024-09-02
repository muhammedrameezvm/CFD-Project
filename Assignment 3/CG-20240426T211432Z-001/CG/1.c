#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void conjugate_gradient(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N]);
void initialisation(int N,int n,double T[N],double b[N],double A_p[N],double A_e[N],double A_w[N],double A_n[N],double A_s[N],double beta_sq,double T_top,double T_bottom,double T_left,double T_right);
double dot_product(int n,double A[n],double B[n]);
void print_final_temp(int N,int n,double T[N],double T_top,double T_bottom,double T_left,double T_right,double dx,double dy);


void main()
{
int N_x,N_y;
double T_top,T_bottom,T_left,T_right,L,H,dx,dy,beta,beta_sq;

//--------------------reading inputs--------------------
FILE *fp1;
fp1 = fopen("input.txt","r");


fscanf(fp1,"Number of grids in x-direction = %d\n",&N_x);
fscanf(fp1,"Number of grids in y-direction = %d\n",&N_y);
fscanf(fp1,"Temperature at top = %lf\n",&T_top);
fscanf(fp1,"Temperature at bottom = %lf\n",&T_bottom);
fscanf(fp1,"Temperature at left = %lf\n",&T_left);
fscanf(fp1,"Temperature at right = %lf\n",&T_right);
fscanf(fp1,"Width of the domain = %lf\n",&L);
fscanf(fp1,"Height of the domain = %lf\n",&H);

//printf("%lf\n",L);

//--------------------Declaring 1-d array size for lexicographic arrangement --------------
int N = N_x * N_y; //--------------All grid points (128 x 128)-------
int N_sq = N * N;	//----------size of matrix A (128 x 128 X 28 x 128)-----------------
double T[N+1],b[N+1],A_p[N+1],A_e[N+1],A_w[N+1],A_n[N+1],A_s[N+1]; //---------------declaring 5 diagonals of A matrix (Ap , Ae , Aw , An , As)

dx = L / N_x;
dy = H / N_y;


beta = dx / dy;
beta_sq = beta * beta;
//printf("%d",N_sq);

initialisation(N+1,N_x,T,b,A_p,A_e,A_w,A_n,A_s,beta_sq,T_top,T_bottom,T_left,T_right);

conjugate_gradient(N+1,N_x,T,b,A_p,A_e,A_w,A_n,A_s);

print_final_temp(N+1,N_x,T,T_top,T_bottom,T_left,T_right,dx,dy);










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
  	//b[i] = 0;
  	
  	
  	//-----------------Left boundary----------------
  	if(i % n == 1 )
  		{
  		A_p[i] += 1;
  		A_w[i] = 0;
  		b[i] = 2*T_left;
  		}
  	//----------------Right boundary----------------
  	if(i % n == 0 )
  		{
  		A_p[i] += 1;
  		A_e[i] = 0;
  		b[i] = 2*T_right;
  		}
  		
  	//----------------Bottom boundary------------
  	if(i >= N-n)
  		{
  		A_p[i] += beta_sq;
  		A_s[i] = 0;
  		b[i] = 2*T_bottom; 
  		}
  	//-----------------Top boundary----------------
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
{double Ax[N],r[N],r_new[N],d[N],Ad[N],alpha,beta,residue;
int iteration = 1;
FILE *fp2;
fp2 = fopen("CG_Iter_vs_residue.dat","w");



    //--------------------calculating AT or Ax , residue , direction------------------
    for(int i =1;i<N;i++)
    {
        Ax[i] = A_n[i]*T[i-n] + A_w[i]*T[i-1] + A_p[i]*T[i] + A_e[i]*T[i+1] + A_s[i]*T[i+n]; 
        r[i] = b[i] - Ax[i];
        d[i] = r[i]; 
    }
    
    
    
//fprintf(fp2,"r\t b \t Ax \n");
//for(int i =1;i<N;i++)fprintf(fp2,"%.3lf\t %.3lf \t %.3lf \n",r[i],b[i],Ax[i]);


do{
//-------------------calculating Ad[i]------------------
for(int i =1;i<N;i++)
    {
        Ad[i] = A_n[i]*d[i-n] + A_w[i]*d[i-1] + A_p[i]*d[i] + A_e[i]*d[i+1] + A_s[i]*d[i+n];
    }

//---------------calculating Alpha---------------
alpha = dot_product(N,r,r) / dot_product(N,d,Ad);
    
//---------------calculating new temp----------------
for(int i =1;i<N;i++)T[i] += alpha*d[i];

if(iteration % 50 == 0){
for(int i =1;i<N;i++)
{
Ax[i] = A_n[i]*T[i-n] + A_w[i]*T[i-1] + A_p[i]*T[i] + A_e[i]*T[i+1] + A_s[i]*T[i+n]; 
r_new[i] = b[i] - Ax[i];
}
}
//------------------calculating new residue---------------
else for(int i =1;i<N;i++)r_new[i] = r[i] - alpha*Ad[i];

//--------------------calculating beta--------------------
beta = dot_product(N,r_new,r_new) / dot_product(N,r,r);

residue = 0;
//----------------------calculating new direction-----------------------
for(int i =1;i<N;i++)
	{
	d[i] = r_new[i] + beta*d[i];
	residue += pow(r_new[i],2);
	r[i] = r_new[i];
	}
residue = sqrt(residue);
if(iteration % 1 == 0)printf("\nIteration = %d ,Residue = %lf ",iteration,residue); 
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
void print_final_temp(int N,int n,double T[N],double T_top,double T_bottom,double T_left,double T_right,double dx,double dy)
{int j=1,k=1;
FILE *fp2;
fp2 = fopen("T.plt","w");
fprintf(fp2,"VARIABLES = \"X\", \"Y\", \"T\"\n");
fprintf(fp2,"ZONE T = \"BLOCK1\", I = 128, J = 128, F = POINT\n\n");
for(int i = 1; i <= n+1 ; i++)fprintf(fp2,"%1f \t %1f \t %1f \n",0*dx,(i-1)*dy,T_top);
fprintf(fp2,"%1f \t %1f \t %1f \n",j*dx,0*dy,T_left);
for(int i = 1; i < N-n ; i++)
	{if (i%n != 0){fprintf(fp2,"%1f \t %1f \t %1f \n",j*dx,k*dy,0.25*(T[i]+T[i+1]+T[i+n]+T[i+n+1]));
	k++;}
	else {fprintf(fp2,"%1f \t %1f \t %1f \n%1f \t %1f \t %1f \n",j*dx,k*dy,T_right,(j+1)*dx,0*dy,T_left);
	j++;
	k = 1;}
	}
for(int i = N-n+1; i <= N ; i++){fprintf(fp2,"%1f \t %1f \t %1f \n",n*dx,(k-1)*dy,T_bottom);
k++;}

}	

