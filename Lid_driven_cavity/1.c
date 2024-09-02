#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main()
{
int i,j,m,n;
printf("\nEnter the number of nodes in x direction\n");
scanf("%d",&m);
printf("\nEnter the number of nodes in y direction\n");
scanf("%d",&n);
float l=1.0,h=1.0,siold[m][n],si[m+1][n+1],Re,w[m][n],wold[m+1][n+1],u[m][n],v[m][n];
printf("\nEnter the value of Reynolds Number\n");
scanf("%f",&Re);
float dx, dy;
dx=l/(m-1);
dy=h/(n-1);
float b = dx/dy;
float dxsq,dysq,k;
dxsq = pow(dx,2);
dysq = pow(dy,2);
k = 2*((1/dxsq)+(1/dysq));
printf("%f",k);
for(i=1;i<=m;i++)
	{for(j=1;j<=n;j++)
	{	if(j == n)u[i][n] = 1.0;
		else u[i][j] = 0.0;
	}
	}
	for(i=1;i<=m;i++)
	{for(j=1;j<=n;j++)
		{v[i][j] = 0.0;
		si[i][j] = 0.0;
		
		if (j == 1)  //LEFT BC
		{
		w[i][j] = 0;
		}
		else if (j == n)  
			{
			
			w[i][n]=(-2*u[i][n])/dy; //Top BC
			}
		else if (i == 1) 
				{
				
				w[i][j] = 0; //Left BC
				}
		else if (i == m)
				{
			
				w[i][j] = 0; //Inner points
				}
		else 
			{
		
			w[i][j]= 0.0;
			}
			
			}
	}
	
	for(i=1;i<=m;i++)printf("%lf \t",w[i][n]);

	int iter = 0;
	float s_err,v_err;
	FILE *file1= fopen("Stream_Error_400.dat","w");
	
	
	for(iter=1;iter<=6000;iter++)
	{
		for(i=1;i<=m;i++)
		{
			for(j=1;j<=n;j++)
			{
				siold[i][j] = si[i][j];
				wold[i][j] = w[i][j];	//storing old points
			}
		}
		s_err=0;
		for(j=2;j<n;j++)
		{
			for(i=2;i<m;i++)
			{	
				
			
				si[i][j] = (((si[i+1][j]+si[i-1][j])/dxsq)+((si[i][j+1]+si[i][j-1])/dysq) + w[i][j])/k;
				s_err = s_err + pow((si[i][j] - siold[i][j]),2);
				 
				
			}
		}
		
		v_err = 0;
		for(j=2;j<n;j++)
		{
			for(i=2;i<n;i++)
			{
			u[i][j] = (si[i][j+1] - si[i][j-1])/(2*dy);
				v[i][j] = (si[i-1][j] - si[i+1][j])/(2*dx);
			
				w[i][j]=(Re/k)*(((w[i+1][j]+w[i-1][j])/(dxsq*Re)) + ((w[i][j+1]+w[i][j-1])/(dysq*Re))
				-((u[i][j]*(w[i+1][j]-w[i-1][j])/(2*dx)) +(v[i][j]*(w[i][j+1]-w[i][j-1])/(2*dy)))); 
			v_err = v_err + pow((w[i][j] - wold[i][j]),2);
			
				
				
				
			}
		}
			
		v_err =sqrt(fabs(v_err))/((m-2)*(n-2));
		s_err =sqrt(fabs(s_err))/((m-2)*(n-2));
	
		
		printf("Iteration %d\t",iter);
		printf("SError %.10f\tVError %.10f\n",s_err,v_err);
		fprintf(file1, "%d\t%.10f\t%.10f\n",iter,s_err,v_err);
		
		
		
		
		for(j=1;j<=n;j++)
	{for(i=1;i<=m;i++)
		{
		
		if (j == n)  
		{
		w[i][j]=-2*(si[i][n-1]-si[i][n]+(u[i][j]*dy))/dysq;
		
		}
		else if (i == m)  
			{
			w[i][j]=-2*(si[m-1][j]-si[m][j])/dxsq;
			
			}
		else if (j == 1) 
				{
				 w[i][j]=-2*(si[i][2]-si[i][0])/dysq;
				
				}
		else if (i == 1)
				{
				w[i][j]=-2*(si[2][j]-si[1][j])/dxsq;
				
				}
			
			
			}
	} 
		
	 if(s_err < pow(10,-6) && v_err < pow(10,-6) )break;
	 }
	fclose(file1);

for(j=2;j<n;j++)
			{
			for(i=2;i<m;i++)
				{
				  u[i][j]=(si[i][j+1]-si[i][j-1])/(2*dy);
        	v[i][j]=(si[i-1][j]-si[i+1][j])/(2*dx);
				}
				}

FILE *fp2;
fp2=fopen("400_Stream_234103312.plt","w");
fprintf(fp2,"VARIABLES = \"X\", \"Y\", \"PHI\"\n");
fprintf(fp2,"ZONE T = \"BLOCK1\", I = 101, J = 101, F = POINT\n\n");
for(int i=1;i<=m;i++){
	for(int j=1;j<=n;j++){
		fprintf(fp2,"%1f \t %1f \t %1f \n",i*dx,j*dy,si[i][j]);}
		}
		
FILE *fp3;
fp3=fopen("400_Vorticity_234103312.plt","w");
fprintf(fp3,"VARIABLES = \"X\", \"Y\", \"W\"\n");
fprintf(fp3,"ZONE T = \"BLOCK1\", I = 101, J = 101, F = POINT\n\n");
for(int i=1;i<=m;i++){
	for(int j=1;j<=n;j++){
		fprintf(fp3,"%1f \t %1f \t %1f \n",i*dx,j*dy,w[i][j]);}
		}
		
FILE *fp5;
fp5=fopen("Velocity_400.dat","w");

for(i=1;i<=m;i++){
	for(j=1;j<=n;j++){
		fprintf(fp5,"%f \t %f \t %f \t %f\n",i*dx,j*dy,u[i][j],v[i][j]);}
		}
		
FILE *fp8;
fp8 = fopen("all_400.dat","w");
for(i=1;i<=m;i++){
	for(j=1;j<=n;j++){
		fprintf(fp8,"%f \t %f \t %f \t %f \t %f \t %f\n",i*dx,j*dy,si[i][j],w[i][j],u[i][j],v[i][j]);}
		}
		
FILE *fp6;
fp6 = fopen("u_centerline_400.dat","w");
for( j = 1;j<=n;j++)fprintf(fp6,"%f \t%f \n",u[(m-1)/2][j],j*dy);

FILE *fp7;
fp7 = fopen("v_centerline_400.dat","w");
for( i = 1;i<=m;i++)fprintf(fp7,"%f \t%f \n",v[i][(n-1)/2],i*dx);     
printf("%f",k);
}
