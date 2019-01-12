#include<stdio.h>
#include<math.h>
#include<stdlib.h>
const double mi=940.0,Isp=315.0*9.8;
double mu;
//total steps
#define N 376
//time step
double t0=0,tf=375.0;
double Td[N-1];
double h;
//initial conditions
double rf=1735100.0,thetaf=0.20944,uf=0.0,vf=0.0;//init1/*ri=1755010.0,thetai=(double)0,ui=(double)1690,vi=(double)0,rmoon=1735000.0,*/
double x0[4]={1755010.0,0.0,1690.0,0.0};//init1
//variable used to store state
double x[N][4];
//control matrices
double ux[N],uy[N];
//error bound 10^-4
double eb=1;//angleconstr
//variables used in weight
double w0[16]={1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0};
double W[N][16];
//to calculate B_s sensitivity matrix
double B[4][2]={{0,0},{0,0},{1,0},{0,1}};
//double R[2][2]={{1,0},{0,1}};
double R[2][2];
double Rinv[2][2];
double A_lambda[4][4];
double B_lambda[4][1];
//vaariable for mass
double newm;
double mass[N];
//variables for thrust
double T[N];
double beta[N];
double B_s1[4][2], Bs[N][8],Blf[4][1], BsT[2][4];
//static double ts[N],tc[N];
int itr;
//functions
void ode4state();
void ode4weight();
void ode4mass();
void r_mat(int t,double r[2][2]);
void matmul(double *a,double *b,double *c,int m,int n,int p);
void retdfdx(double *dfdx,double ri,double ui,double vi);
void matinv(double *a,double *inv);
void main()
{
	double ur[2][1],u_intr[2][1],u_intr2[4][1],u_intr3[2][4];
	double err[4][1];
	double wx[4][4],ksw1[4][2],ksw1t[2][4],ks[2][4],alamb1[4][4];
	double Alambdainv[4][4];
	int i,j,k,m;
	mu=4.90278*(pow(10,12));
	h=(tf-t0)/(N-1);
	T[0]=0.0;//required
	for(j=0;j<N;j++)
		{
			ux[j]=0.0;
			uy[j]=0.0;
			}
	//iterations start--not converging at 4 with float instead of double
	for(itr=0;itr<10;itr++)
	{
	//	printf("itring\n");
		ode4state();
		err[0][0]=rf-x[N-1][0];
		err[1][0]=thetaf-x[N-1][1];
		err[2][0]=uf-x[N-1][2];
		err[3][0]=vf-x[N-1][3];
/*	for(j=0;j<4;j++)
		{
			for(k=0;k<1;k++)
			{
				printf("%lf  ",err[j][k]);
			}
			printf(";%d\n",itr);
		}*/
		if((fabs(err[0][0])<=eb)&&(fabs(err[1][0])<=eb)&&(fabs(err[2][0])<=eb)&&(fabs(err[3][0])<=eb))
		{
				ode4mass();
			for(i=0;i<N;i++)
			{
				T[i]=mass[i]*(sqrt(pow(ux[i],2)+pow(uy[i],2)));
				if(i>0)
				{
				Td[i-1]=(T[i]-T[(i-1)]);
				}
				beta[i]=(atan2(uy[i],ux[i])*180/3.14159);
				if(beta[i]<=0)
				{
				beta[i]+=360.0;
			    }
			//	printf("%lf \n",beta[i]);
			}
			printf("convergence complete at %d \n",itr);
			break;
		}
		ode4weight();
			for(k=0;k<4;k++)//assigning zero
			{
				for(m=0;m<4;m++)
				{
					A_lambda[m][k]=0.0;
				}
				B_lambda[k][0]=0.0;
			}
		for(j=0;j<N;j++)
		{
			for(k=0;k<4;k++)//reshaping
			{
				for(m=0;m<4;m++)
				{
					wx[m][k]=W[j][(4*k)+m];
				}
	    	}
	    	matmul((double*)wx,(double*)B,(double*)B_s1,4,4,2);

		   	Bs[j][0]=B_s1[0][0];
	    	Bs[j][1]=B_s1[1][0];
	    	Bs[j][2]=B_s1[2][0];
	    	Bs[j][3]=B_s1[3][0];
	    	Bs[j][4]=B_s1[0][1];
	    	Bs[j][5]=B_s1[1][1];
	    	Bs[j][6]=B_s1[2][1];
	    	Bs[j][7]=B_s1[3][1];
	    	//B_lf(:,:,iji) =  reshape(B_s2,4,2)*u(iji,:)' ; matlab code
	    	ur[0][0]=ux[j];
	    	ur[1][0]=uy[j];
	    	matmul((double*)B_s1,(double*)ur,(double*)Blf,4,2,1);
	    	r_mat((j+1),R);
	    /*	double ds8=1.0,ds;
		double p1=5.536*pow(10,-8);
	double p2=-8.397*pow(10,-5);
	double p3=0.04809;
	double p4=-12.41;
	double p5=1242.2;
	//ds = p1*tf^4 + p2*tf^3 + p3*tf^2 + p4*tf + p5 ;
	ds=(p1*pow(tf,4))+(p2*pow(tf,3))+(p3*pow(tf,2))+(p4*pow(tf,1))+p5;
	R[0][0]=ds8*exp(ds*x[N-1-j+1][1]);
			R[1][1]=ds8*exp(ds*x[N-1-j+1][1]);
			R[0][1]=0.0;R[1][0]=0.0;*/
	       	for(k=0;k<2;k++)//reshaping
			{
				for(m=0;m<4;m++)
				{
					ksw1[m][k]=Bs[j][(4*k)+m];
					ksw1t[k][m]=Bs[j][(4*k)+m];
				}

			}
			Rinv[0][0]=R[1][1]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[0][1]=-R[0][1]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[1][0]=-R[1][0]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[1][1]=R[0][0]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
	    	matmul((double*)Rinv,(double*)ksw1t,(double*)ks,2,2,4);//define ks
	    	matmul((double*)ksw1,(double*)ks,(double*)alamb1,4,2,4);
	       	for(k=0;k<4;k++)//updating Alambda
			{
				for(m=0;m<4;m++)
				{
				A_lambda[m][k]=A_lambda[m][k]+(alamb1[m][k]*h);	//precision error upto tens place

				}
				B_lambda[k][0]=B_lambda[k][0]+(Blf[k][0]*h);//checked proper precision

	    	}
	   		}
		for(j=0;j<4;j++)
		{
			B_lambda[j][0]=B_lambda[j][0]+err[j][0];//changed it-checked for precison
		}
		matinv((double*)A_lambda,(double*)Alambdainv);
		matmul((double*)Alambdainv,(double*)B_lambda,(double*)u_intr2,4,4,1);
		for(j=0;j<N;j++)
		{
			r_mat((j+1),R);
			Rinv[0][0]=R[1][1]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[0][1]=-R[0][1]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[1][0]=-R[1][0]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			Rinv[1][1]=R[0][0]/((R[0][0]*R[1][1])-(R[0][1]*R[1][0]));
			for(k=0;k<2;k++)//reshaping
			{
				for(m=0;m<4;m++)
				{
					BsT[k][m]=Bs[j][(4*k)+m];
				}
			}
			matmul((double*)Rinv,(double*)BsT,(double*)u_intr3,2,2,4);
			matmul((double*)u_intr3,(double*)u_intr2,(double*)u_intr,2,4,1);
			ux[j]=u_intr[0][0];
			uy[j]=u_intr[1][0];
		}

	}
//	printf("done\n");

}
void matmul(double *a,double *b,double *c,int m,int n,int p)
{
	int i,j,k;
	for(i=0;i<m;i++)
	{
		for(j=0;j<p;j++)
		{
			*(c+(p*i)+j)=0.0;
		}
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<p;j++)
		{
			for(k=0;k<n;k++)
			{
				(*(c+(i*p+j)))+=((*(a+(i*n+k)))*(*(b+(k*p+j))));
			}
		}
	}
}
void ode4state()
{
	double sum[4],k1[4],k2[4],k3[4],k4[4];
	double newx[4];
	int t;
	int j;
	double inter;
	//printf("in ode4state\n");
	//int t;
	//int j;
	//double inter;
	//double sum[4],k1[4],k2[4],k3[4],k4[4];
	//double newx[4];
	memcpy(&newx,&x0,sizeof(x0));
	for (j=0;j<4;j++)
		{
			x[0][j]=newx[j];
		}
		for(t=1;t<N;t++)
	{
		//	printf("stateforloop\n");

		k1[0]=h*(newx[3]);
		k1[1]=h*(newx[2]/newx[0]);
		inter=ux[t];
		k1[2]=h*(-((newx[2]*newx[3])/newx[0])+inter);
		inter=uy[t];
		k1[3]=h*(((newx[2]*newx[2])/newx[0])-(mu/(newx[0]*newx[0]))+inter);
		for (j=0;j<4;j++)
		{
			sum[j]=newx[j]+(k1[j]*0.5);
		}
		k2[0]=h*(sum[3]);
		k2[1]=h*(sum[2]/sum[0]);
		inter=ux[t];
		k2[2]=h*(-((sum[2]*sum[3])/sum[0])+inter);
		inter=uy[t];
		k2[3]=h*(((sum[2]*sum[2])/sum[0])-(mu/(sum[0]*sum[0]))+inter);
		for (j=0;j<4;j++)
		{
			sum[j]=newx[j]+(k2[j]*0.5);
		}
		//CALCULATION of K3
		k3[0]=h*(sum[3]);
		k3[1]=h*(sum[2]/sum[0]);
		inter=ux[t];
		k3[2]=h*(-((sum[2]*sum[3])/sum[0])+inter);
		inter=uy[t];
		k3[3]=h*(((sum[2]*sum[2])/sum[0])-(mu/(sum[0]*sum[0]))+inter);
		for (j=0;j<4;j++)
		{
			sum[j]=newx[j]+(k3[j]);
			}
			//CALCULATION of K4
		k4[0]=h*(sum[3]);
		k4[1]=h*(sum[2]/sum[0]);
		inter=ux[t];
		k4[2]=h*(-((sum[2]*sum[3])/sum[0])+inter);
	    inter=uy[t];
		k4[3]=h*(((sum[2]*sum[2])/sum[0])-(mu/(sum[0]*sum[0]))+inter);
		for(j=0;j<4;j++)
		{
			newx[j]+=((k1[j]+(2.0*k2[j])+(2.0*k3[j])+k4[j])/6.0);
        }
		for (j=0;j<4;j++)
		{
			x[t][j]=newx[j];
		}
	}
}
void ode4weight()
{
	int t,i,j;
	double r1,u1,v1;
	double sum[16],k1[16],k2[16],k3[16],k4[16];
	double neww[16];
	double dw[4][4];
	double dw16[16];
	double w[4][4];
	double dfdx[4][4];
	memcpy(&neww,&w0,sizeof(w0));
/*	for(i=0;i<N;i++)
	{
		r_m[i]=x[i][0];
		u_m[i]=x[i][2];
		v_m[i]=x[i][3];
	}*/
for(i=0;i<4;i++)//reshaping
		{
			for(j=0;j<4;j++)
			{
				w[j][i]=-neww[(4*i)+j];
				W[N-1][(4*i)+j]=neww[(4*i)+j];
			}
	    }
	for(t=N-2;t>=0;t--)
	{
	    r1=x[N-1-t][0];
	    u1=x[N-1-t][2];
	    v1=x[N-1-t][3];
	    retdfdx((double*)dfdx,r1,u1,v1);
		matmul((double*)w,(double*)dfdx,(double*)dw,4,4,4);
		//obtain k1
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				dw16[(4*i)+j]=dw[j][i];
			}
		}
	    for(i=0;i<16;i++)
	    {
	    	k1[i]=-h*dw16[i];//k1 values
	    	sum[i]=neww[i]+(k1[i]*0.5);//w0+(k1/2)
		}
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				w[j][i]=-sum[(4*i)+j];//w now holds -(w0+(k1/2))
			}
	    }
	    retdfdx((double*)dfdx,r1,u1,v1);
		matmul((double*)w,(double*)dfdx,(double*)dw,4,4,4);
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				dw16[(4*i)+j]=dw[j][i];
			}
	    }
	    for(i=0;i<16;i++)
	    {
	    	k2[i]=-h*dw16[i];//k2 values
	    	sum[i]=neww[i]+(k2[i]*0.5);//w0+(k2/2)
		}
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				w[j][i]=-sum[(4*i)+j];//w now holds w0+(k2/2)
			}
	    }
		//obtain k3
		retdfdx((double*)dfdx,r1,u1,v1);
		matmul((double*)w,(double*)dfdx,(double*)dw,4,4,4);
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				dw16[(4*i)+j]=dw[j][i];//dw=f(t+(h/2),w0+(k2/2))
			}
	    }
	    for(i=0;i<16;i++)
	    {
	    	k3[i]=-h*dw16[i];//k3 values
	    	sum[i]=neww[i]+k3[i];//w0+k3
		}
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				w[j][i]=-sum[(4*i)+j];//w now holds w0+k3
			}
	    }
	    //obtain k4
	    retdfdx((double*)dfdx,r1,u1,v1);
	    matmul((double*)w,(double*)dfdx,(double*)dw,4,4,4);
	    for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				dw16[(4*i)+j]=dw[j][i];//dw=f(t+h,w0+k3)
			}
	    }
	    for(i=0;i<16;i++)
	    {
	    	k4[i]=-h*dw16[i];//k4 values
		}
		for(j=0;j<16;j++)
		{
			neww[j]+=((k1[j]+(2.0*k2[j])+(2.0*k3[j])+k4[j])/6.0);
	      }
           	for(i=0;i<4;i++)//reshaping
		{
			for(j=0;j<4;j++)
			{
				w[j][i]=-neww[(4*i)+j];
				W[t][(4*i)+j]=neww[(4*i)+j];
			}
	    }
	}
}
void r_mat(int t,double R[2][2])//R profile
{
	double ds8=1.0,ds;
		double p1=5.536*pow(10,-8);
	double p2=-8.397*pow(10,-5);
	double p3=0.04809;
	double p4=-12.41;
	double p5=1242.2;
	//ds = p1*tf^4 + p2*tf^3 + p3*tf^2 + p4*tf + p5 ;
	ds=(p1*pow(tf,4))+(p2*pow(tf,3))+(p3*pow(tf,2))+(p4*pow(tf,1))+p5;
	R[0][0]=ds8*exp(ds*x[N-1-t][1]);
			R[1][1]=ds8*exp(ds*x[N-1-t][1]);
			R[0][1]=0.0;R[1][0]=0.0;
}
void retdfdx(double *dfdx,double ri,double ui,double vi)
{
	*(dfdx+0)=*(dfdx+1)=*(dfdx+2)=*(dfdx+5)=*(dfdx+7)=*(dfdx+9)=*(dfdx+13)=*(dfdx+15)=0.0;*(dfdx+3)=1.0;
	*(dfdx+4)=-ui/pow(ri,2);
	*(dfdx+6)=1.0/ri;
	*(dfdx+8)=(ui*vi)/pow(ri,2);
	*(dfdx+10)=-vi/ri;
	*(dfdx+11)=-ui/ri;
	*(dfdx+12)=((2*mu)/pow(ri,3))-(pow(ui,2)/pow(ri,2));
	*(dfdx+14)=(2*ui)/ri;
}
void matinv(double *a,double *inv)
{
	double det=((*(a+0))*(*(a+5))*(*(a+10))*(*(a+15))) + ((*(a+0))*(*(a+6))*(*(a+11))*(*(a+13))) + ((*(a+0))*(*(a+7))*(*(a+9))*(*(a+14)))
			 + ((*(a+1))*(*(a+4))*(*(a+11))*(*(a+14))) + ((*(a+1))*(*(a+6))*(*(a+8))*(*(a+15))) + ((*(a+1))*(*(a+7))*(*(a+10))*(*(a+12)))
		     + ((*(a+2))*(*(a+4))*(*(a+9))*(*(a+15))) + ((*(a+2))*(*(a+5))*(*(a+11))*(*(a+12))) + ((*(a+2))*(*(a+7))*(*(a+8))*(*(a+13)))
             + ((*(a+3))*(*(a+4))*(*(a+10))*(*(a+13))) + ((*(a+3))*(*(a+5))*(*(a+8))*(*(a+14))) + ((*(a+3))*(*(a+6))*(*(a+9))*(*(a+12)))
	         - ((*(a+0))*(*(a+5))*(*(a+11))*(*(a+14))) - ((*(a+0))*(*(a+6))*(*(a+9))*(*(a+15))) - ((*(a+0))*(*(a+7))*(*(a+10))*(*(a+13)))
	         - ((*(a+1))*(*(a+4))*(*(a+10))*(*(a+15))) - ((*(a+1))*(*(a+6))*(*(a+11))*(*(a+12))) - ((*(a+1))*(*(a+7))*(*(a+8))*(*(a+14)))
	         - ((*(a+2))*(*(a+4))*(*(a+11))*(*(a+13))) - ((*(a+2))*(*(a+5))*(*(a+8))*(*(a+15))) - ((*(a+2))*(*(a+7))*(*(a+9))*(*(a+12)))
	         - ((*(a+3))*(*(a+4))*(*(a+9))*(*(a+14))) - ((*(a+3))*(*(a+5))*(*(a+10))*(*(a+12))) - ((*(a+3))*(*(a+6))*(*(a+8))*(*(a+13)));

	double b11=((*(a+5))*(*(a+10))*(*(a+15))) + ((*(a+6))*(*(a+11))*(*(a+13))) + ((*(a+7))*(*(a+9))*(*(a+14))) -
	           ((*(a+5))*(*(a+11))*(*(a+14))) - ((*(a+6))*(*(a+9))*(*(a+15))) - ((*(a+7))*(*(a+10))*(*(a+13)));
	double b12=((*(a+1))*(*(a+11))*(*(a+14))) + ((*(a+2))*(*(a+9))*(*(a+15))) + ((*(a+3))*(*(a+10))*(*(a+13))) -
	           ((*(a+1))*(*(a+10))*(*(a+15))) - ((*(a+2))*(*(a+11))*(*(a+13))) - ((*(a+3))*(*(a+9))*(*(a+14)));
	double b13=((*(a+1))*(*(a+6))*(*(a+15))) + ((*(a+2))*(*(a+7))*(*(a+13))) + ((*(a+3))*(*(a+5))*(*(a+14))) -
	           ((*(a+1))*(*(a+7))*(*(a+14))) - ((*(a+2))*(*(a+5))*(*(a+15))) - ((*(a+3))*(*(a+6))*(*(a+13)));
	double b14=((*(a+1))*(*(a+7))*(*(a+10))) + ((*(a+2))*(*(a+5))*(*(a+11))) + ((*(a+3))*(*(a+6))*(*(a+9))) -
	           ((*(a+1))*(*(a+6))*(*(a+11))) - ((*(a+2))*(*(a+7))*(*(a+9))) - ((*(a+3))*(*(a+5))*(*(a+10)));
	double b21=((*(a+4))*(*(a+11))*(*(a+14))) + ((*(a+6))*(*(a+8))*(*(a+15))) + ((*(a+7))*(*(a+10))*(*(a+12))) -
	           ((*(a+4))*(*(a+10))*(*(a+15))) - ((*(a+6))*(*(a+11))*(*(a+12))) - ((*(a+7))*(*(a+8))*(*(a+14)));
	double b22=((*(a+0))*(*(a+10))*(*(a+15))) + ((*(a+2))*(*(a+11))*(*(a+12))) + ((*(a+3))*(*(a+8))*(*(a+14))) -
	           ((*(a+0))*(*(a+11))*(*(a+14))) - ((*(a+2))*(*(a+8))*(*(a+15))) - ((*(a+3))*(*(a+10))*(*(a+12)));
	double b23=((*(a+0))*(*(a+7))*(*(a+14))) + ((*(a+2))*(*(a+4))*(*(a+15))) + ((*(a+3))*(*(a+6))*(*(a+12))) -
	           ((*(a+0))*(*(a+6))*(*(a+15))) - ((*(a+2))*(*(a+7))*(*(a+12))) - ((*(a+3))*(*(a+4))*(*(a+14)));
	double b24=((*(a+0))*(*(a+6))*(*(a+11))) + ((*(a+2))*(*(a+7))*(*(a+8))) + ((*(a+3))*(*(a+4))*(*(a+10))) -
	           ((*(a+0))*(*(a+7))*(*(a+10))) - ((*(a+2))*(*(a+4))*(*(a+11))) - ((*(a+3))*(*(a+6))*(*(a+8)));
    double b31=((*(a+4))*(*(a+9))*(*(a+15))) + ((*(a+5))*(*(a+11))*(*(a+12))) + ((*(a+7))*(*(a+8))*(*(a+13))) -
	           ((*(a+4))*(*(a+11))*(*(a+13))) - ((*(a+5))*(*(a+8))*(*(a+15))) - ((*(a+7))*(*(a+9))*(*(a+12)));
	double b32=((*(a+0))*(*(a+11))*(*(a+13))) + ((*(a+1))*(*(a+8))*(*(a+15))) + ((*(a+3))*(*(a+9))*(*(a+12))) -
	           ((*(a+0))*(*(a+9))*(*(a+15))) - ((*(a+1))*(*(a+11))*(*(a+12))) - ((*(a+3))*(*(a+8))*(*(a+13)));
	double b33=((*(a+0))*(*(a+5))*(*(a+15))) + ((*(a+1))*(*(a+7))*(*(a+12))) + ((*(a+3))*(*(a+4))*(*(a+13))) -
	           ((*(a+0))*(*(a+7))*(*(a+13))) - ((*(a+1))*(*(a+4))*(*(a+15))) - ((*(a+3))*(*(a+5))*(*(a+12)));
	double b34=((*(a+0))*(*(a+7))*(*(a+9))) + ((*(a+1))*(*(a+4))*(*(a+11))) + ((*(a+3))*(*(a+5))*(*(a+8))) -
	           ((*(a+0))*(*(a+5))*(*(a+11))) - ((*(a+1))*(*(a+7))*(*(a+8))) - ((*(a+3))*(*(a+4))*(*(a+9)));
	double b41=((*(a+4))*(*(a+10))*(*(a+13))) + ((*(a+5))*(*(a+8))*(*(a+14))) + ((*(a+6))*(*(a+9))*(*(a+12))) -
	           ((*(a+4))*(*(a+9))*(*(a+14))) - ((*(a+5))*(*(a+10))*(*(a+12))) - ((*(a+6))*(*(a+8))*(*(a+13)));
	double b42=((*(a+0))*(*(a+9))*(*(a+14))) + ((*(a+1))*(*(a+10))*(*(a+12))) + ((*(a+2))*(*(a+8))*(*(a+13))) -
	           ((*(a+0))*(*(a+10))*(*(a+13))) - ((*(a+1))*(*(a+8))*(*(a+14))) - ((*(a+2))*(*(a+9))*(*(a+12)));
	double b43=((*(a+0))*(*(a+6))*(*(a+13))) + ((*(a+1))*(*(a+4))*(*(a+14))) + ((*(a+2))*(*(a+5))*(*(a+12))) -
	           ((*(a+0))*(*(a+5))*(*(a+14))) - ((*(a+1))*(*(a+6))*(*(a+12))) - ((*(a+2))*(*(a+4))*(*(a+13)));
	double b44=((*(a+0))*(*(a+5))*(*(a+10))) + ((*(a+1))*(*(a+6))*(*(a+8))) + ((*(a+2))*(*(a+4))*(*(a+9))) -
	           ((*(a+0))*(*(a+6))*(*(a+9))) - ((*(a+1))*(*(a+4))*(*(a+10))) - ((*(a+2))*(*(a+5))*(*(a+8)));

	*(inv+0)=b11/det;*(inv+1)=b12/det;*(inv+2)=b13/det;*(inv+3)=b14/det;
	*(inv+4)=b21/det;*(inv+5)=b22/det;*(inv+6)=b23/det;*(inv+7)=b24/det;
	*(inv+8)=b31/det;*(inv+9)=b32/det;*(inv+10)=b33/det;*(inv+11)=b34/det;
	*(inv+12)=b41/det;*(inv+13)=b42/det;*(inv+14)=b43/det;*(inv+15)=b44/det;
}
void ode4mass()
{
	int t;
	newm=mi;
	mass[0]=newm;
	double sum,k1,k2,k3,k4;
	for(t=1;t<N;t++)
	{
		k1=-(h*(newm)*sqrt((pow(ux[t],2))+(pow(uy[t],2))))/Isp;
		sum=newm+(k1*0.5);
		k2=-(h*(sum)*sqrt((pow(ux[t],2))+(pow(uy[t],2))))/Isp;
		sum=newm+(k2*0.5);
		k3=-(h*(sum)*sqrt((pow(ux[t],2))+(pow(uy[t],2))))/Isp;
		sum=newm+(k3);
		k4=-(h*(sum)*sqrt((pow(ux[t],2))+(pow(uy[t],2))))/Isp;
		newm+=((k1+(2.0*k2)+(2.0*k3)+k4)/6.0);
		mass[t]=newm;
	}
}

