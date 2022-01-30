double Tau=1.0;
double Itau=1.0;

/*

This file is used with all.cc and field.h. It extends the code in order to deal
with inertial particles.

*/


void RK4step_I(double t,double x,double y,double u, double v,double h,double *newstatus) {
	double xyp[2],field1[4],field2[4],field3[4],field4[4],thalf,tfull,halfh,xtmp,ytmp,utmp,vtmp,oldvalue[4];
	int i;
	//RK step 1
	halfh=h/2.0;
	oldvalue[0]=x;
	oldvalue[1]=y;
	oldvalue[2]=u;
	oldvalue[3]=v;
	
	//printf("t:%lf\n",t);
	pField->xyp(t,x,y,xyp);
	field1[0]=u;
	field1[1]=v;
	field1[2]=Itau*(xyp[0]-u);
	field1[3]=Itau*(xyp[1]-v);
	
	thalf=t+halfh;
	xtmp=x+field1[0]*halfh;
	ytmp=y+field1[1]*halfh;
	utmp=u+field1[2]*halfh;
	vtmp=v+field1[3]*halfh;
		

	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp);
	field2[0]=utmp;
	field2[1]=vtmp;
	field2[2]=Itau*(xyp[0]-utmp);
	field2[3]=Itau*(xyp[1]-vtmp);
	
	xtmp=x+field2[0]*halfh;
	ytmp=y+field2[1]*halfh;
	utmp=u+field2[2]*halfh;
	vtmp=v+field2[3]*halfh;
	
	
	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp);
	field3[0]=utmp;
	field3[1]=vtmp;
	field3[2]=Itau*(xyp[0]-utmp);
	field3[3]=Itau*(xyp[1]-vtmp);
	
	xtmp=x+field3[0]*h;
	ytmp=y+field3[1]*h;
	utmp=u+field3[2]*h;
	vtmp=v+field3[3]*h;
	
	
	tfull=t+h;
	
	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp);
	
	field4[0]=utmp;
	field4[1]=vtmp;
	field4[2]=Itau*(xyp[0]-utmp);
	field4[3]=Itau*(xyp[1]-vtmp);
	
	
	for (i=0;i<4;i++) {
		newstatus[i]=oldvalue[i]+h/6.0*(field1[i]+2.0*(field2[i]+field3[i])+field4[i]);
		//newstatus[i]=oldvalue[i]+h*field1[i];
		
		//printf("oldstatus %d:%lf newstatus %d:%lf\n",i,oldvalue[i],i,newstatus[i]);
		}
	}
	
	
int RK4_I(double *tspan,double *x0,int numx0,int Nstep,double *trj)
{
/*
Like RK4_I_tau, but computing the trajectories. The trajectories are
stored in trj as (x1(t1)...x1(tmax) y1(1)...y1(tmax) ... u1(t1)... u1(tn)... v1(t1)... v1(tn)... ... ... xn(t1)...xn(tmax)...).
NB:
tn=tspan(1)-h This is because Nstep is actually the number
of points over the trajectories; going from tspan(0) to tspan(1) in Nstep
would require Nstep+1 points!

*/
int j1,j2,indx,indy,indu,indv,pos0;
double h,t,xyp1[2],x,y,u,v,xtmp,ytmp,RKout[4];

h=(tspan[1]-tspan[0])/(double)(Nstep);

t=tspan[0];
//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0;j2++) {
	pos0=4*j2*Nstep;
	trj[pos0]=x0[4*j2];
	trj[pos0+Nstep]=x0[4*j2+1];
	trj[pos0+2*Nstep]=x0[4*j2+2];
	trj[pos0+3*Nstep]=x0[4*j2+3];
	}
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0;j2++) {
		indx=j1+4*j2*Nstep;
		indy=indx+Nstep;
		indu=indx+2*Nstep;
		indv=indx+3*Nstep;
		
		x=trj[indx-1];
		y=trj[indy-1];
		u=trj[indu-1];
		v=trj[indv-1];

		RK4step_I(t,x,y,u,v,h,RKout);
		trj[indx]=RKout[0];
		trj[indy]=RKout[1];		
		trj[indu]=RKout[2];
		trj[indv]=RKout[3];		
		
		}
	
	}
return 0;
}


int RK4_I_tau(double *tspan,double *x0,int numx0,int Nstep,double delta,int *elmax,double *tau)
{
/*
FSLEs for inertial particles.
This function integrates a number numx0 of 2d points which initial conditions are
contained in the vector x0 as: (x1 y1 u1 v1 x2 y2 u2 v2 x3 y3... un vn). The idea
is to put in x1 y1 a point, and in the other variables other points around
x1 y1. The funxtion computes the distance of all trajectories from x1 y1 and
stops when one has a distance greater than delta. It then provides the time
(in tau) and the number of the trajectory (in elmax). (I.e., "3" if the trajectory started)
at point x3 y3.) The integrator is a Runge-Kutta of 4th order, with fixed
time step (Nstep number of steps in the time window tspan[1], tspan[2]).
*/
int j1,j2,indx,indy;
double h,t,tfull,x,y,u,v,dist,dx,dy,RKout[4];
double *newstatus,*oldstatus,*swapstatus;

*tau=0.0;
*elmax=-1;

newstatus=new double[numx0*4];
oldstatus=new double[numx0*4];

h=(tspan[1]-tspan[0])/(double)Nstep;


t=tspan[0];
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0*4;j2++) oldstatus[j2]=x0[j2];
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	tfull=t+h;	

	for(j2=0;j2<numx0;j2++) {
		x=oldstatus[j2*4];
		y=oldstatus[j2*4+1];
		u=oldstatus[j2*4+2];
		v=oldstatus[j2*4+3];

		RK4step_I(t,x,y,u,v,h,RKout);

		newstatus[4*j2]=RKout[0];
		newstatus[4*j2+1]=RKout[1];
		newstatus[4*j2+2]=RKout[2];
		newstatus[4*j2+3]=RKout[3];
		
		
		if(j2>0)
			{
			dx=newstatus[4*j2]-newstatus[0];
			dy=newstatus[4*j2+1]-newstatus[1];
			
			dist=sqrt(dx*dx+dy*dy);
	
			if(dist>delta) {
				
				*elmax=j2;
				*tau=tfull-tspan[0]; //ADDED: -tspan[0]
				j2=numx0*4;
				j1=Nstep;
				}
		

			}
	
		
		//printf("x=%lf\t y=%lf\n",trj[j1],trj[j1+Nstep]);
		}

	swapstatus=oldstatus;
	oldstatus=newstatus;
	newstatus=swapstatus;
	}
delete oldstatus;
delete newstatus;
return 0;
}















