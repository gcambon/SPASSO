#include "mdefs.h"
#include "math.h"
//field.h 1.0.1
//improving lon/lat fit in altimetry_field.xyp (num_x->num_x-1, num_y->num_y-1)
//field.h 1.0.2          03/11/2004
//lon->lon-360 (lon2posx, imagen_field) and 0-330, not 0-331
//added laminar flow (useful for testing inertial particles)
//field.h 1.0.3          09/04/2005
//adding ERA40 wind datasets....... OK.
//field.h 1.0.4 	23/04/2005
//changing the integrator for working at the poles....... seems ok
//field.h 1.0.5 25/05/2005
//adding the POMME field... seems ok
//field.h 1.0.6 27/05/2005
//new dataset for the POMME field... seems ok
//field.h 1.0.7 27/05/2005
//adding finite-time with manifold directions...........
//adding the possibility of freezing the the velocity field at a specific time....... seems ok (22may2006: correction)
//field.h 1.0.8 12/04/2006
//adding AVISO Mediterranean..... seems ok
//field.h 1.0.8 27/04/2006
//changing the aviso field to cover a large time frame..... seems ok
//field.h 1.0.9 20/02/2007
//adding the Kerguelen plateau field.........
//field.h 1.0.10 6/11/2009
//adding noise to the velocity field



int scalecount(double *,int,double,int *,int *,int *);


// CONSTANTS and VARIABLES

const double pi=3.14159265358979;
int verb=1;
const double RTerra=6370e5; // in cm!
// GENERAL FUNCTIONS


void boxmuller(double *y1,double *y2) {
float x1, x2, w;
 
         do {
                 x1 = (2.0 * rand())/RAND_MAX - 1.0;
                 x2 = (2.0 * rand())/RAND_MAX - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         *y1 = x1 * w;
         *y2 = x2 * w;
	}


double lon2posx(double lon) {
	double a1=8;
	double a0=42;
	double posx;
	lon=lon;
	if(lon>180.0) lon-=360.0;
	posx=a1*lon+a0;
	return posx;
	}
	
double lat2posy(double lat) {
	double a4=0.000015445556;
	double a3=-0.001187709723;
	double a2=0.070791787735;
	double a1=6.526851079341;
	double a0=-242.193123554951;
	double lat2,posy;
	
	lat2=lat*lat;
	
	posy=a4*lat2*lat2+a3*lat2*lat+a2*lat2+a1*lat+a0;
	return posy;
	}

double posy2lat(double posy) {
	double a4=0.04035068698576;
	double a3=-0.14932066760361;
	double a2=-1.42895092126193;
	double a1=16.73888922849881;
	double a0=30.24146555579375;
	double lat,posy2;
	
	posy=posy/155;
	posy2=posy*posy;
	
	lat=a4*posy2*posy2+a3*posy2*posy+a2*posy2+a1*posy+a0;
	return lat;
	}

double mod(double a,double b) {
	int div;
	double r;
	div=int(a/b);
	r=a-(double)div*b;
	if(r<0.0) r+=b;
	return r;
	}

float intpl(float *x,double lambdax,double lambday,double lambdat) {
float ix1,ix2,ix;

ix1=(1-lambdax)*(1-lambday)*x[0]+lambdax*(1-lambday)*x[1]+lambdax*lambday*x[2]+(1-lambdax)*lambday*x[3];
ix2=(1-lambdax)*(1-lambday)*x[4]+lambdax*(1-lambday)*x[5]+lambdax*lambday*x[6]+(1-lambdax)*lambday*x[7];
ix=(1-lambdat)*ix1+lambdat*ix2;
return ix;
}

class license{
	public:
	license() {
		printf("Lamta 0.3 Copyright (C) 2009 11 6 Francesco d'Ovidio (francesco.dovidio@iscpif.fr).\n\nThis program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License\nalong with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n");
		}
	} License;
	
// THE FIELD FAMILY 

class field{
	public:
	field() {default_par();}
	int datatype; //1 is U and V in deg./sec, 2 for U and V in cm/sec
	int disttype; //1 is Euclidean, 2 is over the sphere.
	int gridtype; //0 is flat, 1 is sphere, regular spacing in deg., 2 is sphere, Mercator grid
	int noise; //0 no noise, not zero noisy velocity field
	double noiselevel;
	int freeze; //freeze<>0, to freeze the velocity field at time frozentime.
	double frozentime;
	void freezewarning() {
		if(freeze) printf("WARNING: time frozen at t=%lf. Use unfreezetime() to release.\n",frozentime);
		else printf("Time is not frozen.\n");
		}
	void noisewarning() {
		if(noise) printf("WARNING: Noisy velocity field (sigma=%lf). Use field_noise(0,0); to remove noise.\n",noiselevel);
		else printf("Time is not frozen.\n");
		}
		
	int xyp(double t,double x,double y, double *xp) {
	
		if(disttype==2) {
			while(x<-180) x+=360;
			while(x>180) x-=360;
			}
		
		if(freeze) xypt(frozentime,x,y,xp);
		else xypt(t,x,y,xp);
		
		return 0;
		}
	
	virtual int xypt(double t,double x,double y, double *xp) {
		xp[0]=0.0;
		xp[1]=0.0;
		return 0;
		}
	virtual int rawdata_out(long indx, double *U,double *V) {
		*U=0;
		*V=0;
		printf("No data in memory! (Field not initialised or field is a function.)");
		return 0;
		}
	
	virtual void print_par() {printf("Empty field is active.\n");};
	virtual void default_par() {datatype=1;//1 is U and V in deg./sec, 2 for U and V in cm/sec
		noise=0;
		freeze=0;
		datatype=1;
		disttype=1;
		gridtype=1;
		
		};
	virtual void set_par(double *) {};
	} Genericfield, *pField=&Genericfield;

class meanderingjet: public field {
	public:
	meanderingjet() {default_par();}
	double B0,c,k,epsi,phi,omega;
	void default_par() {
		printf("Meandering jet created.\n");
		B0=1.2;c=.12;k=2*pi/7.5;epsi=.3;phi=pi/2;omega=.1;
		}
	void print_par() {
		printf("Meandering jet is active:\nB0=%f c=%f k=%f epsi=%f phi=%f omega=%f\n",B0,c,k,epsi,phi,omega);
		}
	int xypt(double t,double x,double y,double *xp)
		{
		double B,BK,s,S,f1,f2,df1x,df2x,invsqrtf2,dfx,expr0,expr1,dpsix,dpsiy,BK2,s2;

		if (xp==NULL) return 1;

		B=B0+epsi*cos(omega*t+phi);
		BK=B*k;
		BK2=BK*BK;
		s=sin(k*x);
		s2=s*s;
		S=cos(k*x);
		f1=y-B*S;
		f2=1+BK2*s2;
		df1x=BK*s;
		df2x=2*k*BK2*s*S;
		invsqrtf2=1/sqrt(f2);
		dfx=-.5*f1/pow(f2,(double)3/2)*df2x+df1x*invsqrtf2;
		expr0=tanh(f1*invsqrtf2);
		expr1=(-1+expr0*expr0);

		dpsix=expr1*dfx;
		dpsiy=expr1*invsqrtf2+c;

		xp[0]=dpsiy;
		xp[1]=-dpsix;

		return 0;
		}
	void set_par(double *pl) {
		B0=pl[0];c=pl[1];k=pl[2];epsi=pl[3];phi=pl[4];omega=pl[5];
		}
	} Meanderingjet;

class homocline: public field {
	public:
	homocline() {default_par();}
	double a,b,c,omega1,omega2;
	void default_par() {
		printf("Homocline created.\n");
		a=1.0;b=.5;c=2.0,omega1=omega2=2*pi/7.5;
		datatype=1;
		disttype=1;
		
		}
	void print_par() {
		printf("Homocline is active:\na=%f b=%f c=%f omega1=%f omega2=%f\n",a,b,c,omega1,omega2);
		}
	int xypt(double t,double x,double y,double *xp)
		{
		double u,v,expn,b1,a1,x1;

		if (xp==NULL) return 1;
		
		b1=b*(2+sin(omega2*t));
		//a1=a*(2+cos(omega2*t));
		a1=a;
		
		expn=exp(-x*x-b1*y*y);
		
		u=-2*a1*b1*y*expn;
		v=2*a1*x*expn-c*(2+sin(omega1*t));
		
		
		
		
		
		xp[0]=u;
		xp[1]=v;

		return 0;
		}
	void set_par(double *pl) {
		a=pl[0];b=pl[1];c=pl[2];omega1=pl[3];omega2=pl[4];
		}
	} Homocline;


class periodicflow: public field {
	//E. Shuckburgh and P. Haynes, Phys. Fluids 15, 3342-4457 (2003). 
	//The field is on a sphere, x and y are lon and lat and should be in radiant.
	public:
	periodicflow() {default_par();}
	double a,b,epsilon1,epsilon2,omega;
	void default_par() {
		printf("Periodicflow (on the sphere) created.\n");
		a=0.25;b=0.05;epsilon1=0.25;epsilon2=0.0125;omega=1.0;disttype=2;
		}
	void print_par() {
		printf("Periodic flow is active:\na=%lf b=%lf epsilon1=%lf epsilon2=%lf omega=%lf\n",a,b,epsilon1,epsilon2,omega);
		}
	int xypt(double t,double x,double y,double *xp)
		{
		double u,v,U,V,mu,dpsidmu,mu2,lambda,phi;
		double cv=pi/180;
		lambda=x*cv;
		phi=y*cv;
		mu=sin(phi);
		mu2=mu*mu;
		
		dpsidmu=b-3*a*mu2-3*mu*sqrt(1-mu2)*(epsilon1*cos(lambda)+epsilon2*cos(lambda+omega*t));
		U=-dpsidmu;
		
		V=-sqrt(1-mu2)*(epsilon1*sin(lambda)+epsilon2*sin(lambda+omega*t));

		if (xp==NULL) return 1;

		u=U*cos(phi);
		v=V*cos(phi);

		xp[0]=u/cv;
		xp[1]=v/cv;

		return 0;
		}
	void set_par(double *pl) {
		a=pl[0];b=pl[1];epsilon1=pl[2];epsilon2=pl[3];omega=pl[4];
		}
	} Periodicflow;




class linearsaddle: public field {
	public:
	linearsaddle() {default_par();}
	double lambda1,lambda2;
	void default_par() {
		printf("Saddle created.\n");
		lambda1=.3;lambda2=-.2;
		}
	void print_par() {
		printf("Linear saddle is active:\n lambda1=%f lambda2=%f\n",lambda1,lambda2);
		}
	int xypt(double t,double x,double y,double *xp) {
		xp[0]=lambda1*x;
		xp[1]=lambda2*y;
		return 0;
		}
	void set_par(double *pl) {
		lambda1=pl[0];lambda2=pl[1];
		}
	} Linearsaddle;

class laminarflow: public field {
	public:
	laminarflow() {default_par();}
	double u,v;
	void default_par() {
		printf("Laminar flow.\n");
		u=.3;v=-.2;
		}
	void print_par() {
		printf("Laminar flow is active:\n u=%f v=%f\n",u,v);
		}
	int xypt(double t,double x,double y,double *xp) {
		xp[0]=u;
		xp[1]=v;
		return 0;
		}
	void set_par(double *pl) {
		u=pl[0];v=pl[1];
		}
	
	} Laminarflow;

class vonkarman2: public field {
	public:
	vonkarman2() {default_par();}
	double a,w,R0,Tc,alpha,y0,L,u0;
	void default_par() {
		printf("Von Karman2 field created.\n");
		a=1.0;w=35.06;R0=.35;Tc=1.0;alpha=2.0;y0=0.3;L=2.0;u0=14.0;
		}
	void print_par() {
		printf("Von Karman2 is active:\n a=%f w=%f R0=%f Tc=%f alpha=%f y0=%f L=%f u0=%f\n",a,w,R0,Tc,alpha,y0,L,u0);
		}
	int xypt(double t,double x,double y,double *xp) {
		double x1,x2,dpsidx,dpsidy;
		
		x1=1 + L*mod(t/Tc,1);
		x2=1 + L*mod((t - Tc/2.)/Tc,1);
		
		dpsidx=(2*a*x*(-1 + Sqrt(Power(x,2) + Power(y,2)))*
      ((1 - Power(E,-(Power(-1 + x,2)/Power(alpha,2)) - 
             Power(y,2)))*u0*y - 
        (w*Abs(Sin((Pi*t)/Tc)))/
         Power(E,R0*(Power(x - x1,2) + 
             Power(alpha,2)*Power(y - y0,2))) + 
        (w*Abs(Sin((Pi*(t - Tc/2.))/Tc)))/
         Power(E,R0*(Power(x - x2,2) + 
             Power(alpha,2)*Power(y + y0,2)))))/
    (Power(E,a*Power(-1 + Sqrt(Power(x,2) + Power(y,2)),
         2))*Sqrt(Power(x,2) + Power(y,2))) + 
   (1 - Power(E,-(a*Power(-1 + 
            Sqrt(Power(x,2) + Power(y,2)),2))))*
    ((2*Power(E,-(Power(-1 + x,2)/Power(alpha,2)) - 
           Power(y,2))*u0*(-1 + x)*y)/Power(alpha,2) + 
      (2*R0*w*(x - x1)*Abs(Sin((Pi*t)/Tc)))/
       Power(E,R0*(Power(x - x1,2) + 
           Power(alpha,2)*Power(y - y0,2))) - 
      (2*R0*w*(x - x2)*Abs(Sin((Pi*(t - Tc/2.))/Tc)))/
       Power(E,R0*(Power(x - x2,2) + 
           Power(alpha,2)*Power(y + y0,2))));
		
		
		dpsidy=(2*a*y*(-1 + Sqrt(Power(x,2) + Power(y,2)))*
      ((1 - Power(E,-(Power(-1 + x,2)/Power(alpha,2)) - 
             Power(y,2)))*u0*y + 
        (w*Abs(Cos((Pi*t)/Tc)))/
         Power(E,R0*(Power(x - x2,2) + 
             Power(alpha,2)*Power(y + y0,2))) - 
        (w*Abs(Sin((Pi*t)/Tc)))/
         Power(E,R0*(Power(x - x1,2) + 
             Power(alpha,2)*Power(y - y0,2)))))/
    (Power(E,a*Power(-1 + Sqrt(Power(x,2) + Power(y,2)),
         2))*Sqrt(Power(x,2) + Power(y,2))) + 
   (1 - Power(E,-(a*Power(-1 + 
            Sqrt(Power(x,2) + Power(y,2)),2))))*
    (u0 - Power(E,-(Power(-1 + x,2)/Power(alpha,2)) - 
         Power(y,2))*u0 + 
      2*Power(E,-(Power(-1 + x,2)/Power(alpha,2)) - 
         Power(y,2))*u0*Power(y,2) - 
      (2*Power(alpha,2)*R0*w*(y + y0)*
         Abs(Cos((Pi*t)/Tc)))/
       Power(E,R0*(Power(x - x2,2) + 
           Power(alpha,2)*Power(y + y0,2))) + 
      (2*Power(alpha,2)*R0*w*(y - y0)*
         Abs(Sin((Pi*t)/Tc)))/
       Power(E,R0*(Power(x - x1,2) + 
           Power(alpha,2)*Power(y - y0,2))));
		
		
		xp[0]=dpsidy;
		xp[1]=-dpsidx;
		
		return 0;
		}
	} Vonkarman2;

class vonkarman: public field {
	public:
	vonkarman() {default_par();}
	double a,w,R0,Tc,alpha,y0,L,u0;
	void default_par() {
		printf("Von Karman field created.\n");
		a=1.0;w=35.06;R0=.35;Tc=1.0;alpha=2.0;y0=0.3;L=2.0;u0=14.0;
		}
	void print_par() {
		printf("Von Karman is active:\n a=%f w=%f R0=%f Tc=%f alpha=%f y0=%f L=%f u0=%f\n",a,w,R0,Tc,alpha,y0,L,u0);
		}
	int xypt(double t,double x,double y,double *xp) {
		
		double xv1,xv2,yv1,yv2,f,g,h,fx,gx,fy,gy,h1,h2,d11,d1x,d1y,d2,d2x,d2y,d3,d3x,d3y,d4,dpsidx,dpsidy;
		double s,sx,sy,g1x,g1y,g2x,g2y,g1,g2,x2,y2,alpha2,rho;
		
		x2=x*x;
		y2=y*y;
		rho=sqrt(x2+y2);
		xv1=1.0+L*mod(t/Tc,1.0);
		xv2=1.0+L*mod((t - Tc/2.0)/Tc,1.0);
		
		yv1=y0;
		yv2=-y0;
		alpha2=alpha*alpha;
		
		h1=fabs(sin(Pi*t/Tc));
		h2=fabs(sin(Pi*(t-Tc/2)/Tc));
		
		d11=(-2*a*(rho-1))/rho;
		d1x=x*d11;
		d1y=y*d11;
		f=1-exp(-a*(rho-1)*(rho-1));
		fx=-exp(-a*(rho-1)*(rho-1))*d1x;
		fy=-exp(-a*(rho-1)*(rho-1))*d1y;
		
		d4=exp(-(x2-2*x+1)/alpha2-y2);
		s=1-d4;
		sx=d4*(1/alpha2*(2*x-2));
		sy=d4*(2*y);
		
		d2x=-2*R0*(x-xv1);
		d2y=-2*R0*alpha2*(y-yv1);
		d3x=-2*R0*(x-xv2);
		d3y=-2*R0*alpha2*(y-yv2);
		
		
		g1=exp(-R0*((x-xv1)*(x-xv1)+alpha2*(y-yv1)*(y-yv1)));
		g2=exp(-R0*((x-xv2)*(x-xv2)+alpha2*(y-yv2)*(y-yv2)));
		g=-w*h1*g1+w*h2*g2+u0*y*s;
		g1x=g1*d2x;
		g1y=g1*d2y;
		g2x=g2*d3x;
		g2y=g2*d3y;
		gx=-w*h1*g1x+w*h2*g2x+u0*y*sx;
		gy=-w*h1*g1y+w*h2*g2y+u0*s+u0*y*sy;
		
		dpsidx=fx*g+gx*f;
		dpsidy=fy*g+gy*f;
		
		if (rho<1.0) {
		dpsidx=0.0;
		dpsidy=0.0;
		}
		
		xp[0]=dpsidy;
		xp[1]=-dpsidx;
		
		return 0;
		}
	} Vonkarman;
class vonkarman_m: public field {
	public:
	vonkarman_m() {default_par();}
	double a,w,R0,Tc,alpha,y0,L,u0,Dc;
	void default_par() {
		printf("Von Karman_m field created.\n");
		a=1.78e-14;w=61e6;R0=6.22e-15;Tc=7e6;alpha=2.0;y0=22.5e5;L=75e6;u0=10.0;Dc=75e5;
		//a=1.0;w=35.06;R0=.35;Tc=1.0;alpha=2.0;y0=0.3;L=2.0;u0=14.0;Dc=1.0;
		}
	void print_par() {
		printf("Von Karman_m (Mediterranean) is active:\n a=%lf w=%lf R0=%lf Tc=%lf alpha=%lf y0=%lf L=%lf u0=%lf Dc=%lf\n",a,w,R0,Tc,alpha,y0,L,u0,Dc);
		}
	int xypt(double t,double x,double y,double *xp) {
		
		double xv1,xv2,yv1,yv2,f,g,h,fx,gx,fy,gy,h1,h2,d11,d1x,d1y,d2,d2x,d2y,d3,d3x,d3y,d4,dpsidx,dpsidy;
		double Dc2,s,sx,sy,g1x,g1y,g2x,g2y,g1,g2,x2,y2,alpha2,rho;
		
		Dc2=Dc*Dc;
		x2=x*x;
		y2=y*y;
		rho=sqrt(x2+y2);
		xv1=Dc+L*mod(t/Tc,1.0);
		xv2=Dc+L*mod((t - Tc/2.0)/Tc,1.0);
		
		yv1=y0;
		yv2=-y0;
		alpha2=alpha*alpha;
		
		
		h1=fabs(sin(Pi*t/Tc));
		h2=fabs(sin(Pi*(t-Tc/2)/Tc));
		
		
		d11=(-2*a*(rho-Dc))/rho;
		d1x=x*d11;
		d1y=y*d11;
		
		f=1-exp(-a*(rho-Dc)*(rho-Dc));
		fx=-exp(-a*(rho-Dc)*(rho-Dc))*d1x;
		fy=-exp(-a*(rho-Dc)*(rho-Dc))*d1y;
		
		d4=exp(-1/Dc2*((x2-2*Dc*x+Dc*Dc)/alpha2+y2));
		s=1-d4;
		sx=d4*(1/(alpha2*Dc2)*(2*x-2*Dc));
		sy=1/Dc2*d4*(2*y);
		
				
		d2x=-2*R0*(x-xv1);
		d2y=-2*R0*alpha2*(y-yv1);
		d3x=-2*R0*(x-xv2);
		d3y=-2*R0*alpha2*(y-yv2);
		
		
		g1=exp(-R0*((x-xv1)*(x-xv1)+alpha2*(y-yv1)*(y-yv1)));
		g2=exp(-R0*((x-xv2)*(x-xv2)+alpha2*(y-yv2)*(y-yv2)));
		g=-w*h1*g1+w*h2*g2+u0*y*s;
		g1x=g1*d2x;
		g1y=g1*d2y;
		g2x=g2*d3x;
		g2y=g2*d3y;
		gx=-w*h1*g1x+w*h2*g2x+u0*y*sx;
		gy=-w*h1*g1y+w*h2*g2y+u0*s+u0*y*sy;
		
		
		dpsidx=fx*g+gx*f;
		dpsidy=fy*g+gy*f;
		
		if (rho<Dc) {
		dpsidx=0.0;
		dpsidy=0.0;
		}
		
		xp[0]=dpsidy;
		xp[1]=-dpsidx;
		
		
		return 0;
		}
	} Vonkarman_m;

class lut_field: public field {
	public:
	lut_field() {default_par();}
	double xi,xf,yi,yf,ti,tf,latLUT[8192],invdlat;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated,numlatLUT;
	void default_par() {
		printf("Lut created.\n");
		xi=xf=yi=yf=ti=tf=0.0;datax=NULL;datay=NULL;memoryisallocated=0;
		numx=numy=numt=0;
		datatype=1;
		disttype=2;
		numlatLUT=8192;
		invdlat=1;
		}
	void print_par() {
		printf("LUT field is active:\n xi=%f xf=%f numx=%d yi=%f yf=%f numy=%d ti=%f tf=%f numt=%d memoryisallocated=%d\n",xi,xf,numx,yi,yf,numy,ti,tf,numt,memoryisallocated);
		freezewarning();
		noisewarning();
		}
	void set_par(double *par) {
		int j;
		double xL,c,ddeg;
		xi=par[0];
		xf=par[1];
		numx=(int)par[2];
		yi=par[3];
		yf=par[4];
		numy=(int)par[5];
		ti=par[6];
		tf=par[7];
		numt=(int)par[8];
		
		if(memoryisallocated) {
			delete datax;
			delete datay;
			printf("U and V freed in LUT field.\n");
			memoryisallocated=0;
			}
		
		ddeg=(xf-xi)/(numx-1);
		c=1/(ddeg/180.0*pi);
		for(j=0;j<numlatLUT;j++) {
			xL=((double)j)/numlatLUT*(yf-yi)+yi;
			xL=xL/180.0*pi;
			latLUT[j]=c*(log(cos(xL/2)+sin(xL/2))-log(cos(xL/2)-sin(xL/2)));
			//printf("j:%d\txL:%lf\t%lf\n",j,xL,latLUT[j]);
			}
		invdlat=1.0/(yf-yi)*(numlatLUT-1);
		}
	
	double clatLUT(double lat) {
		int latp1;
		latp1=(int)((lat-yi)*invdlat);
		//printf("latp1:%d\n",latp1);
		if(latp1<0) latp1=0;
		if(latp1>numlatLUT) latp1=numlatLUT;
		return latLUT[latp1];
		}
	
	
	double mercatorpos(double lat) {
		double Fphi0,Fphii,c,ddeg,ilat;
		Fphi0=clatLUT(yi);
		Fphii=clatLUT(lat);
		ilat=Fphii-Fphi0;
		//printf("Fphi0:%lf Fphii:%lf\n",Fphi0,Fphii);
		return ilat;
		}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		posx=(x-xi)/(xf-xi)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		
		if(gridtype==2) {
			posy=mercatorpos(y);
			//printf("mercatorpos:%lf\n",posy);
			}
		else posy=(y-yi)/(yf-yi)*(numy-1);
		
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=(t-ti)/(tf-ti)*(numt-1);
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			//printf("No data!\n");
			return 1;
			}
			
	/*
		if((ix<0)||(ix>=numx-1)||(iy<0)||(iy>=numy-1)||(it<0)||(it>=numt-1)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("Out of limits of stored data!\n");
			return 2;
			}
	*/	
		
		if(ix<0) ix=0;
		if(iy<0) iy=0;
		if(it<0) it=0;
		if(ix>numx-1) ix=numx-1;
		if(iy>numy-1) iy=numy-1;
		if(it>numt-1) it=numt-1;	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) ixl=numx-1;
		if(iyl>numy-1) iyl=numy-1;
		if(itl>numt-1) itl=numt-1;	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Lut_field;

class imagen_field: public field {
	public:
	imagen_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double Day_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated;
	void default_par() {
		lat_i=30.2415;
		lat_f=45.4424;
		lon_i=-5.25;
		lon_f=36.0;
		Day_in_secs=60*60*24;
		numx=331;numy=156;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		printf("IMAGEN field is created.\n");
		}
	void print_par() {
		printf("IMAGEN field is active:\nnumt=%d memoryisallocated=%d\n",numt,memoryisallocated);
		}
	void set_par(double *par) {}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		posx=lon2posx(x);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=lat2posy(y);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/Day_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) ix=0;
		if(iy<0) iy=0;
		if(it<0) it=0;
		if(ix>numx-1) ix=numx-1;
		if(iy>numy-1) iy=numy-1;
		if(it>numt-1) it=numt-1;	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) ixl=numx-1;
		if(iyl>numy-1) iyl=numy-1;
		if(itl>numt-1) itl=numt-1;	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Imagen_field;
	
class altimetry_field: public field {
	public:
	altimetry_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double TenDays_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated;
	void default_par() {
		lat_i=30;
		lat_f=45.8;
		lon_i=-5;
		lon_f=35.8;
		TenDays_in_secs=60*60*24*10; //10 days!
		numx=205;numy=80;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		printf("Altimetry field is (re?)initialised.\n");
		}
	void print_par() {
		printf("Altimetry field is active:\nnumt=%d memoryisallocated=%d\n",numt,memoryisallocated);
		}
	void set_par(double *par) {}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/TenDays_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;
		//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
		}
		if(iy<0) {iy=0;printf("Warning: lat. out of range (%d), using value at the grid border.\n",iy);}
		if(it<0) {it=0;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}
		if(ix>numx-1) {ix=numx-1;printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);}
		if(iy>numy-1) {iy=numy-1;printf("Warning: lat. out of range (%d, using value at the grid border).\n",iy);}
		if(it>numt-1) {it=numt-1;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) {ixl=numx-1;printf("Warning: long. (ixl) out of range (%d).\n",ixl);}
		if(iyl>numy-1) {iyl=numy-1;printf("Warning: lat. (iyl) out of range (%d).\n",iyl);}
		if(itl>numt-1) {itl=numt-1;printf("Warning: time (itl) out of range (%d).\n",itl);}	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Altimetry_field;

class aviso_Med_field: public field {
	public:
	aviso_Med_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double SevenDays_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated,framesskipped;
	void default_par() {
		lat_i=30;
		lat_f=46;
		lon_i=-5;
		lon_f=36.875;
		SevenDays_in_secs=60*60*24*7; //7 days!
		numx=336;numy=129;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		framesskipped=0;
		printf("Aviso Med. field is (re?)initialised.\n");
		}
	void print_par() {
		printf("Aviso Med. field is active:\nnumt=%d memoryisallocated=%d framesskipped=%d\n",numt,memoryisallocated,framesskipped);
		printf("Time t=0 is second 0 of 1 January 2003. Field code=61.\n");
		
		}
	void set_par(double *par) {}
	
	int allocate_mem(unsigned long sz) {
		int retcode=-2;
		delete_mem();
		printf("Trying to reserve two %ld float arrays...\n",sz);

		datax=new float[sz];
		datay=new float[sz];

		if ((datax==NULL)||(datay==NULL)) {
			printf("Not enough memory! Please try with less days.\n");
			retcode=-1;
	 		}
			
		else {
			memoryisallocated=1;
			retcode=0;
			}
		return retcode;
		}
		
	int delete_mem() {
		int retcode=-2;
		if(memoryisallocated) {
			printf("Deleting previously allocated memory...\n");
			delete datax;
			delete datay;
			memoryisallocated=0;
			retcode=0;
			}
		else {
			printf("Data memory empty.\n");retcode=1;}
		return retcode;
		}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/SevenDays_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;
		//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
		}
		if(iy<0) {iy=0;printf("Warning: lat. out of range (%d), using value at the grid border.\n",iy);}
		if(it<0) {it=0;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}
		if(ix>numx-1) {ix=numx-1;printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);}
		if(iy>numy-1) {iy=numy-1;printf("Warning: lat. out of range (%d, using value at the grid border).\n",iy);}
		if(it>numt-1) {it=numt-1;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) {ixl=numx-1;printf("Warning: long. (ixl) out of range (%d).\n",ixl);}
		if(iyl>numy-1) {iyl=numy-1;printf("Warning: lat. (iyl) out of range (%d).\n",iyl);}
		if(itl>numt-1) {itl=numt-1;printf("Warning: time (itl) out of range (%d).\n",itl);}	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Aviso_Med_field;

class kerguelen_field: public field {
	public:
	kerguelen_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double OneDay_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated,framesskipped;
	void default_par() {
		lat_i=-59.7500;
		lat_f=-35.2500;
		lon_i=57.25;
		lon_f=109.75;
		OneDay_in_secs=60*60*24; //1 day!
		numx=211;numy=99;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		framesskipped=0;
		datatype=1;
		disttype=2;
		printf("Kerguelen field is (re?)initialised.\n");
		}
	void print_par() {
		printf("Kerguelen field is active:\nnumt=%d memoryisallocated=%d framesskipped=%d\n",numt,memoryisallocated,framesskipped);
		printf("Time t=0 is second 0 of 1 July 2003. Field code=1001.\n");
		
		}
	void set_par(double *par) {}
	
	int allocate_mem(unsigned long sz) {
		int retcode=-2;
		delete_mem();
		printf("Trying to reserve two %ld float arrays...\n",sz);

		datax=new float[sz];
		datay=new float[sz];

		if ((datax==NULL)||(datay==NULL)) {
			printf("Not enough memory! Please try with less days.\n");
			retcode=-1;
	 		}
			
		else {
			memoryisallocated=1;
			retcode=0;
			}
		return retcode;
		}
		
	int delete_mem() {
		int retcode=-2;
		if(memoryisallocated) {
			printf("Deleting previously allocated memory...\n");
			delete datax;
			delete datay;
			memoryisallocated=0;
			retcode=0;
			}
		else {
			printf("Data memory empty.\n");retcode=1;}
		return retcode;
		}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/OneDay_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;
			//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
			}
		if(iy<0) {iy=0;
			//printf("Warning: lat. out of range (%d), using value at the grid border.\n",iy);
			}
		if(it<0) {it=0;
			//printf("Warning: time out of range (%d), using value at the grid border.\n",it);
			}
		if(ix>numx-1) {ix=numx-1;
			//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
			}
		if(iy>numy-1) {iy=numy-1;
			//printf("Warning: lat. out of range (%d, using value at the grid border).\n",iy);
			}
		if(it>numt-1) {it=numt-1;
			//printf("Warning: time out of range (%d), using value at the grid border.\n",it);
			}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) {ixl=numx-1;printf("Warning: long. (ixl) out of range (%d).\n",ixl);}
		if(iyl>numy-1) {iyl=numy-1;printf("Warning: lat. (iyl) out of range (%d).\n",iyl);}
		if(itl>numt-1) {itl=numt-1;printf("Warning: time (itl) out of range (%d).\n",itl);}	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Kerguelen_field;

class windERA40_field: public field {
	public:
	windERA40_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double SixHours_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated;
	void default_par() {
		lat_i=90;
		lat_f=-90;
		lon_i=-180;
		lon_f=179;
		SixHours_in_secs=6*60*60; //data every six hours!
		numx=360;numy=181;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		datatype=2;
		disttype=2;
		printf("WindERA40 field is (re?)initialised.\n");
		}
	void print_par() {
		printf("WindERA40 field is active:\nnumt=%d memoryisallocated=%d datatype=%d disttype=%d\n",numt,memoryisallocated,datatype,disttype);
		}
	void set_par(double *par) {datatype=(int)par[0];disttype=(int)par[1];}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		//printf("x:%lf y:%lf\n",x,y);
		while(y>90) y=180-y;
		while(y<-90) y=-180-y;
		while(x>179) x=x-360;
		while(x<-180) x=x+360;
		
		 
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/SixHours_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;printf("Error (windERA40): long. out of range (%d, x=%lf)! This should never happen!\n",ix,x);}
		if(ix>numx-1) {ix=numx-1;printf("Error (windERA40): long. out of range (%d, x=%lf)! This should never happen!\n",ix,x);}
		if(iy<0) {iy=0;printf("Error (windERA40): lat. out of range (%d, y=%lf)! This should never happen!\n",iy,y);}
		if(iy>numy-1) {iy=numy-1;printf("Error (windERA40): lat. out of range (%d, y=%lf)! This should never happen!\n",iy,y);}
		
		if(it<0) {it=0;printf("Warning (windERA40): time out of range (%d, t=%lf), using value at the grid border.\n",it,t);}
		if(it>numt-1) {it=numt-1;printf("Warning (windERA40): time out of range (%d, t=%lf), using value at the grid border.\n",it,t);}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) ixl=0;
		if(iyl>numy-1) iyl=0;
		if(itl>numt-1) itl=0;
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		//printf("xp1:%lf yp1:%lf\n",xp1,yp1);

		return 0;
		}
	int rawdata_out(long indx, double *U,double *V) {
		if(memoryisallocated) {
			*U=datax[indx];
			*V=datay[indx];
			}
		else {
			*U=0;
			*V=0;
			printf("WindERA40_field not initialised.\n");
			}
		
		return 0;
		} 
	} WindERA40_field;

class pomme_field: public field {
	public:
	pomme_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double SevenDays_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated;
	void default_par() {
		lat_i=-9.5465;
		lat_f=59.6440;
		lon_i=-79.5;
		lon_f=-.5;
		SevenDays_in_secs=60*60*24*7; //7 days!
		numx=238;numy=420;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		datatype=1;
		disttype=1;
		printf("Pomme field is (re?)initialised.\n");
		}
	void print_par() {
		printf("Pomme field is active:\nnumt=%d memoryisallocated=%d\n",numt,memoryisallocated);
		}
	void set_par(double *par) {}
	
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		while(y>90) y=180-y;
		while(y<-90) y=-180-y;
		while(x>179) x=x-360;
		while(x<-180) x=x+360;
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/SevenDays_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;
		//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
		}
		if(iy<0) {iy=0;printf("Warning: lat. out of range (%d), using value at the grid border.\n",iy);}
		if(it<0) {it=0;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}
		if(ix>numx-1) {ix=numx-1;printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);}
		if(iy>numy-1) {iy=numy-1;printf("Warning: lat. out of range (%d, using value at the grid border).\n",iy);}
		if(it>numt-1) {it=numt-1;printf("Warning: time out of range (%d), using value at the grid border.\n",it);}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) {ixl=numx-1;printf("Warning: long. (ixl) out of range (%d).\n",ixl);}
		if(iyl>numy-1) {iyl=numy-1;printf("Warning: lat. (iyl) out of range (%d).\n",iyl);}
		if(itl>numt-1) {itl=numt-1;printf("Warning: time (itl) out of range (%d).\n",itl);}	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Pomme_field;


class pco2_field: public field {
	public:
	pco2_field() {default_par();}
	double lat_i,lat_f,lon_i,lon_f;
	double OneDay_in_secs;
	float *datax,*datay;
	int numx,numy,numt,memoryisallocated;
	void default_par() {
		lat_i=37.3;
		lat_f=45.8;
		lon_i=-22;
		lon_f=-14.5;
		OneDay_in_secs=60*60*24; //1 day!
		numx=151;numy=171;numt=0;
		datax=NULL;datay=NULL;memoryisallocated=0;
		datatype=1;
		disttype=1;
		printf("Pco2 field is (re?)initialised.\n");
		}
	void print_par() {
		printf("Pco2 field is active:\nnumt=%d memoryisallocated=%d\n",numt,memoryisallocated);
		}
	void set_par(double *par) {}
	
	int allocate_mem(unsigned long sz) {
		int retcode=-2;
		delete_mem();
		printf("Trying to reserve two %ld float arrays...\n",sz);

		datax=new float[sz];
		datay=new float[sz];

		if ((datax==NULL)||(datay==NULL)) {
			printf("Not enough memory! Please try with less days.\n");
			retcode=-1;
	 		}
			
		else {
			memoryisallocated=1;
			retcode=0;
			}
		return retcode;
		}
		
	int delete_mem() {
		int retcode=-2;
		if(memoryisallocated) {
			printf("Deleting previously allocated memory...\n");
			delete datax;
			delete datay;
			memoryisallocated=0;
			retcode=0;
			}
		else {
			printf("Data memory empty.\n");retcode=1;}
		return retcode;
		}
	int xypt(double t,double x,double y,double *xyp) {
		double posx,posy,post,lambdax,lambday,lambdat;
		float xpi[8],ypi[8],xp1,yp1;
		int ix,iy,it,ixl,iyl,itl;
		
		while(y>90) y=180-y;
		while(y<-90) y=-180-y;
		while(x>179) x=x-360;
		while(x<-180) x=x+360;
		
		posx=(x-lon_i)/(lon_f-lon_i)*(numx-1);
		ix=(int)posx;
		lambdax=posx-(double)ix;
			
		posy=(y-lat_i)/(lat_f-lat_i)*(numy-1);
		iy=(int)posy;
		lambday=posy-(double)iy;
			
		post=t/OneDay_in_secs;
		it=(int)post;
		lambdat=post-(double)it;
		
		
		if ((datax==NULL)||(datay==NULL)) {
			xyp[0]=0.0;
			xyp[1]=0.0;
			printf("No data!\n");
			return 1;
			}
			
	
		
		if(ix<0) {ix=0;
		//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
		}
		if(iy<0) {iy=0;
		//printf("Warning: lat. out of range (%d), using value at the grid border.\n",iy);
		}
		if(it<0) {it=0;
		//printf("Warning: time out of range (%d), using value at the grid border.\n",it);
		}
		if(ix>numx-1) {ix=numx-1;
		//printf("Warning: long. out of range (%d), using value at the grid border.\n",ix);
		}
		if(iy>numy-1) {iy=numy-1;printf("Warning: lat. out of range (%d, using value at the grid border).\n",iy);
		}
		if(it>numt-1) {it=numt-1;
		//printf("Warning: time out of range (%d), using value at the grid border.\n",it);
		}	
		
		ixl=ix+1;
		iyl=iy+1;
		itl=it+1;
		
		if(ixl>numx-1) {ixl=numx-1;
		//printf("Warning: long. (ixl) out of range (%d).\n",ixl);
		}
		if(iyl>numy-1) {iyl=numy-1;
		//printf("Warning: lat. (iyl) out of range (%d).\n",iyl);
		}
		if(itl>numt-1) {itl=numt-1;
		//printf("Warning: time (itl) out of range (%d).\n",itl);
		}	
		
//printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
//printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
				
		xpi[0]=datax[ix+iy*numx+it*numx*numy];
		xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
		xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
		xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
		
		xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
		xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
		xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
		xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
		
		ypi[0]=datay[ix+iy*numx+it*numx*numy];
		ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
		ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
		ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
		
		ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
		ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
		ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
		ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
		
		xp1=intpl(xpi,lambdax,lambday,lambdat);
		yp1=intpl(ypi,lambdax,lambday,lambdat);
		
		xyp[0]=(double)xp1;
		xyp[1]=(double)yp1;
		return 0;
		}
	} Pco2_field;



// FUNCTIONS FOR HANDLING THE OBJECTS, INTEGRATING ETC.



int select_sys(int sys_code) {
	int ret=0;
	switch(sys_code) {
		case 1: pField=&Meanderingjet; break;
		case 11: pField=&Periodicflow; break;
		case 2: pField=&Linearsaddle; break;
		case 3: pField=&Vonkarman; break;
		case 31: pField=&Vonkarman2; break;
		case 32: pField=&Vonkarman_m; break;
		case 4: pField=&Lut_field; break;
		case 5: pField=&Imagen_field; break;
		case 6: pField=&Altimetry_field; break;
		case 61: pField=&Aviso_Med_field; break;
		case 7: pField=&Laminarflow; break;
		case 8: pField=&WindERA40_field; break;
		case 9: pField=&Pomme_field; break;
		case 1001: pField=&Kerguelen_field; break;
		case 1011: pField=&Pco2_field; break;
		case 2001: pField=&Homocline; break;
		default: pField=&Genericfield; break;
		}
	return ret;
	}

void pol2cart(double lon,double lat,double *x,double *y,double *z) {
	*z=RTerra*sin(lat/360.0*2*pi);
	*x=RTerra*cos(lat/360.0*2*pi)*cos(lon/360.0*2*pi);
	*y=RTerra*cos(lat/360.0*2*pi)*sin(lon/360.0*2*pi);
	}

void cart2pol(double x,double y,double z,double *lon,double *lat) {
	/*
	if(z>RTerra) z=RTerra;
	else if(z<-RTerra) z=-RTerra;
	*lat=asin(z/RTerra)*360.0/(2*pi);
	printf("cart2pol:\n lat:%lf lat2:%lf\n",*lat,atan2(z,sqrt(x*x+y*y))*360/(2*pi));
	*/
	*lat=atan2(z,sqrt(x*x+y*y))*360.0/(2*pi);
	*lon=atan2(y,x)*360.0/(2*pi);
	//printf("cart2pol:\nx:%lf y:%lf z:%lf lon:%lf lat:%lf test:%lf\n",x,y,z,*lon,*lat,atan2(y,x)*360.0/(2*pi));
	}
	
void Tdist(double x1,double y1,double x2,double y2,double *d) {
	int distswitch;
	double dx,dy,dz,dist1,dist,cx1,cx2,cy1,cy2,cz1,cz2,lat1,lon1,lat2,lon2;
	//printf("Tdist\n");
	if(pField->gridtype==0)
		distswitch=0;
	else
		distswitch=1;
		
		
	switch(distswitch) {
		case 1: {
			lon1=x1/360.0*2*pi;
			lat1=y1/360.0*2*pi;
			lon2=x2/360.0*2*pi;
			lat2=y2/360.0*2*pi;
			dist=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*360.0/(2.0*pi);
			} break;
		
		
		/*
		case 2: {
			//printf("Tdist case 2 (Euclidean on the sphere)\n");
			pol2cart(x1,y1,&cx1,&cy1,&cz1);
			pol2cart(x2,y2,&cx2,&cy2,&cz2);
			
			dx=cx2-cx1;
			dy=cy2-cy1;
			dz=cz2-cz1;
			dist1=sqrt(dx*dx+dy*dy+dz*dz);
			dist=2*asin(dist1/(2*RTerra))*360.0/(2*pi);
			
		
			} break;
		*/
		default:
			{
			dx=x2-x1;
			dy=y2-y1;
			dist=sqrt(dx*dx+dy*dy);
			} break;
		}
	*d=dist;
	return;
	}
	

void RK1cartesianstep(double lon,double lat,double h,double u,double v,double *lonn,double *latn)
	{
	double x,y,z,xn,yn,zn;
	pol2cart(lon,lat,&x,&y,&z);
	lat=lat/360*2*pi;
	lon=lon/360*2*pi;
	//NB: It assumes that u and v are in cm/sec, not deg./sec!
	xn=x+(-u*sin(lon)-v*cos(lon)*sin(lat))*h;
	yn=y+(u*cos(lon)-v*sin(lat)*sin(lon))*h;
	zn=z+(v*cos(lat))*h;
	cart2pol(xn,yn,zn,lonn,latn);
	//printf("RK1cartesianstep:\nxn:%lf yn:%lf zn:%lf lonn:%lf latn:%lf\n",xn,yn,zn,*lonn,*latn);
	}
	
void RK4cartesianstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;
	//RK step 1
	halfh=h/2.0;
	
	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	//xtmp=x+xyp1[0]*halfh;
	//ytmp=y+xyp1[1]*halfh;
	RK1cartesianstep(x,y,halfh,xyp1[0],xyp1[1],&xtmp,&ytmp);
		
	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	//xtmp=x+xyp2[0]*halfh;
	//ytmp=y+xyp2[1]*halfh;
	RK1cartesianstep(x,y,halfh,xyp2[0],xyp2[1],&xtmp,&ytmp);
		
	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	//xtmp=x+xyp3[0]*h;
	//ytmp=y+xyp3[1]*h;
	RK1cartesianstep(x,y,h,xyp3[0],xyp3[1],&xtmp,&ytmp);
		
	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp4);
	
	//newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	//newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
	RK1cartesianstep(x,y,h/6.0,xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0],xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1],newstatus,newstatus+1);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Euler method, next line only for testing! 
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//RK1cartesianstep(x,y,h,xyp1[0],xyp1[1],newstatus,newstatus+1);
	
	//done in RK4degstep
	//while(y>90) y=180-y;
	//while(y<-90) y=-180-y;
	//while(x>179) x=x-360;
	//while(x<-180) x=x+360;
		
	
	}
	

	
void RK4degstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;
	
	//RK step 1
	halfh=h/2.0;
	
	////TEST
	
	while(x>180.0) x-=360.0;
	while(x<-180.0) x+=360.0;
	while(y>90.0) y=180.0-y;
	while(y<-90.0) y=-180.0-y;
	
	
	//////////////
	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	xtmp=x+xyp1[0]*halfh;
	ytmp=y+xyp1[1]*halfh;
	//printf("xyp1[0]=%lf xyp1[1]=%lf\n",xyp1[0],xyp1[1]);	
	
	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	xtmp=x+xyp2[0]*halfh;
	ytmp=y+xyp2[1]*halfh;
		
	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	xtmp=x+xyp3[0]*h;
	ytmp=y+xyp3[1]*h;
		
	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp4);
	
	newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
	
	
	}
	
void RK4flatstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;
	
	//RK step 1
	halfh=h/2.0;
	
	
	//////////////
	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	xtmp=x+xyp1[0]*halfh;
	ytmp=y+xyp1[1]*halfh;
	//printf("xyp1[0]=%lf xyp1[1]=%lf\n",xyp1[0],xyp1[1]);	
	
	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	xtmp=x+xyp2[0]*halfh;
	ytmp=y+xyp2[1]*halfh;
		
	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	xtmp=x+xyp3[0]*h;
	ytmp=y+xyp3[1]*h;
		
	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp4);
	
	newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
	
	
	}
void RK4step(double t,double x,double y,double h,double *newstatus) {
	double noise1,noise2;
	if (pField->gridtype==0) RK4flatstep(t,x,y,h,newstatus);
	else
		if (pField->datatype==2) RK4cartesianstep(t,x,y,h,newstatus);
		else RK4degstep(t,x,y,h,newstatus);
	if(pField->noise) {
			boxmuller(&noise1,&noise2);
			if(pField->datatype==1) noise1=noise1/cos(y/180*pi);
			newstatus[0]+=noise1*(pField->noiselevel)*sqrt(fabs(h));
			newstatus[1]+=noise2*(pField->noiselevel)*sqrt(fabs(h));
			}
	}	

int RK4(double *tspan,double *x0,int numx0,int Nstep,double *trj)
{
/*
Like RK4_tau, but computing the trajectories. The trajectories are
stored in trj as (x1(t1)...x1(tmax) y1(1)...y1(tmax) ... xn(t1)...xn(tmax)).
*/
int j1,j2,indx,indy;
double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,RKout[2],noise1,noise2;

h=(tspan[1]-tspan[0])/(double)Nstep;
halfh=h/2.0;

t=tspan[0];
//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0;j2++) {
	trj[2*j2*Nstep]=x0[2*j2];
	trj[(2*j2+1)*Nstep]=x0[2*j2+1];
	}
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	//printf("j1=%d Nstep=%d\n",j1,Nstep);
	for(j2=0;j2<numx0;j2++) {
		indx=j1+2*j2*Nstep;
		indy=j1+(2*j2+1)*Nstep;
		x=trj[indx-1];
		y=trj[indy-1];

		RK4step(t,x,y,h,RKout);
		
		//printf("RKout[0]=%lf RKout[1]=%lf\n",RKout[0],RKout[1]);
		trj[indx]=RKout[0];
		trj[indy]=RKout[1];		
		}
	
	}
return 0;
}




int RK4_tau(double *tspan,double *x0,int numx0,int Nstep,double delta,int *elmax,double *tau)
{
/*
This function integrates a number numx0 of 2d points which initial conditions are
contained in the vector x0 as: (x1 y1 x2 y2 x3 y3... xn yn). The idea
is to put in x1 y1 a point, and in the other variables other points around
x1 y1. The function computes the distance of all trajectories from x1 y1 and
stops when one has a distance greater than delta. It then provides the time
(in tau) and the number of the trajectory (in elmax). (I.e., "3" if the trajectory started)
at point x3 y3.) The integrator is a Runge-Kutta of 4th order, with fixed
time step (Nstep number of steps in the time window tspan[1], tspan[2]).
*/
int j1,j2,indx,indy;
double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,RKout[2];
double *newstatus,*oldstatus,*swapstatus;

*tau=0.0;
*elmax=-1;

newstatus=new double[numx0*2];
oldstatus=new double[numx0*2];

h=(tspan[1]-tspan[0])/(double)Nstep;
halfh=h/2.0;

t=tspan[0];
//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	tfull=t+h;
	

	for(j2=0;j2<numx0;j2++) {
		x=oldstatus[j2*2];
		y=oldstatus[j2*2+1];

		RK4step(t,x,y,h,RKout);

		newstatus[2*j2]=RKout[0];
		newstatus[2*j2+1]=RKout[1];
		
		if(j2>0)
			{
			//dx=newstatus[2*j2]-newstatus[0];
			//dy=newstatus[2*j2+1]-newstatus[1];
			//dist=sqrt(dx*dx+dy*dy);
			Tdist(newstatus[0],newstatus[1],newstatus[2*j2],newstatus[2*j2+1],&dist);
			if(dist>delta) {
				
				*elmax=j2;
				*tau=tfull-tspan[0]; //ADDED: -tspan[0]
				j2=numx0*2;
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

int RK4_delta(double *tspan,double *x0,int numx0,int Nstep,double *direxp,double *delta,double *dirx0)
{
/*
Like RK4_tau, but returning the maximum distance between the trajectories
starting in x2 y2... xn yn and the one starting in x1 y1.
Used for the Finite Time Lyapunov Exponent.
*/
int j1,j2,indx,indy;
double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,lastf[2];
double *newstatus,*oldstatus,*swapstatus;



newstatus=new double[numx0*2];
oldstatus=new double[numx0*2];

h=(tspan[1]-tspan[0])/(double)Nstep;
halfh=h/2.0;

t=tspan[0];
//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	//printf("t=%f\n",t);
	for(j2=0;j2<numx0;j2++) {
		x=oldstatus[j2*2];
		y=oldstatus[j2*2+1];

		RK4step(t,x,y,h,newstatus+2*j2);
		
		}

		
		

	swapstatus=oldstatus;
	oldstatus=newstatus;
	newstatus=swapstatus;
	}

*delta=0.0;
*direxp=-1.0;
pField->xyp(t,newstatus[0],newstatus[1],lastf);
*dirx0=2*pi+atan(lastf[1]/lastf[0]);
for(j2=1;j2<numx0;j2++) {
	dx=newstatus[2*j2]-newstatus[0];
	dy=newstatus[2*j2+1]-newstatus[1];
	//dist=sqrt(dx*dx+dy*dy);
	
	Tdist(newstatus[0],newstatus[1],newstatus[2*j2],newstatus[2*j2+1],&dist);
	if(dist>(*delta)) {
		*delta=dist;
		*direxp=2*pi+atan(dy/dx);
		}
	}
		

delete oldstatus;
delete newstatus;
if (dist==0) {
	printf("RK4_delta:\nError! dist=0 dx:%lf dy:%lf u:%lf v:%lf \n",dx,dy,lastf[0],lastf[1]);
	printf("x0_x%lf x0_y%lf\n",x0[0],x0[1]);
	}
	
	
return 0;
}

int RK4_mix(double *tspan,double *x0,int numx0,int Nstep,double sc,int *flagct,int *ndisk,int *quality)
{
/*
Like RK4_delta, but returning the number of discs of radius sc needed to cover the trajectories
starting in x1 y1... xn yn.
*/
int j1,j2,indx,indy;
double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,lastf[2];
double *newstatus,*oldstatus,*swapstatus;



newstatus=new double[numx0*2];
oldstatus=new double[numx0*2];

h=(tspan[1]-tspan[0])/(double)Nstep;
halfh=h/2.0;

t=tspan[0];
//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
//printf("numx0=%d\n",numx0);
for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];
	
for (j1=1;j1<Nstep;j1++) {
	
	t=tspan[0]+(j1-1)*h;
	//printf("t=%f\n",t);
	for(j2=0;j2<numx0;j2++) {
		x=oldstatus[j2*2];
		y=oldstatus[j2*2+1];

		RK4step(t,x,y,h,newstatus+2*j2);
		
		}

		
		

	swapstatus=oldstatus;
	oldstatus=newstatus;
	newstatus=swapstatus;
	}

*ndisk=0;

scalecount(oldstatus,numx0,sc,flagct,ndisk,quality);	

delete oldstatus;
delete newstatus;
if ((*ndisk<=0)||(*ndisk>numx0)) {
	printf("RK4_mix:\nError!: wrong disk number (ndisk=%d)!\n",*ndisk);
	//printf("dx:%lf dy:%lf u:%lf v:%lf \n",dx,dy,lastf[0],lastf[1]);
	//printf("x0_x%lf x0_y%lf\n",x0[0],x0[1]);
	}
	
	
return 0;
}



int scalecount(double *v,int N,double sc,int *flagct,int *ndisk,int *quality) {

int ct,ct1,ct2,already_covered;
double dist;
*ndisk=0;
flagct[0]=1;
for(ct=1;ct<N;ct++) flagct[ct]=0; 

for(ct1=1;ct1<N;ct1++) {
	//printf("scalecount:\nct1:%d\n",ct1);
	already_covered=0;
	for(ct2=0;ct2<ct1;ct2++) if(flagct[ct2]){
		Tdist(v[2*ct1],v[2*ct1+1],v[2*ct2],v[2*ct2+1],&dist);
		//printf("scalecount:\nct1:%d ct2:%d dist:%lf\n",ct1,ct2,dist);
		if(dist<sc) {
			already_covered=1;
			ct2=ct1;
			}
		}
	if(already_covered==0) {flagct[ct1]=1;*ndisk++;}
	}
//quality: not yet implemented!
return 0;
}
	


// OTHER GLOBAL VARIABLES


