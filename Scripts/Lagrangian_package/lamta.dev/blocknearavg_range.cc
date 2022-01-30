#include "oct.h"



DEFUN_DLD (blocknearavg_range, args, ,
  "Octave file out=blocknearavg_range(M,n,r1,r2).")
{

int i,j,i2,j2,n,x,y,ctavg;
double avg,r1,r2;
Matrix M=args(0).matrix_value();
y=args(0).columns();
x=args(0).rows();

//printf("x:%d y:%d\n",x,y);
n=(int)args(1).scalar_value();
r1=(double)args(2).scalar_value();
r2=(double)args(3).scalar_value();

Matrix out(x,y);

for(i=0;i<x;i++) for(j=0;j<y;j++) {
		avg=0;
		ctavg=0;
		//printf("i:%d j:%d\n",i,j);
		for(i2=i-n;i2<=(i+n);i2++) for(j2=j-n;j2<=(j+n);j2++) 
			if((i2>=0)&&(i2<x)&&(j2>=0)&&(j2<y)&&(M(i2,j2)>=r1)&&(M(i2,j2)<=r2)) {
				//printf("i2:%d j2:%d\n",i2,j2);
				avg+=M(i2,j2);
				ctavg++;
				}
		if(ctavg) avg=avg/ctavg;
		else avg=M(i,j);
		
		out(i,j)=avg;
		}
		
return octave_value(out);	
		
}		
