//#define hypot gsl_hypot

#include <stdio.h>
#include <iostream>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "Poly.h"

 double testfunc(double x)
 {
//	 return cos(x);
//	 return 1/(1+x);
	 return sqrt(1+x);

 }
 
 void outputPoly(FILE *f, int size,int idx,gsl_vector* data)
 {
	fprintf(f,"%10.8f",gsl_vector_get(data,idx));

	if(idx+1<size)
	{
		fprintf(f,"+x*(");
		outputPoly(f,size,idx+1,data);
		fprintf(f,")");
	}
 }

 int fit(int p)
     {
       int i, n;
       double xi, yi, ei, chisq;
       gsl_matrix *X, *cov;
       gsl_vector *y, *w, *c;
     
       n = 15;
     
       X = gsl_matrix_alloc (n, p);
       y = gsl_vector_alloc (n);
       w = gsl_vector_alloc (n);
     
       c = gsl_vector_alloc (p);
       cov = gsl_matrix_alloc (p, p);

	   // range of x
       double xs=-0.5;
       double xe=0.5;

       for (i = 0; i < n; i++)
         {
//		   xi=3.1415926/4.0*i;
			 double step=(xe-xs)/(n-1);
			 xi=xs+step*i;
			yi=testfunc(xi);
			ei=0.01;
     
           printf ("%6.3f %6.3f \n", xi, yi);

		   double v=1.0;
		   for(int j=0;j<p;j++) {
	           gsl_matrix_set (X, i, j, v);
			   v*=xi;
		   }
/*
           gsl_matrix_set (X, i, 0, 1.0);
           gsl_matrix_set (X, i, 1, xi);
           gsl_matrix_set (X, i, 2, xi*xi);
*/           
           gsl_vector_set (y, i, yi);
           gsl_vector_set (w, i, 1.0/(ei*ei));
         }
     
       {
         gsl_multifit_linear_workspace * work 
           = gsl_multifit_linear_alloc (n, p);
//         gsl_multifit_linear (X, y, c, cov,&chisq, work);
         gsl_multifit_wlinear (X, w, y, c, cov,&chisq, work);
         gsl_multifit_linear_free (work);
       }
     
     #define C(i) (gsl_vector_get(c,(i)))
     #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
     
       {
         printf ("# best fit: Y = %g + %g X + %g X^2\n", C(0), C(1), C(2));
		   printf("# best fit: Y=");
		   for(int i=0;i<p;i++) {
			   if(C(i)>=0) 
				   printf("+");
			   printf("%g",C(i));
			   int j=i;
			   while(j-->0)
				   printf("*x");

		   }
         printf ("\n# chisq = %g\n", chisq);
       }
     
	   char fn[200];
	   sprintf(fn,"..\\\\graph\\\\graph-p%d.plt",p);
		FILE * f=fopen(fn,"w");
		if(f!=NULL) {
			fprintf(f,"set key left box\nset samples 150\n");
			fprintf(f,"plot [%g:%g]",xs,xe);
//			fprintf(f," 1000*(1/(1+x)-(");
			fprintf(f," 1000*(sqrt(1+x)-(");
			outputPoly(f,p,0,c);
			fprintf(f," ))");

			fprintf(f,"\n");
			fclose(f);
		}

       gsl_matrix_free (X);
       gsl_vector_free (y);
       gsl_vector_free (w);
       gsl_vector_free (c);
       gsl_matrix_free (cov);
     
       return 0;
     }




// (x'y''-x''y')^2 = (x'^2+y'^2)^3 *{[(xe-2x)*Fy-(ye-2y)Fx)]/(2EI)}^2

     /* Paraboloid centered on (p[0],p[1]), with  
        scale factors (p[2],p[3]) and minimum p[4] */
     
     double
     my_f (const gsl_vector *v, void *params)
     {
       double x, y;
       double *p = (double *)params;
       
       x = gsl_vector_get(v, 0);
       y = gsl_vector_get(v, 1);
      
       return p[2] * (x - p[0]) * (x - p[0]) +
                p[3] * (y - p[1]) * (y - p[1]) + p[4]; 
     }

void test_poly()
{
	Poly a= Poly(3);
	Poly b= Poly(3);
	a.m_c[0]=1;
	b.m_c[1]=2;
	a.m_c[2]=1;
	b.m_c[2]=5;
	cout<<"a=" << a.toString()<<endl;
	cout<<"b=" << b.toString()<<endl;
	Poly c=b;
	cout<<"c=" << c.toString()<<endl;
	c=a+b;
	cout<<"a+b=" << c.toString()<<endl;
	c=a-b;
	cout<<"a-b=" << c.toString()<<endl;
	c=a*b;
	cout<<"a*b=" << c.toString()<<endl;

	Poly d=c.derivative();
	cout<<"d=derivative of c" << d.toString()<<endl;

	Poly e=d.integral();
	cout<<"e=integral of d" << e.toString()<<endl;

	Poly xx= Poly(3);
	xx.m_c[1]=1;
	xx.m_c[2]=1;
//	xx.m_c[3]=0;

	Poly inv= xx.inverse(60);
	cout<<"inv(x)" << inv.toString()<<endl;
	double x0=-0.3;
	double v1=1/(1+x0+x0*x0);
	double v2=inv.value(x0);
	cout<<"v1=" << v1<<"; v2=" << v2<<endl;
}

//swing

void swing()
{
	gsl_matrix_float * A = gsl_matrix_float_calloc(2,2);

	gsl_matrix_float_set(A,0,0,2.0);
	gsl_matrix_float_set(A,0,1,1.0);

	gsl_matrix_float_set(A,1,0,1.0);
	gsl_matrix_float_set(A,1,1,1.0);

	gsl_vector_float * v=gsl_vector_float_calloc(2);
	gsl_vector_float_set(v,0,10.0);
	gsl_vector_float_set(v,1,20.0);

	gsl_blas_strsv(CblasUpper,CblasNoTrans,CblasNonUnit,A,v);
//	gsl_blas_strsv(CblasLower,CblasNoTrans,CblasNonUnit,A,v);

	cout<<"v[0]="<<gsl_vector_float_get(v,0)<<"v[1]="<<gsl_vector_float_get(v,1)<<endl;
}



typedef struct _SwingData {
	double a1,a2,a3,m1,m2,m3,s1,s2,s3;
	double r1,r2,r3;
	double F1,F2,F3,P1,P2;
} SwingData;

typedef struct _SwingVars {
	double f1,f2,f3,ac1,ac2;
} SwingVars;

int swing_step(SwingData &d,SwingVars &vars)
{
	double rss1s=d.r1*d.s1*d.s1*sin(d.a1);
	double rss1c=d.r1*d.s1*d.s1*cos(d.a1);
	double rss2s=d.r2*d.s2*d.s2*sin(d.a2);
	double rss2c=d.r2*d.s2*d.s2*cos(d.a2);
	double rss3c=d.r3*d.s3*d.s3*cos(d.a3);
/*
	double a_data[] = { 
		-sin(d.a1)/d.m1,	sin(d.a2)/d.m1,		-sin(d.a3)/d.m1,	-d.r1*cos(d.a1),	0,
		-cos(d.a1)/d.m1,	-cos(d.a2)/d.m1,	cos(d.a3)/d.m1,		d.r1*sin(d.a1),		0,
		0,					-sin(d.a2)/d.m2,	0,					-d.r1*cos(d.a1),	-d.r2*cos(d.a2),
		0,					cos(d.a2)/d.m2,	0,						d.r1*sin(d.a1),		-d.r2*sin(d.a2),
		0,					0,					-cos(d.a3)/d.m3,	d.r1*sin(d.a1),		d.r3*sin(d.a3),
	};

	double b_data[] = { 
		-d.P1/d.m1-rss1s,
		d.F1/d.m1-rss1c,
		-d.P2/d.m2-rss1s+rss2c,
		d.F2/d.m2-rss1c+rss2c,
		-d.F3/d.m3-rss1c-rss3c
	};
*/
	double res=60.0;

	double a_data[] = { 
		sin(d.a1)/d.m1,	-sin(d.a2)/d.m1,	sin(d.a3)/d.m1,		d.r1*cos(d.a1),		0,
		cos(d.a1)/d.m1,	cos(d.a2)/d.m1,		-cos(d.a3)/d.m1,	-d.r1*sin(d.a1),	0,
		0,				sin(d.a2)/d.m2,		0,					d.r1*cos(d.a1),		d.r2*cos(d.a2),
		0,				-cos(d.a2)/d.m2,	0,					-d.r1*sin(d.a1),	d.r2*sin(d.a2),
		0,				0,					-cos(d.a3)/d.m3,	d.r1*sin(d.a1),		d.r3*sin(d.a3),
	};

	double b_data[] = { 
		d.P1/d.m1-rss1s,
		-d.F1/d.m1-rss1c,
		d.P2/d.m2-rss1s-rss2s-res*(d.r1*d.s1*cos(d.a1)+d.r2*d.s2*cos(d.a2)),
		-d.F2/d.m2-rss1c+rss2c,
		-d.F3/d.m3-rss1c-rss3c
	};


	gsl_matrix_view m = gsl_matrix_view_array (a_data, 5, 5);
	gsl_vector_view b = gsl_vector_view_array (b_data, 5);
	gsl_vector *x = gsl_vector_alloc (5);
	int s;
/*
	printf ("A=\n");
	for(int i=0;i<5;i++)  {
		for(int j=0;j<5;j++) {
			printf ("%8.4f, ",gsl_matrix_get(&m.matrix,i,j));
		}
		printf ("\n");
	}
	printf ("b=\n");
	for(int j=0;j<5;j++) {
		printf ("%d: %8.4f\n",j,gsl_vector_get(&b.vector,j));
	}
*/
	gsl_permutation * p = gsl_permutation_alloc (5);
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	vars.f1=gsl_vector_get(x,0);
	vars.f2=gsl_vector_get(x,1);
	vars.f3=gsl_vector_get(x,2);
	vars.ac1=gsl_vector_get(x,3);
	vars.ac2=gsl_vector_get(x,4);

//	printf ("x = \n");
//	gsl_vector_fprintf (stdout, x, "%8.5f");

	gsl_permutation_free (p);
	gsl_vector_free (x);

	return 0;
}

 int main(void)
{

//	test_poly();

//	for(int p=9;p<14;p++)
//		fit(p);

//	swing();

	double m1=0.005;
	double m2=0.01;
	double m3=0.018;
	double dt=0.0001;
	double rps=10; // round per second
	double rspeed=2*3.141596*10; // angle speed
	double g=9.8;
	double r=0.1;
	double gr=r*rspeed*rspeed;

	SwingData d = {
//	double a1,a2,a3,m1,m2,m3,s1,s2,s3;
		0.0,0.0,0.0,
		m1,m2,m3,
		0.0,0.0,0.0,

//	double r1,r2,r3;
		0.008,0.01,0.05,
//	double F1,F2,F3,P1,P2;
		gr*d.m1,gr*d.m2,gr*d.m3,
		m1,m2
	};

	SwingVars vars;

	char fn[200];
	sprintf(fn,"..\\\\graph\\\\data\\\\swing.dat");
	FILE * f=fopen(fn,"w");

	for(int i=0;i<20000;i++) {
		double t=i*dt;
		d.a3=asin(d.r1*sin(d.a1)/d.r3);
		d.s3=d.r1*d.a1*cos(d.a1)/(d.r3*cos(d.a3));
		double r=cos(rspeed*t);
		d.P1=g*d.m1*r;
		d.P2=g*d.m2*r;
		printf("%3d: P=(%5.3f, %5.3f) a=(%5.3f, %5.3f, %5.3f)  s=(%5.2f, %5.2f, %5.2f)  ",i,d.P1,d.P2, d.a1,d.a2,d.a3,d.s1,d.s2,d.s3);
		// print P1, P2, a1, a2, s1,s1, ac1, ac2
		fprintf(f,"%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f ",i,d.P1,d.P2, d.a1,d.a2, d.s1,d.s2 );
		swing_step(d,vars);
		printf("f=(%6.4f, %6.4f, %6.4f)  ac=(%6.2f, %6.2f)\n",vars.f1,vars.f2,vars.f3,vars.ac1,vars.ac2);
		fprintf(f,"  %8.5f %8.5f\n",vars.ac1,vars.ac2);
		d.a1+=d.s1*dt;
		d.a2+=d.s2*dt;
		d.s1+=vars.ac1*dt;
		d.s2+=vars.ac2*dt;
	}

	fclose(f);
	char ch=getchar();
//	char in;
//	cin>>in;
	exit(0);
//==============================================
    double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
     
    const gsl_multimin_fminimizer_type *T = 
         gsl_multimin_fminimizer_nmsimplex;
       gsl_multimin_fminimizer *s = NULL;
       gsl_vector *ss, *x;
       gsl_multimin_function minex_func;
     
       size_t iter = 0;
       int status;
       double size;
     
       /* Starting point */
       x = gsl_vector_alloc (2);
       gsl_vector_set (x, 0, 5.0);
       gsl_vector_set (x, 1, 7.0);
     
       /* Set initial step sizes to 1 */
       ss = gsl_vector_alloc (2);
       gsl_vector_set_all (ss, 1.0);
     
       /* Initialize method and iterate */
       minex_func.n = 2;
       minex_func.f = my_f;
       minex_func.params = par;
     
       s = gsl_multimin_fminimizer_alloc (T, 2);
       gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
     
       do
         {
           iter++;
           status = gsl_multimin_fminimizer_iterate(s);
           
           if (status) 
             break;
     
           size = gsl_multimin_fminimizer_size (s);
           status = gsl_multimin_test_size (size, 1e-2);
     
           if (status == GSL_SUCCESS)
             {
               printf ("converged to minimum at\n");
             }
     
           printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
                   iter,
                   gsl_vector_get (s->x, 0), 
                   gsl_vector_get (s->x, 1), 
                   s->fval, size);
         }
       while (status == GSL_CONTINUE && iter < 100);
       
       gsl_vector_free(x);
       gsl_vector_free(ss);
       gsl_multimin_fminimizer_free (s);
     
       return status;
     }
