#include <stdio.h>
#include <iostream>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

#include "Poly.h"
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
     
 int main(void)
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
