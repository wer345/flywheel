#include "Poly.h"
#include <stdlib.h>
#include <memory.h>

Poly::Poly(int size)
{
	m_size=size;
	m_c=(double*) malloc(size*sizeof(double));
	for(int i=0;i<size;i++)
		m_c[i]=0.0;
}

Poly& Poly::operator=(const Poly &v) {

	if(m_size!=v.m_size) {
		free(m_c);
		m_size=v.m_size;
		m_c=(double*) malloc(v.m_size*sizeof(double));
	}
	for(int i=0;i<m_size;i++)
		m_c[i]=v.m_c[i];

    return *this;  
/*
	Poly result = Poly(m_size);
	for(int i=0;i<m_size;i++)
		result.m_c[i]=v.m_c[i];
	return result;
*/
  }

void Poly::resize(int size)
{
	if(m_size!=size) {
		m_size=size;
		double * m_c1=(double*) malloc(size*sizeof(double));
		for(int i=0;i<size;i++) 
			m_c1[i]=m_c[i];
		free(m_c);
	}
}

Poly Poly::operator + (const Poly &other) 
	{
	int size=m_size;
	if(size<other.m_size) 
		size=other.m_size;
	Poly result = Poly(size);

	for(int i=0;i<m_size;i++) 
		result.m_c[i]=m_c[i]+other.m_c[i];

    return result; 
  }

Poly Poly::operator - (const Poly &other) 
	{
	int size=m_size;
	if(size<other.m_size) 
		size=other.m_size;
	Poly result = Poly(size);

	for(int i=0;i<m_size;i++) 
		result.m_c[i]=m_c[i]-other.m_c[i];

    return result; 
  }


Poly Poly::operator * (const Poly &other) 
{
	int size=m_size+other.m_size-1;
	Poly result = Poly(size);

	for(int i=0;i<m_size;i++) 
		for(int j=0;j<other.m_size;j++) 
		{
			result.m_c[i+j] +=m_c[i]*other.m_c[j];
		}

    return result; 
  }

Poly Poly::operator * (const double &other) 
{
	Poly result = Poly(m_size);

	for(int i=0;i<m_size;i++) 
			result.m_c[i] =m_c[i]*other;
    return result; 
  }

Poly Poly::derivative()
{
	if(m_size<=1)
		return Poly(1);

	Poly result= Poly(m_size-1);
	for(int i=0;i<m_size-1;i++)
		result.m_c[i]=(i+1)*m_c[i+1];
	return result;
}

Poly Poly::integral()
{
	Poly result= Poly(m_size+1);
	result.m_c[0]=0.0;
	for(int i=1;i<m_size+1;i++)
		result.m_c[i]=m_c[i-1]/i;
	return result;
}

double Poly::value(const double x)
{
	if(m_size<=0)
		return 0.0;
	double v=m_c[m_size-1];

	for(int i=m_size-2;i>=0;i--)
		v=m_c[i]+x*v;
	return v;
}
/*
  return 1/(1+this) 
  level is the number of interation
*/
Poly Poly::inverse(int level)
{
	Poly result=Poly(1);
	result.m_c[0]=1.0;

	while(level-->0) {
//		cout<<"result 0:" << result.toString()<<endl;
		Poly r=((*this)*result)*(-1);
		r.m_c[0]=1.0;
//		cout<<"result 1:" << r.toString()<<endl;
		result=r;

	}
	return result;
}
