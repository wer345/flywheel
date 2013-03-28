#pragma once

#include <iostream>
#include <sstream>
#include <string>
using namespace std;
class Poly 
{
public:
	Poly(int size);
	int m_size;
	double * m_c;
	void resize(int size);

	Poly & operator = (const Poly &v);
	Poly operator + (const Poly &other);
	Poly operator - (const Poly &other);
	Poly operator * (const Poly &other);
	Poly Poly::operator * (const double &other) ;

	Poly derivative();
	Poly integral();
	Poly inverse(int level);
	double value(const double x);

	const std::string toString()
	{
		std::stringstream ss;
		ss << "[" << m_size << "]";
		for (int i=0;i<m_size;i++) {
			if(i==0)
				ss<<": ";
			else 
				ss<<", ";
			ss<<m_c[i];
		}
		return ss.str();
	}

};


