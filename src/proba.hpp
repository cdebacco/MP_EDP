/*
   Copyright 2009-2010 Alfredo Braunstein and Riccardo Zecchina

   This file is part of MSGSTEINER (Max Sum for generalized steiner problems on graphs).

   MSGSTEINER is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   MSGSTEINER is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MSGSTEINER; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */




#ifndef PROBA_H
#define PROBA_H

#include <iostream>
#include <assert.h>
#include <math.h>

//double const inf = std::numeric_limits<double>::max();
double const inf = 100000.;// 1e10;

class Proba {
public:	
  Proba(int depth) : p_(new double[depth]), depth(depth) {
		assert(depth>=0);
		for (int d = 0; d<depth;d++ )
			p_[d] = 0.;
	}
  
  Proba(int depth, double value) : p_(new double[depth]), depth(depth) {
		assert(depth>=0);
		for (int d = 0; d<depth;d++ )
			p_[d] = value;
	}

	Proba(Proba const & other) : p_(new double[other.depth]), depth(other.depth) {
		*this = other;
	}
	Proba & operator=(Proba const & other) {
		for (int d = 0; d<depth;d++ )
			p_[d] = other[d];
		return *this;
	}
	~Proba() { delete[] p_; }
	double & operator[](int d) {
		return p_[d];
	}
	double const & operator[](int d) const{
		return p_[d];
	}
	friend std::ostream & operator<<(std::ostream & ost, Proba const & p) {
		for (int d = 0; d < p.depth; ++d)
			ost << p[d] << " ";
		return ost;
	}
	Proba & operator+=(Proba const & b) {
		for (int d = 0; d < depth; ++d)
			p_[d] += b[d];
		return *this;
	}

	Proba & operator*=(double b) {
		for (int d = 0; d<depth;d++ )
			p_[d] *= b;
		return *this;
	}

	friend double l8dist(Proba const & a, Proba const & b) {
		double n = 0.;
		for (int d = 0; d < a.depth; ++d)
			n = std::max(n, fabs(a[d] - b[d]));
			//n += fabs(a[d] - b[d]);
		return n;
	}

	friend void swap(Proba & A, Proba & B) {
		std::swap(A.depth, B.depth);
		std::swap(A.p_, B.p_);
	}
 	
private:
	double * p_;
   int depth;
};


double l8dist(Proba const & a,  Proba const & b);



#endif

