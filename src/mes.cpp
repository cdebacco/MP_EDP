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




#include "mes.hpp"


std::ostream & operator<<(std::ostream & ost, Mes const & m)
{
	ost << "E: " << m.E << std::endl;
	return ost;
}

double l8dist(Mes const & a,  Mes const & b)  
{
	double l8 = 0.;
//	if (a.B < -inf/2 || b.B < -inf/2)
//		return inf;
//	if (a.D < -inf/2 || b.D < -inf/2)
//		return inf;
	
	l8 = std::max(l8, l8dist(a.E, b.E));
	//l8=l8dist(a.E, b.E);
	return l8;
}


Mes operator+(Mes const & a, Mes const & b) 
{
	Mes out = a;
	out.E += b.E;
	return out;
}


Mes operator*(double cc, Mes const & b) 
{
	Mes out = b;
	out.E *= cc;
	return out;
}
/*
void init_rand(Mes & m) {
  int bound=m.depth();
  for (int d=0; d<bound; d++) {
    m.E[d] =1.;
  }
}
*/

