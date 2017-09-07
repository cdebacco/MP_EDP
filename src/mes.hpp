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




#ifndef MES_H
#define MES_H



#include "proba.hpp"

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>


struct Mes {
  Mes(int depth) : E(depth), depth_(depth) {}
  Mes(int depth, double value) : E(depth,value), depth_(depth) {}
  int depth() const { return depth_; }

  Proba E;

  friend std::ostream & operator<<(std::ostream & ost, Mes const & m);

  double minimum() const {
    double m = inf;
    for (int d = 0; d <depth_ ; ++d)
      m = std::min(m,E[d]);
    return m;
  }

  Mes & operator=(Mes const & other) {
    E = other.E;
    depth_ = other.depth_;
    return *this;
  }
	
  void reduce() {
    //double m = minimum();
	  double m= E[0];
    // assert(m >=0.);

    for (int d = depth_; d--; ) {
      E[d] -= m;
    }
  }
  friend void swap(Mes & u, Mes & v)
  {
    swap(u.E, v.E);
  }


private:
  int depth_;
};


double l8dist(Mes const & a,  Mes const & b);

Mes operator+(Mes const & a, Mes const & b);

Mes operator*(double c, Mes const & b);


template<class R>
void init_rand(Mes & m, R & mes_real01) {
  int bound=m.depth();
  for (int d = bound; d--; ) {
    m.E[d] = 0*mes_real01();
	}
}


#endif

