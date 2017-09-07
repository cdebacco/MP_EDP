/*
   Copyright 2015 Caterina De Bacco

   This file is part of Edge-Noverlap  (Cavity equations for not-overlapping communications in problems on graphs with edge constraints).

- MATCHING ROUTINE
- DIRECT EDGE (S,R)
- CALCULATE MAX M ACCOMODATED
  
  If you use this software, please cite:

  Altarelli, F., Braunstein, A., Dall’Asta, L., De Bacco, C., & Franz, S. (2015). The edge-disjoint path problem on random graphs by message-passing. PloS one, 10(12), e0145222.

  MP_EDP is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  MP_EDP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Noverlap; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <lemon/matching.h>
#include <lemon/full_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/math.h>

#include <boost/config.hpp>

#include "mes.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>

#include <fstream>
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <math.h>
#include <iomanip>
#include <boost/limits.hpp>
#include <queue>
#include <time.h>
#include <ctime> 

using namespace boost;
using namespace std;

/* ----------------GRAPH VARIABLES ----------------*/
int V=500;
int M=2;                        // number of communications
int N_real=500;                 // number of realizations
int depth= 2*M+1;               // 
int k=3;  

/* ----------------CONVERGENCE VARIABLES ----------------*/
bool tree = false;
int decision = 15;      //  program converges after this # repeats of the decision variables
double beta = 0.;       // reinforcement parameter
double rein = 0.000;    // rein= beta*it
int maxit = 30;         // max iteration time
int T_eq=10;            // check convergence only after this time
double tolerance = 0.01;      //  if err< tolerance than coincide++
double bias_intensity=0.;     // communication bias
double partition=0;
unsigned int seed=0;          // gen_int seed --> instance seed
unsigned int rseed=0;         // double gen seed and mes seed --> randomizes bias and initial messages
double inf2=100; // to be used for a fake direct edge S,R
double inf1=1e10; // subs inf so that higher gap btw inf2 and inf1
double err_max=1e300; // max err in the converge() routine

/* ----------------RESULTS VARIABLES ----------------*/
double Elink=0.,Enode=0.;      // energy per link and per node
vector<double> E_n;
double tot_lenght=0;          //  sum_i L_i    , i=1,...,M and L_i is the length of communication i
double single_tot_length=0;   //  same as tot_lenght but considering shortest path lengths 
int direct_path=0;  // number of direct paths S,R
bool flag_SR=0;     // if flag_SR==1 then input pairs (S,R) from file; else generate them randomly
bool flag_output_messages=0;  // if 1 then output_messages at convergence
/* ----------------I/O FILES ----------------*/
std::ofstream out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,outd;
std::ifstream in1, in2,in3;
std::string folder="steiner16 ";
std::string output_folder="/output/";
std::string com="1";
std::string input_folder="/input/";
std::string instance="50ins1.dat";


//--------------------------------------------------------------------
/* ---------------- HISTOGRAM ROUTINES: needed to collect optional stats about the final path lenghts -----------------------------*/
//--------------------------------------------------------------------

void add_histo(int path_size, vector< pair<int,int> >&  histo){
  bool deja_ajoute = false;
  pair<int,int> paire= make_pair(path_size,1);
  for(unsigned int i=0;i<histo.size();i++){
    if(histo[i].first==path_size){histo[i].second++; i=histo.size(); deja_ajoute= true;}
  }
  if(!deja_ajoute){histo.push_back(paire);}
}

//To sort the histogram
struct sort_pred {
  bool operator()(const pair<int,int> &left, const pair<int,int> &right) {
    return left.first < right.first;
  }
};
vector< pair<int,int> >  histo;int Z=0;
vector< pair<int,int> >  histo_d;int Z_d=0;
vector< pair<int,int> >  histo_k; int Z_k=0;


//----------------RANDOM NUMBER GENERATORS----------------------------------------------------

boost::mt19937 gen,mes_gen,gen_int, gen_int2;

uniform_real<double> const uni_dist(0,1);
variate_generator<boost::mt19937 &, uniform_real<> > real01(gen, uni_dist);
variate_generator<boost::mt19937 &, uniform_real<> > mes_real01(mes_gen, uni_dist);

// Pick an int random number btw [min,max] ;  used to assign random pairs (S,R)
int roll_die(int min, int max) {          
  boost::uniform_int<> dist(min, max );
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(gen_int, dist);
  return die();
  //return (int) rand()%(max-min)+min;
}
// Same as roll_die but with a different seed --> used to purmute updated verticesin iterate()
int roll_die2(int min, int max) {
  boost::uniform_int<> dist(min, max );
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(gen_int2, dist);
  return die();
  //return (int) rand()%(max-min)+min;
}

//--------------------------------------------------------------------
/* ----------------COMPLETE GRAPH FOR MATCHING TOOLS -----------------------------*/
//--------------------------------------------------------------------

typedef double real_t;
typedef lemon::FullGraph graph_t;
typedef graph_t::EdgeMap<real_t> weightmap_t;

//cancels out neighbors of i in w, resets everything when going out of context
template<class G, class M, class S>
class TmpEraser {
public:
  TmpEraser(G const & g, M & w, S & save, int i) : g(g), w(w),
                     save(save), i(i)
  {
    for (int k = 0, n = g.nodeNum(); k < n; ++k) if (k != i) {
   save[k] = w[g.edge(g(i), g(k))];
   w[g.edge(g(i), g(k))] = 0;
      }

  }
  ~TmpEraser() {
    for (int k = 0, n = g.nodeNum(); k < n; ++k) if (k != i)
                     w[g.edge(g(i), g(k))] = save[k];
  }
private:
  G const & g;
  M & w;
  S & save;
  int i;
};




//--------------------------------------------------------------------
/* ---------------- GRAPH TOOLS -----------------------------*/

void randomize(Mes & m) {
  int bound=m.depth();
  for (int d = bound; d--; ) {
    m.E[d] =0.*real01();
  }
}


//-------------------------------------------------------------------------------


struct EdgeProperties {
  EdgeProperties() : ij(depth), ji(depth), c(1.), q(0),x(0), bias(2*M+1,0.), cij(2*M+1,0.),cji(2*M+1,0.) {
    randomize(ij);
    randomize(ji);
    //    c=1.+tolerance*0.1*real01();
    c=1.+0.001*real01();
    bias[0]=0.;
    for(unsigned int i=1; i<=M;i++){
      bias[i]=tolerance*bias_intensity*mes_real01();
      bias[i+M]=bias[i];
      cij[i]=bias[i]+c;
      cij[i+M]=bias[i+M]+c;
      cji[i]=bias[i]+c;
      cji[i+M]=bias[i+M]+c;
    }
    cij[0]=bias[0]+c;
    cji[0]=bias[0]+c;
  }
  Mes ij;
  Mes ji;
  double c; // edge cost
  int q;   // # of communications passing through it
  int x;
  vector<double> bias; // bias[mu]  constant since the beginning
  vector<double> cij,cji; // cost edge ij [mu], update with the reinforcement
};

//--------------------------------------------------------------------

struct VertexProperties  {
  VertexProperties() : G(0), c(0){}
  VertexProperties(string const & name) : name(name), type(normal), G(0),c(0){}
  string name;
  enum typeEnum {
    root,
    terminal,
    normal
  };
  typeEnum type;
  bool isroot() const { return  c!=0; }

  int G;    //  communication number such that Lambda_mu !=0
  int c;    //  Lambda=0,+1,-1
};
//--------------------------------------------------------------------

typedef adjacency_list<setS, vecS, undirectedS,
             VertexProperties, EdgeProperties> GraphBase;

typedef graph_traits<GraphBase>::vertex_iterator vertex_iterator;
typedef graph_traits<GraphBase>::out_edge_iterator edge_iterator;
typedef graph_traits<GraphBase>::edge_iterator graph_edge_iterator;
typedef graph_traits<GraphBase>::edge_descriptor Edge;
typedef graph_traits<GraphBase>::vertex_descriptor Vertex;

struct Graph : public GraphBase {
  Vertex rootid;
};

Graph g,g_f;

//--------------------------------------------------------------------

struct VertexProperties1  {
  VertexProperties1() : c(0),d(0), prm(0) {}
  VertexProperties1(string const & name) : name(name), type(normal),c(0), d(0), prm(0){}
  string name;
  enum typeEnum {
    root,
    terminal,
    normal
  };
  typeEnum type;

  int c;      // colour
  int d;      // distance from the source
  int prm;    // permanent colour (instead of removing the vertex)
};

struct EdgeProperties1 {
  EdgeProperties1() :  q(0) {
    q=0;
}
  int q;   // edge colour if a communication passes through it
};


//--------------------------------------------------------------------

typedef adjacency_list<setS,vecS, undirectedS,VertexProperties1,EdgeProperties1> SimpleGraph;

typedef graph_traits<SimpleGraph>::edge_descriptor Edge1;
typedef graph_traits<SimpleGraph>::vertex_descriptor Vertex1;
typedef graph_traits<SimpleGraph>::vertex_iterator vertex_iterator1;
typedef graph_traits<SimpleGraph>::out_edge_iterator edge_iterator1;
typedef graph_traits<SimpleGraph>::edge_iterator graph_edge_iterator1;

SimpleGraph g1;       //  used to calculate shortest paths via BFS

//--------------------------------------------------------------------
/* ---------------- FUNCTIONS -----------------------------*/
//--------------------------------------------------------------------
 
//------------------------------------------------------------------------------------------

// Used to assign random pairs (S,R) if flag_SR==0
void shuffle(int a[],int N)
{
  for (int i=N-1;i>=1;i--)
    {
      int j=(int)roll_die(0,i-1); 
      int ai=a[i];
      int aj=a[j];
      a[i]=aj;
      a[j]=ai;
    }
}

// Used to purmute updated verticesin iterate()
void shuffle2(int a[],int N)
{
  for (int i=N-1;i>=1;i--)
    {
      int j=(int)roll_die2(0,i-1); 
      int ai=a[i];
      int aj=a[j];
      a[i]=aj;
      a[j]=ai;
    }
}

//------------------------------------------------------------------------------
/*--------------- BREADTH-FRIST-SEARCH ALGORITHM WITH PREDECESSOR -------------*/

// Doesn't modify the graph, used to calculate single_paths...

int bfs_single(Vertex1 &  s, Vertex1 &  r, SimpleGraph & g1) 
{
  assert(s!=r);
  vertex_iterator1 vit,vend;
  //erase colours and distances of g1
  for (tie(vit,vend)=vertices(g1); vit!=vend; vit++) {
    g1[*vit].c=0;g1[*vit].d=0;
  }
  
  vector<Vertex1> queue;
  vector< pair<Vertex1,Vertex1> > parent;
  queue.push_back(s);
  parent.push_back(make_pair(s,s));
  g1[s].c=1;

  while (queue.size()!=0 && g1[r].c==0) {
    Vertex1 p=queue[0];
    edge_iterator1 eit,eend;
    for (tie(eit,eend)=out_edges(p,g1); eit!=eend; eit++) {
      Vertex1 u=(target(*eit,g1)==p) ? source(*eit,g1) : target(*eit,g1);
      if (g1[u].c==0) {
   g1[u].c=1;
   g1[u].d=g1[p].d+1; 
   queue.push_back(u);
   parent.push_back(make_pair(u,p));
      }
    }
    //FIFO
    queue.erase(queue.begin());
  } 
  if (g1[r].d!=0) {
    Vertex1 follower=r;
    int length=0;
    while (follower!=s) {
      //build the path
      for (int i=parent.size()-length-1; i>=0; i--) {
   if (parent[i].first==follower) {
     follower=parent[i].second;
     length++;
     break;
   } 
      }
    }

    return g1[r].d;
  }// end if(g1[r].d!=0)

  return g1[r].d;
  
  
}

//-------------------------------------------------------------------------------
inline Mes & getMes(Edge e) 
{
  return source(e, g) < target(e, g) ? g[e].ji : g[e].ij;
}


inline Mes & getMesInv(Edge e) 
{
  return source(e, g) > target(e, g) ? g[e].ji : g[e].ij;
}

inline vector<double>  & getCost(Edge e) 
{
  return source(e, g) < target(e, g) ? g[e].cji : g[e].cij;
}

inline double E_120( Mes const & a, Mes const & b, int mu )
{
  return a.E[mu]+ b.E[mu+M]-a.E[0]- b.E[0];
}


//-----------------OUTPUT MESSAGES--------------------------------------------------------------

void output_messages(ofstream & cout1)
{
  
  vertex_iterator vit, vend;
  for (tie(vit, vend) = vertices(g); vit != vend; ++vit) {

    edge_iterator eit, eend;
    tie(eit, eend) = out_edges(*vit, g);
    
    Vertex v = *vit;
   
    // if(g[v].G==2){
 cout1 << "----------------------------------" << endl 
     << "Vertex " << g[*vit].name << endl << "c: " << g[v].c << endl 
     << "G: " << g[v].G << endl;
    for (; eit != eend; ++eit) {
      Vertex w = target(*eit, g);
      cout1 << g[v].name << " -> " <<  g[w].name << endl;
      Mes m = getMesInv(*eit);
      cout1 << m << endl;
      cout1 << "bias ij: "<<g[*eit].c-1.<<endl;
      unsigned int size=g[*eit].bias.size();
      cout1<<"bias[mu]: ";
      for(unsigned int i=0;i<size;i++){cout1<<g[*eit].bias[i]<<" ";}    
      cout1<<endl;
}
  
  }  
}

//---------------------------------------------------------------------

//------------------------------------------------------------------------------

/* ---------------MATCHING  UPDATE ROUTINE: this is the main Message-Passing routine --------------------------------------------*/
//-------------------------------------------------------------------------------

real_t updateMS(Vertex i){
  int k,a,source;
  if (g[i].isroot()){
    source=(g[i].c==1)? 1:0;    //source=1 for a sender , source=0 for a receiver  
    k=out_degree(i,g)+1;        // in the kth position we put edge ai
    a=g[i].G;                   // name/color/id of the root communication
  }else {k=out_degree(i,g); a=0; }

  vector< vector<double> > cost(k,vector<double>(2*M+1)); // costs of the links i..
  edge_iterator eit,eend;
  int l=0;
  for( tie(eit,eend)=out_edges(i,g);eit!=eend;++eit, l++){
    for(int mu=0;mu<2*M+1;mu++){cost[l][mu]=(i>target(*eit,g)) ? g[*eit].cij[mu] : g[*eit].cji[mu] ;}
  }
  if(g[i].isroot()){
    for(int mu=1;mu<2*M+1;mu++)if(mu!=a){cost[k-1][mu]=inf1;}
    cost[k-1][a]=-inf1; cost[k-1][a+M]=-inf1; 

}

  graph_t const complete_g(k);  
  weightmap_t w(complete_g);
  vector<Mes> h;
  //  h.reserve(k);
  // inserts costs

  for( tie(eit,eend)=out_edges(i,g);eit!=eend;++eit){
     Mes const & in =getMes(*eit);
     h.push_back(in);
  }
  // insert messages (energy of fake edge a
  if(g[i].isroot()){
    Mes in_a(depth, inf);
    in_a.E[a]=(source==1) ? -inf1 : inf1;
    in_a.E[a+M]=(source==1) ? inf1 : -inf1;
    h.push_back(in_a);
 } // end if root

  for(int l=0;l<k;l++){
    for(int nu=0;nu<2*M+1;nu++){h[l].E[nu]=-h[l].E[nu];}     //reverse sign so to pass to a max instead of a min
    for(int j=0;j<l;j++){
      real_t x=h[l].E[0]+h[j].E[0]; //case of nothing passing through li and ij
      for(int nu=1;nu<=M;nu++){
         double m=std::max( h[l].E[nu]+h[j].E[nu+M], h[l].E[nu+M]+h[j].E[nu] );    
         x=std::max(x,m);
      } // end cycle for over nu
      w[complete_g.edge(complete_g(l),complete_g(j))]=x;
    }// end cycle for over j
  }// end cycle over l
  
  lemon::MaxWeightedMatching<graph_t,weightmap_t> mm(complete_g,w);
  vector<real_t>savel(k),savej(k);
  vector<Mes> U(k,Mes(2*M+1,-inf1));
  for(int l=0;l<k;++l){
    TmpEraser<graph_t, weightmap_t,vector<real_t> > rl(complete_g, w, savel, l);
    mm.run();
    real_t const mw=mm.matchingWeight();
    U[l].E[0]=mw;
    for(int j=0;j<l;++j){
      TmpEraser<graph_t, weightmap_t,vector<real_t> > rj(complete_g, w, savej, j);
      mm.run();
      real_t const mw=mm.matchingWeight();
      for(int nu=1;nu<2*M+1;nu++){
         U[l].E[nu]=std::max(U[l].E[nu],-cost[l][nu] + mw + h[j].E[nu]);
         U[j].E[nu]=std::max(U[j].E[nu],-cost[j][nu] + mw + h[l].E[nu]);
      } // end cycle for over nu
    } // end cycle for over j
  } // end cycle for over l

  /*---------------  UPDATE OUT-MESSAGES OF NODE i -------------- */
  l=0;
  real_t eps=0.;
  static Mes old(depth);
  for(tie(eit,eend)=out_edges(i,g); eit!=eend;++eit, ++l){
    Mes & out =getMesInv(*eit);
    swap(old,out);
    for(int nu=0;nu<2*M+1;nu++){out.E[nu]=-U[l].E[nu];}
    out.reduce();
    eps = max(eps, l8dist(out, old));   // calculate error  
  }  // end cycle over eit

  /*---------------  REINFORCEMENT UPDATE -------------- */

  for(tie(eit,eend)=out_edges(i,g); eit!=eend;++eit){
    Mes const & in=getMes(*eit);
    Mes const & out=getMesInv(*eit);
    for(int nu=1;nu<=M;nu++){
      double Hij=out.E[nu]+in.E[nu+M]-(g[*eit].c+g[*eit].bias[nu]);
      double HijInv=out.E[nu+M]+in.E[nu]-(g[*eit].c+g[*eit].bias[nu+M]);
      vector<double>  & c_out = getCost(*eit);
      c_out[nu]=g[*eit].c+g[*eit].bias[nu]+rein*Hij;
      c_out[nu+M]=g[*eit].c+g[*eit].bias[nu+M]+rein*HijInv;
    }
    g[*eit].cij[0]=g[*eit].c+g[*eit].bias[0]+rein*(out.E[0]+in.E[0]-(g[*eit].c+g[*eit].bias[0]));
  } // end for out edges i

  return eps;

}

/* --------------- END UPDATE ROUTINE --------------------------------------------*/
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
double iterate()
{  
  int* permutation=new int[num_vertices(g)];
  double eps=0.; 

  for (unsigned int j = 0; j < num_vertices(g); ++j) permutation[j] = j;
  shuffle2(permutation, num_vertices(g));

  for (unsigned j = 0; j < num_vertices(g); ++j) eps = max(eps, updateMS(permutation[j]));

  delete[] permutation;

  return eps;
}


//-------------------------------------------------------------------------------
int idx(string const & id, Graph & g)
{
  static map <string, int> idx_map;
  map<string, int>::iterator mit = idx_map.find(id);
  if (mit == idx_map.end()) 
    return idx_map[id] = add_vertex(VertexProperties(id), g);
  return mit->second;
}

int idx1(string const & id, SimpleGraph & g1)
{
  static map <string, int> idx_map1;
  map<string, int>::iterator mit1 = idx_map1.find(id);
  if (mit1 == idx_map1.end()) 
    return idx_map1[id] =add_vertex(VertexProperties1(id), g1);
  return mit1->second;
}

//-------------------------------------------------------------------------------
Vertex vertex_g(Vertex1 v)
{
  int i=idx1(g1[v].name,g1);
  return vertex(i,g);
}

//-------------------------------------------------------------------------------
Vertex1 vertex_g1(Vertex v)
{
  int i=idx(g[v].name,g);
  return vertex(i,g1);
}


//-------------------------------------------------------------------------------

void read_graph(istream & file, Graph & a, SimpleGraph & b)
{  string tok;
  while (file >> tok){
    if (tok == "E") {
      string i, j;
      file >> i >> j ;
      Edge e = add_edge(vertex(idx(i,a), a), vertex(idx(j,a), a), a).first;
      Edge1 e1 = add_edge(vertex(idx1(i,b), b), vertex(idx1(j,b), b), b).first;

    }
  }

}

//-------------------------------------------------------------------------------

/*------------------------  ASSIGN PAIRS (S,R) RANDOMLY AND PERFORMS BFS TO CALCULATE SHORTEST PATHS  ---------------------------- */

int assign_couples( vector<int> & s_path, vector< pair<int,int> > & l_s_path, Graph & g, SimpleGraph & g1, int m, vector< pair<int,int> > & st)
{
  s_path.clear();
  l_s_path.clear();
  int length=0;
  int bfs_l=0;
  unsigned int size=num_vertices(g);
  int* node=new int[size];
  for (unsigned int i=0; i<size; i++) {node[i]=i;}
  shuffle(node,size);

  for (int j=0, i=0; j<2*m; i++,j+=2) {
    
    int s=node[j];
    int t=node[j+1];
    if (s!=t) {
      l_s_path.push_back(make_pair(s,t));
      g[s].G = i+1;
      g[s].c = 1; 
      g[t].G = i+1;
      g[t].c = -1;
      Vertex1 s1=vertex_g1(s);
      Vertex1 t1=vertex_g1(t);
          
      bfs_l=bfs_single(s1,t1,g1);
  
      // assign direct edge btw S,R
      Edge e=add_edge(vertex(s,g),vertex(t,g),g).first;
      
      if(bfs_l>1){
        g[e].c=inf2;
        vector<double> extra_cost(2*M+1,inf2);
        g[e].cij=extra_cost;g[e].cji=extra_cost;
        }  // if (S,R) edge already exists does not assign infty cost to it
     
      length+=bfs_l;
      s_path.push_back(bfs_l);
      st.push_back(make_pair(s,t)); // build up a vector so to insert in g1 s,r edges after having calculated the single paths
    }      //  end if(s!=t)
    
  }
  assert(l_s_path.size()==(unsigned int) m);
  cout << " Single path total length="<<length<<endl;
  
  delete[] node;
  
  return length;
}

/*    ----------------     ASSIGNS PAIRS (S,R) FROM FILE --> ONLY IF flag_SR==1    ----------------------------------------*/
int assign_couples_from_file( vector<int> & s_path, vector< pair<int,int> > & l_s_path, Graph & g, SimpleGraph & g1, int m, vector< pair<int,int> > & st,istream & in)
{
  s_path.clear();
  l_s_path.clear();
  int length=0;
  int bfs_l=0;
  unsigned int size=2*m;
  int* node=new int[size];
  int h=0;
  while(in.good() && h<2*m){
    string ids,idt;
    in>>ids>>idt;
    int s=vertex(idx(ids,g),g);
    int t=vertex(idx(idt,g),g);
    node[h]=s;node[h+1]=t;h+=2;
  }

  for (int j=0, i=0; j<2*m; i++,j+=2) {
    int s=node[j];
    int t=node[j+1];
    if (s!=t) {
      l_s_path.push_back(make_pair(s,t));
      g[s].G = i+1;
      g[s].c = 1; 
      g[t].G = i+1;
      g[t].c = -1;
      Vertex1 s1=vertex_g1(s);
      Vertex1 t1=vertex_g1(t);
          
      bfs_l=bfs_single(s1,t1,g1);
      
      // assign direct edge btw S,R
      Edge e=add_edge(vertex(s,g),vertex(t,g),g).first;
      if(bfs_l>1){
        g[e].c=inf2;
        vector<double> extra_cost(2*M+1,inf2);
        g[e].cij=extra_cost;g[e].cji=extra_cost;
        }  // if (S,R) edge already exists does not assign infty cost to it
            // g[e].x=1;  // sign it so that you'l not updat it
      
      length+=bfs_l;
      s_path.push_back(bfs_l);
      st.push_back(make_pair(s,t)); // build up a vector so to insert in g1 s,r edges after having calculated the single paths
    }
    
  }
  assert(l_s_path.size()==(unsigned int) m);
  cout << " Single path total length="<<length<<endl;
  
  delete[] node;
  
  return length;
}
//-------------------------------------------------------------------------------


//----------- Mathematica-style--------------------------------------------------------------------
void out_graph(ofstream & file, SimpleGraph gr)
{
  graph_edge_iterator1 eit, eend;
  for (tie(eit, eend) = edges(gr); eit != eend; ++eit) {
    Edge e=*eit;
    Vertex1 i = target(e, gr); 
    Vertex1 j = source(e, gr);
    file<<"{ ";  
    file << gr[i].name << " -> " << gr[j].name <<" , 0 "<<" } "<<","<< endl;
  } 
}
/*--------------------------------------------------------------------
  -------------------------------------------------------------------*/

//---------------------E LINK---------------------------------------------------------------------------//

double marginal_Energy_link(Edge eit) 
{  
  double mintot=inf;

  Mes const & out = getMesInv(eit);
  Mes const & in = getMes(eit);
  double mintot0 =out.E[0]+in.E[0] ;
  for(int mu=1;mu<=M;++mu){ 
    double m=std::min(out.E[mu]+in.E[mu+M]-g[eit].c,out.E[mu+M]+in.E[mu]-g[eit].c);
    if (m<= mintot) { mintot =m;}
  } 
  
  if(mintot<mintot0){
    Elink+=mintot/(double)M;  
    return mintot;
  }
  else  {
    Elink+=mintot0/(double)M;  
    return mintot0;
  }
  
} 

//---------------------E NODE ---------------------------------------------------------------------------//

double marginal_Energy_node( Vertex v) 
{  
  return E_n[v];
} 

//---------------------E TOTAL ---------------------------------------------------------------------------//

pair<double,double> energy(Graph g){
  double e_link=0.,e_node=0.;
  vertex_iterator vit, vend;
  graph_edge_iterator geit, geend;
  for (tie(vit, vend) = vertices(g); vit != vend; ++vit) {e_node+=marginal_Energy_node(*vit)/(double)num_vertices(g);}
  for (tie(geit, geend) = edges(g); geit != geend; ++geit) {e_link+=marginal_Energy_link(*geit)/(double)num_edges(g);}
  return make_pair(e_node,e_link);
}
//-------------------------------------------------------------------------------

//-----------MARGINAL EDGE---------------------------------------------------------

int marginal(Edge eit) 
{  
  int nu=0;
  
  Mes const & out = getMesInv(eit);
  Mes const & in = getMes(eit);
  double mintot =out.E[0]+in.E[0] ; assert((int)mintot==0);
  for(int mu=1;mu<=M;++mu){ 
    double m=std::min(out.E[mu+M]+in.E[mu]-g[eit].c,out.E[mu]+in.E[mu+M]-g[eit].c);
    if (m < mintot) {
      mintot = m;
      nu=mu;
    }
  }
  return nu;
} 

/* --------------- L_tot --------------------------------------------*/

// calculate total length inside the convergence, does not keep track of the single length of the paths, only the total one.
// Uses marginal()
int L_tot()
{
  int total_delta=0;   // difference in each single edge current w.r.t. the previous step
  //  int l=0;
  graph_edge_iterator geit, geend;
  for (tie(geit, geend) = edges(g); geit != geend; ++geit) {    
    int d=0;
    int nu=0;   // colour of the path passing through eit
    nu=marginal(*geit);
    d=(g[*geit].x-nu==0)? 0:1;  // count 1 if the currents at that edge are different at t and t+1    
    g[*geit].x=nu;    // update current
    total_delta+=d;
//    if (nu!=0) {l++;}
  }
  return total_delta;
}

//------------------COVERGENCE-------------------------------------------------------------

pair<int,double> converge( ofstream & out)
{
  int Tmax=(beta==0)? maxit : (int ( 1./beta))*10;
  int it = 0;
  double err;
  int coincide=0, coincide2=0;
  //  int length=0,oldlength=0;
  int delta=0;
  cout<<"T max="<<Tmax<<"  beta="<<beta<<endl;
  do {
    rein=beta*it;
    err = iterate();

    //    oldlength=length; // length at previous step
    delta=L_tot();   // calculate total length difference edge by edge 
    out<<it<<" "<<delta<<" "<<err<<" "<<V<<" "<<M<<" "<<maxit<<" "<<tolerance<< endl;
    // cout<<it<<" "<<delta<<" "<<err<<" "<<V<<" "<<M<<" "<<maxit<<" "<<tolerance<< endl;

    if(it>T_eq) { if (/*err<tolerance*/(delta==0)/*&&(length>=single_tot_length)*/) {
      coincide++;
      } else if(err<tolerance){coincide2++;} else {coincide = 0;coincide2=0;}
    }// end if it>t_eq
  } while (coincide < decision && ++it < Tmax && coincide2< decision && err< err_max);
  

  return make_pair(it,delta);
}

//---------------EDGE PATH-----------------------------------------------------

//calculate the path from source to receiver using E_edge. Output= Total length and failure=T/F
pair<int,bool> edge_path(int & violations, vector <int> single_path, vector< pair<int,int> > l_s_path)
{
  bool failure=false;
  int tot_lenght=0;
  vector< int > h_lenghts;  // hipotetical length of the paths (discarded if reailzation fails)
  vector< int > h_delta_l;
  int length_s;
  
  for (int m=0; m<M; m++) {
    Vertex s=l_s_path[m].first;
    Vertex r=l_s_path[m].second ;

    Vertex t=s;
    Vertex previous=s;
    int length=0;
    int mu=g[s].G;
    while(t!=r){
      if (failure==true) {break;}
      int no_path=0;   // check if no edges find a path
      edge_iterator eit, eend;
      for (tie(eit,eend)=out_edges(t,g); eit!=eend; ++eit) {
   // check to not pick the already picked edge
   if ( (target(*eit,g)!=previous && source(*eit,g)!=previous) || (previous==s && t==s) ){
     Edge e=edge(source(*eit,g),target(*eit,g),g).first;
     int out_m=0;  // checks if there are more coms exiting t
     int a = marginal(*eit);  
     if (a==mu) {
       out6<<"{ ";
       out6 << g[t].name << " -> " ;
       no_path++;
       length++;
       out_m++;  
       previous=t;
       t=(target(*eit,g)==t) ? source(*eit,g) : target(*eit,g);
     
       g[e].q++;  // colour edge
      
       if (g[e].q>1) {failure=true;violations++; cout<<" g[*eit].q="<<g[e].q<<" "<<V<<" "<<M<<endl;}
       out6 << g[t].name << " ,  " << mu <<" } "<<" , "<< endl;
       break;
     }
    
   } // end *eit!=e_prev

   if (t==r) {
     failure=false;break;
   }
      } // end cycle over out_edges(target)
      if(t==s && no_path==0){
   t=r;length=1;
       out6<<"{ ";
       out6 << g[s].name << " -> " ;
       out6 << g[r].name << " ,  " << mu <<" } "<<" , "<< endl;
      }// add a direct path S,R      
      if(t!=s && no_path==0) {failure=true; break;}
    }  // end while
    if (t!=r) {failure=true;out6<<"mu= "<<mu<<" failed"<<endl; violations++;break;}
    else {
      tot_lenght+=length;
      length_s=single_path[m];
      h_lenghts.push_back(length); 
      if(length==1 && length_s>1){direct_path++;}
      if(length>1){      int delta_l=length-length_s;
   h_delta_l.push_back(delta_l);}
    }

  } // end cycle over m


  if (failure==false) {
    //    assert(h_lenghts.size()==h_delta_l.size());
    for (unsigned int i=0; i< h_lenghts.size(); i++) {
      add_histo(h_lenghts[i],histo);Z++;
      if(h_delta_l.size()!=0) {add_histo(h_delta_l[i], histo_d);Z_d++;}
    }
    h_lenghts.clear();
    h_delta_l.clear();
    return make_pair(tot_lenght,failure);
  }
  else {return make_pair(0,failure);}
}

//---------OUTPUT TREE -----------------------------------------------------------
// Consider that inside edge_path() one already outputs the paths ... maybe this other function is redundant
void output_tree(ofstream & out)
{
  graph_edge_iterator eit, eend;
  for (tie(eit, eend) = edges(g); eit != eend; ++eit) {
    int a = marginal(*eit);  
    if (a>=0) {
      Vertex i = target(*eit, g); 
      Vertex j = source(*eit, g);  
      out<<"{ ";
      out << g[i].name << " -> " << g[j].name << " ,  " << a <<" } "<<" , "<< endl;
            
    } 
  }
} 


/*---------------------AVERAGE DEGREE-----------------------------------------------*/

pair<double, int> average_k(SimpleGraph & g1, vector< pair<int,int> > & histo_k, int & Z)
{ 
  int ave_degree=0.;
  int num_vert=0;
  pair<double, int> result;
  vertex_iterator1 vit,vend;
  for (tie(vit,vend)=vertices(g1); vit!=vend; vit++){
      int k=0;
      edge_iterator1 eit,eend;
      for (tie(eit,eend)=out_edges(*vit,g1); eit!=eend; eit++) {
   Edge e=edge(*vit, target(*eit,g1),g).first;
   if(g[e].q==0){ave_degree++;k++; }
}  // end cycle over neighbours of *vit
      if(k>0) {num_vert++;g1[*vit].prm=k;}else{g1[*vit].prm=0;} //denotes the degree of node vit with .prm
      add_histo(k, histo_k);Z++;
  
  } // end for vertices
  if(num_vert==0)cout<<"V'= "<<num_vert<<endl;
  result= make_pair( ((double)ave_degree/(double)num_vert), num_vert );
  return result;
}
 
/*---------------------GRAPH CONNECTIVITY-----------------------------------------------*/

pair<int, bool> connectivity(SimpleGraph & g1, int num_vert)
{ 
  int max_size=0;  // stores the value of the max cluster size 
  int n_of_search=20; //number of times we search for a giant component
  bool connection=false;
  int coloured_vertices; // monitors if we picked v0 in the max cluster or not
  int j=0;
  //pick a starting node for the search
  do{
    Vertex1 v0=roll_die2(0,num_vertices(g1)-1);
    while (g1[v0].prm==0) {v0=roll_die2(0,num_vertices(g1)-1);}
    // fix all colour to 0
    vertex_iterator1 vit,vend;
    for (tie(vit,vend)=vertices(g1); vit!=vend; vit++) {
      g1[*vit].c=0;/*g1[*vit].d=0;g1[*vit].prm=0;*/
    }
    coloured_vertices=1;  // counts how many vertices we visit
    vector<Vertex1> queue;
    vector< pair<Vertex1,Vertex1> > parent;
    queue.push_back(v0);
    parent.push_back(make_pair(v0,v0));
    g1[v0].c=1;
    
    while (queue.size()!=0 ) {
      Vertex1 p=queue[0];
      edge_iterator1 eit,eend;
      for (tie(eit,eend)=out_edges(p,g1); eit!=eend; eit++) {
   Edge e=edge(source(*eit,g1), target(*eit,g1),g).first;
   if(g[e].q==0){   // if the edge is empty  
     Vertex1 u=(target(*eit,g1)==p) ? source(*eit,g1) : target(*eit,g1);
     if (g1[u].prm>1) { //if u has at least one ohter free edge
       if (g1[u].c==0) {
         g1[u].c=1;
         coloured_vertices++;
         queue.push_back(u);
         parent.push_back(make_pair(u,p));
       }  // end if u not coloured yet
     } //end if u has at least 2 free edges
   } //end if eit.q==0
      }  // end  cycle for over edges eit
      //FIFO
      queue.erase(queue.begin());
      }  // end while queue not empty
    j++; max_size=std::max(coloured_vertices,max_size);
    if (coloured_vertices==num_vert) {connection=true;}
  }while( (max_size<(int)((double)V*0.5)) && (j<n_of_search) );
  return make_pair(max_size,connection);
}

/* -------------------------------------------------------------------*/


namespace po = boost::program_options;

po::variables_map parse_command_line(int ac, char ** av)
{
  po::options_description desc("Usage: " + string(av[0]) + " <option> ... \n\twhere <option> is one or more of");
  desc.add_options()
    ("help", "produce help message")
    ("folder,f", po::value(&folder)->default_value("blrand"), "set dataset folder name")
    ("input_folder,S", po::value(&input_folder)->default_value("/input/"), "set name of input folder")
    ("output_folder,E", po::value(&output_folder)->default_value("/output/"), "set output folder name")
    ("instance,h", po::value(&instance)->default_value("50ins1.dat"), "file specific instance (S,R)")
    ("com,c", po::value(&com)->default_value("1"), "set communication number, beginning of filename")
    ("flag_SR,d", po::value(&flag_SR)->default_value(0), "flag for fixing or not specific instance (S,R): if flag_SR==1 then input pairs from file")
    ("M,m", po::value(&M)->default_value(1), "set number of communications")
    ("N_real,N", po::value(&N_real)->default_value(1), "set number of realizations")
    ("degree,k", po::value(&k)->default_value(3), "set average degree")
    ("maxit,t", po::value(&maxit)->default_value(500), "set maximum number of iterations")
    ("tolerance,e", po::value(&tolerance)->default_value(0.01), "set convergence tolerance")
    ("bias_intensity,b", po::value<double>(&bias_intensity)->default_value(0.), "set color bias")
    ("seed,s", po::value<unsigned>(), "sets instance seed")
    ("rseed,z", po::value<unsigned>(), "sets biases seed")
    ("messages,M", po::value(&flag_output_messages)->default_value(0) ,"if 1 then output messages on convergence")
    ("beta,g", po::value<double>(&beta)->default_value(0.000), "sets reinforcement parameter beta")
    ("decision,y", po::value<int>(&decision)->default_value(20), "program converges after this # repeats of the decision variables");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("seed")) {
    unsigned s = vm["seed"].as<unsigned>();
    gen_int.seed(s);
  }
  if (vm.count("rseed")) {
    unsigned s = vm["rseed"].as<unsigned>();
    gen.seed(s);mes_gen.seed(s);gen_int2.seed(s);
  }
  
  if (vm.count("help")) {
    cerr << desc << "\n";
    exit(1);
  }

  return vm;
}

/* ---------------- MAIN ---------------------------------------------*/
//--------------------------------------------------------------------

int main(int ac, char** av)
{
  
  cout.setf(ios_base::fixed, ios_base::floatfield);
  po::variables_map vm = parse_command_line(ac, av);
  
  // FILE
  string file="../data/"+folder+output_folder+com;

  // Variables
  SimpleGraph g0_1;
  Graph g0;
  
  vector <int> single_path;  // list of the lenghts of the paths in the single case
  vector< pair<int,int> > l_s_path;  //list of the s,r for all the communications

  clock_t t=0.;//time to process one m (N_realizations)
  int fail_conv=0;// number of time convergence fails
  int fail_constr_edge=0;     // number of times (over the converged ones) that the overlap constraint is violated calculating E_edge
  
  // COMMUNICATIONS
  depth=2*M+1;
  double fm=(double)M/(double)V;
  
  /*-------------  FILE I/O ------------------------------*/

  //in1.open(("/Users/Caterina/Dropbox/Lavoro/Noverlap/Data/"+folder+"/k3/Single/adjacency.dat").c_str());
  in1.open(("../data/"+folder+"/input/adjacency.dat").c_str());

  out1.open((file+"result.dat").c_str(), ios::app);
  out2.open((file+"f_histo.dat").c_str(), ios::app);
  out3.open((file+"f_histo_d.dat").c_str(), ios::app);
  out4.open((file+"histo_k.dat").c_str(), ios::app);
  out6.open((file+"tree.dat").c_str());
  out8.open((file+"energy.dat").c_str(),ios::app);

  out10.open((file+"error.dat").c_str());
  out11.open((file+"SRpairs.dat").c_str());
  out12.open((file+"graph_mathe.dat").c_str()); 
  
  read_graph(in1,g0,g0_1);
  V=(int)num_vertices(g0);
  int E=(int)num_edges(g0);
  out12<<"#  V="<<V<<" E="<<E<<" <k>="<<((double)2*E)/(double)V<<endl;
  out_graph(out12,g0_1);

  /*--------------------------------------------------*/
  
  for (int r=0; r<N_real; r++) {  // start cycle over realizations
    
    t=clock();  //CHRONO starts
    
    int tot_edge_path=0;// length calculated by EDGE of the total length for the single realization
    int tot_single_path=0;// length calculated by greedy of the SINGLE length for the single realization
    int time=0;// convergence time for the real considered
    double err=0.;
    bool connection=false;// flags that stops the iteration if constraints are violated (for the single real)
    int loop_stuck=0;// # of paths converged but that find a long length due to loops 
    bool edge_overlap=false;
    int N_vertices=0;// number of vertices after removal
    double degree=0.;           // average degree after removal
    double cluster_size=0.;
    vector< pair<int,int> > st;  //list of s,t for g1 so to add direct links
    
    direct_path=0;            //  counts # of paths of length=1 connecting S and R
    Elink=0.;
    E_n.resize(V);

    g.clear(); g1.clear();
    g=g0;g1=g0_1;
  
    if(flag_SR==1)in2.open(("../data/"+folder+"/input/"+instance).c_str());    
    //pick randomly M coms and store the average path length single GREEDY
    if(flag_SR==0)tot_single_path=assign_couples(single_path,l_s_path,g, g1,M, st);
    else tot_single_path=assign_couples_from_file(single_path,l_s_path,g, g1,M, st,in2); 

    g1.clear();g1=g0_1;
    out11<<" N_real= "<<r<<"  ---------------------"<<endl;
    for(unsigned int i=0, n=st.size(); i<n; i++){
      int s=st[i].first;
      int t=st[i].second;
      Edge e1=add_edge(vertex(s,g1), vertex(t,g1),g1).first;
      out11<<g[s].name<<" "<<g[t].name<<" "<<i+1<<endl;
    }
    out11<<"     ---------------------"<<endl;
    out10<<"# N_real="<<r<<"  Single path="<< tot_single_path<<endl;
    out6<<"# N_real="<<r<<"  Single path="<< tot_single_path<<endl;
    // if(r<10) out7<<"# N_real="<<r<<"  Single path="<< tot_single_path<<endl;
    tie(time,err)=converge(out10);
   
    // if convergence fails:
    if (time==maxit){fail_conv=1;}else{fail_conv=0;}

    //CALCULATE ENERGY_LINK/NODE
    tie(Enode,Elink)=energy(g);

    // CALCULATE PATHS BY EGDE
    tie(tot_edge_path,edge_overlap)=edge_path(loop_stuck, single_path, l_s_path);
    // IF PATHS EXIST
    if(edge_overlap==false){
      cout << "Total not overlapping path= "<<tot_edge_path<<endl;
      cout<< "# of direct paths="<<direct_path<<endl;
      tie(degree,N_vertices)= average_k(g1,histo_k,Z_k);
      cout<<"<k>= "<<degree<<endl;      
      tie(cluster_size,connection)= connectivity(g1,N_vertices);
      cluster_size/=(double)N_vertices;
      fail_constr_edge=0;
      cout<<"cluster size= "<<cluster_size<<endl;
    }
    else{fail_constr_edge=1;}
    out6<<" -----  failure= "<<fail_constr_edge<<endl;
    t = clock() - t;
    
    out1<<fm<<" "<<tot_single_path<<" "<<tot_edge_path<<" ";
    out1<<tot_single_path/(double)M<<" "<<tot_edge_path/(double)M<<" ";
    out1<<loop_stuck<<" "<<loop_stuck/(double)M<<" ";
    out1<<time<<" ";
    out1<<fail_conv<<" "<<fail_constr_edge<<" ";
    out1<<degree<<" "<<N_vertices<<" ";
    out1<<cluster_size<<" "<<connection<<" ";
    out1<< V<<" "<<M<<" "<<((float)t)/CLOCKS_PER_SEC<<" ";
    out1<<tolerance<<" "<<decision<<" "<<maxit<<" ";
    out1<<err<<" "<<r<<" "<<seed<<" "<<rseed<<" "<<direct_path<<" "<<((double)direct_path/(double)M)<<" ";
    out1<<beta<<endl;
    
    out8<<fail_constr_edge<<" "<<endl;
    out8<<fm<<" "<<(Enode-(double)k*0.5*Elink)/fm<<" "<<Enode/fm<<" "<<Elink/fm<<" ";
    out8<<Enode-(double)k*0.5*Elink<<" "<<Enode<<" "<<Elink<<" ";
    out8<<tot_edge_path/(double)M<<" "<<V<<" "<<M<<" "<<r<<" "<<seed<<" "<<rseed<<endl;

    // ------------------HISTOGRAM PLOT -----------------------------------
    
    out2<<"#  M="<<M<<"  fm="<<fm<<"  ------------------------------"<<endl;
    sort(histo.begin(), histo.end(), sort_pred());
    for(unsigned int i=0;i<histo.size();i++) out2 << histo[i].first << " " << (double) histo[i].second/(double)Z << endl;
    out2<<endl<<"#-------------------------------------------------------------"<<endl<<endl;
     
    out3<<"#  M="<<M<<"  fm="<<fm<<"  ------------------------------"<<endl;
    sort(histo_d.begin(), histo_d.end(), sort_pred());
    for(unsigned int i=0;i<histo_d.size();i++)out3 << histo_d[i].first << " " << (double) histo_d[i].second/(double)Z_d << endl;
    out3<<endl<<"#-------------------------------------------------------------"<<endl<<endl; 
    
    sort(histo_k.begin(), histo_k.end(), sort_pred());
    out4<<fm<<" ";
    for(unsigned int i=0;i<histo_k.size();i++){out4 <<histo_k[i].second/(double)Z_k <<" ";out4 << histo_k[i].first;}
    out4<<degree<<" "<<N_vertices<<" ";
    out4<<cluster_size<<" ";
    out4<<connection<<" "<<V<<" "<<M<<endl;
    
    if(flag_output_messages==1 && r==0){
      out9.open((file+"messages.dat").c_str());
      output_messages(out9);
      out9.close();
    }
    // if(r<10){ output_tree(out7);out7<<" -----  failure= "<<fail_constr_edge<<endl;    }

    if(flag_SR==1){
      seed++;rseed++; 
      gen.seed(rseed);mes_gen.seed(rseed);gen_int2.seed(rseed);gen_int.seed(rseed);
      in2.close();
    }
  }  // end cycle over realizations
  
  //--------------------------------------------------------------------------------
  
  out1.close();
  in1.close();

  out2.close();
  out3.close();
  out4.close();
  out6.close();
  out8.close();
  out10.close();
  out11.close();
  
  return 1;
}


