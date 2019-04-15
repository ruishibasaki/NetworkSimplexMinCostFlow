//
//  networksimplex.h
//  NetworkSimplexMCF
//
//  Created by Rui Shibasaki on 31/05/16.
//  Copyright (c) 2016 Rui Shibasaki. All rights reserved.
//

#ifndef __NetworkSimplexMCF__networksimplex__
#define __NetworkSimplexMCF__networksimplex__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "estruturas.hpp"

bool operator ==(const Arc & a, const Arc & b);


struct spanningTree{
    
    int nNodes, nArcs, root;
    std::vector< std::vector<int> > T;
    std::vector< std::vector<int> > L;
    std::vector< std::vector<int> > U;
    std::vector< std::vector<double> > x;
    std::vector<double> phi;  //tamanho nNodes
    std::vector<int> pred;  //tamanho nNodes
    std::vector<int> depth; //tamanho nNodes
    std::vector<int> thread;    //tamanho nNodes
    std::vector<Arc> elegibles;
    Arc lastChecked;
    
    bool existsIn(int, int, const std::vector< std::vector<int> >  &);

    //organiza a estrutura essa estrutura (spanningTree)
    void make_tree_struct();
    
    //usada dentro do make_tree_struct() para estabelcer o thread(i);
    void depth_fisrt_search();
    
   //calcule custos reduzidos e pi.
    void computePhi(const Graph &);
    
    //compute reduced costs
    double computeRC(const Graph &, int , int );
    
    //compute flows
    void computeFlow(const Graph &);
    
    //inicializa adicionando nó 0 e arcos
    void initialize(Graph &);
    
    //check optimality; if not optimal add eligible arcs until a
    //number int max of elegibles arcs is reached
    bool checkOptimality(const Graph &, int);
    
    Arc enteringArc(const Graph &);
    
    Arc leavingArc(const Graph &, const Arc &);
    
    double computeDelta(const Graph & , const Arc & , int, int );
    
};

void updateTree(spanningTree &, const Graph &, Arc , Arc );

//Main algorithm
double NetworkSimplex(spanningTree &, Graph &, int, int);

void readInstance(Graph &, std::string);


#endif /* defined(__NetworkSimplexMCF__networksimplex__) */
