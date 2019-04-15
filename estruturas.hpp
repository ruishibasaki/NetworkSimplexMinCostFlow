//
//  estruturas.hpp
//  NetworkSimplexMCF
//
//  Created by Rui Shibasaki on 17/06/16.
//  Copyright (c) 2016 Rui Shibasaki. All rights reserved.
//

#ifndef NetworkSimplexMCF_estruturas_hpp
#define NetworkSimplexMCF_estruturas_hpp


struct Graph{
    
    int nArcs, nNodes;
    //-1 if does not exists; 0 added arcs; 1 original arcs
    std::vector< std::vector<int> >arcs;
    std::vector< std::vector<double> >capacities;
    std::vector< std::vector<double> > costs;
    std::vector< double > b;
    std::vector< std::vector<int> > neighbors;
    
    Graph & operator = (const Graph & g){
        this->nArcs = g.nArcs;
        this->nNodes = g.nNodes;
        this->capacities = g.capacities;
        this->costs = g.costs;
        for(int i= 0; i<g.nArcs; ++i){
            this->arcs[i] = g.arcs[i];
        }
        for(int i= 0; i<g.nNodes; ++i){
            this->neighbors[i] = g.neighbors[i];
        }
        return *this;
    }
};

struct Arc{
    int i, j;
    
    Arc(){};
    Arc(int i, int j){
        this->i = i;
        this->j = j;
    }
    
    Arc & operator=(Arc a){
        
        std::swap(i, a.i);
        std::swap(j, a.j);
        
        return *this;
    }
    
};

#endif
