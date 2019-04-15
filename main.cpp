//
//  main.cpp
//  NetworkSimplexMCF
//
//  Created by Rui Shibasaki on 31/05/16.
//  Copyright (c) 2016 Rui Shibasaki. All rights reserved.
//

#include "networksimplex.hpp"



int main(int argc, const char * argv[]) {
    
    //grafo de 6 nós
  /*  Graph graph;
    graph.nNodes = 6;
    graph.arcs.resize(graph.nNodes+1);
    graph.capacities.resize(graph.nNodes+1);
    graph.b.resize(graph.nNodes+1, 0);
    graph.costs.resize(graph.nNodes+1);
    for (int i=0; i<=graph.nNodes; ++i) {
        graph.arcs[i].resize(graph.nNodes+1);
        graph.capacities[i].resize(graph.nNodes+1);
        graph.costs[i].resize(graph.nNodes+1);
        for (int j=0; j<=graph.nNodes; ++j) {
            graph.arcs[i][j] = -1;
        }
    }
    
    graph.nArcs = 9;
    graph.arcs[1][2] = 1;
    graph.capacities[1][2] = 8;
    graph.costs[1][2] = 3;
    graph.arcs[1][3] = 1;
    graph.capacities[1][3] = 3;
    graph.costs[1][3] = 2;
    graph.arcs[2][3] = 1;
    graph.capacities[2][3] = 3;
    graph.costs[2][3] = 2;
    graph.arcs[2][4] = 1;
    graph.capacities[2][4] = 7;
    graph.costs[2][4] = 5;
    graph.arcs[3][5] = 1;
    graph.capacities[3][5] = 3;
    graph.costs[3][5] = 4;
    graph.arcs[2][5] = 1;
    graph.capacities[2][5] = 2;
    graph.costs[2][5] = 2;
    graph.arcs[5][4] = 1;
    graph.capacities[5][4] = 4;
    graph.costs[5][4] = 5;
    graph.arcs[4][6] = 1;
    graph.capacities[4][6] = 5;
    graph.costs[4][6] = 3;
    graph.arcs[5][6] = 1;
    graph.capacities[5][6] = 6;
    graph.costs[5][6] = 4;
    
   // b > 0  oferta b < 0  demanda

    graph.b[1] = 5;
    graph.b[2] = -1;
    graph.b[3] = -1;
    graph.b[4] = -1;
    graph.b[5] = -1;
    graph.b[6] = -1;
*/
 //GRAFO DE 10 NÓS
/*    Graph graph;
    graph.nNodes = 10;
    graph.arcs.resize(graph.nNodes+1);
    graph.capacities.resize(graph.nNodes+1);
    graph.b.resize(graph.nNodes+1, 0);
    graph.costs.resize(graph.nNodes+1);
    for (int i=0; i<=graph.nNodes; ++i) {
        graph.arcs[i].resize(graph.nNodes+1);
        graph.capacities[i].resize(graph.nNodes+1);
        graph.costs[i].resize(graph.nNodes+1);
        for (int j=0; j<=graph.nNodes; ++j) {
            graph.arcs[i][j] = -1;
        }
    }

  graph.nArcs = 21;
  graph.arcs[1][2] = 1;
  graph.capacities[1][2] = 7;
  graph.costs[1][2] = 2;
  graph.arcs[1][3] = 1;
  graph.capacities[1][3] = 5;
  graph.costs[1][3] = 3;
  graph.arcs[2][5] = 1;
  graph.capacities[2][5] = 3;
  graph.costs[2][5] = 2;
  graph.arcs[2][4] = 1;
  graph.capacities[2][4] = 6;
  graph.costs[2][4] = 5;
  graph.arcs[3][5] = 1;
  graph.capacities[3][5] = 7;
  graph.costs[3][5] = 5;
  graph.arcs[3][4] = 1;
  graph.capacities[3][4] = 4;
  graph.costs[3][4] = 2;
  graph.arcs[4][6] = 1;
  graph.capacities[4][6] = 5;
  graph.costs[4][6] = 4;
  graph.arcs[5][7] = 1;
  graph.capacities[5][7] = 3;
  graph.costs[5][7] = 1;
  graph.arcs[6][8] = 1;
  graph.capacities[6][8] = 4;
  graph.costs[6][8] = 3;
  graph.arcs[6][9] = 1;
  graph.capacities[6][9] = 8;
  graph.costs[6][9] = 6;
  graph.arcs[7][6] = 1;
  graph.capacities[7][6] = 3;
  graph.costs[7][6] = 1;
  graph.arcs[7][9] = 1;
  graph.capacities[7][9] = 3;
  graph.costs[7][9] = 4;
  graph.arcs[8][10] = 1;
  graph.capacities[8][10] = 5;
  graph.costs[8][10] = 5;
  graph.arcs[9][8] = 1;
  graph.capacities[9][8] = 2;
  graph.costs[9][8] = 6;
  graph.arcs[9][10] = 1;
  graph.capacities[9][10] = 6;
  graph.costs[9][10] = 4;
  graph.arcs[8][4] = 1;
  graph.capacities[8][4] = 5;
  graph.costs[8][4] = 4;
  graph.arcs[9][5] = 1;
  graph.capacities[9][5] = 2;
  graph.costs[9][5] = 3;
  graph.arcs[6][3] = 1;
  graph.capacities[6][3] = 7;
  graph.costs[6][3] = 5;
  graph.arcs[2][1] = 1;
  graph.capacities[2][1] = 6;
  graph.costs[2][1] = 5;
  graph.arcs[10][1] = 1;
  graph.capacities[10][1] = 3;
  graph.costs[10][1] = 4;
  graph.arcs[5][3] = 1;
  graph.capacities[5][3] = 6;
  graph.costs[5][3] = 5;
  */
    
  // b > 0  oferta b < 0  demanda
/*  graph.b[1] = 0;
  graph.b[2] = 0;
  graph.b[3] = 0;
  graph.b[4] = -2;
  graph.b[5] = 0;
  graph.b[6] = 0;
  graph.b[7] = 0;
  graph.b[8] = 0;
  graph.b[9] = 2;
  graph.b[10] = 0;
  */
    
    //escolher max de iterações principais e max de iterações menores
    int nmajor = 10;
    int nminor = 10;
    double mincost;
    spanningTree bSTTsol;
    
    //PARA O CASO DE LEITURA DO GRAFO A PARTIR DE UM ARQUIVO
    //DESCOMENTAR AS LINHAS ABAIXO E COMENTAR A DEFINIÇÃO DO GRAFO ACIMA
    //compilar: g++ main.cpp networksimplex.cpp -o NetworkSimplex
    //linha de commando:  NetworkSimplex   instancia    max_main_iter    max_minor_iter
    //    Graph graph;
    //readInstance(graph, argv[1]);
  //  nmajor = atof(argv[2]);
   // nminor = atof(argv[3]);
    

    
    mincost = NetworkSimplex(bSTTsol, graph, nmajor, nminor);
    

    
    if (mincost >= 0) {
        std::cout<<"optimal cost: "<<mincost<<std::endl;
    
        for (int i=0; i<=graph.nNodes; ++i) {
            for (int j=0; j<=graph.nNodes; ++j) {
                if (graph.arcs[i][j]>0)
                        std::cout<<"X ("<<i<<","<<j<<") = "<<bSTTsol.x[i][j]<<std::endl;
            }
        }
    }
    else std::cout<<"Infeasible !! "<<mincost<<std::endl;
    
    return 0;
}




