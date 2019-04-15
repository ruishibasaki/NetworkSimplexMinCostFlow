//
//  networksimplex.cpp
//  NetworkSimplexMCF
//
//  Created by Rui Shibasaki on 31/05/16.
//  Copyright (c) 2016 Rui Shibasaki. All rights reserved.
//

#include "networksimplex.hpp"

void spanningTree::initialize(Graph & graph){
    
    
    nNodes = graph.nNodes;
    nArcs = graph.nArcs;
    phi.resize(nNodes+1);
    x.resize(nNodes+1);
    T.resize(nNodes+1);
    U.resize(nNodes+1);
    L.resize(nNodes+1);
    for (int i=0; i<=nNodes; ++i){
        x[i].resize(nNodes+1,0);
        T[i].resize(nNodes+1,0);
        U[i].resize(nNodes+1,0);
        L[i].resize(nNodes+1,0);
        
    }
    
    
    for (int i=1; i<=nNodes; ++i) {
        if (graph.b[i] >= 0) {
            T[i][0] = 1;
            if (graph.arcs[i][0]<0){
                graph.costs[i][0] = 1e10;
                graph.arcs[i][0] = 0;
                graph.capacities[i][0] = graph.b[i];
                graph.nArcs++;
            }
        }
        else{
            T[0][i] = 1;
            if (graph.arcs[0][i]<0){
                graph.costs[0][i] = 1e10;
                graph.arcs[0][i] = 0;
                graph.capacities[0][i] = -graph.b[i];
                graph.nArcs++;
            }
        }
    }
    
    
    for (int i=0; i<=nNodes; ++i) {
        for(int j=0; j<=nNodes; ++j){
            if (graph.arcs[i][j]>0)
                     L[i][j] = 1;
        }
    }
    
    make_tree_struct();
    
    
    computeFlow(graph);
    computePhi(graph);
    lastChecked.i = 0;
    lastChecked.j = 0;
  
}


void spanningTree::make_tree_struct() {

    root = 0;
    thread.assign(nNodes+1,0);
    depth.assign(nNodes+1,0);
    pred.assign(nNodes+1, -1);
    std::vector<bool> marked(nNodes+1, false);
    std::vector<int>recover;
    
    
    depth_fisrt_search();
    
    depth[root]= 0;
    for (int p=1; p<=nNodes; ++p){
        int cont= 0;
        int i = p;
        while ( pred[i] > (-1) ) {
            ++cont;
            i = pred[i];
        }
            
        depth[p]=cont;
        
    }
    
}

void spanningTree::depth_fisrt_search(){
    
    std::vector<bool> marked(nNodes+1, false);
    std::vector<int> stack;
    std::vector<int> order(nNodes+1);
    std::vector<int>::iterator it;
    
    marked[root] = true;
    order[0]=root;
    stack.push_back(root);
    pred[root] = -1;
    
    int i,j, next=1;
    
    while (stack.size()>0){
        it = stack.end()-1;
        i = (*it);
        j=0;
        for (int p=0; p<=nNodes; p++)
            if((T[i][p]|| T[p][i]) && marked[p] == false){
                j = p;
                break;
            }
        if(j>0){
            order[next] = j;
            pred[j] = i;
            ++next;
            marked[j] = true;
            stack.push_back(j);
        }else{
            stack.pop_back();
        }
    }
   
    for (int p=0; p<nNodes; ++p)
        thread[order[p]] = order[p+1];
    
    thread[order[nNodes]] = 0;
    
}


void spanningTree::computePhi(const Graph & graph){
    phi[root] = 0;
    int j = thread[root];
    int i;
    while (j != root) {
        i = pred[j];
        if (existsIn(i, j, T))
            phi[j] = phi[i] - graph.costs[i][j];
        if (existsIn(j, i, T))
            phi[j] = phi[i] + graph.costs[j][i];
        
        j = thread[j];
    }
}

double spanningTree::computeRC(const Graph & graph,int i, int j){
    double rc;
    rc =  graph.costs[i][j] - phi[i] +phi[j] ;
    return rc;
    
}

bool operator ==(const Arc & a, const Arc & b){
    
    if (a.i == b.i && a.j == b.j) {
        return true;
    }
    
    return false;
}


bool spanningTree::existsIn(int i, int j, const std::vector< std::vector<int> >  & T ){
    
    if (T[i][j]) {
        return true;
    }
    return false;
}


void spanningTree::computeFlow(const Graph & graph){
    std::vector<double> balance(graph.b);
    std::vector<int> Tline(nNodes+1);
    Tline[0] = 0;
    for (int p=1; p<=nNodes; ++p) {
        Tline[p] = thread[Tline[p-1]];
    }

    
    for (int i=0; i<=nNodes; ++i){
        for (int j=0; j<=nNodes; ++j) {
            if (U[i][j]) {
                x[i][j] = graph.capacities[i][j];
                balance[i] -= graph.capacities[i][j];
                balance[j] += graph.capacities[i][j];
            }
        }
    }
    for (int i=0; i<=nNodes; ++i){
        for (int j=0; j<=nNodes; ++j) {
            if (L[i][j])
                x[i][j] = 0;
        }
    }
    int i,j;
    while (*(Tline.end()-1) != 0) {
        j = *(Tline.end()-1);
        i = pred[j];
        if (existsIn(i,j, T))
            x[i][j] = -balance[j];
        else if (existsIn(j,i, T))
                x[j][i] = balance[j];
        balance[i] += balance[j];
        Tline.pop_back();
    }
    
    
}

bool spanningTree::checkOptimality(const Graph & graph, int nMajor){
    
    Arc aux;
    int cont=0;
    while (cont < graph.nArcs && elegibles.size()< nMajor) {
        
        if (lastChecked.i == nNodes && lastChecked.j == nNodes) {
            lastChecked.i = 0;
            lastChecked.j = 0;
            //std::cout<<"zerou"<<std::endl;
        }

        
        for (int i=lastChecked.i; i<=nNodes; ++i) {
            for (int j=lastChecked.j; j<=nNodes; ++j) {
            
                if (elegibles.size()>= nMajor)
                    break;
            
            
                if(graph.arcs[i][j]>=0){
                    ++cont;
                    if(computeRC(graph, i, j) > 0 && existsIn(i,j, U)){
                        Arc e(i,j);
                        if (std::find(elegibles.begin(), elegibles.end(), e) == elegibles.end()) {
                            elegibles.push_back(e);

                        }
                    }
                    if(computeRC(graph, i, j) < 0 && existsIn(i,j, L)){
                        Arc e(i,j);
                        if (std::find(elegibles.begin(), elegibles.end(), e) == elegibles.end()) {
                            elegibles.push_back(e);
                        }
                    }
                }
            
                aux.i = i;
                aux.j = j;
            }
            if (elegibles.size()>= nMajor)
                break;
        }
        lastChecked = aux;
    }
    
    return (elegibles.size() == 0);
    
}

Arc spanningTree::enteringArc(const Graph & graph){
    double max = 0;
    double aux;
    Arc enteringArc (-1,-1);
    std::vector<Arc>::iterator it = elegibles.begin();
    std::vector<Arc>::iterator maxIt = elegibles.end();
    
    while (it != elegibles.end()) {
        aux = computeRC(graph, (*it).i, (*it).j);
       // std::cout<<"current! : "<<(*it).i<<" "<<(*it).j<<std::endl;
        
        if (aux >= 0 && existsIn((*it).i,(*it).j, L)) {
         //   std::cout<<"erase! : "<<(*it).i<<" "<<(*it).j<<std::endl;
            elegibles.erase(it);
            --it;
        }
        else if (aux <= 0 && existsIn((*it).i,(*it).j, U)) {
          //      std::cout<<"erase! : "<<(*it).i<<" "<<(*it).j<<std::endl;
                elegibles.erase(it);
                --it;
            
            }
            else if (max< std::abs(aux)){
                    max = std::abs(aux);
                    maxIt = it;
                }
        ++it;
    }
    if (maxIt != elegibles.end() && elegibles.size() >0) {
        //std::cout<<"max! : "<<(*maxIt).i<<" "<<(*maxIt).j<<std::endl;
        enteringArc = *maxIt;
        elegibles.erase(maxIt);
       // std::cout<<"size"<<elegibles.size()<<std::endl;
    }

    return enteringArc;
}

Arc spanningTree::leavingArc(const Graph & graph, const Arc & enterArcKL){
    double minDelta=0;
    double delta;
    Arc leavingArc(-1,-1);
    int apex;
    int i,j;
    
    //identify mindelta
    i = enterArcKL.i;
    j = enterArcKL.j;
    
    if (existsIn(enterArcKL.i, enterArcKL.j, L)) {
         minDelta = graph.capacities[i][j];
    }else{
        minDelta = x[i][j];
    }
    
    while (i != j) {
        if (depth[i] > depth[j]){
            delta = computeDelta(graph, enterArcKL,pred[i], i);
            i = pred[i];
            
            if (minDelta > delta) {
                minDelta = delta;
            }
        }
        else if(depth[j] > depth[i]){
                delta = computeDelta(graph, enterArcKL,j, pred[j]);
                j = pred[j];
            
                if (minDelta > delta) {
                    minDelta = delta;
                }
            }
            else{
                    delta = computeDelta(graph, enterArcKL, pred[i], i);
                    i = pred[i];
    
                    if (minDelta > delta) {
                        minDelta = delta;
                    }
                    delta = computeDelta(graph, enterArcKL,j , pred[j]);
                    j = pred[j];
                
                    if (minDelta > delta) {
                        minDelta = delta;
                    }

            }
    }
    apex = i;
    
    //Check leaving last blocking arc to leave and push flows
    std::cout<<" minDelta = " <<minDelta<<std::endl;
    i = enterArcKL.i;
    j = enterArcKL.j;

    Arc lastblocking(-1,-1);
    if (existsIn(enterArcKL.i, enterArcKL.j, L))
        x[i][j] += minDelta;
    else x[i][j] -= minDelta;
    
    
    while (i != apex) {
            
        if (minDelta == computeDelta(graph, enterArcKL, pred[i], i)) {
            if (leavingArc.i ==-1 && leavingArc.j==-1) {
                if(existsIn(pred[i], i, T)){leavingArc.i = pred[i]; leavingArc.j = i;}
                if(existsIn(i, pred[i], T)){leavingArc.i = i; leavingArc.j = pred[i];}
            }
        }
        
        if (existsIn(enterArcKL.i, enterArcKL.j, L)) {
            if (existsIn(pred[i],i, T))
                x[pred[i]][i] += minDelta;
            if (existsIn(i, pred[i], T))
                x[i][pred[i]] -= minDelta;
        }else{
            if (existsIn(pred[i],i, T))
                x[pred[i]][i] -= minDelta;
            if (existsIn(i, pred[i], T))
                x[i][pred[i]] += minDelta;
        }
        
        i = pred[i];
        
    }
    while(j != apex){
        
        if (minDelta == computeDelta(graph,enterArcKL, j, pred[j])) {
                if(existsIn(j, pred[j], T)){lastblocking.i = j; lastblocking.j = pred[j];}
                if(existsIn(pred[j], j, T)){lastblocking.i = pred[j]; lastblocking.j = j;}
        }
        
        if (existsIn(enterArcKL.i, enterArcKL.j, L)) {
            if (existsIn(pred[j],j, T))
                x[pred[j]][j] -= minDelta;
            if (existsIn(j, pred[j], T))
                x[j][pred[j]] += minDelta;
        }else{
            if (existsIn(pred[j],j, T))
                x[pred[j]][j] += minDelta;
            if (existsIn(j, pred[j], T))
                x[j][pred[j]] -= minDelta;
        }
        
        j = pred[j];
    
    }
    
    //update leaving arc
    
    if (lastblocking.i !=-1 && lastblocking.j!=-1){
        leavingArc = lastblocking;
    }else if (leavingArc.i == -1 && leavingArc.j == -1)
        leavingArc = enterArcKL;

    return leavingArc;
}

double spanningTree::computeDelta(const Graph & graph, const Arc & enterArcKL, int i, int j){
    double delta = DBL_MAX;
    
    if (existsIn(enterArcKL.i, enterArcKL.j, L)){
        if (existsIn(i, j, T)){
            delta = graph.capacities[i][j] - x[i][j];
        }else if (existsIn(j, i, T))
                delta = x[j][i];
    }else{
        if (existsIn(i, j, T))
            delta = x[i][j];
        else if (existsIn(j, i, T))
            delta = graph.capacities[j][i] - x[j][i];
    }

    
    return delta;
}

void updateTree(spanningTree & bSTTsol, const Graph & graph, Arc entering, Arc leaving){
    
    if (entering == leaving) {
        if (bSTTsol.existsIn(entering.i, entering.j, bSTTsol.L)) {
            if (bSTTsol.x[leaving.i][leaving.j] > 0 ) {
                bSTTsol.L[leaving.i][leaving.j] = 0;
                bSTTsol.U[leaving.i][leaving.j] = 1;
            }
        }
        else if (bSTTsol.existsIn(entering.i, entering.j, bSTTsol.U)){
                if (bSTTsol.x[leaving.i][leaving.j] == 0 ) {
                    bSTTsol.L[leaving.i][leaving.j] = 1;
                    bSTTsol.U[leaving.i][leaving.j] = 0;
                }
        }
    }
    else{
        // if entering != leaving update sets
        bSTTsol.T[leaving.i][leaving.j] = 0;
        bSTTsol.T[entering.i][entering.j] = 1;
        
        if (bSTTsol.U[entering.i][entering.j]) {
            bSTTsol.U[entering.i][entering.j] = 0;
        }else if(bSTTsol.L[entering.i][entering.j]) {
                bSTTsol.L[entering.i][entering.j] = 0;
              }
        
        if (bSTTsol.x[leaving.i][leaving.j] == 0) {
            bSTTsol.L[leaving.i][leaving.j] = 1;
        }
        else if(bSTTsol.x[leaving.i][leaving.j] == graph.capacities[leaving.i][leaving.j]){
                bSTTsol.U[leaving.i][leaving.j] = 1;
             }
        
        //update potentials
        int y, z;
        double change;
        if (bSTTsol.depth[leaving.i]< bSTTsol.depth[leaving.j]) {
            y=leaving.j;
        }else y = leaving.i;
        
        
        int k =entering.i;
        for (int p=0; p<bSTTsol.pred.size(); ++p) {
            if (k == bSTTsol.root) {
                change = -bSTTsol.computeRC(graph, entering.i, entering.j);
                break;
            }
            if (k == y){
                change = +bSTTsol.computeRC(graph, entering.i, entering.j);
                break;
            }
            k = bSTTsol.pred[k];
        }
        
        bSTTsol.phi[y] += change;
        z = bSTTsol.thread[y];
        while (bSTTsol.depth[z] > bSTTsol.depth[y]) {
            bSTTsol.phi[z] += change;
            z = bSTTsol.thread[z];
        }
        
        bSTTsol.make_tree_struct();
        
    }
    
}

double NetworkSimplex(spanningTree & bSTTsol, Graph & graph, int nMajor, int nMinor){
    double mincost = 0;
    bSTTsol.initialize(graph);
    Arc enteringArc(-1,-1);
    Arc leavingArc(-1,-1);
    bool isOptimal = false;
    int minorIteration;
    //candidate list pivot rule
    //major iteration
    
    while (!isOptimal) {
        isOptimal = bSTTsol.checkOptimality(graph, nMajor);
        minorIteration = 0;
        
        while (minorIteration < nMinor) {
            //minor iteration
            enteringArc = bSTTsol.enteringArc(graph);
           
            if (enteringArc.i == -1 || enteringArc.j == -1)
                break;
           // std::cout<<"entrar arco: ("<<enteringArc.i<<", "<<enteringArc.j<<")"<<std::endl;
            
            leavingArc = bSTTsol.leavingArc(graph, enteringArc);
           // std::cout<<"sair arco: ("<<leavingArc.i<<","<<leavingArc.j<<")"<<std::endl;
            
            updateTree(bSTTsol, graph, enteringArc, leavingArc);
            ++minorIteration;
            
            if (bSTTsol.elegibles.size() == 0) {
                break;
            }
        }
       
    }
    bool feasible = true;
    for (int i=0; i<=graph.nNodes; ++i) {
        for (int j=0; j<=graph.nNodes; ++j) {
            mincost += bSTTsol.x[i][j] * graph.costs[i][j];
            if (graph.arcs[i][j]<1 && bSTTsol.x[i][j]>0) {
                feasible = false;
            }
        }
    }
    if (!feasible) {
        mincost = -1000;
    }
    
    return mincost;
    
}

void readInstance(Graph & g, std::string input){
    std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open ufl datafile: %s\n "<<fname;
        abort();
    }
    
    std::string s;
    std::istringstream ss;
    
    getline(file,s);
    getline(file,s);
    ss.str(s);
    
    // read number of locations and number of customers
    
    ss>>g.nNodes;
    ss>>g.nArcs;
    
    ss.clear();
    
    d_k.resize(ndemands);
    arcs.resize(narcs);
    
    for(int i=0; i<narcs; ++i){
        getline(file,s);
        ss.str(s);
        ss>>arcs[i].i;
        ss>>arcs[i].j;
        ss>>arcs[i].c;
        ss>>arcs[i].capa;
        ss>>arcs[i].f;
        ss.clear();
    }
    
    for(int i=0; i<ndemands; ++i){
        getline(file,s);
        ss.str(s);
        ss>> d_k[i].O >> d_k[i].D >> d_k[i].quantity;
        ss.clear();
    }
    
    bijk.resize(ndemands);
    for(int k=0;k<ndemands;++k){
        bijk[k].resize(narcs);
        for(int e=0;e<narcs;++e){
            bijk[k][e] = d_k[k].quantity;
            if(bijk[k][e] > arcs[e].capa)
                bijk[k][e] = arcs[e].capa;
        }
    }
    

    g.arcs.resize(nnodes+1);
    g.capacities.resize(nnodes+1);
    g.costs.resize(nnodes+1);
    g.b.resize(nnodes+1, 0);
    for(int i=0; i<=nnodes;++i){
        g.arcs[i].resize(nnodes+1, -1);
        g.capacities[i].resize(nnodes+1,0);
        g.costs[i].resize(nnodes+1,0);
    }
    
    for(int e=0; e<narcs; ++e){
        g.arcs[ arcs[e].i ][ arcs[e].j ] = 1;
    }
    
    file.close();
}




