//
//  greedy.h
//  TSP
//
//  Created by Sora Sugiyama on 9/26/24.
//

#ifndef greedy_h
#define greedy_h

#include "disjoint_set.h"
#include <algorithm>
#include <vector>
#include <iostream>

namespace TSP{

double greedy(int rout[],problem &X){
    // Store the edge weight alongside the endpoints so the sort comparator is a
    // plain double compare instead of two 2-D lookups (X.W[a][b]) per call,
    // which dominated the O(n^2 log n) sort.
    struct Edge{double w;int a,b;};
    std::vector<Edge>edges;
    edges.resize(X.dimension*(X.dimension-1)>>1);
    disjoint_set ds;
    ds.resize(X.dimension);

    int i,j,idx=0;
    for(i=0;i<X.dimension;i++){
        for(j=0;j<i;j++){
            edges[idx++]={X.W[j][i],j,i};
        }
    }

    std::sort(edges.begin(),edges.end(),[](const Edge&a,const Edge&b){
        return a.w<b.w;
    });

    idx=0;
    std::vector<std::vector<int> >G;
    G.resize(X.dimension);
    
    for(auto &e:edges){
        if(G[e.a].size()>=2||G[e.b].size()>=2)continue;
        if(ds.Disjoint(e.a,e.b)){
            ds.Union(e.a,e.b);
            G[e.a].push_back(e.b);
            G[e.b].push_back(e.a);
        }
    }
    
    int u=0;
    for(i=0;i<X.dimension;i++){
        if(G[i].size()==1){u=i;break;}
    }
    
    std::vector<bool>vis=std::vector<bool>(X.dimension,false);
    vis[u]=true;
    while(1){
        bool stop=true;
        rout[idx++]=u;
        for(int &v:G[u]){
            if(!vis[v]){
                vis[v]=true;
                stop=false;
                u=v;
                break;
            }
        }
        if(stop)break;
    }
    
    if(X.dimension!=idx){
        std::cout<<"FAIL"<<std::endl;
    }
    return loss(rout,X);
}

}

#endif /* greedy_h */
