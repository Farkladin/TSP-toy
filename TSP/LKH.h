//
//  LKH.h
//  TSP
//
//  Created by Sora Sugiyama on 10/22/24.
//

#ifndef LKH_h
#define LKH_h

#include <algorithm>
#include <vector>
#include <iostream>
#include <random>
#include "TSP.h"

namespace TSP{

class Helsagon{
    // The sequential search only ever reads the first KNN candidates of each
    // node's neighbor list (see the i<5 loops below), so there is no reason to
    // build and fully sort the whole n x n matrix. Keep just the KNN nearest
    // neighbors per node via partial_sort: O(n^2 log KNN) time and O(n*KNN)
    // memory instead of O(n^2 log n) / O(n^2). The pair (weight,index) ordering
    // is a strict total order, so the retained prefix is identical to a full
    // sort's prefix.
    static const int KNN=5;
    void NNM(std::vector<std::vector<int> >&para, problem &X){
        const int n=X.dimension;
        const int k=std::min(n,KNN);
        para.assign(n,std::vector<int>(k));
        std::vector<std::pair<double,int> >vp(n);
        int i,j;
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                vp[j]={X.W[i][j]+(i==j)*1e18,j};
            }
            std::partial_sort(vp.begin(),vp.begin()+k,vp.end());

            for(j=0;j<k;j++)para[i][j]=vp[j].second;
        }
    }
    
    std::vector<bool>use,vis;
    std::vector<int>nxt,pre,pos;
    std::vector<std::pair<int,int> >adj;
    std::vector<std::vector<int> >nnm;
    
public:
    
    void HelsagonNNM(problem &X){
        nnm.clear();
        NNM(nnm,X);
        adj.clear();
        adj.resize(X.dimension);
    }
    
    void edge(const int &u,const int &v,bool flag){
        if(flag){
            if(adj[v].first==-1)adj[v].first=u;
            else adj[v].second=u;
            
            if(adj[u].first==-1)adj[u].first=v;
            else adj[u].second=v;
            
        }else{
            adj[v]={(adj[v].first==u?-1:adj[v].first),(adj[v].second==u?-1:adj[v].second)};
            adj[u]={(adj[u].first==v?-1:adj[u].first),(adj[u].second==v?-1:adj[u].second)};
        }
    }
    
    // if construct==false, function do only the Hamiltonicity check.
    bool constructNewTour(int rout[],const int &n,bool construct=true){
        // Hot path: called very frequently during backtracking. Reuse the
        // buffer instead of reallocating a fresh vector<bool> every call.
        if((int)vis.size()!=n)vis.assign(n,false);
        else std::fill(vis.begin(),vis.end(),false);

        // Hamiltonicity Check
        int u=0,i;
        vis[u]=true;
        for(i=1;i<n;i++){
            if(adj[u].first<0||adj[u].second<0)return false;
            if(!vis[adj[u].first]){
                vis[adj[u].first]=true;
                u=adj[u].first;
            }else if(!vis[adj[u].second]){
                vis[adj[u].second]=true;
                u=adj[u].second;
            }else{
                // Tour can not be close
                return false;
            }
        }

        if(adj[u].first<0||adj[u].second<0)return false;
        if(adj[u].first!=0&&adj[u].second!=0)return false;
        if(!construct)return true;
        
        fill(vis.begin(),vis.end(),false);
        rout[0]=0;
        vis[0]=true;
        for(i=1;i<n;i++){
            u=rout[i-1];
            if(!vis[adj[u].first]){
                vis[adj[u].first]=true;
                rout[i]=adj[u].first;
            }else{
                vis[adj[u].second]=true;
                rout[i]=adj[u].second;
            }
        }
        return true;
    }

    double LKHbacktracking(int rout[],problem &X,int x,int y,double g1,double g2,int k){
        if(k>5)return -1;
        if(k>2){
            if(constructNewTour(rout, X.dimension,false))return -1;
        }
        int z1,z2,i,j;
        double found=-1;
        const int kn=(int)nnm[y].size();
        for(i=0;i<kn;i++){
            z1=nnm[y][i];
            
            if(use[z1]||adj[z1].first==y||adj[z1].second==y)continue;
            use[z1]=true;
            for(j=0;j<2;j++){
                z2=j==0?pre[z1]:nxt[z1];
                
                if(use[z2]||adj[z2].first==x||adj[z2].second==x)continue;
                use[z2]=true;
                
                g1+=X.W[z1][z2];
                edge(z1,z2,false);
                
                g2+=X.W[y][z1];
                edge(y,z1,true);
                if(g1>g2+X.W[z2][x]){
                    edge(z2,x,true);
                    if(constructNewTour(rout,X.dimension)){
                        return loss(rout,X);
                    }
                    edge(z2,x,false);
                }
                found=LKHbacktracking(rout, X, x, z2, g1, g2, k+1);
                if(found>=0)return found;
                
                edge(y,z1,false);
                edge(z1,z2,true);
                use[z2]=false;
                g2-=X.W[y][z1];
                g1-=X.W[z1][z2];
            }
            use[z1]=false;
        }
        
        // fail
        return -1;
    }
  
    // Non-recursive equivalent of LKHbacktracking: the same depth-bounded
    // sequential edge-exchange search, unrolled onto an explicit stack of
    // call frames instead of the call stack. x is constant across all
    // frames (as in the recursive version); each frame resumes its (i,j)
    // loop position after its simulated recursive call returns.
    double LKHbacktrackingIterative(int rout[],problem &X,int x,int y,double g1,double g2,int k){
        struct Frame{
            int y;
            double g1,g2;
            int k;
            int i,j;
            int z1,z2;
        };
        std::vector<Frame> st;
        st.push_back({y,g1,g2,k,0,0,-1,-1});

        double childResult=-1;
        bool haveChildResult=false;

        while(!st.empty()){
            Frame &f=st.back();

            if(haveChildResult){
                haveChildResult=false;
                if(childResult>=0)return childResult;
                edge(f.y,f.z1,false);
                edge(f.z1,f.z2,true);
                use[f.z2]=false;
                f.g2-=X.W[f.y][f.z1];
                f.g1-=X.W[f.z1][f.z2];
                f.j++;
            }else{
                if(f.k>5){
                    childResult=-1;haveChildResult=true;st.pop_back();continue;
                }
                if(f.k>2&&constructNewTour(rout,X.dimension,false)){
                    childResult=-1;haveChildResult=true;st.pop_back();continue;
                }
            }

            bool pushedChild=false;
            const int kn=(int)nnm[f.y].size();
            for(;f.i<kn;f.i++){
                if(f.j==0){
                    f.z1=nnm[f.y][f.i];
                    if(use[f.z1]||adj[f.z1].first==f.y||adj[f.z1].second==f.y)continue;
                    use[f.z1]=true;
                }
                for(;f.j<2;f.j++){
                    f.z2=f.j==0?pre[f.z1]:nxt[f.z1];
                    if(use[f.z2]||adj[f.z2].first==x||adj[f.z2].second==x)continue;
                    use[f.z2]=true;

                    f.g1+=X.W[f.z1][f.z2];
                    edge(f.z1,f.z2,false);
                    f.g2+=X.W[f.y][f.z1];
                    edge(f.y,f.z1,true);

                    if(f.g1>f.g2+X.W[f.z2][x]){
                        edge(f.z2,x,true);
                        if(constructNewTour(rout,X.dimension)){
                            return loss(rout,X);
                        }
                        edge(f.z2,x,false);
                    }

                    st.push_back({f.z2,f.g1,f.g2,f.k+1,0,0,-1,-1});
                    pushedChild=true;
                    break;
                }
                if(pushedChild)break;
                use[f.z1]=false;
                f.j=0;
            }
            if(pushedChild)continue;

            childResult=-1;
            haveChildResult=true;
            st.pop_back();
        }
        return childResult;
    }

    double LKHstepIterative(int rout[], problem &X){
        const int n=X.dimension;
        int u,v,w;
        // Reuse buffers across steps: pos/pre/nxt are fully overwritten below,
        // only use[] must be cleared to false.
        if((int)use.size()!=n){use.assign(n,false);pos.resize(n);pre.resize(n);nxt.resize(n);}
        else std::fill(use.begin(),use.end(),false);

        int x1,x2,i;
        for(int i=0;i<n;i++){
            w=rout[i];u=rout[(i+1)%n];v=rout[(i-1+n)%n];
            nxt[rout[i]]=u;
            pre[rout[i]]=v;
            pos[rout[i]]=i;
            adj[w]={u,v};
        }

        double newTour=-1;
        for(i=0;i<n;i++){
            x1=rout[i];
            x2=nxt[x1];
            use[x1]=use[x2]=true;
            edge(x1,x2,false);
            newTour=LKHbacktrackingIterative(rout, X, x1, x2, X.W[x1][x2], 0.0, 2);
            if(newTour>0)return newTour;

            edge(x1,x2,true);
            use[x1]=use[x2]=false;
        }

        return -1;
    }

    double LKHstep(int rout[], problem &X){
        const int n=X.dimension;
        int u,v,w;
        // Reuse buffers across steps: pos/pre/nxt are fully overwritten below,
        // only use[] must be cleared to false.
        if((int)use.size()!=n){use.assign(n,false);pos.resize(n);pre.resize(n);nxt.resize(n);}
        else std::fill(use.begin(),use.end(),false);

        int x1,x2,i;
        for(int i=0;i<n;i++){
            w=rout[i];u=rout[(i+1)%n];v=rout[(i-1+n)%n];
            nxt[rout[i]]=u;
            pre[rout[i]]=v;
            pos[rout[i]]=i;
            adj[w]={u,v};
        }

        // find k-opt
        double newTour=-1;
        for(i=0;i<n;i++){
            x1=rout[i];
            x2=nxt[x1];
            use[x1]=use[x2]=true;
            edge(x1,x2,false);
            newTour=LKHbacktracking(rout, X, x1, x2, X.W[x1][x2], 0.0, 2);
            if(newTour>0)return newTour;
            
            edge(x1,x2,true);
            use[x1]=use[x2]=false;
        }
        
        // fail
        return -1;
    }
    
};

void LKH(int rout[],problem &X,int step=100000,int pr=1,bool verbose=true){
    if(verbose)std::cout<<"*--------------------------------------"<<std::endl;
    Helsagon solver;
    solver.HelsagonNNM(X);
    double res=-1;
    for(int i=1;i<=step;i++){
        res=solver.LKHstep(rout,X);
        if(res<0){
            if(verbose){
                std::cout<<"Fail to find better tour.\n"<<std::endl;
                std::cout<<"Final Step : "<<i<<"\n"<<"Best tour : "<<loss(rout,X)<<std::endl;
                std::cout<<"---------------------------------------*"<<std::endl;
            }
            return;
        }
        if(verbose&&i%pr==0){
            std::cout<<"Step "<<i<<"\n"<<"Best tour : "<<res<<std::endl;
            std::cout<<"---------------------------------------"<<std::endl;
        }
    }
    if(verbose){
        if(res<0)std::cout<<"Fail to find better tour."<<std::endl;
        std::cout<<"Best tour : "<<loss(rout,X)<<std::endl;
        std::cout<<"--------------------------------------*"<<std::endl;
    }
}

// Perturbs the tour with a random double bridge (4-opt) move, the standard
// kick used by iterated/chained Lin-Kernighan since it cannot be undone by
// a single sequential 2-opt/3-opt move.
void doubleBridge(std::vector<int>&rout,std::mt19937 &rng){
    const int n=(int)rout.size();
    if(n<8)return;
    std::uniform_int_distribution<int> dist(1,n-1);
    std::vector<int> cuts;
    while((int)cuts.size()<3){
        int c=dist(rng);
        if(std::find(cuts.begin(),cuts.end(),c)==cuts.end())cuts.push_back(c);
    }
    std::sort(cuts.begin(),cuts.end());
    const int a=cuts[0],b=cuts[1],c=cuts[2];
    std::vector<int> next;
    next.reserve(n);
    next.insert(next.end(),rout.begin(),rout.begin()+a);
    next.insert(next.end(),rout.begin()+b,rout.begin()+c);
    next.insert(next.end(),rout.begin()+a,rout.begin()+b);
    next.insert(next.end(),rout.begin()+c,rout.end());
    rout=next;
}

// Iterative LKH: same descent to a local optimum as LKH(), but backed by
// LKHbacktrackingIterative's explicit-stack search instead of the recursive
// LKHbacktracking. Produces identical results to LKH(), without recursion.
void IterativeLKH(int rout[],problem &X,int step=100000,int pr=1,bool verbose=true){
    if(verbose)std::cout<<"*--------------------------------------"<<std::endl;
    Helsagon solver;
    solver.HelsagonNNM(X);
    double res=-1;
    for(int i=1;i<=step;i++){
        res=solver.LKHstepIterative(rout,X);
        if(res<0){
            if(verbose){
                std::cout<<"Fail to find better tour.\n"<<std::endl;
                std::cout<<"Final Step : "<<i<<"\n"<<"Best tour : "<<loss(rout,X)<<std::endl;
                std::cout<<"---------------------------------------*"<<std::endl;
            }
            return;
        }
        if(verbose&&i%pr==0){
            std::cout<<"Step "<<i<<"\n"<<"Best tour : "<<res<<std::endl;
            std::cout<<"---------------------------------------"<<std::endl;
        }
    }
    if(verbose){
        if(res<0)std::cout<<"Fail to find better tour."<<std::endl;
        std::cout<<"Best tour : "<<loss(rout,X)<<std::endl;
        std::cout<<"--------------------------------------*"<<std::endl;
    }
}

// Chained LKH (Martin/Otto/Felten, as used in Applegate/Cook/Rohe's CLK):
// descend to a local optimum, then repeatedly perturb with a random double
// bridge kick and re-descend, keeping the kicked tour only if it is not
// worse than the current one, tracking the best tour seen overall.
void ChainedLKH(int rout[],problem &X,int kicks=50,int step=100000,int pr=1){
    const int n=X.dimension;
    std::mt19937 rng(std::random_device{}());

    IterativeLKH(rout,X,step,pr,false);
    std::vector<int> cur(rout,rout+n);
    std::vector<int> best=cur;
    double curLoss=loss(cur,X),bestLoss=curLoss;
    std::cout<<"[ChainedLKH] kick 0 loss="<<curLoss<<" best="<<bestLoss<<std::endl;

    for(int k=1;k<=kicks;k++){
        std::vector<int> kicked=cur;
        doubleBridge(kicked,rng);
        IterativeLKH(kicked.data(),X,step,pr,false);
        double l=loss(kicked,X);
        if(l<=curLoss){
            cur=kicked;
            curLoss=l;
            if(l<bestLoss){
                bestLoss=l;
                best=kicked;
            }
        }
        std::cout<<"[ChainedLKH] kick "<<k<<" loss="<<l<<" current="<<curLoss<<" best="<<bestLoss<<std::endl;
    }
    std::copy(best.begin(),best.end(),rout);
}

}

#endif /* LKH_h */
