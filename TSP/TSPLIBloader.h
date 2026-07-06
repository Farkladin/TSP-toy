//
//  TSPLIBloader.h
//  TSP
//
//  Created by Sora Sugiyama on 9/26/24.
//

#ifndef TSPLIBloader_h
#define TSPLIBloader_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unistd.h>
#include "TSP.h"
#include "weight.h"

namespace TSP{

std::string PATH="";

void SetPATH(){
    std::cout<<"new path : "<<std::flush;
    std::cin>>PATH;
}

void makeFile(){
    while(access((PATH+"ALL_tsp.tar.gz").c_str(),F_OK)==-1){
        std::cout<<"File path reset required"<<std::endl;
        SetPATH();
    }

    std::system(("mkdir -p "+PATH+"ALL_tsp").c_str());
    std::system(("tar -zxvf "+PATH+"ALL_tsp.tar.gz -C "+PATH+"ALL_tsp").c_str());
    std::system(("find "+PATH+"ALL_tsp -name '*.gz' -exec gzip -d {} \\;").c_str());
}

void ReadFile(const std::string fileName,problem &X){
    std::ifstream FILE(PATH+"ALL_tsp/"+fileName+".tsp");
    
    if(FILE.fail()){
        std::cout<<"File does not exist."<<std::endl;
        FILE.close();
        exit(1);
        return;
    }

    std::string junk1,junk2;
    getline(FILE,junk1);
    getline(FILE,junk1);
    getline(FILE,junk1);

    FILE>>junk1>>junk2;
    if(junk2==":")FILE>>X.dimension;
    else X.dimension=std::stod(junk2);
    
    FILE>>junk1>>junk2;
    if(junk2==":")FILE>>X.WeightType;
    else X.WeightType=junk2;
    X.W=std::vector<std::vector<double> >(X.dimension,std::vector<double>(X.dimension,1e18));
    if(X.WeightType=="EXPLICIT"){
        std::string FORMAT;
        int i,j;
        
        FILE>>junk1>>junk2;
        if(junk2==":")FILE>>FORMAT;
        else FORMAT=junk2;
        
        while(1){
            getline(FILE,junk1);
            if(junk1=="EDGE_WEIGHT_SECTION")break;
        }
        
        if(FORMAT=="FULL_MATRIX"){
            for(i=0;i<X.dimension;i++){
                for(j=0;j<X.dimension;j++)FILE>>X.W[i][j];
            }
        }else if(FORMAT=="UPPER_ROW"){
            for(i=0;i<X.dimension;i++){
                for(j=i+1;j<X.dimension;j++){
                    FILE>>X.W[i][j];
                    X.W[j][i]=X.W[i][j];
                }
            }
        }else if(FORMAT=="LOWER_ROW"){
            for(i=0;i<X.dimension;i++){
                for(j=0;j<i;j++){
                    FILE>>X.W[i][j];
                    X.W[j][i]=X.W[i][j];
                }
            }
        }else if(FORMAT=="UPPER_DIAG_ROW"){
            for(i=0;i<X.dimension;i++){
                for(j=i;j<X.dimension;j++){
                    FILE>>X.W[i][j];
                    X.W[j][i]=X.W[i][j];
                }
            }
        }else if(FORMAT=="LOWER_DIAG_ROW"){
            for(i=0;i<X.dimension;i++){
                for(j=0;j<i+1;j++){
                    FILE>>X.W[i][j];
                    X.W[j][i]=X.W[i][j];
                }
            }
        }else if(FORMAT=="UPPER_COL"){
            for(i=0;i<X.dimension;i++){
                for(j=0;j<i;j++){
                    FILE>>X.W[j][i];
                    X.W[i][j]=X.W[j][i];
                }
            }
        }else if(FORMAT=="LOWER_COL"){
            for(i=0;i<X.dimension;i++){
                for(j=i+1;j<X.dimension;j++){
                    FILE>>X.W[j][i];
                    X.W[i][j]=X.W[j][i];
                }
            }
        }else if(FORMAT=="UPPER_DIAG_COL"){
            for(i=0;i<X.dimension;i++){
                for(j=0;j<i+1;j++){
                    FILE>>X.W[j][i];
                    X.W[i][j]=X.W[j][i];
                }
            }
        }else if(FORMAT=="LOWER_DIAG_COL"){
            for(i=0;i<X.dimension;i++){
                for(j=i;j<X.dimension;j++){
                    FILE>>X.W[j][i];
                    X.W[i][j]=X.W[j][i];
                }
            }
        }else{
            std::cout<<"WRONG WEIGHT FORMAT"<<std::endl;
            return;
        }
    }else{
        while(1){
            getline(FILE,junk1);
            if(junk1=="NODE_COORD_SECTION")break;
        }
        
        X.coord.resize(X.dimension);
        for(auto &it:X.coord){
            if(X.WeightType=="EUC_3D"||X.WeightType=="MAX_3D"||X.WeightType=="MAN_3D"){
                it.resize(3);
                FILE>>junk1>>it[0]>>it[1]>>it[2];
            }else{
                it.resize(2);
                FILE>>junk1>>it[0]>>it[1];
            }
        }
        
        TSP::makeWeightMatrix(X);
    }
    
    FILE.close();
    
    std::ifstream optFile(PATH+"ALL_tsp/"+fileName+".opt.tour");
    if(optFile.fail()){
        std::cout<<"Opt tour (file) not exist"<<std::endl;
        optFile.close();
        return;
    }
    while(1){
        getline(optFile,junk1);
        if(junk1=="TOUR_SECTION")break;
    }
    
    X.optTour.resize(X.dimension);
    for(auto &it:X.optTour){
        optFile>>it;
        it--;
    }
    
    optFile.close();
}

}

#endif /* TSPLIBloader_h */
