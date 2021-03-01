#include "../RpcL2MuonSA/RPC_FCBM.h"
#include <iostream>
#include <string>

using namespace std;

//void RPC_FCBM::setFrameRegion(const char* station, float& Zmin, float& Rmin, float& Zmax, float& Rmax)
//{
//  double ZERO_LIMIT = 1e-5;
//  double deltaR = 0.3;
//  double deltaZ = 0.3;
//  int nPlauBI = 0; 
//  int nPlauBM = 0; 
//  int nPlauBO = 0; 
//  int nClus_fitInn = (SARPCCluster_fitInnSlope->at(0)).size();
//  int nClus_fitMid = (SARPCCluster_fitMidSlope->at(0)).size();
//  int nClus_fitOut = (SARPCCluster_fitOutSlope->at(0)).size();
//  for(unsigned int iClus_fit = 0; iClus_fit < nClus_fitInn; iClus_fit++){
//    if(isOld){
//      if(!SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit)) continue;
//    }
//    nPlauBI++;
//  }
//  for(unsigned int iClus_fit = 0; iClus_fit < nClus_fitMid; iClus_fit++){
//    if(isOld){
//      if(!SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit)) continue;
//    }
//    nPlauBM++;
//  }
//  for(unsigned int iClus_fit = 0; iClus_fit < nClus_fitOut; iClus_fit++){
//    if(isOld){
//      if(!SARPCCluster_isPlausibleFitOut->at(0).at(iClus_fit)) continue;
//    }
//      nPlauBO++;
//  }
//   
//  std::vector<float> vecZ, vecR;
//  float Z[4] = {0, 0, 0, 0};
//  float R[4] = {0, 0, 0, 0};
//  if(strcmp(station,"inner") == 0){
//    for(int iPlau=0; iPlau<nPlauBI; iPlau++){
//      setMdtAllRegion(iPlau, "inner", Z, R);
//      for(int i = 0; i < 4; i++){
//        vecZ.push_back(Z[i]);
//        vecR.push_back(R[i]);
//      }
//    }
//  }
//  if(strcmp(station,"middle") == 0){
//    for(int iPlau=0; iPlau<nPlauBM; iPlau++){
//      setMdtAllRegion(iPlau, "middle", Z, R);
//      for(int i = 0; i < 4; i++){
//        vecZ.push_back(Z[i]);
//        vecR.push_back(R[i]);
//      }
//    }
//  }
//  if(strcmp(station,"outer") == 0){
//    for(int iPlau=0; iPlau<nPlauBO; iPlau++){
//      setMdtAllRegion(iPlau, "outer", Z, R);
//      for(int i = 0; i < 4; i++){
//        vecZ.push_back(Z[i]);
//        vecR.push_back(R[i]);
//      }
//    }
//  }
//
//  Zmin = *std::min_element(vecZ.begin(), vecZ.end());
//  Zmax = *std::max_element(vecZ.begin(), vecZ.end());
//  Rmin = *std::min_element(vecR.begin(), vecR.end());
//  Rmax = *std::max_element(vecR.begin(), vecR.end());
//  
//  Zmin -= (fabs(Zmin) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 2 : deltaZ;
//  Zmax += (fabs(Zmax) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 2 : deltaZ;
//  Rmin -= (fabs(Rmin) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 2 : deltaR;
//  Rmax += (fabs(Rmax) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 2 : deltaR;
//  cout << "setFrameRegion--DEBUG::Zmin/Zmax/Rmin/Rmax = " << Zmin << "/" << Zmax << "/" << Rmin << "/" << Rmax << endl;
//  
//}

void RPC_FCBM::setFrameRegion(const char* station, float& Zmin, float& Rmin, float& Zmax, float& Rmax)
{
  float ZERO_LIMIT = 1e-5;
  float deltaR = 0.3;
  float deltaZ = 0.3;
  std::vector<float> vecZ, vecR;
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    float Z[4] = {0, 0, 0, 0};
    float R[4] = {0, 0, 0, 0};
    if(strcmp(station,"inner") == 0) setMdtRegion(iMuon, "inner", Z, R);
    else if(strcmp(station,"middle") == 0) setMdtRegion(iMuon, "middle", Z, R);
    else if(strcmp(station,"outer") == 0) setMdtRegion(iMuon, "outer", Z, R);
    
    for(int i = 0; i < 4; i++){
      vecZ.push_back(Z[i]);
      vecR.push_back(R[i]);
    }
  }//iMuon loop end
  Zmin = *std::min_element(vecZ.begin(), vecZ.end());
  Zmax = *std::max_element(vecZ.begin(), vecZ.end());
  Rmin = *std::min_element(vecR.begin(), vecR.end());
  Rmax = *std::max_element(vecR.begin(), vecR.end());
  
  Zmin -= (fabs(Zmin) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 4 : deltaZ;
  Zmax += (fabs(Zmax) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 4 : deltaZ;
  Rmin -= (fabs(Rmin) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 3 : deltaR;
  Rmax += (fabs(Rmax) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 3 : deltaR;
//  Zmin -= (fabs(Zmin) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 5 : deltaZ;
//  Zmax += (fabs(Zmax) >= ZERO_LIMIT) ? fabs(Zmax - Zmin) / 5 : deltaZ;
//  Rmin -= (fabs(Rmin) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 3 : deltaR;
//  Rmax += (fabs(Rmax) >= ZERO_LIMIT) ? fabs(Rmax - Rmin) / 3 : deltaR;
  cout << "setFrameRegion--DEBUG::Zmin/Zmax/Rmin/Rmax = " << Zmin << "/" << Zmax << "/" << Rmin << "/" << Rmax << endl;
  
}

void RPC_FCBM::setMdtRegion(unsigned int iMuon, const char* station, float* Z, float* R)
{
  double SlopeMin = 0;
  double SlopeMax = 0;
  if(strcmp(station,"inner") == 0){
    cout << "setMdtRegion--DEBUG::SAEtaMin->at(" << iMuon << ") = " << SAEtaMin->at(iMuon)[1] << endl; 
    cout << "setMdtRegion--DEBUG::SAEtaMax->at(" << iMuon << ") = " << SAEtaMax->at(iMuon)[1] << endl; 
    SlopeMin = tan(2*atan(exp(- SAEtaMin -> at(iMuon) [0])));
    SlopeMax = tan(2*atan(exp(- SAEtaMax -> at(iMuon) [0])));
    Z[0] = (SArMax -> at(iMuon)[0] / 1000.)/SlopeMin;
    Z[1] = (SArMax -> at(iMuon)[0] / 1000.)/SlopeMax;
    Z[2] = (SArMin -> at(iMuon)[0] / 1000.)/SlopeMax;
    Z[3] = (SArMin -> at(iMuon)[0] / 1000.)/SlopeMin;
    R[0] = (SArMax -> at(iMuon)[0] / 1000.);
    R[1] = (SArMax -> at(iMuon)[0] / 1000.);
    R[2] = (SArMin -> at(iMuon)[0] / 1000.);
    R[3] = (SArMin -> at(iMuon)[0] / 1000.);
  }
  else if(strcmp(station,"middle") == 0){
    cout << "setMdtRegion--DEBUG::SAEtaMin->at(" << iMuon << ") = " << SAEtaMin->at(iMuon)[1] << endl; 
    cout << "setMdtRegion--DEBUG::SAEtaMax->at(" << iMuon << ") = " << SAEtaMax->at(iMuon)[1] << endl; 
    SlopeMin = tan(2*atan(exp(- SAEtaMin -> at(iMuon) [1])));
    SlopeMax = tan(2*atan(exp(- SAEtaMax -> at(iMuon) [1])));
    Z[0] = (SArMax -> at(iMuon)[1] / 1000.)/SlopeMin;
    Z[1] = (SArMax -> at(iMuon)[1] / 1000.)/SlopeMax;
    Z[2] = (SArMin -> at(iMuon)[1] / 1000.)/SlopeMax;
    Z[3] = (SArMin -> at(iMuon)[1] / 1000.)/SlopeMin;
    R[0] = (SArMax -> at(iMuon)[1] / 1000.);
    R[1] = (SArMax -> at(iMuon)[1] / 1000.);
    R[2] = (SArMin -> at(iMuon)[1] / 1000.);
    R[3] = (SArMin -> at(iMuon)[1] / 1000.);
  }
  else if(strcmp(station, "outer")== 0){
    SlopeMin = tan(2*atan(exp(- SAEtaMin -> at(iMuon) [2])));
    SlopeMax = tan(2*atan(exp(- SAEtaMax -> at(iMuon) [2])));
    Z[0] = (SArMax -> at(iMuon)[2] / 1000.)/SlopeMin;
    Z[1] = (SArMax -> at(iMuon)[2] / 1000.)/SlopeMax;
    Z[2] = (SArMin -> at(iMuon)[2] / 1000.)/SlopeMax;
    Z[3] = (SArMin -> at(iMuon)[2] / 1000.)/SlopeMin;
    R[0] = (SArMax -> at(iMuon)[2] / 1000.);
    R[1] = (SArMax -> at(iMuon)[2] / 1000.);
    R[2] = (SArMin -> at(iMuon)[2] / 1000.);
    R[3] = (SArMin -> at(iMuon)[2] / 1000.);
  }
  cout << "setMdtRegion--DEBUG::" << "SlopeMin/Max = " << SlopeMin << "/" << SlopeMax << endl;
}

void RPC_FCBM::setMdtAllRegion(int iClus, const char* station, float* Z, float* R)
{
  double SlopeMin = 0;
  double SlopeMax = 0;
  std::cout << "MdtHit Region for clusterRoad" << std::endl;
  std::cout << SARPCCluster_rMin->at(0).at(1).size() << "/" << SARPCCluster_rMin->at(0).at(2).size() << std::endl;
  if(strcmp(station,"inner") == 0){
    SlopeMin = tan(2*atan(exp(- SARPCCluster_zMin -> at(0).at(0).at(iClus))));
    SlopeMax = tan(2*atan(exp(- SARPCCluster_zMax -> at(0).at(0).at(iClus))));
    Z[0] = (SARPCCluster_rMax -> at(0).at(0)[iClus] / 1000.)/SlopeMin;
    Z[1] = (SARPCCluster_rMax -> at(0).at(0)[iClus] / 1000.)/SlopeMax;
    Z[2] = (SARPCCluster_rMin -> at(0).at(0)[iClus] / 1000.)/SlopeMax;
    Z[3] = (SARPCCluster_rMin -> at(0).at(0)[iClus] / 1000.)/SlopeMin;
    R[0] = (SARPCCluster_rMax -> at(0).at(0)[iClus] / 1000.);
    R[1] = (SARPCCluster_rMax -> at(0).at(0)[iClus] / 1000.);
    R[2] = (SARPCCluster_rMin -> at(0).at(0)[iClus] / 1000.);
    R[3] = (SARPCCluster_rMin -> at(0).at(0)[iClus] / 1000.);
  }
  else if(strcmp(station,"middle") == 0){
    std::cout << "middle station" << endl;
    SlopeMin = tan(2*atan(exp(- SARPCCluster_zMin -> at(0).at(1).at(iClus))));
    SlopeMax = tan(2*atan(exp(- SARPCCluster_zMax -> at(0).at(1).at(iClus))));
    Z[0] = (SARPCCluster_rMax -> at(0).at(1).at(iClus) / 1000.)/SlopeMin;
    Z[1] = (SARPCCluster_rMax -> at(0).at(1).at(iClus) / 1000.)/SlopeMax;
    Z[2] = (SARPCCluster_rMin -> at(0).at(1).at(iClus) / 1000.)/SlopeMax;
    Z[3] = (SARPCCluster_rMin -> at(0).at(1).at(iClus) / 1000.)/SlopeMin;
    R[0] = (SARPCCluster_rMax -> at(0).at(1).at(iClus) / 1000.);
    R[1] = (SARPCCluster_rMax -> at(0).at(1).at(iClus) / 1000.);
    R[2] = (SARPCCluster_rMin -> at(0).at(1).at(iClus) / 1000.);
    R[3] = (SARPCCluster_rMin -> at(0).at(1).at(iClus) / 1000.);
    std::cout << "zMin/zMax/rMin/rMax = " << Z[2] <<"/"<< Z[0] <<"/"<< R[2] <<"/"<< R[0] << std::endl;
  }
  else if(strcmp(station, "outer")== 0){
    SlopeMin = tan(2*atan(exp(- SARPCCluster_zMin -> at(0).at(2).at(iClus))));
    SlopeMax = tan(2*atan(exp(- SARPCCluster_zMax -> at(0).at(2).at(iClus))));
    Z[0] = (SARPCCluster_rMax -> at(0).at(2)[iClus] / 1000.)/SlopeMin;
    Z[1] = (SARPCCluster_rMax -> at(0).at(2)[iClus] / 1000.)/SlopeMax;
    Z[2] = (SARPCCluster_rMin -> at(0).at(2)[iClus] / 1000.)/SlopeMax;
    Z[3] = (SARPCCluster_rMin -> at(0).at(2)[iClus] / 1000.)/SlopeMin;
    R[0] = (SARPCCluster_rMax -> at(0).at(2)[iClus] / 1000.);
    R[1] = (SARPCCluster_rMax -> at(0).at(2)[iClus] / 1000.);
    R[2] = (SARPCCluster_rMin -> at(0).at(2)[iClus] / 1000.);
    R[3] = (SARPCCluster_rMin -> at(0).at(2)[iClus] / 1000.);
  }
}
