#include "../RpcL2MuonSA/RPC_FCBM.h"

using namespace std;

void RPC_FCBM::setRPCFitEtaR(){
//calculate rpcfitR and Eta at each RPC plane
  SARPCFitInnR.clear();
  SARPCFitInnEta.clear();
  SARPCFitMidR.clear();
  SARPCFitMidEta.clear();
  SARPCFitOutR.clear();
  SARPCFitOutEta.clear();
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    if (SAsAddress->at(iMuon) == 0){
      //Large Sector
      SARPCFitInnR.push_back(6800.);
      SARPCFitMidR.push_back(7478.);
      SARPCFitOutR.push_back(9832.);

    }
    else if (SAsAddress->at(iMuon) == 1){
      //Small Sector
      SARPCFitInnR.push_back(7820.);
      SARPCFitMidR.push_back(8365.);
      SARPCFitOutR.push_back(10229.);
    }
    else if (SAsAddress->at(iMuon) == 2){
      //Large Special Sector (trial)
      SARPCFitInnR.push_back(7740.);
      SARPCFitMidR.push_back(8100.);
      SARPCFitOutR.push_back(11000.);

    }
    else if (SAsAddress->at(iMuon) == 3){
      //Small Special Sector (trial)
      SARPCFitInnR.push_back(7850.);
      SARPCFitMidR.push_back(8400.);
      SARPCFitOutR.push_back(10400.);
    }
    else{
      SARPCFitInnR.push_back(-9999.);
      SARPCFitMidR.push_back(-9999.);
      SARPCFitOutR.push_back(-9999.);
    }
  }
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    float rpcFitInnEta = -9999.;
    float rpcFitInnSlope = SARPCFitInnSlope->at(iMuon);
    float rpcFitInnOffset = SARPCFitInnOffset->at(iMuon);
    float rpcFitInnR = SARPCFitInnR.at(iMuon);
    
    float rpcFitMidEta = -9999.;
    float rpcFitMidSlope = SARPCFitMidSlope->at(iMuon);
    float rpcFitMidOffset = SARPCFitMidOffset->at(iMuon);
    float rpcFitMidR = SARPCFitMidR.at(iMuon);
    
    float rpcFitOutEta = -9999.;
    float rpcFitOutSlope = SARPCFitOutSlope->at(iMuon);
    float rpcFitOutOffset = SARPCFitOutOffset->at(iMuon);
    float rpcFitOutR = SARPCFitOutR.at(iMuon);
    if(rpcFitInnR > 0){
      float theta = atan(rpcFitMidSlope*(rpcFitInnR)/(rpcFitInnR - rpcFitMidOffset));
      float abs_tan = fabs(tan(theta/2));
      int sign_tan = tan(theta/2)/abs_tan;
      rpcFitInnEta = -sign_tan*log(abs_tan);
    }
    else{
    }
    if(rpcFitMidR > 0){
      float theta = atan(rpcFitMidSlope*(rpcFitMidR)/(rpcFitMidR - rpcFitMidOffset));
      float abs_tan = fabs(tan(theta/2));
      int sign_tan = tan(theta/2)/abs_tan;
      rpcFitMidEta = -sign_tan*log(abs_tan);
    }
    else{
    }
    if(rpcFitOutR > 0){
      float theta = atan(rpcFitOutSlope*(rpcFitOutR)/(rpcFitOutR - rpcFitOutOffset));
      float abs_tan = fabs(tan(theta/2));
      int sign_tan = tan(theta/2)/abs_tan;
      rpcFitOutEta = -sign_tan*log(abs_tan);
    }
    else{
    }
    SARPCFitInnEta.push_back(rpcFitInnEta);
    SARPCFitMidEta.push_back(rpcFitMidEta);
    SARPCFitOutEta.push_back(rpcFitOutEta);
  }

}
