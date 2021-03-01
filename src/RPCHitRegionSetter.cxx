#include "../RpcL2MuonSA/RPC_FCBM.h"

using namespace std;
void RPC_FCBM::setRPCHitRegion(unsigned int nMuon, vector<vector<double>>*& SARPCHitEta, vector<vector<double>>*& SARPCHitPhi, double* rpcHitRegion){
  unsigned int nCol = SARPCHitEta->size();
  double RegionEtaMin = 999.;
  double RegionEtaMax = -999.;
  double RegionPhiMin = 999.;
  double RegionPhiMax = -999.;
  for(unsigned int iCol = 0; iCol < nCol; iCol++){
    unsigned int nEta = (SARPCHitEta->at(iCol)).size();
    for(unsigned int iEta = 0; iEta < nEta; iEta++){
      if((SARPCHitEta->at(iCol))[iEta] < RegionEtaMin) RegionEtaMin = (SARPCHitEta->at(iCol))[iEta];
      if((SARPCHitEta->at(iCol))[iEta] > RegionEtaMax) RegionEtaMax = (SARPCHitEta->at(iCol))[iEta];
    }
    unsigned int nPhi = (SARPCHitPhi->at(iCol)).size();
    for(unsigned int iPhi = 0; iPhi < nPhi; iPhi++){
      if((SARPCHitPhi->at(iCol))[iPhi] < RegionPhiMin) RegionPhiMin = (SARPCHitPhi->at(iCol))[iPhi];
      if((SARPCHitPhi->at(iCol))[iPhi] > RegionPhiMax) RegionPhiMax = (SARPCHitPhi->at(iCol))[iPhi];
    }
  }


  rpcHitRegion[0] = RegionEtaMin-0.05;
  rpcHitRegion[1] = RegionPhiMin-0.1;
  rpcHitRegion[2] = RegionEtaMax+0.05;
  rpcHitRegion[3] = RegionPhiMax+0.1;

}

