#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TVector3.h"

using namespace std;

bool RPC_FCBM::is2muIn1RoI()
{
  bool barrelFlag = true;
  bool rpcFitFlag = true;
  bool OffsegFlag = true;
  bool isOffsegMidFlag = true;
  bool isOffsegOutFlag = true;
  bool OfflinePtFlag = true;
  bool OfflineZ0Flag = true;
  bool isClusFitMidFlag = true;
  bool isSeparateMuon = true;
  bool isRpcSuccess = true;
  bool isSingleRoI = true;
  //    bool isEFpass = false;
  bool isEFpass = true;
  vector<pair<int, int>> MuonPair;
  vector<vector<double>> OffsegMid;
  vector<double> forClear;
  forClear.clear();
  OffsegMid.assign(nMuon, forClear);

  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    if(SAsAddress->at(iMuon) == -1 || SAsAddress->at(iMuon) > 4) barrelFlag == false;
    if(EFPass->at(iMuon)!=1) isEFpass = false;
    //      if(EFPass->at(iMuon)==1) isEFpass = true;
    int countOffsegMid = 0;
    int countOffsegOut = 0;
    if(L1nRoI->at(iMuon) != 1 ) isSingleRoI = false;
    if(fabs(OfflineExtEta->at(iMuon)) >= 1.05) barrelFlag = false;
    //if(fabs(SARPCFitMidSlope->at(iMuon)) <= ZERO_LIMIT) rpcFitFlag = false;
    int nSeg = OfflineNumSegment->at(iMuon);
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      TVector3 vecOffseg;
      vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
      if(fabs(vecOffseg.Eta()) >= 1.05) OffsegFlag = false;
      if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
        countOffsegMid++;
        OffsegMid.at(iMuon).push_back(vecOffseg.Perp());
      }
      if(vecOffseg.Perp() >= 8500.) countOffsegOut++;

    }
    if(countOffsegMid < 1) isOffsegMidFlag = false; 
    if(countOffsegOut < 1) isOffsegOutFlag = false; 
    if(OfflinePt->at(iMuon) < 4) OfflinePtFlag=false;
    if(fabs(OfflineZ0->at(iMuon)) > 50)OfflineZ0Flag = false;

    //if((SARPCCluster_fitMidSlope->at(iMuon)).size() == 0) isClusFitMidFlag = false;
    if(nMuon > 1){
      for(unsigned int jMuon = iMuon+1; jMuon < nMuon; jMuon++){
        float dextEta = sqrt(pow(OfflineExtEta->at(iMuon), 2) + pow(OfflineExtEta->at(jMuon), 2));
        float dextPhi = sqrt(pow(OfflineExtPhi->at(iMuon), 2) + pow(OfflineExtPhi->at(jMuon), 2));
        float dR = sqrt(pow(dextPhi, 2)+ pow(dextEta, 2));

        if(dR >= 0.08) isSeparateMuon = false;
      }
    }
    if(SAisRPCFailure->at(iMuon)) isRpcSuccess = false;
  }
  bool OffsegQualityMid = false;
  if(nMuon == 2 && isOffsegMidFlag){
    int nOff_0 = OffsegMid.at(0).size();
    int nOff_1 = OffsegMid.at(1).size();
    for(int iOff_0 = 0; iOff_0 < nOff_0; iOff_0++){
      for(int iOff_1 = 0; iOff_1 < nOff_1; iOff_1++){
        if(OffsegMid.at(0).at(iOff_0) != OffsegMid.at(1).at(iOff_1)) OffsegQualityMid = true;
      }
    }
  }
  //    if(!isEFpass) continue;
  //    if(isEFpass) continue;
  if(!barrelFlag) return false;
  //if(!rpcFitFlag) continue;
  if(!OffsegFlag) return false;
  //if(!isOffsegMidFlag) continue;
  if(!OffsegQualityMid) return false;
  if(!isOffsegOutFlag) return false;
  if(!OfflinePtFlag) return false;
  if(!OfflineZ0Flag) return false;
  //if(!isSeparateMuon) continue;  // <--- never use this when you study close by event!!
  if(!isRpcSuccess) return false;
  if(!(isCloseByCut(nMuon, MuonPair))) return false;
  if(!(isSingleRoI)) return false;
    //=======================================================
  return true;
}

bool RPC_FCBM::is2muInBarrel()
{
  bool isBarrel = true;
  if(nMuon != 2) return false;
  for(int iMuon=0; iMuon < nMuon; iMuon++){
    if(std::fabs(OfflineExtEta->at(iMuon)) >= 1.05) isBarrel = false;
    if(std::fabs(L1Eta->at(iMuon)) >= 1.05) isBarrel = false; 
  }
  if(!isBarrel) return false;
  return true;
}
