#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TVector3.h"

using namespace std;
void RPC_FCBM::RpcClusterSetter(TGraphErrors* grRPC,
    vector<float>& clusgX,
    vector<float>& clusgY,
    vector<float>& clusgZ,
    vector<int>& clusLayer,
    vector<int>& clusMeasPhi )
{
    for(unsigned int iClus = 0; iClus < (clusgX).size(); iClus++){
      
      TVector3 vecClus;
      vecClus.SetXYZ((clusgX).at(iClus), (clusgY).at(iClus), (clusgZ).at(iClus));

      bool isPhiFlag = false;
      if((clusMeasPhi).at(iClus)){
        isPhiFlag = true;
      }
      int layID = (clusLayer).at(iClus);
      if(isPhiFlag){
        grRPC[layID].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
        grRPC[layID].SetPointError(iClus, 0.05, 0.003);

      }
      else{
        grRPC[layID].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
        grRPC[layID].SetPointError(iClus, 0.003, 0.05);
      }

    /*  if(isPhiFlag){
        if((clusLayer->at(iMuon)).at(iClus) == 0 || (clusLayer->at(iMuon)).at(iClus) == 1){
          grRPC1[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC1[iMuon].SetPointError(iClus, 0.05, 0.003);

          if((clusLayer->at(iMuon)).at(iClus) == 1){
            isOutLay[iMuon].push_back(1);
          }
          else{
            isOutLay[iMuon].push_back(0);
          }
        }
        if((clusLayer->at(iMuon)).at(iClus) == 2 || (clusLayer->at(iMuon)).at(iClus) == 3) {
          grRPC2[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC2[iMuon].SetPointError(iClus, 0.05, 0.003);
          if((clusLayer->at(iMuon)).at(iClus) == 3){
            isOutLay2[iMuon].push_back(1);
          }
          else{
            isOutLay2[iMuon].push_back(0);
          }
        }
        if((clusLayer->at(iMuon)).at(iClus) == 4 || (clusLayer->at(iMuon)).at(iClus) == 5){
          grRPC3[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC3[iMuon].SetPointError(iClus, 0.05, 0.003);
          if((clusLayer->at(iMuon)).at(iClus) >= 5){
            isOutLay3[iMuon].push_back(1);
          }
          else{
            isOutLay3[iMuon].push_back(0);
          }
        }
      }
      else{
        if((clusLayer->at(iMuon)).at(iClus) == 0 || (clusLayer->at(iMuon)).at(iClus) == 1){
          grRPC1[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC1[iMuon].SetPointError(iClus, 0.001, 0.05);
          if((clusLayer->at(iMuon)).at(iClus) == 1){
            isOutLay[iMuon].push_back(1);
          }
          else{
            isOutLay[iMuon].push_back(0);
          }
        }
        if((clusLayer->at(iMuon)).at(iClus) == 2 || (clusLayer->at(iMuon)).at(iClus) == 3) {
          grRPC2[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC2[iMuon].SetPointError(iClus, 0.001, 0.05);
          if((clusLayer->at(iMuon)).at(iClus) == 3){
            isOutLay2[iMuon].push_back(1);
          }
          else{
            isOutLay2[iMuon].push_back(0);
          }
        }
        if((clusLayer->at(iMuon)).at(iClus) == 4 || (clusLayer->at(iMuon)).at(iClus) == 5){
          grRPC3[iMuon].SetPoint(iClus, vecClus.PseudoRapidity(), vecClus.Phi());
          grRPC3[iMuon].SetPointError(iClus, 0.001, 0.05);
          if((clusLayer->at(iMuon)).at(iClus) >= 5){
            isOutLay3[iMuon].push_back(1);
          }
          else{
            isOutLay3[iMuon].push_back(0);
          }
        }

      }*/
    
  }
}
