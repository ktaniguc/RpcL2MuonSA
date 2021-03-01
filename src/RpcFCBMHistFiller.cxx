#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TVector3.h"
#include <math.h>
#include "../RpcL2MuonSA/RpcDispPreparator.h"
using namespace std;


void RPC_FCBM::FillNoCut(){
  vector<vector<double>> distEachRoIrpcEta; 
  distEachRoIrpcEta.clear();
  vector<vector<double>> distEachRoIrpcPhi; 
  distEachRoIrpcPhi.clear();
  getDistToEachRoI(SARoIEta, SARPCHitEta, distEachRoIrpcPhi, SARPCHitMeasuresPhi, false);
  if(nMuon == 2){
    h_OffEtaPhi->Fill(OfflineEta->at(0)-OfflineEta->at(1), OfflinePhi->at(0)-OfflinePhi->at(1)); 
    double deta = OfflineEta->at(0)-OfflineEta->at(1);
    double dphi = OfflinePhi->at(0)-OfflinePhi->at(1);
    double dR = sqrt(pow(deta, 2) + pow(dphi, 2));
    h_OffdR->Fill(dR);
  }

  if(!distEachRoIrpcPhi.empty()){
    vector<vector<double>>::iterator itrDistPhi = distEachRoIrpcPhi.begin();
    vector<vector<double>>::iterator itrDistPhiEnd = distEachRoIrpcPhi.end();
    for(itrDistPhi; itrDistPhi != itrDistPhiEnd; itrDistPhi++){
      if(!(*itrDistPhi).empty()){
        vector<double>::iterator itrDistPhiParam = (*itrDistPhi).begin();
        vector<double>::iterator itrDistPhiParamEnd = (*itrDistPhi).end();
        for(itrDistPhiParam; itrDistPhiParam != itrDistPhiParamEnd; itrDistPhiParam++){
          h_distnoCutEachRoIrpcPhi->Fill(*itrDistPhiParam);
        }
      }
    }
  }

  getDistToEachRoI(SARoIEta, SARPCHitEta, distEachRoIrpcEta, SARPCHitMeasuresPhi, true);
  if(!distEachRoIrpcEta.empty()){
    vector<vector<double>>::iterator itrDistEta = distEachRoIrpcEta.begin();
    vector<vector<double>>::iterator itrDistEtaEnd = distEachRoIrpcEta.end();
    for(; itrDistEta != itrDistEtaEnd; itrDistEta++){
      if(!(*itrDistEta).empty()){
        vector<double>::iterator itrDistEtaParam = (*itrDistEta).begin();
        vector<double>::iterator itrDistEtaParamEnd = (*itrDistEta).end();
        for(; itrDistEtaParam != itrDistEtaParamEnd; itrDistEtaParam++){
          h_distnoCutEachRoIrpcEta->Fill(*itrDistEtaParam);
        }
      }
    }
  }
  ////////===================================================//////////
  //To know what is Offline Segment 
  //=========================================================//////////
  //
/*
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    int nSeg = OfflineNumSegment->at(iMuon);
    int nClus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      TVector3 vecOffseg;
      vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
      TVector3 vecOffsegDir;
      vecOffsegDir.SetXYZ((OfflineSegmentPx->at(iMuon)).at(iSeg), (OfflineSegmentPy->at(iMuon)).at(iSeg), (OfflineSegmentPz->at(iMuon)).at(iSeg));
      int etaLayID = (int)(OfflineSegmentNTrigEtaLayers->at(iMuon)).at(iSeg);
      int phiLayID = (int)(OfflineSegmentNPhiLayers->at(iMuon)).at(iSeg);
      h_OffSegZR->Fill((OfflineSegmentZ->at(iMuon)).at(iSeg), vecOffseg.Perp());
      h_OffSegXY->Fill((OfflineSegmentX->at(iMuon)).at(iSeg)/1000., (OfflineSegmentY->at(iMuon)).at(iSeg)/1000.);

      h_OffSegEtaPhi->Fill(vecOffseg.PseudoRapidity(), vecOffseg.Phi());
      //h_OffSegNEtaPhiLayers->Fill((OfflineSegmentNTrigEtaLayers->at(iMuon)).at(iSeg), (OfflineSegmentNPhiLayers->at(iMuon)).at(iSeg));
      if(phiLayID == 0){

        h_OffSegXY_etaLayer[etaLayID]->Fill((OfflineSegmentX->at(iMuon)).at(iSeg)/1000., (OfflineSegmentY->at(iMuon)).at(iSeg)/1000.);
      }

      //==========calculate the minimum distance between OffsegMid and clusterRoad
      if(6500. >= vecOffseg.Perp() || vecOffseg.Perp() >= 8500.) continue;
      float min_eta = 999;
      for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
        RpcDispPreparator clustersFit;
        int sAddress = SAsAddress->at(iMuon);
        float rpcFitInn[3] = {-9999., -9999., -9999.};    //{Phi, Slope, Offset}
        float rpcFitMid[3] = {-9999., -9999., -9999.};
        float rpcFitOut[3] = {-9999., -9999., -9999.};
        float rpcFitInnEta = 999;
        float rpcFitMidEta = 999;
        float rpcFitOutEta = 999;
        cout << "FillnoCut()--DEBUG::clusterFitMidSlope=" << 
	    (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit) << endl;
	cout << "FillnoCut()--DEBUG::clusterFitMidOffset=" <<  (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit) << endl;
        if(!(SARPCCluster_isUsingMidCluster->at(iMuon)).at(iClus_fit)) continue;

        clustersFit.SetParameter_rpcFit(
	    1, 
	    rpcFitInn, 
	    rpcFitMid, 
	    rpcFitOut,
	    (SARPCCluster_fitMidPhi->at(iMuon)).at(iClus_fit),  
	    (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit),
	    (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit),
	    (SARPCCluster_fitOutPhi->at(iMuon)).at(iClus_fit),
	    (SARPCCluster_fitOutSlope->at(iMuon)).at(iClus_fit),
	    (SARPCCluster_fitOutOffset->at(iMuon)).at(iClus_fit)
	    );
        clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);
        float delta_eta = fabs(rpcFitMidEta - vecOffseg.PseudoRapidity());
        if(min_eta > delta_eta) {
          min_eta = delta_eta;
        }
      }//nClus_fit loop
      h_distMidOffClusfitEtaMin->Fill(min_eta);
    }//nSeg loop
    cout << "FillnoCut()--DEBUG::OfflineSegmentSector = ";
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      cout << (OfflineSegmentSector->at(iMuon)).at(iSeg) << "/";
    }
    cout << endl;
    cout << "FillnoCut()--DEBUG::OfflineSegmentChamberIndex = ";
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      cout << (OfflineSegmentChamberIndex->at(iMuon)).at(iSeg) << "/";
    }
    cout << endl;
    cout << "FillnoCut()--DEBUG::OfflineSegmentEtaIndex = ";
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      cout << (OfflineSegmentEtaIndex->at(iMuon)).at(iSeg) << "/";
    }
    cout << endl;
    cout << "FillnoCut()--DEBUG::OfflineSegmentNPhiLayers = ";
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      cout << (OfflineSegmentNPhiLayers->at(iMuon)).at(iSeg) << "/";
    }
    cout << endl;
    cout << "FillnoCut()--DEBUG::OfflineSegmentNTrigEtaLayers = ";
    for(int iSeg = 0; iSeg < nSeg; iSeg++){
      cout << (OfflineSegmentNTrigEtaLayers->at(iMuon)).at(iSeg) << "/";
    }
    cout << endl;
  }//Offsegstudy iMuon loop
*/
    //=========================================================//////////

  //its histogram for checking the dEta and dPhi from RoI to Offline

  setRPCFitEtaR();
for(int iPt = 0; iPt<(int)SAptclus->at(0).size(); iPt++){
  m_h_ptclus->Fill(SAptclus->at(0).at(iPt));
}
double ptlead = -999999;
double ptsub = 999999;
if(EFPass->at(0) == 1)m_h_SApt->Fill(SAPt->at(0));

for(int iMuon = 0; iMuon< nMuon; iMuon++){
  if(ptlead < OfflinePt->at(iMuon)) ptlead = OfflinePt->at(iMuon);
  if(ptsub > OfflinePt->at(iMuon)) ptsub = OfflinePt->at(iMuon);
}
      
h_Offpt_lead->Fill(ptlead);
h_Offpt_sublead->Fill(ptsub);
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
//    if(SAOvRmPass->at(iMuon)){ 
//    if(EFPass->at(iMuon)!=1){ 
      h_Offpt->Fill(OfflinePt->at(iMuon));
//    }
    h_SARPCFitInnEta->Fill(SARPCFitInnEta.at(iMuon));
    h_SARPCFitInnR->Fill(SARPCFitInnR.at(iMuon));
    h_SARPCFitMidEta->Fill(SARPCFitMidEta.at(iMuon));
    h_SARPCFitMidR->Fill(SARPCFitMidR.at(iMuon));
    h_SARPCFitOutEta->Fill(SARPCFitOutEta.at(iMuon));
    h_SARPCFitOutR->Fill(SARPCFitOutR.at(iMuon));

    //Fill rpcHit position(x-y) for each layer ID
    unsigned int nHits = (SARPCHitLayer->at(iMuon)).size();
    for(unsigned int iHits = 0; iHits < nHits; iHits++){
      int LayID = (SARPCHitLayer->at(iMuon)).at(iHits);
      h_rpcHitXYLay[LayID]->Fill((SARPCHitX->at(iMuon)).at(iHits)/1000, (SARPCHitY->at(iMuon)).at(iHits)/1000);

    }
    //======Fill the distance between Offline segment and rpcCluster==========//
    vector<float> OffsegMidEta, OffsegMidPhi;
    OffsegMidEta.clear();
    OffsegMidPhi.clear();
    vector<float> OffsegMidSlope;
    OffsegMidSlope.clear();

   
    //calc OfflineSegment eta, phi and select only the information of middle OfflineSegment
    int countMidOffSeg = 0;
    for(int iSeg = 0; iSeg < OfflineNumSegment->at(iMuon); iSeg++){
      TVector3 vecOffseg;
      vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
      TVector3 vecOffsegDir;
      vecOffsegDir.SetXYZ((OfflineSegmentPx->at(iMuon)).at(iSeg), (OfflineSegmentPy->at(iMuon)).at(iSeg), (OfflineSegmentPz->at(iMuon)).at(iSeg));

      if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
	OffsegMidEta.push_back(vecOffseg.PseudoRapidity());
	OffsegMidPhi.push_back(vecOffseg.Phi());
	float Slope = vecOffsegDir.Perp()/(OfflineSegmentPz->at(iMuon)).at(iSeg);
	cout << "OffsegMidSlope = " << Slope << endl;
	cout << "OffsegMidEta/Phi = " << vecOffseg.PseudoRapidity() << "/" << vecOffseg.Phi() << endl;
	OffsegMidSlope.push_back(Slope);
        countMidOffSeg++;
      }
    }//iSeg Loop 
    h_nOffSegMid->Fill(countMidOffSeg);  //middle offline segment number
    int nClus = (SARPCCluster_clusterMeasPhi->at(iMuon)).size();
    int nClus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
    //    if(!(OffsegMidEta.empty())){
    int nOff = OffsegMidEta.size();
/*    for(int iOff = 0; iOff < nOff; iOff++){

      //search the minimum value of distance between Cluster to Offline segment
      float distPhiMin = 999;
      float distEtaMin = 999;
      //calc Cluster eta, phi
      for(int iClus = 0; iClus < nClus; iClus++){
        TVector3 RPCClusterPos;
      RPCClusterPos.SetXYZ((SARPCCluster_gX->at(iMuon)).at(iClus), (SARPCCluster_gY->at(iMuon)).at(iClus), (SARPCCluster_gZ->at(iMuon)).at(iClus));
        //cut the clusters belonging Outer station
        if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 4){
          float distPhi = -(RPCClusterPos.Phi()-OffsegMidPhi.at(iOff));
          float distEta = -(RPCClusterPos.Eta()-OffsegMidEta.at(iOff));
          h_distMidOffClus->Fill(distEta, distPhi);
          if((SARPCCluster_clusterMeasPhi->at(iMuon)).at(iClus)){
            h_distMidOffClusPhi->Fill(fabs(distPhi));
            if(fabs(distPhiMin) > fabs(distPhi)) distPhiMin = distPhi;
          }
          else{
            h_distMidOffClusEta->Fill(fabs(distEta));
            if(fabs(distEtaMin) > fabs(distEta)) distEtaMin = distEta;
          }

          //calc the distance between cluster center and Offline segment
          float deltaR = sqrt(pow(distPhi, 2) + pow(distEta, 2)); 
          h_distMidOffClusR->Fill(deltaR);
        }
      }//iClus loop
      for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
        RpcDispPreparator clustersFit;
        int sAddress = SAsAddress->at(iMuon);
        float rpcFitInn[3] = {-9999., -9999., -9999.};    //{Phi, Slope, Offset}
      float rpcFitMid[3] = {-9999., -9999., -9999.};
      float rpcFitOut[3] = {-9999., -9999., -9999.};
      float rpcFitInnEta = -999;
      float rpcFitMidEta = -999;
      float rpcFitOutEta = -999;
      if(!(SARPCCluster_isUsingMidCluster->at(iMuon)).at(iClus_fit)) continue;
      clustersFit.SetParameter_rpcFit(
          1, 
          rpcFitInn, 
          rpcFitMid, 
          rpcFitOut,
          (SARPCCluster_fitMidPhi->at(iMuon)).at(iClus_fit),  
          (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutPhi->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutSlope->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutOffset->at(iMuon)).at(iClus_fit)
          );
      clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);
      float distOffFit = fabs(rpcFitMidEta - OffsegMidEta.at(iOff));
      h_distMidOffFit->Fill(distOffFit);
      float slopeOffFit = fabs((SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit)-OffsegMidSlope.at(iOff));
      h_slopeMidOffFit->Fill(slopeOffFit);
      }//iClus_fit loop


      h_distMidOffClusPhiMin->Fill(fabs(distPhiMin));
      h_distMidOffClusEtaMin->Fill(fabs(distEtaMin));
      h_distMidOffClusMin->Fill(distEtaMin, distPhiMin);
    }//iOff loop
    //counting the number of fakeRoad
    //
    //
    h_nClusFitMid->Fill(nClus_fit);
    int countFakeRoad = 0;
    int countTruthRoad = 0;
    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
      RpcDispPreparator clustersFit;
      int sAddress = SAsAddress->at(iMuon);
      float rpcFitInn[3] = {-9999., -9999., -9999.}; 
      float rpcFitMid[3] = {-9999., -9999., -9999.};
      float rpcFitOut[3] = {-9999., -9999., -9999.};
      float rpcFitInnEta = -999;
      float rpcFitMidEta = -999;
      float rpcFitOutEta = -999;
      clustersFit.SetParameter_rpcFit(
          1, 
          rpcFitInn, 
          rpcFitMid, 
          rpcFitOut,
          (SARPCCluster_fitMidPhi->at(iMuon)).at(iClus_fit),  
          (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutPhi->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutSlope->at(iMuon)).at(iClus_fit),
          (SARPCCluster_fitOutOffset->at(iMuon)).at(iClus_fit)
          );
      clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);
      h_ClusFitMidEta->Fill(rpcFitMidEta);
      h_ClusFitOutEta->Fill(rpcFitOutEta);
      bool isTruelikeRoad = false;
      for(int iOff = 0; iOff < nOff; iOff++){
        if(fabs(rpcFitMidEta-OffsegMidEta.at(iOff)) <= 0.01 || 
            fabs((SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit)-OffsegMidSlope.at(iOff)) <= 0.3 ) isTruelikeRoad = true;
        //fabs((SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit)-OffsegMidSlope.at(iOff)) <= 0.1*i ) isTruelikeRoad = true;
      }//iOff loop
      if(isTruelikeRoad){ countTruthRoad++; }
      else { countFakeRoad++; }

    }//iClus_fit loop
    h_fakeRoadMid->Fill(countFakeRoad);
    h_truthRoadMid->Fill(countTruthRoad);
    h_faketruthRoadMid->Fill(countFakeRoad, countTruthRoad);
    //  }//!empty()

    ////////////////////////////////////////////////////////////////////////////
    //count the number of clusterRoad with (without) using middle cluster information
    ////////////////////////////////////////////////////////////////////////////
    int countWithMidInfo = 0; //the number of clusterRoads with calculated by middle station's cluster
    int countWithOutInfo = 0;
    if((SARPCCluster_isUsingMidCluster->at(iMuon)).at(iClus_fit)){ countWithMidInfo++;}
    else{ countWithOutInfo++;}
*/

    //======Fill the distance between Offline segment and rpcCluster end======//

    //======Fill the distance between RoI and rpcCluster==========//

    int countLay2Phi = 0;
    int countLay3Phi = 0;
    int countLay2Eta = 0;
    int countLay3Eta = 0;
    int countHit2Phi = 0;
    int countHit3Phi = 0;
    int countHit2Eta = 0;
    int countHit3Eta = 0;
    for(int iClus = 0; iClus < (SARPCCluster_clusterSize->at(iMuon)).size(); iClus++){
      h_clusterSize->Fill((SARPCCluster_clusterSize->at(iMuon)).at(iClus));
      if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) >= 0 && (SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 2){
        h_clusterSizeLay01->Fill((SARPCCluster_clusterSize->at(iMuon)).at(iClus));
      }
      else if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) > 1 && (SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 4){
        h_clusterSizeLay23->Fill((SARPCCluster_clusterSize->at(iMuon)).at(iClus));

      }
      else if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) >= 4 && (SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 8){
        h_clusterSizeLayOut->Fill((SARPCCluster_clusterSize->at(iMuon)).at(iClus));

      }
    }
    for(int iClus = 0; iClus < nClus; iClus++){
      //use only RPC2 info
      if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) > 1 && (SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 4){
        TVector3 RPCClusterPos;
        RPCClusterPos.SetXYZ((SARPCCluster_gX->at(iMuon)).at(iClus), (SARPCCluster_gY->at(iMuon)).at(iClus), (SARPCCluster_gZ->at(iMuon)).at(iClus));
        float distRoIPhi = -(RPCClusterPos.Phi()-SARoIPhi->at(iMuon));
        float distRoIEta = -(RPCClusterPos.Eta()-SARoIEta->at(iMuon));
        float deltaRoIR = sqrt(pow(distRoIPhi, 2) + pow(distRoIEta, 2)); 
        h_distRoIClus->Fill(distRoIEta, distRoIPhi);
        h_distRoIClusR->Fill(deltaRoIR);


        //fill cluster size

        //count the amount of clusters for each layers (2,3) which are in deltaR < 0.1
        if(deltaRoIR < 0.1){
          if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) == 2){
            if((SARPCCluster_clusterMeasPhi->at(iMuon)).at(iClus)){
              countLay2Phi++;
            }
            else{
              countLay2Eta++;
            }
          }
          else if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) == 3){
            if((SARPCCluster_clusterMeasPhi->at(iMuon)).at(iClus)){
              countLay3Phi++;
            }
            else{
              countLay3Eta++;
            }
          }
        }//deltaR < 0.1
      }//layer if
    }//iClus loop


    
    h_nClusLay2Phi->Fill(countLay2Phi);
    h_nClusLay2Eta->Fill(countLay2Eta);
    h_nClusLay3Phi->Fill(countLay3Phi);
    h_nClusLay3Eta->Fill(countLay3Eta);
    h_nClusLay23Phi->Fill(countLay2Phi + countLay3Phi);
    h_nClusLay23Eta->Fill(countLay2Eta + countLay3Eta);
    h_nClusLay23EtaPhi->Fill(countLay2Phi + countLay3Phi + countLay2Eta + countLay3Eta);
    h_nClusHit23Phi->Fill(countHit2Phi + countHit3Phi);
    h_nClusHit23Eta->Fill(countHit2Eta + countHit3Eta);

    //======Fill the distance between RoI and rpcCluster end======//

    //check size of rpcHit(Eta, Phi) for each muons
    for(int jMuon = iMuon + 1; jMuon < nMuon; jMuon++){
      cout << "check size of rpcHit(Eta, Phi) for each muons" << endl;
      cout << "SARPCHitEta->at(iMuon).size() = " << (SARPCHitEta->at(iMuon)).size() << endl;
      cout << "SARPCHitEta->at(jMuon).size() = " << (SARPCHitEta->at(jMuon)).size() << endl;
      
      if((SARPCHitEta->at(iMuon)).size() == (SARPCHitEta->at(jMuon)).size()){
        h_isSameSize_rpcHitEta->Fill(0);
      }
      else{
        h_isSameSize_rpcHitEta->Fill(1);
      }
      if((SARPCHitPhi->at(iMuon)).size() == (SARPCHitPhi->at(jMuon)).size()){
        h_isSameSize_rpcHitPhi->Fill(0);
      }
      else{
        h_isSameSize_rpcHitPhi->Fill(1);
      }
    }



    for(int iHits = 0; iHits < (SARPCHitEta->at(iMuon)).size(); iHits++){

      if((SARPCHitLayer->at(iMuon)).at(iHits) > 1 && (SARPCHitLayer->at(iMuon)).at(iHits) < 4){
        float distRoIPhi = -((SARPCHitPhi->at(iMuon)).at(iHits)-SARoIPhi->at(iMuon));
        float distRoIEta = -(SARPCHitEta->at(iMuon).at(iHits)-SARoIEta->at(iMuon));
        float deltaRoIR = sqrt(pow(distRoIPhi, 2) + pow(distRoIEta, 2)); 
        if(deltaRoIR < 0.1){
          if((SARPCHitLayer->at(iMuon)).at(iHits) == 2){
            if((SARPCHitMeasuresPhi->at(iMuon)).at(iHits)){
              countHit2Phi++;
            }
            else{
              countHit2Eta++;
            }
          }
          else if((SARPCHitLayer->at(iMuon)).at(iHits) == 3){
            if((SARPCHitMeasuresPhi->at(iMuon)).at(iHits)){
              countHit3Phi++;
            }
            else{
              countHit3Eta++;
            }
          }
        
        }//deltaR < 0.1
      }//layer if


    }//iHits loop

  }//iMuon loop

  //check rpcHitMid position is same or not for each muon pair
  if(nMuon == 1){
    h_isSameRpcFitMid->Fill(2);
  }
  else if(nMuon == 2){
    if( SARPCFitMidPhi->at(0) == SARPCFitMidPhi->at(1) && SARPCFitMidEta.at(0) == SARPCFitMidEta.at(1)){
      h_isSameRpcFitMid->Fill(0);
    }
    else{
      h_isSameRpcFitMid->Fill(1);
    }
  }

  else if(nMuon == 3){
    if( SARPCFitMidPhi->at(0) == SARPCFitMidPhi->at(1) && SARPCFitMidEta.at(0) == SARPCFitMidEta.at(1)){
      h_isSameRpcFitMid->Fill(0);
    }
    else{
      h_isSameRpcFitMid->Fill(1);
      if( SARPCFitMidPhi->at(1) == SARPCFitMidPhi->at(2) && SARPCFitMidEta.at(1) == SARPCFitMidEta.at(2)){
        h_isSameRpcFitMid->Fill(0);
      }
      else{
        h_isSameRpcFitMid->Fill(1);
    
        if( SARPCFitMidPhi->at(2) == SARPCFitMidPhi->at(0) && SARPCFitMidEta.at(2) == SARPCFitMidEta.at(0)){
          h_isSameRpcFitMid->Fill(0);
        }
        else{
          h_isSameRpcFitMid->Fill(1);
        }
      }
    }
  }
}

void RPC_FCBM::FillBranchHist(){
/*  size_t count = 2;
  auto tplBranchElmnt = tie(nMuon, OfflinePt);
  for(size_t t = 0; t < count; t++){
    size_t st = tuple_size<tplBranchElmnt>::get<t>(tplBranchElmnt);
    for(int i = 0;i < st; i++){
      cout << get<t>(tplBranchElmnt).at(i) << endl;
    }
  }
*/
  setFillBranch(h_Branch[0], nMuon); 
  setFillBranch_v1(h_Branch[1], OfflinePt); 
  setFillBranch_v1(h_Branch[2], OfflineEta); 
  setFillBranch_v1(h_Branch[3], OfflineZ0); 
  setFillBranch_v1(h_Branch[4], OfflinePhi); 
  setFillBranch_v1(h_Branch[5], OfflineNumSegment); 
  setFillBranch_v2(h_Branch[6], OfflineSegmentX); 
  setFillBranch_v2(h_Branch[7], OfflineSegmentY); 
  setFillBranch_v2(h_Branch[8], OfflineSegmentZ); 
  setFillBranch_v2(h_Branch[9],  OfflineSegmentPx);               
  setFillBranch_v2(h_Branch[10], OfflineSegmentPy);               
  setFillBranch_v2(h_Branch[11], OfflineSegmentPz);               
  setFillBranch_v2(h_Branch[12], OfflineSegmentSector);           
  setFillBranch_v1(h_Branch[13], L1Pt);                            
  setFillBranch_v1(h_Branch[14], L1Eta);                           
  setFillBranch_v1(h_Branch[15], L1Phi);                           
  setFillBranch_v1(h_Branch[16], L1RoINumber);                     
  setFillBranch_v1(h_Branch[17], L1RoISector);                     
  setFillBranch_v1(h_Branch[18], SAPt);                            
  setFillBranch_v1(h_Branch[19], SAEta);                           
  setFillBranch_v1(h_Branch[20], SAPhi);                           
  setFillBranch_v1(h_Branch[21], SAEtaMS);                         
  setFillBranch_v1(h_Branch[22], SAEtaBE);                         
  setFillBranch_v1(h_Branch[23], SAPhiMS);                         
  setFillBranch_v1(h_Branch[24], SAPhiBE);                         
  setFillBranch_v1(h_Branch[25], SABarrelRadius); 
  setFillBranch_v1(h_Branch[26], SAsAddress);                      
  setFillBranch_v1(h_Branch[27], SARoIEta);                        
  setFillBranch_v1(h_Branch[28], SARoIPhi);                        
  setFillBranch_v1(h_Branch[29], SARoINumber);                    
  setFillBranch_v1(h_Branch[30], SARoISector);                    
  setFillBranch_v2(h_Branch[31], SARPCHitX);                      
  setFillBranch_v2(h_Branch[32], SARPCHitY);                      
  setFillBranch_v2(h_Branch[33], SARPCHitZ);                      
  setFillBranch_v2(h_Branch[34], SARPCHitR);                      
  setFillBranch_v2(h_Branch[35], SARPCHitEta);                    
  setFillBranch_v2(h_Branch[36], SARPCHitPhi);                    
  setFillBranch_v2(h_Branch[37], SARPCHitMeasuresPhi);            
 // setFillBranch_v2(h_Branch[38], SARPCHitStationName);            
  setFillBranch_v2(h_Branch[39], SARPCHitStationNumber);          
  setFillBranch_v1(h_Branch[40], SARPCFitInnPhi);                 
  setFillBranch_v1(h_Branch[41], SARPCFitInnSlope);               
  setFillBranch_v1(h_Branch[42], SARPCFitInnOffset);              
  setFillBranch_v1(h_Branch[43], SARPCFitMidPhi);                 
  setFillBranch_v1(h_Branch[44], SARPCFitMidSlope);               
  setFillBranch_v1(h_Branch[45], SARPCFitMidOffset);              
  setFillBranch_v1(h_Branch[46], SARPCFitOutPhi);                 
  setFillBranch_v1(h_Branch[47], SARPCFitOutSlope);               
  setFillBranch_v1(h_Branch[48], SARPCFitOutOffset);              
  setFillBranch_v2(h_Branch[49], SARoadAw);                       
  setFillBranch_v2(h_Branch[50], SARoadBw);                       
  setFillBranch_v1(h_Branch[52], OfflinePt); 
  setFillBranch_v1(h_Branch[53], SAPt); 
  setFillBranch_v2(h_Branch[58], SARPCCluster_fitMidSlope);                       
  setFillBranch_v2(h_Branch[59], SARPCCluster_fitMidOffset);                       
  setFillBranch_v2(h_Branch[60], SARPCCluster_fitMidPhi);                       

  for(int iMuon = 0; iMuon < nMuon; iMuon++)
  {
    int nSegment = OfflineNumSegment->at(iMuon);
    for(int iSeg = 0; iSeg < nSegment; iSeg++)
    {
      TVector3 OfflineSegmentR;
      OfflineSegmentR.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
      h_Branch[51]->Fill(OfflineSegmentR.Perp());
    }

  }


  h_Branch[54]->Fill((int)(SARPCFitMidSlope->size()));
  if(nMuon == 2)
  {
    if(SARPCFitMidSlope->at(0) == SARPCFitMidSlope->at(1)){
      h_Branch[55]->Fill(1);
    }
    else{
      h_Branch[55]->Fill(0);
    }
  }
  else if(nMuon == 3){
    if(SARPCFitMidSlope->at(0) == SARPCFitMidSlope->at(1)){
      h_Branch[55]->Fill(1);
    }
    else if(SARPCFitMidSlope->at(1) == SARPCFitMidSlope->at(2)){
      h_Branch[55]->Fill(1);
    }
    else if(SARPCFitMidSlope->at(2) == SARPCFitMidSlope->at(0)){
      h_Branch[55]->Fill(1);
    }
    else{
      h_Branch[55]->Fill(0);
    }
  }
  else{
    h_Branch[55]->Fill(0);
  }
  setFillBranch_v2(h_Branch[56], SARPCHitLayer);           
  if(nMuon == 2){
    h_OffEtaPhi->Fill(OfflineEta->at(0)-OfflineEta->at(1), OfflinePhi->at(0)-OfflinePhi->at(1)); 
  }
  setFillBranch_v2(h_Branch[57], SARPCHitTime);                      
}



template <class Type> void RPC_FCBM::setFillBranch_v2(TH1D*& h_Branch, vector<vector<Type>>*& BranchName){
  for(auto& itrBranch : *BranchName)
  {
    if(!itrBranch.empty()){
      for(auto& itrBranchElmnt : itrBranch){
        h_Branch->Fill(itrBranchElmnt);
      }
    }
  }
}

template <class Type> void RPC_FCBM::setFillBranch_v1(TH1D*& h_Branch, vector<Type>*& BranchName){
  for(auto& itrBranch : *BranchName)
  {
    h_Branch->Fill(itrBranch);
    
  }
}

template <class Type> void RPC_FCBM::setFillBranch(TH1D*& h_Branch, Type& BranchName){
  h_Branch->Fill(BranchName);
}



void RPC_FCBM::FillFCBM(vector<pair<int, int>>& MUpair){
  //=====fill the number of RpcClusters at each layers ==========// 
  if(nMuon == 2){
    int nClus[8] = {};
    int nClusMeasPhi[8] = {};
    for(unsigned int iClus = 0; iClus < (SARPCCluster_clusterSize->at(0)).size(); iClus++){

      int iLayer = (SARPCCluster_clusterLayer->at(0)).at(iClus);
        nClus[iLayer]++;
      if((SARPCCluster_clusterMeasPhi->at(0)).at(iClus)){
        nClusMeasPhi[iLayer]++;
      }
    }
    for(int i = 0; i < 8; i++){
      h_rpcClusPhiSizeLay_CloseBy[i]->Fill(nClusMeasPhi[i]);
      h_rpcClusEtaSizeLay_CloseBy[i]->Fill(nClus[i]-nClusMeasPhi[i]);
    }
  }

  //=====fill info about rpcCluster end =================//  
  
  
  vector < vector < double > > rpcHitEta;
  vector < vector < double > > rpcHitPhi;
  vector < vector < unsigned int > > rpcHitMeasuresPhi;
  vector < float >  roiEta;
  vector < float >  roiPhi;
  getCloseBy(SARPCHitEta, rpcHitEta);
  getCloseBy(SARPCHitPhi, rpcHitPhi);
  getCloseBy(SARPCHitMeasuresPhi, rpcHitMeasuresPhi);
  getCloseBy(SARoIEta, roiEta);
  getCloseBy(SARoIPhi, roiPhi);

  int nMu = roiEta.size();
  for(int iMu = 0; iMu < nMu; iMu++){
    for(int iElmnt = 0; iElmnt < rpcHitEta[iMu].size(); iElmnt++){
      float distRoIrpcEta = 0;
      float distRoIrpcPhi = 0;
      if(rpcHitMeasuresPhi[iMu].at(iElmnt)){
        distRoIrpcPhi = fabs(roiPhi.at(iMu)-rpcHitPhi[iMu].at(iElmnt));
        h_distRoIrpcPhi->Fill(distRoIrpcPhi);
      }
      else
      {
        distRoIrpcEta = fabs(roiEta.at(iMu)-rpcHitEta[iMu].at(iElmnt));
        h_distRoIrpcEta->Fill(distRoIrpcEta);
      }
    }
  }

  //RPC Cluster num at RPC2
  /*
  if(nMuon == 2){
    int countEta = 0;
    int countPhi = 0;
    vector<bool> ClusMeasuresPhi;
    ClusMeasuresPhi.clear();
    for(int i = 0; i < (SARPCHitMeasuresPhi->at(0)).size(); i++){
      if((SARPCHitMeasuresPhi->at(0)[i])){
        cout << "ktaniguc::CHECK:SARPCHitMeasuresPhiDependency: RPCHitEta=" << (SARPCHitEta->at(0))[i] << endl;
      }
      else{
        cout << "ktaniguc::CHECK:SARPCHitMeasuresPhiDependency: RPCHitPhi=" << (SARPCHitPhi->at(0))[i] << endl;
      }
    }
    for(unsigned int j = 0; j < (SARPCCluster_gX->at(0)).size(); j++){
      TVector3 CheckClus;
      CheckClus.SetXYZ((SARPCCluster_gX->at(0)).at(j), (SARPCCluster_gY->at(0)).at(j), (SARPCCluster_gZ->at(0)).at(j));
      cout << "ktaniguc::CHECK:SARPCHitMeasuresPhiDependency: RPCClusterPhi=" << CheckClus.Phi() << endl;
      cout << "ktaniguc::CHECK:SARPCHitMeasuresEtaDependency: RPCClusterEta=" << CheckClus.PseudoRapidity() << endl;
      bool isPhi = false;

      for(int i = 0; i < (SARPCHitMeasuresPhi->at(0)).size(); i++){
        if(!(SARPCHitMeasuresPhi->at(0)[i])){
          cout << "===========================================================" << endl;
          cout << "deltaPhi = " << (SARPCHitPhi->at(0))[i] - CheckClus.Phi() << endl;
          cout << "===========================================================" << endl;

        }
        else{
          cout << "deltaEta = " << (SARPCHitEta->at(0))[i] - CheckClus.PseudoRapidity() << endl;
          cout << "===========================================================" << endl;
          if(fabs((SARPCHitPhi->at(0))[i] - CheckClus.Phi()) < 0.0009 ){
            isPhi = true;
          }
        }
      }
      ClusMeasuresPhi.push_back(isPhi);
    }
    for(int i = 0; i < ClusMeasuresPhi.size(); i++){
      cout << "ktaninguc::CHECK: ClusMeasuresPhi=" << ClusMeasuresPhi.at(i) << endl;;
    }

    cout << "size :: ClusMeasuresPhi/Cluster_gX = " << ClusMeasuresPhi.size() << "/" << (SARPCCluster_gZ->at(0)).size() << endl;
//=============================MeasuresPhi end=================================================

    for(unsigned int iClus = 0; iClus < (SARPCCluster_gX->at(0)).size(); iClus++){
      TVector3 vecClus;
      if(SAsAddress->at(0) == 0){
        vecClus.SetXYZ((SARPCCluster_gX->at(0)).at(iClus), (SARPCCluster_gY->at(0)).at(iClus), (SARPCCluster_gZ->at(0)).at(iClus));
        if(7200. < vecClus.Perp() && vecClus.Perp() <= 8500.){
          if(ClusMeasuresPhi.at(iClus)){
            countPhi++;
          }else{
            countEta++;
          }
        }
      }
      else{
        vecClus.SetXYZ((SARPCCluster_gX->at(0)).at(iClus), (SARPCCluster_gY->at(0)).at(iClus), (SARPCCluster_gZ->at(0)).at(iClus));
        if(8000. < vecClus.Perp() && vecClus.Perp() <= 9000.){
          if(ClusMeasuresPhi.at(iClus)){
            countPhi++;
          }else{
            countEta++;
          }
        }
      }
      TVector3 vecCluster;
      vecCluster.SetXYZ((SARPCCluster_gX->at(0)).at(iClus), (SARPCCluster_gY->at(0)).at(iClus), (SARPCCluster_gZ->at(0)).at(iClus));
      if(!(ClusMeasuresPhi.at(iClus))){
        h_rpcClusEta_CloseBy->Fill(vecCluster.PseudoRapidity());
      }else{
        h_rpcClusPhi_CloseBy->Fill(vecCluster.Phi());
      }
    }  
    h_rpc2ClusEtaNum_CloseBy->Fill(countEta);
    h_rpc2ClusPhiNum_CloseBy->Fill(countPhi);
    

  }
  */
  //isSAMERpcFitMidPosition
  for(int iPair = 0; iPair < MUpair.size(); iPair++){
    if( SARPCFitMidPhi->at(MUpair[iPair].first) == SARPCFitMidPhi->at(MUpair[iPair].second) && SARPCFitMidEta.at(MUpair[iPair].first) == SARPCFitMidEta.at(MUpair[iPair].second)){
      h_isSameRpcFitMid_CloseBy->Fill(0);
    }
    else{
      h_isSameRpcFitMid_CloseBy->Fill(1);
    }
  }

  // check the histgrams of distance from RoI to OfflinePosition
  
  for(int iPair = 0; iPair < MUpair.size(); iPair++){
    int muID1 = MUpair[iPair].first;
    int muID2 = MUpair[iPair].second;
    //check the distribution of delta R (RoI - cluster)
    int nClus = (SARPCCluster_gX->at(muID1)).size();
    for(int iClus = 0; iClus < nClus; iClus++){
      //use only RPC2 info
      if((SARPCCluster_clusterLayer->at(muID1)).at(iClus) > 1 && (SARPCCluster_clusterLayer->at(muID1)).at(iClus) < 4){
        TVector3 RPCClusterPos;
        RPCClusterPos.SetXYZ((SARPCCluster_gX->at(muID1)).at(iClus), (SARPCCluster_gY->at(muID1)).at(iClus), (SARPCCluster_gZ->at(muID1)).at(iClus));
        float distRoIPhi = -(RPCClusterPos.Phi()-SARoIPhi->at(muID1));
        float distRoIEta = -(RPCClusterPos.Eta()-SARoIEta->at(muID1));
        float deltaRoIR = sqrt(pow(distRoIPhi, 2) + pow(distRoIEta, 2)); 
        h_distRoIClusR_CloseBy->Fill(deltaRoIR);
      }
    }
    
    int nClus2 = (SARPCCluster_gX->at(muID2)).size();
    for(int iClus = 0; iClus < nClus2; iClus++){
      //use only RPC2 info
      if((SARPCCluster_clusterLayer->at(muID2)).at(iClus) > 1 && (SARPCCluster_clusterLayer->at(muID2)).at(iClus) < 4){
        TVector3 RPCClusterPos;
        RPCClusterPos.SetXYZ((SARPCCluster_gX->at(muID2)).at(iClus), (SARPCCluster_gY->at(muID2)).at(iClus), (SARPCCluster_gZ->at(muID2)).at(iClus));
        float distRoIPhi = -(RPCClusterPos.Phi()-SARoIPhi->at(muID2));
        float distRoIEta = -(RPCClusterPos.Eta()-SARoIEta->at(muID2));
        float deltaRoIR = sqrt(pow(distRoIPhi, 2) + pow(distRoIEta, 2)); 
        h_distRoIClusR_CloseBy->Fill(deltaRoIR);
      }
    }
    
    double distOffEta = OfflineEta->at(muID1)-OfflineEta->at(muID2);
    double distOffPhi = OfflinePhi->at(muID1)-OfflinePhi->at(muID2);
    h_distOffEta_CloseBy->Fill(distOffEta);
    h_distOffPhi_CloseBy->Fill(distOffPhi);
    
    //rpcfit inner
    double distRpcFitInnEta = SARPCFitInnEta.at(muID1)-SARPCFitInnEta.at(muID2);
    double distRpcFitInnPhi = SARPCFitInnPhi->at(muID1)-SARPCFitInnPhi->at(muID2);
    h_distRpcFitInnEta_CloseBy->Fill(distRpcFitInnEta);
    h_distRpcFitInnPhi_CloseBy->Fill(distRpcFitInnPhi);
    
    //rpcfit middle
    double distRpcFitMidEta = SARPCFitMidEta.at(muID1)-SARPCFitMidEta.at(muID2);
    double distRpcFitMidPhi = SARPCFitMidPhi->at(muID1)-SARPCFitMidPhi->at(muID2);
    h_distRpcFitMidEta_CloseBy->Fill(distRpcFitMidEta);
    h_distRpcFitMidPhi_CloseBy->Fill(distRpcFitMidPhi);
    
    //rpcfit outer
    double distRpcFitOutEta = SARPCFitOutEta.at(muID1)-SARPCFitOutEta.at(muID2);
    double distRpcFitOutPhi = SARPCFitOutPhi->at(muID1)-SARPCFitOutPhi->at(muID2);
    h_distRpcFitOutEta_CloseBy->Fill(distRpcFitOutEta);
    h_distRpcFitOutPhi_CloseBy->Fill(distRpcFitOutPhi);
     
    h_RoInumDelEta->Fill( distOffEta, 0);
    h_RoInumDelPhi->Fill( distOffPhi, 0);
  
    //RPCHit distance from rpc hit layer to rpcfit 
    vector<double> DeltaRPCOutPhi;
    vector<double> DeltaRPCOutEta;
    DeltaRPCOutPhi.clear();
    DeltaRPCOutEta.clear();
    vector<double> DeltaRPCMidPhi;
    vector<double> DeltaRPCMidEta;
    DeltaRPCMidPhi.clear();
    DeltaRPCMidEta.clear();

    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(muID1)).size(); iRPCHit++){
      unsigned int layID = (SARPCHitLayer->at(muID1)).at(iRPCHit);
      if(layID == 0 || layID == 1){
        if((SARPCHitMeasuresPhi->at(muID1)).at(iRPCHit)){
          h_distRpcHitToFitInnPhi_CloseBy->Fill(SARPCFitMidPhi->at(muID1)-(SARPCHitPhi->at(muID1)).at(iRPCHit));
        }
        else{
          h_distRpcHitToFitInnEta_CloseBy->Fill(SARPCFitInnEta.at(muID1)-(SARPCHitEta->at(muID1)).at(iRPCHit));
        }
      }
      else if(layID == 2 || layID == 3){
        if((SARPCHitMeasuresPhi->at(muID1)).at(iRPCHit)){
          h_distRpcHitToFitMidPhi_CloseBy->Fill(SARPCFitMidPhi->at(muID1)-(SARPCHitPhi->at(muID1)).at(iRPCHit));
          DeltaRPCMidPhi.push_back(fabs(SARPCFitMidPhi->at(muID1)-(SARPCHitPhi->at(muID1)).at(iRPCHit)));
        }
        else{
          h_distRpcHitToFitMidEta_CloseBy->Fill(SARPCFitMidEta.at(muID1)-(SARPCHitEta->at(muID1)).at(iRPCHit));
          DeltaRPCMidEta.push_back(fabs(SARPCFitMidEta.at(muID1)-(SARPCHitEta->at(muID1)).at(iRPCHit)));
        }
      }
      else if(layID == 4 || layID == 5){
        if((SARPCHitMeasuresPhi->at(muID1)).at(iRPCHit)){
          h_distRpcHitToFitOutPhi_CloseBy->Fill(SARPCFitOutPhi->at(muID1)-(SARPCHitPhi->at(muID1)).at(iRPCHit));
          DeltaRPCOutPhi.push_back(fabs(SARPCFitOutPhi->at(muID1)-(SARPCHitPhi->at(muID1)).at(iRPCHit)));
        }
        else{
          h_distRpcHitToFitOutEta_CloseBy->Fill(SARPCFitOutEta.at(muID1)-(SARPCHitEta->at(muID1)).at(iRPCHit));
          DeltaRPCOutEta.push_back(fabs(SARPCFitOutEta.at(muID1)-(SARPCHitEta->at(muID1)).at(iRPCHit)));
        }
      }
    }
    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(muID2)).size(); iRPCHit++){
      unsigned int layID = (SARPCHitLayer->at(muID2)).at(iRPCHit);
      if(layID == 0 || layID == 1){
        if((SARPCHitMeasuresPhi->at(muID2)).at(iRPCHit)){
        h_distRpcHitToFitInnPhi_CloseBy->Fill(SARPCFitMidPhi->at(muID2)-(SARPCHitPhi->at(muID2)).at(iRPCHit));
        }
        else{
        h_distRpcHitToFitInnEta_CloseBy->Fill(SARPCFitInnEta.at(muID2)-(SARPCHitEta->at(muID2)).at(iRPCHit));
        }
      }
      else if(layID == 2 || layID == 3){
        if((SARPCHitMeasuresPhi->at(muID2)).at(iRPCHit)){
          h_distRpcHitToFitMidPhi_CloseBy->Fill(SARPCFitMidPhi->at(muID2)-(SARPCHitPhi->at(muID2)).at(iRPCHit));
          DeltaRPCMidPhi.push_back(fabs(SARPCFitMidPhi->at(muID2)-(SARPCHitPhi->at(muID2)).at(iRPCHit)));
        }
        else{
          h_distRpcHitToFitMidEta_CloseBy->Fill(SARPCFitMidEta.at(muID2)-(SARPCHitEta->at(muID2)).at(iRPCHit));
          DeltaRPCMidEta.push_back(fabs(SARPCFitMidEta.at(muID2)-(SARPCHitEta->at(muID2)).at(iRPCHit)));
        }
      }
      else if(layID == 4 || layID == 5){
        if((SARPCHitMeasuresPhi->at(muID2)).at(iRPCHit)){
          DeltaRPCOutPhi.push_back(fabs(SARPCFitOutPhi->at(muID2)-(SARPCHitPhi->at(muID2)).at(iRPCHit)));
          h_distRpcHitToFitOutPhi_CloseBy->Fill(SARPCFitOutPhi->at(muID2)-(SARPCHitPhi->at(muID2)).at(iRPCHit));
        }
        else{
          h_distRpcHitToFitOutEta_CloseBy->Fill(SARPCFitOutEta.at(muID2)-(SARPCHitEta->at(muID2)).at(iRPCHit));
          DeltaRPCOutEta.push_back(fabs(SARPCFitOutEta.at(muID2)-(SARPCHitEta->at(muID2)).at(iRPCHit)));
        }
      }
    }
    vector<double>::iterator itrDelrpcOutPhi = min_element( DeltaRPCOutPhi.begin(), DeltaRPCOutPhi.end() );
    vector<double>::iterator itrDelrpcOutEta = min_element( DeltaRPCOutEta.begin(), DeltaRPCOutEta.end() );
    if(!DeltaRPCOutPhi.empty()){
      double DelrpcOutPhi = *itrDelrpcOutPhi;
      cout << DelrpcOutPhi << endl;

      h_DeltaRpcOutMinPhi_CloseBy->Fill(DelrpcOutPhi);
    }
    if(!DeltaRPCOutEta.empty()){
      double DelrpcOutEta = *itrDelrpcOutEta;
      h_DeltaRpcOutMinEta_CloseBy->Fill(DelrpcOutEta);
      cout << DelrpcOutEta << endl;
    }
   
    vector<double>::iterator itrDelrpcMidPhi = min_element( DeltaRPCMidPhi.begin(), DeltaRPCMidPhi.end() );
    vector<double>::iterator itrDelrpcMidEta = min_element( DeltaRPCMidEta.begin(), DeltaRPCMidEta.end() );
    if(!DeltaRPCMidPhi.empty()){
      double DelrpcMidPhi = *itrDelrpcMidPhi;

      h_DeltaRpcMidMinPhi_CloseBy->Fill(DelrpcMidPhi);
    }
    if(!DeltaRPCMidEta.empty()){
      double DelrpcMidEta = *itrDelrpcMidEta;
      h_DeltaRpcMidMinEta_CloseBy->Fill(DelrpcMidEta);
    }
    
    //deltapt (compare the same muon's info)
    h_DeltaPt_OffSA_CloseBy->Fill(OfflinePt->at(muID1) - SAPt->at(muID1));
    h_DeltaPt_OffSA_CloseBy->Fill(OfflinePt->at(muID2) - SAPt->at(muID2));
    int nSARPCHit1 = (SARPCHitEta->at(muID1)).size(); 
    //int nSARPCHit2 = (SARPCHitEta->at(muID2)).size();
    h_rpcHitSize_CloseBy->Fill(nSARPCHit1 );
    //h_rpcHitSize_CloseBy->Fill( nSARPCHit2 );
    
    vector<vector<double>> rpcHitLayID1(8, vector<double>(0));
    rpcHitLayID1.clear();
    vector<vector<double>> rpcHitLayID2(8, vector<double>(0));
    rpcHitLayID2.clear();
    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(muID1)).size(); iRPCHit++){
      unsigned int layID1 = (SARPCHitLayer->at(muID1)).at(iRPCHit);
      rpcHitLayID1[layID1].push_back((SARPCHitEta->at(muID1)).at(iRPCHit));
    }
    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(muID2)).size(); iRPCHit++){
      unsigned int layID2 = (SARPCHitLayer->at(muID2)).at(iRPCHit);
      rpcHitLayID2[layID2].push_back((SARPCHitEta->at(muID2)).at(iRPCHit));
    }
    for(int iLay = 0; iLay < 8; iLay++){
      h_rpcHitSizeLay_CloseBy[iLay]->Fill(rpcHitLayID1[iLay].size());
      h_rpcHitSizeLay_CloseBy[iLay]->Fill(rpcHitLayID2[iLay].size());
    }
  }

}

void RPC_FCBM::FillnoFCBM(){
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    h_DeltaPt_OffSA->Fill(OfflinePt->at(iMuon) - SAPt->at(iMuon));
    h_rpcHitSize->Fill( (SARPCHitEta->at(iMuon)).size() );
    vector<vector<double>> rpcHitLayID(8, vector<double>(0));
    rpcHitLayID.clear();
      //check the distribution of delta R (RoI - cluster)
    int nClus = (SARPCCluster_gX->at(iMuon)).size();
    for(int iClus = 0; iClus < nClus; iClus++){
      //use only RPC2 info
      if((SARPCCluster_clusterLayer->at(iMuon)).at(iClus) > 1 && (SARPCCluster_clusterLayer->at(iMuon)).at(iClus) < 4){
        TVector3 RPCClusterPos;
        RPCClusterPos.SetXYZ((SARPCCluster_gX->at(iMuon)).at(iClus), (SARPCCluster_gY->at(iMuon)).at(iClus), (SARPCCluster_gZ->at(iMuon)).at(iClus));
        float distRoIPhi = -(RPCClusterPos.Phi()-SARoIPhi->at(iMuon));
        float distRoIEta = -(RPCClusterPos.Eta()-SARoIEta->at(iMuon));
        float deltaRoIR = sqrt(pow(distRoIPhi, 2) + pow(distRoIEta, 2)); 
        h_distRoIClusR_noCloseBy->Fill(deltaRoIR);
      }
    }
    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(iMuon)).size(); iRPCHit++){
      unsigned int layID = (SARPCHitLayer->at(iMuon)).at(iRPCHit);
      rpcHitLayID[layID].push_back((SARPCHitEta->at(iMuon)).at(iRPCHit));
    }
    for(int iLay = 0; iLay < 8; iLay++){
      h_rpcHitSizeLay[iLay]->Fill(rpcHitLayID[iLay].size());
    }

    
    for(int iRPCHit = 0; iRPCHit < (SARPCHitLayer->at(iMuon)).size(); iRPCHit++){
      unsigned int layID = (SARPCHitLayer->at(iMuon)).at(iRPCHit);
      if(layID == 0 || layID == 1){
        if((SARPCHitMeasuresPhi->at(iMuon)).at(iRPCHit)){
          h_distRpcHitToFitInnPhi->Fill(SARPCFitMidPhi->at(iMuon)-(SARPCHitPhi->at(iMuon)).at(iRPCHit));
        }
        else{
          h_distRpcHitToFitInnEta->Fill(SARPCFitInnEta.at(iMuon)-(SARPCHitEta->at(iMuon)).at(iRPCHit));
        }
      }
      else if(layID == 2 || layID == 3){
        if((SARPCHitMeasuresPhi->at(iMuon)).at(iRPCHit)){
          h_distRpcHitToFitMidPhi->Fill(SARPCFitMidPhi->at(iMuon)-(SARPCHitPhi->at(iMuon)).at(iRPCHit));
        }
        else{
          h_distRpcHitToFitMidEta->Fill(SARPCFitMidEta.at(iMuon)-(SARPCHitEta->at(iMuon)).at(iRPCHit));
        }
      }
      else if(layID == 4 || layID == 5){
        if((SARPCHitMeasuresPhi->at(iMuon)).at(iRPCHit)){
          h_distRpcHitToFitOutPhi->Fill(SARPCFitOutPhi->at(iMuon)-(SARPCHitPhi->at(iMuon)).at(iRPCHit));
        }
        else{
          h_distRpcHitToFitOutEta->Fill(SARPCFitOutEta.at(iMuon)-(SARPCHitEta->at(iMuon)).at(iRPCHit));
        }
      }
    }
  }
  for(int iMuon = 0; iMuon < nMuon-1; iMuon++){
    for(int jMuon = iMuon+1; jMuon < nMuon; jMuon++){
      double distOffEta = OfflineEta->at(iMuon)-OfflineEta->at(jMuon);
      double distOffPhi = OfflinePhi->at(iMuon)-OfflinePhi->at(jMuon);
      h_distOffEta->Fill(distOffEta);
      h_distOffPhi->Fill(distOffPhi);
      h_RoInumDelEta->Fill( distOffEta, 1);
      h_RoInumDelPhi->Fill( distOffPhi,1);
      double distRpcFitInnEta = SARPCFitInnEta.at(iMuon)-SARPCFitInnEta.at(jMuon);
      double distRpcFitInnPhi = SARPCFitInnPhi->at(iMuon)-SARPCFitInnPhi->at(jMuon);
      h_distRpcFitInnEta->Fill(distRpcFitInnEta);
      h_distRpcFitInnPhi->Fill(distRpcFitInnPhi);
      double distRpcFitMidEta = SARPCFitMidEta.at(iMuon)-SARPCFitMidEta.at(jMuon);
      double distRpcFitMidPhi = SARPCFitMidPhi->at(iMuon)-SARPCFitMidPhi->at(jMuon);
      h_distRpcFitMidEta->Fill(distRpcFitMidEta);
      h_distRpcFitMidPhi->Fill(distRpcFitMidPhi);
      double distRpcFitOutEta = SARPCFitOutEta.at(iMuon)-SARPCFitOutEta.at(jMuon);
      double distRpcFitOutPhi = SARPCFitOutPhi->at(iMuon)-SARPCFitOutPhi->at(jMuon);
      h_distRpcFitOutEta->Fill(distRpcFitOutEta);
      h_distRpcFitOutPhi->Fill(distRpcFitOutPhi);
    }
  }
}


void RPC_FCBM::getCloseBy(auto& param, auto& selectedParam){
  for(UInt_t iMuon = 0; iMuon < nMuon; iMuon++){
    if(isCloseBy(iMuon)){
      selectedParam.push_back(param->at(iMuon));
    }
  }
}
