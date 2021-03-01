#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TVector3.h"
#include <math.h>
#include "../RpcL2MuonSA/RpcDispPreparator.h"
using namespace std;

void RPC_FCBM::FillClusterInfo(){
  setRPCFitEtaR();
  vector<vector<float>> fitMidSlope;
  vector<vector<float>> fitMidOffset;
  fitMidSlope.assign(nMuon, vector<float>(1, 0));
  fitMidOffset.assign(nMuon, vector<float>(1, 0));
  fitMidSlope.clear();
  fitMidOffset.clear();
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    //=======prepare for the loop===================== 
    //set the eta-phi info of OfflineSegment
    vector<float> OffsegMidEta, OffsegMidPhi, OffsegOutEta, OffsegOutPhi;
    vector<float> OffsegMidR, OffsegMidZ, OffsegOutR, OffsegOutZ;
    OffsegMidEta.clear();
    OffsegMidPhi.clear();
    OffsegOutEta.clear();
    OffsegOutPhi.clear();
    OffsegMidR.clear();
    OffsegMidZ.clear();
    OffsegOutR.clear();
    OffsegOutZ.clear();
    vector<float> OffsegMidSlope, OffsegOutSlope;
    OffsegMidSlope.clear();
    OffsegOutSlope.clear();

    vector<float> clusFitMidSlope_plau, clusFitMidOffset_plau;
    vector<float> clusFitOutSlope_plau, clusFitOutOffset_plau;
    clusFitMidSlope_plau.clear();
    clusFitMidOffset_plau.clear();
    clusFitOutSlope_plau.clear();
    clusFitOutOffset_plau.clear();

    vector<float> clusFitMidSlope_best, clusFitMidOffset_best;
    vector<float> clusFitOutSlope_best, clusFitOutOffset_best;
    clusFitMidSlope_best.clear();
    clusFitMidOffset_best.clear();
    clusFitOutSlope_best.clear();
    clusFitOutOffset_best.clear();

    vector<int> n_inSetsMid_plau;
    n_inSetsMid_plau.clear();
    vector<int> n_inSetsOut_plau;
    n_inSetsOut_plau.clear();

    int nClus = (SARPCCluster_clusterMeasPhi->at(iMuon)).size();
    int nClus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
    h_nClusFitMid->Fill(nClus_fit);
    vector<float> clusFitMidEta, clusFitOutEta;
    clusFitMidEta.clear();
    clusFitOutEta.clear();

    vector<vector<int>> n_found;
    n_found.clear();
    n_found = SARPCCluster_n_foundClusters->at(iMuon);
    vector<vector<int>> id_clus;
    id_clus.clear();
    id_clus = SARPCCluster_id_clustersInSets->at(iMuon);
    vector<int> assign;
    assign.clear();
    vector<vector<int>> id_clus_plau;
    id_clus_plau.assign(8, assign);
    vector<vector<int>> n_found_plau;
    n_found_plau.assign(8, assign);

    //isPlausibleFitInnMid
    int countClusterRoadMid = 0;
    int countClusterRoadOut = 0;
    cout << "================check the IDs============" << endl;
    cout << "EventNumber=" << EventNumber << endl; 
    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
      cout << "id.at(" << iClus_fit << ") = ";
      for(int iLay = 0; iLay < 8; iLay++) cout << id_clus.at(iLay).at(iClus_fit) << ",";
      cout << endl;
    }
    cout << "isPlausibleFitInnMid = ";
    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++)cout << SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iClus_fit) << ",";
    cout << endl;
    cout << "isPlausibleFitOut = ";
    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++)cout << SARPCCluster_isPlausibleFitOut->at(iMuon).at(iClus_fit) << ",";
    cout << endl;
    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
      int Nclus = 0;
      for(int ilayer = 0; ilayer < 8; ilayer++){
        if(id_clus.at(ilayer).at(iClus_fit) > -1) Nclus++;
      }
      //if(Nclus < 3) continue;
      if(SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iClus_fit)){
        countClusterRoadMid++;
        clusFitMidSlope_plau.push_back(SARPCCluster_fitMidSlope->at(iMuon).at(iClus_fit));
        clusFitMidOffset_plau.push_back(SARPCCluster_fitMidOffset->at(iMuon).at(iClus_fit));
        for(int i = 0; i < 8; i++) id_clus_plau.at(i).push_back(id_clus.at(i).at(iClus_fit));
        for(int i = 0; i < 8; i++) n_found_plau.at(i).push_back(n_found.at(i).at(iClus_fit));
      }
      if(SARPCCluster_isPlausibleFitOut->at(iMuon).at(iClus_fit)){
        countClusterRoadOut++;
        clusFitOutSlope_plau.push_back(SARPCCluster_fitOutSlope->at(iMuon).at(iClus_fit));
        clusFitOutOffset_plau.push_back(SARPCCluster_fitOutOffset->at(iMuon).at(iClus_fit));
      }
    }//iClus_fit loop
    int N_spinn = 0;
    int N_spmid = 0;
    int N_spout = 0;
    for(int iSPinn=0; iSPinn<SASPCZ_BI->at(0).size(); iSPinn++){
      if( std::fabs(SASPCZ_BI->at(0).at(iSPinn)) < 1e-05 && std::fabs(SASPCR_BI->at(0).at(iSPinn)) < 1e-05)continue;
      N_spinn++;
    }
    for(int iSPmid=0; iSPmid<SASPCZ_BM->at(0).size(); iSPmid++){
      if( std::fabs(SASPCZ_BM->at(0).at(iSPmid)) < 1e-05 && std::fabs(SASPCR_BM->at(0).at(iSPmid)) < 1e-05)continue;
      N_spmid++;
    }
    for(int iSPout=0; iSPout<SASPCZ_BO->at(0).size(); iSPout++){
      if( std::fabs(SASPCZ_BO->at(0).at(iSPout)) < 1e-05 && std::fabs(SASPCR_BO->at(0).at(iSPout)) < 1e-05)continue;
      N_spout++;
    }

    for(int i_spinn = 0; i_spinn<(int)SASPCZ_BI->at(0).size(); i_spinn++){
      double spZ_inn = SASPCZ_BI->at(0).at(i_spinn);
      double spR_inn = SASPCR_BI->at(0).at(i_spinn);
      double spSlope = spZ_inn/spR_inn;

      double dSlope_min = 999;
      for(int i_spmid = 0; i_spmid<(int)SASPCZ_BM->at(0).size(); i_spmid++){
        double spZ_mid = SASPCZ_BM->at(0).at(i_spmid);
        double spR_mid = SASPCR_BM->at(0).at(i_spmid);
        double spSlope_mid = spZ_mid/spR_mid;
        double dSlope = std::fabs(spSlope - spSlope_mid);
        if(dSlope < dSlope_min) dSlope_min = dSlope;
      }
      if(std::fabs(dSlope_min - 999) > 1e-5) m_h_dSlope_innmid->Fill(dSlope_min);
    }

    m_h_countPtclus->Fill(SAptclus->at(0).size());
    h_countSPInn->Fill(N_spinn);
    h_countSPMid->Fill(N_spmid);
    h_countSPOut->Fill(N_spout);

    h_countClusterRoadMid->Fill(countClusterRoadMid);
    h_countClusterRoadOut->Fill(countClusterRoadOut);
    hh_num_clusRoadMidvsOut->Fill(countClusterRoadMid, countClusterRoadOut);
    hh_n_CRvsSPCMid->Fill(countClusterRoadMid, N_spmid);
    if(countClusterRoadMid == 0 && countClusterRoadOut == 0){
      cout << "=========== No cluterRoad are generated event ========" << endl;
      cout << "Offline Ext eta, phi = " << OfflineExtEta->at(iMuon) << "," << OfflineExtPhi->at(iMuon) << endl;
      cout << "======================================================" << endl;
    }
    //calc OfflineSegment eta, phi and select only the information of middle OfflineSegment
    int countMidOffSeg = 0;
    int countOutOffSeg = 0;
    for(int iSeg = 0; iSeg < OfflineNumSegment->at(iMuon); iSeg++){
      TVector3 vecOffseg;
      vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
      TVector3 vecOffsegDir;
      vecOffsegDir.SetXYZ((OfflineSegmentPx->at(iMuon)).at(iSeg), (OfflineSegmentPy->at(iMuon)).at(iSeg), (OfflineSegmentPz->at(iMuon)).at(iSeg));
      h_OffsegEta->Fill(vecOffseg.Eta());

      if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
        OffsegMidEta.push_back(vecOffseg.Eta());
        OffsegMidPhi.push_back(vecOffseg.Phi());
        OffsegMidR.push_back(vecOffseg.Perp());
        OffsegMidZ.push_back(vecOffseg.Z());
        float Slope = vecOffsegDir.Perp()/(OfflineSegmentPz->at(iMuon)).at(iSeg);
        //cout << "OffsegMidSlope = " << Slope << endl;
        //cout << "OffsegMidEta/Phi = " << vecOffseg.Eta() << "/" << vecOffseg.Phi() << endl;
        OffsegMidSlope.push_back(Slope);
        countMidOffSeg++;
      }
      else if(vecOffseg.Perp() >= 8500.){
        OffsegOutEta.push_back(vecOffseg.Eta());
        OffsegOutPhi.push_back(vecOffseg.Phi());
        OffsegOutR.push_back(vecOffseg.Perp());
        OffsegOutZ.push_back(vecOffseg.Z());
        float Slope = vecOffsegDir.Perp()/(OfflineSegmentPz->at(iMuon)).at(iSeg);
        //cout << "OffsegOutSlope = " << Slope << endl;
        //cout << "OffsegOutEta/Phi = " << vecOffseg.Eta() << "/" << vecOffseg.Phi() << endl;
        OffsegOutSlope.push_back(Slope);
        countOutOffSeg++;
      }

    }//iSeg Loop 
    h_nOffSegMid->Fill(countMidOffSeg);  //middle offline segment number
    h_nOffSegOut->Fill(countOutOffSeg);  // the number of outer offline segment

    //the distribution of dSlope(offsegMid - clusterRoadMid), the condition is crossing the clusterRoads at RPC2
    if(countMidOffSeg == 2 && clusFitMidSlope_plau.size() == 2){
      bool crossEvent = true;
      for(int iLayer = 0; iLayer < 4; iLayer++){
        if(iLayer < 2){
          if(id_clus_plau.at(iLayer).at(0) == id_clus_plau.at(iLayer).at(1)){
            crossEvent = false;
          }
        }
        else if(iLayer > 1){
          if(id_clus_plau.at(iLayer).at(0) != id_clus_plau.at(iLayer).at(1)){
            crossEvent = false;
          }
        }
      }
      if(!crossEvent){
        double dR_min = 999;
        double id_off = -1;
        double id_clus = -1;
        //matching with minimum dSLope pair of Offline and clusterRoad
        for(int iOff = 0; iOff < 2; iOff++){
          for(int iClus = 0; iClus < 2; iClus++){
            double dR = fabs(atan(OffsegMidSlope.at(iOff)) - atan(clusFitMidSlope_plau.at(iClus)));
            if(dR < dR_min){
              dR_min = dR;
              id_off = iOff;
              id_clus = iClus;
            }
          }
        }
        cout << "this event is crossing clusterRoad at RPC2 plane, dR_min/id_off/id_clus = " << 
        dR_min << " / " << id_off << " / " << id_clus << endl;
        int id_suboff = (id_off == 0) ? 1:0;
        int id_subclus = (id_clus == 0) ? 1:0;

        double dR_sub = atan(OffsegMidSlope.at(id_suboff)) - atan(clusFitMidSlope_plau.at(id_subclus));
        double dR = atan(OffsegMidSlope.at(id_off)) - atan(clusFitMidSlope_plau.at(id_clus));
        cout << "the other dSlope = " << dR_sub << endl;
        h_dtheta_sub->Fill(dR_sub);
        h_dtheta->Fill(dR);
        h_dtheta_total->Fill(dR_sub);
        h_dtheta_total->Fill(dR);
      }
    }
          

//    for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
//      RpcDispPreparator clustersFit;
//      int sAddress = SAsAddress->at(iMuon);
//      float rpcFitInn[3] = {-9999., -9999., -9999.};    //{Phi, Slope, Offset}
//    float rpcFitMid[3] = {-9999., -9999., -9999.};
//    float rpcFitOut[3] = {-9999., -9999., -9999.};
//    float rpcFitInnEta = -999;
//    float rpcFitMidEta = -999;
//    float rpcFitOutEta = -999;
//    clustersFit.SetParameter_rpcFit(
//        1, 
//        rpcFitInn, 
//        rpcFitMid, 
//        rpcFitOut,
//        (SARPCCluster_fitMidPhi->at(iMuon)).at(0),  
//        (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit),
//        (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit),
//        (SARPCCluster_fitOutPhi->at(iMuon)).at(0),
//        (SARPCCluster_fitOutSlope->at(iMuon)).at(iClus_fit),
//        (SARPCCluster_fitOutOffset->at(iMuon)).at(iClus_fit)
//        );
//    clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);
//
//    h_ClusFitMidEta->Fill(rpcFitMidEta);
//    h_ClusFitOutEta->Fill(rpcFitOutEta);
//    clusFitMidEta.push_back(rpcFitMidEta); 
//    clusFitOutEta.push_back(rpcFitOutEta);
//
//    }//iClus_fit loop end
    h_n_clusRoadMid_total->Fill(nClus_fit);

    //=======middle station===========
    int nOff = OffsegMidEta.size();
    int nOff_out = OffsegOutEta.size();
    for(int iOff = 0; iOff < nOff; iOff++){

      //======about the clusterPosition info=======
      //===========================================
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
    }//iOff loop end


  }// iMuon loop end
}

