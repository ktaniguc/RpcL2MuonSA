#define RPC_FCBM_cxx
//#include "../RpcL2MuonSA/RPC.h"
#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "../RpcL2MuonSA/RpcDispPreparator.h"
#include "../RpcL2MuonSA/RpcClusteringTool.h"
#include "../RpcL2MuonSA/HistData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TFile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <typeinfo>

using namespace std;
const int ZERO_LIMIT = 1e-5;
const double pi = 3.14159265358979323846;

double Res(double param1, double param2);
bool GRLlist(int LumiBlock);
//==================================================================
//main function (argv[1]: PDF_LABEL, argv[2]: INPUT_NTUPLE, argv[3]: IS_DRAW, arg[4]: IS_EVENTDISPLAY, arg[5]: BEGIN_ENTRY, arg[6]: LIMIT_ENTRY, arg[7]: TAP_TYPE, arg[8]: trig_chain, argv[9]: IS_CloseByMuon)
//==================================================================

void RPC_FCBM::Loop_FCBM( int Nevents, RpcL2MuonSA::HistData histData)
{
  if (fChain == 0) return;

  Long64_t nLoop;
  if (Nevents == -1) {
    nLoop = fChain -> GetEntries();
  } else {
    nLoop = Nevents;
  }

  //Long64_t nentries = fChain->GetEntriesFast();
  double entries = fChain->GetEntries();
  cout << "[INFO]: Nentries:" << entries << endl;
  //const int N50 = 14;


  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nLoop;jentry++) {
    int  ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (ientry < 0) break;
    //==================================================================
    //Analysis code for a entry
    //==================================================================

    // Check GRL
    if (RunNumber == 349533){
      if (GRLlist(LumiBlock)){
        continue;
      }
    }
    if(is2muIn1RoI()){
      FillEfficiency();
    }
  } // end of each entry

} // end of RPC::Loop()


void RPC_FCBM::Display_FCBM(int tap_type, int trig_chain, Long64_t begin_entry, Long64_t limit_entry, TString pdf){
  // Prepare Canvas
  TCanvas *c2 = new TCanvas("c2", "c2", 10, 10, 1050, 700);
  c2->SetGrid();
  c2->SetRightMargin(0.20);
  c2->SetLeftMargin(0.23);
  c2->SetBottomMargin(0.20);
  c2->Print(pdf + "[", "pdf");


  if (limit_entry == -1) {
    begin_entry = 0;
    limit_entry = fChain -> GetEntries();
    cout << "-------" << endl;
  }
  // Prepare Loop
  if (fChain == 0) return;
  int nLoop = fChain -> GetEntries();
  //Long64_t nentries = fChain->GetEntriesFast();
  double entries = fChain->GetEntries();
  cout << "[INFO]: Nentries:" << entries << endl;
  Long64_t nbytes = 0, nb = 0;

  int current_entry = 0;
  // Start Loop
  for (Long64_t jentry = begin_entry; jentry<nLoop;jentry++) {
    int  ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (ientry < 0) break;
    // Check GRL
    if (RunNumber == 349533){
      if (GRLlist(LumiBlock)){
        continue;
      }
    }

    //ktaniguc arrange
    bool barrelFlag = true;
    bool rpcFitFlag = true;
    bool OffsegFlag = true;
    bool isOffsegMidFlag = true;
    bool OfflinePtFlag = true;
    bool OfflineZ0Flag = true;
    bool isClusFitMidFlag = true;
    bool TruthLikeRoadFlag = true;
    bool isSingleRoI = true;
    bool ptresFlag = false;

    vector<vector<double>> OffsegMid;
    vector<double> forClear;
    forClear.clear();
    OffsegMid.assign(nMuon, forClear);

    for(int iMuon = 0; iMuon < nMuon; iMuon++){
      if(L1nRoI->at(iMuon) != 1 ) isSingleRoI = false;
      //if(SAsAddress->at(iMuon) == -1 || SAsAddress->at(iMuon) > 4) barrelFlag == false;
      if(fabs(SARoIEta->at(iMuon)) >= 1.05) barrelFlag = false;
      if(fabs(SARPCFitMidSlope->at(iMuon)) <= ZERO_LIMIT) rpcFitFlag = false;
      
      int nSeg = OfflineNumSegment->at(iMuon);
      int countOffsegMid = 0;
      vector<float> OffsegMidEta, OffsegMidSlope;
      OffsegMidEta.clear();
      OffsegMidSlope.clear();
      for(int iSeg = 0; iSeg < nSeg; iSeg++){
        TVector3 vecOffseg;
        vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
        TVector3 vecOffsegDir;
        vecOffsegDir.SetXYZ((OfflineSegmentPx->at(iMuon)).at(iSeg), (OfflineSegmentPy->at(iMuon)).at(iSeg), (OfflineSegmentPz->at(iMuon)).at(iSeg));
        if(fabs(vecOffseg.Eta()) >= 1.05) OffsegFlag = false;
        if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
          countOffsegMid++;
          OffsegMidEta.push_back(vecOffseg.Eta());
          OffsegMid.at(iMuon).push_back(vecOffseg.Perp());
          float Slope = vecOffsegDir.Perp()/(OfflineSegmentPz->at(iMuon)).at(iSeg);
//          cout << "OffsegMidSlope = " << Slope << endl;
//          cout << "OffsegMidEta/Phi = " << vecOffseg.Eta() << "/" << vecOffseg.Phi() << endl;
          OffsegMidSlope.push_back(Slope);
        }

      }
      vector<float> clusFitMidEta;
      int nclus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
      clusFitMidEta.clear();
      for(int iClus_fit = 0; iClus_fit < nclus_fit; iClus_fit++){
        RpcDispPreparator clustersFit;
        int sAddress = SAsAddress->at(iMuon);
        float rpcFitInn[3] = {-9999., -9999., -9999.};    //{Phi, Slope, Offset}
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
            0, //(SARPCCluster_fitMidPhi->at(iMuon)).at(iClus_fit),  
            (SARPCCluster_fitMidSlope->at(iMuon)).at(iClus_fit),
            (SARPCCluster_fitMidOffset->at(iMuon)).at(iClus_fit),
            0, //(SARPCCluster_fitOutPhi->at(iMuon)).at(iClus_fit),
            (SARPCCluster_fitOutSlope->at(iMuon)).at(iClus_fit),
            (SARPCCluster_fitOutOffset->at(iMuon)).at(iClus_fit)
            );
        clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);

        clusFitMidEta.push_back(rpcFitMidEta); 

      }//iClus_fit loop end
      if(countOffsegMid < 1) isOffsegMidFlag = false; 
      if(OfflinePt->at(iMuon) < 4) OfflinePtFlag=false;
      if(fabs(OfflineZ0->at(iMuon)) > 50)OfflineZ0Flag = false;

      if((SARPCCluster_fitMidSlope->at(iMuon)).size() == 0) isClusFitMidFlag = false;
      //////////////////////////////////////////////////////////////
      //TruthLike-selection
      //////////////////////////////////////////////////////
      int nClus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
      int nOff = OffsegMidEta.size();
    }
    vector< pair<int, int> > mupair; // this value is used only for close by cut
    if(!(nMuon == 2)) continue;
    double ptlead = -999999;
    double ptsub = 999999;

    for(int iMuon = 0; iMuon< nMuon; iMuon++){
      if(ptlead < OfflinePt->at(iMuon)) ptlead = OfflinePt->at(iMuon);
      if(ptsub > OfflinePt->at(iMuon)) ptsub = OfflinePt->at(iMuon);
    }
    
    if(SAptclus->at(0).size() == 2){
      double saptlead = -999999;
      double saptsub = 9999999;
      for(int isa =0 ; isa < 2; isa++){
        if(saptlead < fabs(SAptclus->at(0).at(isa))) saptlead = fabs(SAptclus->at(0).at(isa));
        if(saptsub > fabs(SAptclus->at(0).at(isa))) saptsub = fabs(SAptclus->at(0).at(isa));
      }
      double ptres_leading = 1 - saptlead/ptlead;
      double ptres_sub = 1 - saptsub/ptsub;
      if(ptres_leading > 0.7 || ptres_sub > 0.7) ptresFlag = true;
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
//    if(ptlead > 30) continue;
//
    //if(!ptresFlag) continue;
    if(!isCloseByCut(nMuon, mupair)) continue;
    if(!barrelFlag) continue;
    if(!rpcFitFlag) continue;
    if(!OffsegFlag) continue;
//    if(!isOffsegMidFlag) continue;
    if(!OffsegQualityMid) continue;
    if(!OfflinePtFlag)continue;
    if(!OfflineZ0Flag) continue;
    if(!isSingleRoI) continue;
    cout << "EventNumber = " << EventNumber << endl;
    double deltaExtEta = OfflineExtEta->at(0) - OfflineExtEta->at(1);
    double deltaExtPhi = OfflineExtPhi->at(0) - OfflineExtPhi->at(1);
    double deltaExtR_off = sqrt(pow(deltaExtEta, 2) + pow(deltaExtPhi, 2));
    unsigned int nClus_fit = (SARPCCluster_fitMidSlope->at(1)).size();
    int N_clusRoadInn = 0;
    int N_clusRoadMid = 0;
    int N_clusRoadOut = 0;
    if(isOld){
      for(unsigned int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
        if(SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit))N_clusRoadInn++;
        if(SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit))N_clusRoadMid++;
        if(SARPCCluster_isPlausibleFitOut->at(0).at(iClus_fit))N_clusRoadOut++;
      }
    } else {
      N_clusRoadInn = SARPCCluster_fitInnSlope->at(0).size();
      N_clusRoadMid = SARPCCluster_fitMidSlope->at(0).size();
      N_clusRoadOut = SARPCCluster_fitOutSlope->at(0).size();
    }

    for(int iMuon = 0; iMuon < nMuon; iMuon++){
      if(isOld){
        cout << "the plausible clusterRoad flag = {";
        for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
          cout << SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iClus_fit) << ",";

        }
        cout << "}" << endl;
      }
      for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
        cout << "id.at(" << iClus_fit << ") = ";
        for(int iLay = 0; iLay < 8; iLay++) cout << SARPCCluster_id_clustersInSets->at(iMuon).at(iLay).at(iClus_fit) << ",";
        cout << endl;
      }
      int nClus = SARPCCluster_gX->at(iMuon).size();
      for(int i_layer = 0; i_layer < 8; i_layer++){
        for(int iClus = 0; iClus < nClus; iClus++){
          if(SARPCCluster_clusterMeasPhi->at(iMuon).at(iClus)) continue;
          if(i_layer != SARPCCluster_clusterLayer->at(iMuon).at(iClus))  continue;
          TVector3 vecClusPos;
          std::cout << "layer " << i_layer << "...";
          vecClusPos.SetXYZ(SARPCCluster_gX->at(iMuon).at(iClus), SARPCCluster_gY->at(iMuon).at(iClus), SARPCCluster_gZ->at(iMuon).at(iClus));
          cout << "{" << vecClusPos.Perp() << "," << SARPCCluster_gZ->at(iMuon).at(iClus) << "} , ";
        }
        std::cout << std::endl;
      }

    }
    //if(!isClusFitMidFlag) continue;
    //if(!TruthLikeRoadFlag) continue;
    //=======================================================
    if(isOld){
      vector<vector<int>> n_found;
      n_found = SARPCCluster_n_foundClusters->at(0);
    }
    vector<vector<int>> id_clus;
    id_clus = SARPCCluster_id_clustersInSets->at(0);
    setRPCFitEtaR();
    //Offline segment
    TGraph gr_OffsegMid[nMuon] = TGraph(0);
    TGraph gr_OffsegOut[nMuon] = TGraph(0);
    TGraphErrors gr_rpcHitEtaPhi_LayBL[nMuon][8];
    TGraphErrors gr_rpcHitEtaPhi_LayBS[nMuon][8];
    TGraphErrors gr_rpcHitEtaPhi_LayBSP[nMuon][8];
    TGraph gr_rpc1_FitEtaPhi[nMuon] = TGraph();
    TGraph gr_rpc2_FitEtaPhi[nMuon] = TGraph();
    TGraph gr_rpc3_FitEtaPhi[nMuon] = TGraph();
    TGraph gr_roi[nMuon] = TGraph(0);
    TGraph gr_clusFit_rpc1[nMuon] = TGraph(0);
    TGraph gr_clusFit_rpc2[nMuon] = TGraph(0);
    TGraph gr_clusFit_rpc3[nMuon] = TGraph(0);

    TString label_for_sector[nMuon];

    //RPC cluster
    TGraphErrors gr_ClusRPC[nMuon][8];      
    /*   TGraphErrors gr_ClusRPC2[nMuon] = TGraphErrors(0);      
	 TGraphErrors gr_ClusRPC3[nMuon] = TGraphErrors(0);     */

    for(int iMuon = 0; iMuon < nMuon; iMuon++){
      if (SAsAddress->at(iMuon) == 0){
	label_for_sector[iMuon] = "Large";
      } else if (SAsAddress->at(iMuon) == 1){
	label_for_sector[iMuon] = "Small";
      } else if (SAsAddress->at(iMuon) == 2){
	label_for_sector[iMuon] = "Large Special";
      } else if (SAsAddress->at(iMuon) == 3){
	label_for_sector[iMuon] = "Small Special";
      } else if (SAsAddress->at(iMuon) == -1){
	label_for_sector[iMuon] = "Endcap";
      }

      for(int iSeg = 0; iSeg < OfflineNumSegment->at(iMuon); iSeg++){
	TVector3 vecOffseg;
	vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));

	if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
	  gr_OffsegMid[iMuon].SetPoint(iSeg, vecOffseg.Eta(), vecOffseg.Phi());
	}
	if(9500. < vecOffseg.Perp() && vecOffseg.Perp() < 10500.){
	  gr_OffsegOut[iMuon].SetPoint(iSeg, vecOffseg.Eta(), vecOffseg.Phi());
	}
      }


      //rpc hit
      RpcClusteringTool MUcluster;
      MUcluster.sortStrip(gr_rpcHitEtaPhi_LayBL[iMuon], gr_rpcHitEtaPhi_LayBS[iMuon], gr_rpcHitEtaPhi_LayBSP[iMuon], SARPCHitStationNumber->at(iMuon), SARPCHitMeasuresPhi->at(iMuon), SARPCHitEta->at(iMuon), SARPCHitPhi->at(iMuon), SARPCHitLayer->at(iMuon) );

      RpcClusterSetter(gr_ClusRPC[iMuon], 
	  SARPCCluster_gX->at(iMuon), 
	  SARPCCluster_gY->at(iMuon), 
	  SARPCCluster_gZ->at(iMuon), 
	  SARPCCluster_clusterLayer->at(iMuon), 
	  SARPCCluster_clusterMeasPhi->at(iMuon) );
      std::cout << "hoge, after MUcluster.sortStrip" << std::endl;
      //rpcFit
      if(SAsAddress->at(iMuon) == 0 || SAsAddress->at(iMuon) == 1 || SAsAddress->at(iMuon) == 2 || SAsAddress->at(iMuon) == 3){
	gr_rpc1_FitEtaPhi[iMuon].SetPoint(0, SARPCFitInnEta.at(iMuon), SARPCFitMidPhi->at(iMuon));
	gr_rpc2_FitEtaPhi[iMuon].SetPoint(0, SARPCFitMidEta.at(iMuon), SARPCFitMidPhi->at(iMuon));
	gr_rpc3_FitEtaPhi[iMuon].SetPoint(0, SARPCFitOutEta.at(iMuon), SARPCFitOutPhi->at(iMuon));
	//rpcClusterFit
	if(!(SARPCCluster_fitInnSlope)->empty()&&!(SARPCCluster_fitMidPhi)->empty()){
	  int nclus = (SARPCCluster_fitInnSlope->at(iMuon)).size();
	  for(int iclus = 0; iclus < nclus; iclus++){
	    //if((SARPCCluster_isSuccess->at(iMuon)).at(iclus)) continue;
            if(isOld){
              if( !SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iclus)) continue;
            }
            RpcDispPreparator clustersFit;
	    int sAddress = SAsAddress->at(iMuon);
	    float rpcFitInn[3] = {-9999., -9999., -9999.};    //{Phi, Slope, Offset}
	    float rpcFitMid[3] = {-9999., -9999., -9999.};
	    float rpcFitOut[3] = {-9999., -9999., -9999.};
	    float rpcFitInnEta = 0;
	    float rpcFitMidEta = 0;
	    float rpcFitOutEta = 0;
  //          if(!(SARPCCluster_isUsingMidCluster->at(iMuon)).at(iclus)) continue;
	    clustersFit.SetParameter_rpcFit(//isTag,  
	      1, 
	      rpcFitInn, 
	      rpcFitMid, 
	      rpcFitOut,
	      0, //(SARPCCluster_fitMidPhi->at(iMuon)).at(0),  
	      (SARPCCluster_fitMidSlope->at(iMuon)).at(0),
	      (SARPCCluster_fitMidOffset->at(iMuon)).at(0),
	      0, //(SARPCCluster_fitOutPhi->at(iMuon)).at(iclus),
	      (SARPCCluster_fitOutSlope->at(iMuon)).at(iclus),
	      (SARPCCluster_fitOutOffset->at(iMuon)).at(iclus)
	      );
	    clustersFit.SetEta_rpcFit(sAddress, rpcFitInn, rpcFitMid, rpcFitOut, rpcFitInnEta, rpcFitMidEta, rpcFitOutEta);
	    cout << "test-clusterFitInfo" << endl;
	    cout << "MidEta/Phi=" << rpcFitInnEta << "/" << rpcFitMid[0] << endl;
	    cout << "MidEta/Phi=" << rpcFitMidEta << "/" << rpcFitMid[0] << endl;
	    cout << "OutEta/Phi=" << rpcFitOutEta << "/" << rpcFitOut[0] << endl;
	    gr_clusFit_rpc1[iMuon].SetPoint(iclus, rpcFitInnEta, SARPCFitMidPhi->at(iMuon));
	    gr_clusFit_rpc2[iMuon].SetPoint(iclus, rpcFitMidEta, SARPCFitMidPhi->at(iMuon));
	    gr_clusFit_rpc3[iMuon].SetPoint(iclus, rpcFitOutEta, SARPCFitOutPhi->at(iMuon));

	  }
	}

      }
      //ROI
      gr_roi[iMuon].SetPoint(0, SARoIEta->at(iMuon), SARoIPhi->at(iMuon));
    }//iMuon loop end

    double rpcHitRegion[4] = {0., 0., 0., 0.};  //{EtaMin, PhiMin, EtaMax, PhiMax}
  setRPCHitRegion(nMuon, SARPCHitEta, SARPCHitPhi, rpcHitRegion);

////==================create (R,Z) display==========//  
//==========middle station==============//
//
//  double Zmin_BI = 0;
//  double Rmin_BI = 0;
//  double Zmax_BI = 0;
//  double Rmax_BI = 0;
//  setFrameRegion("inner", Zmin_BI, Rmin_BI, Zmax_BI, Rmax_BI);
//  double Zmin_BM = 0;
//  double Rmin_BM = 0;
//  double Zmax_BM = 0;
//  double Rmax_BM = 0;
//  setFrameRegion("middle", Zmin_BM, Rmin_BM, Zmax_BM, Rmax_BM);
//  double Zmin_BO = 0;
//  double Rmin_BO = 0;
//  double Zmax_BO = 0;
//  double Rmax_BO = 0;
//  setFrameRegion("outer", Zmin_BO, Rmin_BO, Zmax_BO, Rmax_BO);
  float Zmin_BI = 0;
  float Rmin_BI = 0;
  float Zmax_BI = 0;
  float Rmax_BI = 0;
  setFrameRegion("inner", Zmin_BI, Rmin_BI, Zmax_BI, Rmax_BI);
  float Zmin_BM = 0;
  float Rmin_BM = 0;
  float Zmax_BM = 0;
  float Rmax_BM = 0;
  setFrameRegion("middle", Zmin_BM, Rmin_BM, Zmax_BM, Rmax_BM);
  float Zmin_BO = 0;
  float Rmin_BO = 0;
  float Zmax_BO = 0;
  float Rmax_BO = 0;
  setFrameRegion("outer", Zmin_BO, Rmin_BO, Zmax_BO, Rmax_BO);

  //======information for each muons==//
  TGraph MdtRegion_BM[nMuon];
  TGraph MdtRegion_BO[nMuon];
  TF1 f_roi[nMuon];
  TF1 f_roi_minus[nMuon];
  TF1 f_roi_plus[nMuon];
  TGraph gr_segment_RZ[nMuon];
  TGraph gr_segment_RZ_BO[nMuon];
  TGraph gr_RpcHit[nMuon];
  TGraph gr_RpcCluster[nMuon];
  TGraph gr_SPC_BI;
  TGraph gr_SPC_BM;
  TGraph gr_SPC_BO;
  TGraph gr_SPC_BI_2;
  TGraph gr_SPC_BM_2;
  TGraph gr_SPC_BO_2;
  TF1 f_road_BI[nMuon];
  TF1 f_road_BM[nMuon];
  TF1 f_road_BO[nMuon];
  int nMdt = SAMDTHitAllR->at(0).size();
  int nMdt_hit = SAMDTHitR->at(0).size();
  cout << "mdt size rpchit/rpcclus = " << nMdt_hit << "/" << nMdt << endl;
  cout << "mdt size rpchit/rpcclus index1 = " << SAMDTHitR->at(1).size() << "/" << SAMDTHitAllR->at(1).size() << endl;
  TGraph gr_mdt;
  TGraph gr_mdt_outlier;
  TGraph gr_mdt_outlier5;
  TGraph gr_mdt_hit;
  TGraph gr_mdt_hit_outlier;
  for(int iMdt = 0; iMdt < nMdt_hit; iMdt++){
    if(SAMDTHitisOutlier->at(0).at(iMdt) == 0){
      gr_mdt_hit.SetPoint(iMdt, SAMDTHitZ->at(0).at(iMdt)/1000, SAMDTHitR->at(0).at(iMdt)/1000);
    }
    else if(SAMDTHitisOutlier->at(0).at(iMdt) != 0 && SAMDTHitisOutlier->at(0).at(iMdt) != 5){
      gr_mdt_hit_outlier.SetPoint(iMdt, SAMDTHitZ->at(0).at(iMdt)/1000, SAMDTHitR->at(0).at(iMdt)/1000);
    }
  }
  for(int iMdt = 0; iMdt < nMdt; iMdt++){
    if(isOld){
      if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0){
        gr_mdt.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
      }
      else if(SAMDTHitAllisOutlier->at(0).at(iMdt) != 0 && SAMDTHitAllisOutlier->at(0).at(iMdt) != 5){
        gr_mdt_outlier.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
      }
      else if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 5){
        gr_mdt_outlier5.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
      }
    } else {
//        if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0 && SAMDTHitAllclusRoadID->at(0).at(iMdt) == 0) gr_mdt.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
//        else if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0 && SAMDTHitAllclusRoadID->at(0).at(iMdt) == 1) gr_mdt_outlier.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
//        else if(SAMDTHitAllisOutlier->at(0).at(iMdt) > -1) gr_mdt_outlier5.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
        if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0) gr_mdt.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
        else if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 1) gr_mdt_outlier.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
        else if(SAMDTHitAllisOutlier->at(0).at(iMdt) > -1) gr_mdt_outlier5.SetPoint(iMdt, SAMDTHitAllZ->at(0).at(iMdt)/1000, SAMDTHitAllR->at(0).at(iMdt)/1000);
    }
  }
  if(!SASPCZ_BI->at(0).empty()){
    int nSPinn = SASPCZ_BI->at(0).size();
    for(int iSP = 0; iSP < nSPinn; iSP++){
      gr_SPC_BI.SetPoint(iSP, SASPCZ_BI->at(0).at(iSP)/1000, SASPCR_BI->at(0).at(iSP)/1000);
      if(iSP > 0)gr_SPC_BI_2.SetPoint(iSP, SASPCZ_BI->at(0).at(iSP)/1000, SASPCR_BI->at(0).at(iSP)/1000);
    }
  }
  if(!SASPCZ_BM->at(0).empty()){
    int nSPmid = SASPCZ_BM->at(0).size();
    for(int iSP = 0; iSP < nSPmid; iSP++){
      gr_SPC_BM.SetPoint(iSP, SASPCZ_BM->at(0).at(iSP)/1000, SASPCR_BM->at(0).at(iSP)/1000);
      if(iSP > 0)gr_SPC_BM_2.SetPoint(iSP, SASPCZ_BM->at(0).at(iSP)/1000, SASPCR_BM->at(0).at(iSP)/1000);
    }
  }
  if(!SASPCZ_BO->at(0).empty()){
    int nSPout = SASPCZ_BO->at(0).size();
    for(int iSP = 0; iSP < nSPout; iSP++){
      gr_SPC_BO.SetPoint(iSP, SASPCZ_BO->at(0).at(iSP)/1000, SASPCR_BO->at(0).at(iSP)/1000);
      if(iSP > 0)gr_SPC_BO_2.SetPoint(iSP, SASPCZ_BO->at(0).at(iSP)/1000, SASPCR_BO->at(0).at(iSP)/1000);
    }
  }

  vector<vector<double>> OffsegSlope, OffsegOffset, OffsegMin, OffsegMax;
  vector<double> forclear;
  forclear.clear();
  int NMu = nMuon;
  OffsegSlope.assign(NMu, forclear);
  OffsegOffset.assign(NMu, forclear);

  OffsegMin.assign(NMu, forclear);
  OffsegMax.assign(NMu, forclear);

  float r_road[4] = {10, 10, 0, 0};
  float roadAw = SARoadAw->at(0)[1];
  float roadBw = SARoadBw->at(0)[1]/1000;
  float rWidth = sqrt(roadAw*roadAw+1)*0.2;
  float z_road_0 = (10-roadBw+rWidth)/roadAw;
  float z_road_1 = (10-roadBw-rWidth)/roadAw;
  float z_road_2 = (0-roadBw-rWidth)/roadAw;
  float z_road_3 = (0-roadBw+rWidth)/roadAw;
  float z_road[4] = {z_road_0, z_road_1, z_road_2, z_road_3};
  std::cout << "z_road = {" << z_road_0 << "/" << z_road_1 << "/" << z_road_2 << "/" << z_road_3 << "}" << std::endl;
  TGraph roadRegion;
  roadRegion = TGraph(4, z_road, r_road);
  roadRegion.SetFillStyle(3002);
  roadRegion.SetFillColorAlpha(9, 0.30);
  
  for(unsigned int iMuon = 0; iMuon < nMuon; iMuon++){
    float Z_mdt[4] = {0, 0, 0, 0};
    float R_mdt[4] = {0, 0, 0, 0};
    float Z_mdt_BO[4] = {0, 0, 0, 0};
    float R_mdt_BO[4] = {0, 0, 0, 0};
    setMdtRegion(iMuon, "middle", Z_mdt, R_mdt);
    MdtRegion_BM[iMuon] = TGraph(4, Z_mdt, R_mdt);
    MdtRegion_BM[iMuon].SetFillStyle(3002);
    MdtRegion_BM[iMuon].SetFillColorAlpha(13, 0.40);
    setMdtRegion(iMuon, "outer", Z_mdt_BO, R_mdt_BO);
    MdtRegion_BO[iMuon] = TGraph(4, Z_mdt_BO, R_mdt_BO);
    MdtRegion_BO[iMuon].SetFillStyle(3002);
    MdtRegion_BO[iMuon].SetFillColorAlpha(13, 0.40);
    
    f_roi[iMuon] = TF1("f_roi", "[0]*x", -20, 20);
    f_roi[iMuon].SetTitle(";Z [m];R [m]");
    f_roi[iMuon].SetParameter(0, tan((2*atan(exp(-SARoIEta->at(iMuon))))));
    f_roi[iMuon].SetLineColor(kYellow+3-iMuon);
    f_roi[iMuon].SetLineWidth(1);
    f_roi[iMuon].SetLineStyle(9);
    f_roi_minus[iMuon] = TF1("f_roi_minus", "[0]*x", -20, 20);
    f_roi_minus[iMuon].SetTitle(";Z [m];R [m]");
    f_roi_minus[iMuon].SetParameter(0, tan((2*atan(exp(-SARoIEta->at(iMuon)+0.05)))));
    f_roi_minus[iMuon].SetLineColor(kYellow+3-iMuon);
    f_roi_minus[iMuon].SetLineWidth(1);
    f_roi_minus[iMuon].SetLineStyle(2);
    f_roi_plus[iMuon] = TF1("f_roi_plus", "[0]*x", -20, 20);
    f_roi_plus[iMuon].SetTitle(";Z [m];R [m]");
    f_roi_plus[iMuon].SetParameter(0, tan((2*atan(exp(-SARoIEta->at(iMuon)-0.05)))));
    f_roi_plus[iMuon].SetLineColor(kYellow+3-iMuon);
    f_roi_plus[iMuon].SetLineWidth(1);
    f_roi_plus[iMuon].SetLineStyle(2);
    
    gr_segment_RZ[iMuon] = TGraph(0);
    gr_segment_RZ_BO[iMuon] = TGraph(0);
    for(unsigned int iSeg = 0; iSeg < OfflineNumSegment->at(iMuon); iSeg++){
      double OffsegR = TMath::Sqrt(pow((OfflineSegmentX->at(iMuon)).at(iSeg), 2) + pow((OfflineSegmentY->at(iMuon)).at(iSeg), 2))/1000.;
      TVector3 vecOffsegDir;
      vecOffsegDir.SetXYZ((OfflineSegmentPx->at(iMuon)).at(iSeg), (OfflineSegmentPy->at(iMuon)).at(iSeg), (OfflineSegmentPz->at(iMuon)).at(iSeg));
      double slope = vecOffsegDir.Perp()/OfflineSegmentPz->at(iMuon).at(iSeg);
      OffsegSlope.at(iMuon).push_back(slope);

      gr_segment_RZ[iMuon].SetPoint(iSeg, (OfflineSegmentZ->at(iMuon)).at(iSeg)/1000., OffsegR);
      gr_segment_RZ_BO[iMuon].SetPoint(iSeg, (OfflineSegmentZ->at(iMuon)).at(iSeg)/1000., OffsegR);
      double offset = OffsegR*1000 - slope*(OfflineSegmentZ->at(iMuon).at(iSeg));
      OffsegOffset.at(iMuon).push_back(offset);
      OffsegMin.at(iMuon).push_back(OfflineSegmentZ->at(iMuon).at(iSeg)/1000 -0.3);
      OffsegMax.at(iMuon).push_back(OfflineSegmentZ->at(iMuon).at(iSeg)/1000 +0.3);
    }
    
    gr_segment_RZ[iMuon].SetMarkerStyle(21);
    gr_segment_RZ[iMuon].SetMarkerSize(2-iMuon*0.5);
    gr_segment_RZ[iMuon].SetMarkerColor(kPink+iMuon);
    gr_segment_RZ[iMuon].GetXaxis()->SetLimits(Zmin_BM,Zmax_BM);
    gr_segment_RZ_BO[iMuon].SetTitle(";Z [m];R [m]");
    gr_segment_RZ_BO[iMuon].SetMarkerStyle(21);
    gr_segment_RZ_BO[iMuon].SetMarkerSize(2-iMuon*0.5);
    gr_segment_RZ_BO[iMuon].SetMarkerColor(kPink+iMuon);
    gr_segment_RZ_BO[iMuon].GetXaxis()->SetLimits(Zmin_BO,Zmax_BO);
    gr_segment_RZ_BO[iMuon].SetTitle(";Z [m];R [m]");

    int nHit = (SARPCHitR->at(iMuon)).size();
    gr_RpcHit[iMuon] = TGraph(nHit);
    for(int iHit = 0; iHit < nHit; iHit++){
      if(SARPCHitStationNumber->at(iMuon).at(iHit) > 0){
        if(SARPCHitMeasuresPhi->at(iMuon).at(iHit)) continue;
        gr_RpcHit[iMuon].SetPoint(iHit, SARPCHitZ->at(iMuon).at(iHit)/1000, SARPCHitR->at(iMuon).at(iHit)/1000);
      }
    }//iHit loop end
    gr_RpcHit[iMuon].SetMarkerStyle(8);
    gr_RpcHit[iMuon].SetMarkerSize(0.85);
    gr_RpcHit[iMuon].SetMarkerColor(kGreen+1);
    //gr_RpcHit[iMuon].SetMarkerColor(kRed);
    
    int nClus = (SARPCCluster_gX->at(iMuon)).size();
    gr_RpcCluster[iMuon] = TGraph(nClus);
    for(int iClus = 0; iClus < nClus; iClus++){
      TVector3 vecClus;
      if((SARPCCluster_clusterMeasPhi->at(iMuon)).at(iClus)) continue;
      vecClus.SetXYZ(SARPCCluster_gX->at(iMuon).at(iClus), SARPCCluster_gY->at(iMuon).at(iClus), SARPCCluster_gZ->at(iMuon).at(iClus) );
      gr_RpcCluster[iMuon].SetPoint(iClus, SARPCCluster_gZ->at(iMuon).at(iClus)/1000, vecClus.Perp()/1000);
    }
    gr_RpcCluster[iMuon].SetMarkerStyle(8);
    gr_RpcCluster[iMuon].SetMarkerSize(1);
    //gr_RpcCluster[iMuon].SetMarkerSize(0.75);
    gr_RpcCluster[iMuon].SetMarkerColor(kBlack);


    f_road_BI[iMuon] = TF1("f_road_BI", "[0]*x+[1]", -9999, 9999);
    f_road_BI[iMuon].SetTitle(";Z [m];R [m]");
    f_road_BI[iMuon].SetParameter(0,SARoadAw -> at(iMuon)[0]);
    f_road_BI[iMuon].SetParameter(1,SARoadBw -> at(iMuon)[0]/1000.);
    f_road_BI[iMuon].SetLineColor(kBlue+1);
    f_road_BI[iMuon].SetLineWidth(1);
    f_road_BI[iMuon].SetLineStyle(2);
    f_road_BM[iMuon] = TF1("f_road_BM", "[0]*x+[1]", -9999, 9999);
    f_road_BM[iMuon].SetTitle(";Z [m];R [m]");
    f_road_BM[iMuon].SetParameter(0,SARoadAw -> at(iMuon)[1]);
    f_road_BM[iMuon].SetParameter(1,SARoadBw -> at(iMuon)[1]/1000.);
    f_road_BM[iMuon].SetLineColor(kBlue+1);
    f_road_BM[iMuon].SetLineWidth(1);
    f_road_BM[iMuon].SetLineStyle(2);
    f_road_BO[iMuon] = TF1("f_road_BO", "[0]*x+[1]", -9999, 9999);
    f_road_BO[iMuon].SetTitle(";Z [m];R [m]");
    f_road_BO[iMuon].SetParameter(0,SARoadAw -> at(iMuon)[2]);
    f_road_BO[iMuon].SetParameter(1,SARoadBw -> at(iMuon)[2]/1000.);
    f_road_BO[iMuon].SetLineColor(kBlue+1);
    f_road_BO[iMuon].SetLineWidth(1);
    f_road_BO[iMuon].SetLineStyle(2);
      

  } //iMuon loop end 
  vector<vector<float>> fitInnSlope, fitInnOffset;
  vector<vector<float>> fitMidSlope, fitMidOffset;
  vector<vector<float>> fitOutSlope, fitOutOffset;
  
  vector<float> forClearFloat;
  forClearFloat.clear();
  fitInnSlope.assign(nMuon, forClearFloat);
  fitInnOffset.assign(nMuon, forClearFloat);
  fitMidSlope.assign(nMuon, forClearFloat);
  fitMidOffset.assign(nMuon, forClearFloat);
  fitOutSlope.assign(nMuon, forClearFloat);
  fitOutOffset.assign(nMuon, forClearFloat);
 
  for(unsigned int iMuon = 0; iMuon < nMuon; iMuon++){ 
    int nClus_fitInn = (SARPCCluster_fitInnSlope->at(iMuon)).size();
    for(unsigned int iClus_fit = 0; iClus_fit < nClus_fitInn; iClus_fit++){
      if(isOld){
        if(!SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iClus_fit) ) continue;
      }
      fitInnSlope.at(iMuon).push_back(SARPCCluster_fitInnSlope->at(iMuon).at(iClus_fit));
      cout << "plausibleFitInnSlope = " << SARPCCluster_fitInnSlope->at(iMuon).at(iClus_fit) << endl;
      fitInnOffset.at(iMuon).push_back(SARPCCluster_fitInnOffset->at(iMuon).at(iClus_fit));
      cout << "plausibleFitInnOffset = " << SARPCCluster_fitInnOffset->at(iMuon).at(iClus_fit) << endl;
    }
    int nClus_fit = (SARPCCluster_fitMidSlope->at(iMuon)).size();
    for(unsigned int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
      if(isOld){
        if(!SARPCCluster_isPlausibleFitInnMid->at(iMuon).at(iClus_fit)) continue;
      }
      fitMidSlope.at(iMuon).push_back(SARPCCluster_fitMidSlope->at(iMuon).at(iClus_fit));
      cout << "plausibleFitMidSlope = " << SARPCCluster_fitMidSlope->at(iMuon).at(iClus_fit) << endl;
      fitMidOffset.at(iMuon).push_back(SARPCCluster_fitMidOffset->at(iMuon).at(iClus_fit));
      cout << "plausibleFitMidOffset = " << SARPCCluster_fitMidOffset->at(iMuon).at(iClus_fit) << endl;
    }
    int nClus_fitOut = (SARPCCluster_fitOutSlope->at(iMuon)).size();
    for(unsigned int iClus_fit = 0; iClus_fit < nClus_fitOut; iClus_fit++){
      if(isOld){
        if(!SARPCCluster_isPlausibleFitOut->at(iMuon).at(iClus_fit)) continue; 
      }
      fitOutSlope.at(iMuon).push_back(SARPCCluster_fitOutSlope->at(iMuon).at(iClus_fit));
      cout << "plausibleFitOutSlope = " << SARPCCluster_fitOutSlope->at(iMuon).at(iClus_fit) << endl;
      fitOutOffset.at(iMuon).push_back(SARPCCluster_fitOutOffset->at(iMuon).at(iClus_fit));
      cout << "plausibleFitOutOffset = " << SARPCCluster_fitOutOffset->at(iMuon).at(iClus_fit) << endl;

    }
  }
  int nPlauBI = fitMidSlope.at(0).size(); 
  int nPlauBM = fitMidSlope.at(0).size(); 
  int nPlauBO = fitOutSlope.at(0).size(); 
  int nOff1 = OffsegSlope.at(0).size();
  int nOff2 = OffsegSlope.at(1).size();
  TF1 f_offseg1[nOff1];
  TF1 f_offseg2[nOff2];
  for(int iOff = 0; iOff < nOff1; iOff++){
    f_offseg1[iOff] = TF1("f_offseg1", "[0]*x+[1]", OffsegMin.at(0).at(iOff), OffsegMax.at(0).at(iOff));
    f_offseg1[iOff].SetTitle(";Z [m];R [m]");
    f_offseg1[iOff].SetParameter(0,OffsegSlope.at(0).at(iOff));
    f_offseg1[iOff].SetParameter(1,OffsegOffset.at(0).at(iOff)/1000.);
    f_offseg1[iOff].SetLineColor(kRed);
    f_offseg1[iOff].SetLineWidth(2);
    f_offseg1[iOff].SetLineStyle(1);
  }
  for(int iOff = 0; iOff < nOff2; iOff++){
    f_offseg2[iOff] = TF1("f_offseg2", "[0]*x+[1]", OffsegMin.at(1).at(iOff), OffsegMax.at(1).at(iOff));
    f_offseg2[iOff].SetTitle(";Z [m];R [m]");
    f_offseg2[iOff].SetParameter(0,OffsegSlope.at(1).at(iOff));
    f_offseg2[iOff].SetParameter(1,OffsegOffset.at(1).at(iOff)/1000.);
    f_offseg2[iOff].SetLineColor(kRed);
    f_offseg2[iOff].SetLineWidth(2);
    f_offseg2[iOff].SetLineStyle(1);
  }

  int nClus_set = fitMidSlope.at(0).size();

  //additional cut for 1road analysis
  int nClus = SARPCCluster_clusterLayer->at(0).size();
  bool isSkip = false;
  for(int iClus = 0; iClus < nClus; iClus++){
    if(SARPCCluster_clusterLayer->at(0).at(iClus) == 0) isSkip = true;
  }
  vector<vector<int>> clus_layer, clus_outlier;
  vector<int> intclear;
  intclear.clear();
  clus_layer.assign(4, intclear);
  clus_outlier.assign(4, intclear);
  for(int iClus = 0; iClus < nClus; iClus++){
    if(SARPCCluster_clusterMeasPhi->at(0).at(iClus)) continue;
    //======check where the test cluster is located======
    float r = sqrt(SARPCCluster_gX->at(0).at(iClus)*SARPCCluster_gX->at(0).at(iClus)+SARPCCluster_gY->at(0).at(iClus)*SARPCCluster_gY->at(0).at(iClus));
    float phi = atan(SARPCCluster_gY->at(0).at(iClus)/SARPCCluster_gX->at(0).at(iClus));
    if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)>0) phi = phi + pi;
    if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)<0) phi = phi - pi;
    float l = sqrt(SARPCCluster_gZ->at(0).at(iClus)*SARPCCluster_gZ->at(0).at(iClus)+r*r);
    float tan = sqrt( (l-SARPCCluster_gZ->at(0).at(iClus))/(l+SARPCCluster_gZ->at(0).at(iClus)) );
    float eta = -log(tan);
    float deta = fabs(SARoIEta->at(0) - eta);
    float dphi = acos(cos( SARoIPhi->at(0) - phi ) );
    std::cout << "deta/dphi = " << deta << "/" << dphi;
    //====================================================
    int i_layer = SARPCCluster_clusterLayer->at(0).at(iClus);
    std::cout << "  layer = " << i_layer << std::endl;
    if((i_layer < 4) && deta < 0.1 && dphi < 0.1) clus_layer.at(i_layer).push_back(1);
    if((i_layer < 4) && deta >= 0.1 && dphi >= 0.1) clus_outlier.at(i_layer).push_back(1);
  }//iClus loop
  int countClus = 0;
  int countOutlier = 0;
  for(int iLay = 0; iLay < 4; iLay++){
    if(clus_layer.at(iLay).size() > 1) countClus++;
    if(clus_outlier.at(iLay).size() > 0) countOutlier++;
  }
//  if(countClus != 1) continue;   //motomoto countClus < 2
  int countPhiClus = 0;
  if(isOld){
    for(unsigned int iClus_phi = 0; iClus_phi < SARPCCluster_isPlausiblePhiInnMid->at(0).size(); iClus_phi++){
      if(SARPCCluster_isPlausiblePhiInnMid->at(0).at(iClus_phi)) countPhiClus++;
    }
  }
//  if(countPhiClus !=2) continue;
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
  int N_pt = SAptclus->at(0).size();
  if(isOld){
    for(int i=0; i< SAptclus->at(0).size(); i++){
      std::cout << "pTclus/SPCsetID : ";
      std::cout << SAptclus->at(0).at(i) << "/" << SAspcSetID->at(0).at(i).at(0) << ","<< SAspcSetID->at(0).at(i).at(1) << "," << SAspcSetID->at(0).at(i).at(2) << std::endl;
    }
  }
  bool isDraw = false;
  if(N_spmid == 1 && N_clusRoadInn == 2) isDraw = true;
//  if(!isDraw)continue;

  TText eventInfo;
  setLegend_eventInfo(&eventInfo);
  TH1* frame_BM = c2->DrawFrame(Zmin_BM,Rmin_BM,Zmax_BM,Rmax_BM);
  //TH1* frame_BM = c2->DrawFrame(Zmin_BO,6.6,Zmax_BM,Rmax_BO);
  frame_BM->GetXaxis()->SetTitle("Z [m]");
  frame_BM->GetYaxis()->SetTitle("R [m]");
  for(int iOff = 0; iOff < nOff1; iOff++){
    f_offseg1[iOff].Draw("f, same");
  }
  for(int iOff = 0; iOff < nOff2; iOff++){
    f_offseg2[iOff].Draw("f, same");
  }
  roadRegion.Draw("f");
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    f_roi[iMuon].Draw("same");
    f_roi_minus[iMuon].Draw("same");
    f_roi_plus[iMuon].Draw("same");
    //gr_segment_RZ[iMuon].Draw("P, same");
    gr_RpcHit[iMuon].Draw("P, same");
    f_road_BM[iMuon].Draw("f, same");
  }//iMuon loop end
//  gr_mdt_hit.Draw("P");
//  gr_mdt_hit_outlier.Draw("P");
//  gr_mdt.SetMarkerStyle(24);
//  gr_mdt.SetMarkerSize(1);
//  gr_mdt.SetMarkerColor(kGreen+2);
//  gr_mdt_outlier.SetMarkerStyle(24);
//  gr_mdt_outlier.SetMarkerSize(1);
//  gr_mdt_outlier.SetMarkerColor(kRed);
//  gr_mdt_outlier5.Draw("P");
//  gr_mdt.Draw("P");
//  gr_mdt_outlier.Draw("P");
//  gr_SPC_BM.SetMarkerSize(1.5);
//  gr_SPC_BM.SetMarkerColor(kBlue);
//  gr_SPC_BM.Draw("P");
//  gr_SPC_BM_2.SetMarkerSize(1.5);
//  gr_SPC_BM_2.SetMarkerColor(kRed);
//  gr_SPC_BM_2.Draw("P");

  ATLASLabel(0.25, 0.9, "Work In Progress");
  myText(0.25, 0.86, 1, "Simulation #sqrt{s}=13TeV");
  myText(0.25, 0.82, 1, "J/#psi#rightarrow#mu#mu");
  eventInfo.Draw();
  c2->Print(pdf, "pdf");
  c2->RedrawAxis();
  delete frame_BM;
  
  current_entry += 1;
  cout << "===" << begin_entry << ": " << current_entry << ": " << limit_entry << endl;
  if (current_entry > limit_entry) {
    cout << "END!!" << endl;
    break;
  }
  c2->Clear();
  } // end of each entry

  c2->Print(pdf + "]", "pdf");
  delete c2;

}



bool RPC_FCBM::isCloseByCut(UInt_t nMuon, vector< pair<int, int>>& MUpair){ 

  bool isSameOffPt = false;
  if(nMuon != 2){
    //   h_CloseByCount->Fill(3);
    return false;
  }
  else{
    for(int iMuon = 0; iMuon < nMuon -1; iMuon++){
      for(int jMuon = iMuon + 1; jMuon < nMuon; jMuon++){
        double deltaNum = SARoINumber->at(iMuon) - SARoINumber->at(jMuon);
        double deltaSec = SARoISector->at(iMuon) - SARoISector->at(jMuon);
        double charge = OfflineCharge->at(iMuon) + OfflineCharge->at(jMuon);
        if(deltaSec == 0 && deltaNum == 0 && charge == 0){
          MUpair.push_back(make_pair(iMuon, jMuon));
        }
        if(OfflinePt->at(iMuon) == OfflinePt->at(jMuon))isSameOffPt = true;
      }
/*      for(int jMuon = iMuon + 1; jMuon < nMuon; jMuon++){
        double distRoIiSAEta = fabs(SARoIEta->at(iMuon) - OfflineEta->at(iMuon));
        double distRoIiSAPhi = fabs(SARoIPhi->at(iMuon) - OfflinePhi->at(iMuon));
        double distRoIjSAEta = fabs(SARoIEta->at(iMuon) - OfflineEta->at(jMuon));
        double distRoIjSAPhi = fabs(SARoIPhi->at(iMuon) - OfflinePhi->at(jMuon));
        double deltaR_i = sqrt(distRoIiSAEta*distRoIiSAEta + distRoIiSAPhi*distRoIiSAPhi);
        double deltaR_j = sqrt(distRoIjSAEta*distRoIjSAEta + distRoIjSAPhi*distRoIjSAPhi);
        if(deltaR_i < 0.2 && deltaR_j < 0.2){
          MUpair.push_back(make_pair(iMuon, jMuon));
        }
      }
    for(int iMuon = 0; iMuon < nMuon -1; iMuon++){
      for(int jMuon = iMuon + 1; jMuon < nMuon; jMuon++){
        double deltaExtEta = OfflineExtEta->at(iMuon) - OfflineExtEta->at(jMuon);
        double deltaExtPhi = OfflineExtPhi->at(iMuon) - OfflineExtPhi->at(jMuon);
        double deltaR = sqrt(pow(deltaExtEta, 2) + pow(deltaExtPhi, 2));
        if(deltaR <= 0.2){
          MUpair.push_back(make_pair(iMuon, jMuon));
        }
      }*/
    }

    if(!MUpair.empty() && !isSameOffPt){
      return true;
    }
    else{
      return false;
    }
  }

}

bool RPC_FCBM::isCloseBy(int iMuon){ 
  double distRoISAEta = fabs(SARoIEta->at(iMuon) - OfflineEta->at(iMuon));
  double distRoISAPhi = fabs(SARoIPhi->at(iMuon) - OfflinePhi->at(iMuon));
  if(distRoISAEta < 0.1 && distRoISAPhi < 0.1){
    return true;
  }
  else return false;
}


void RPC_FCBM::getDistToEachRoI(vector<float>* &roiPos, vector< vector<double>>* &paramPos, vector<vector<double>>& getParamPos, vector<vector<UInt_t>>* &MesPhi, bool isPhi){
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    for(int iParam = 0; iParam < nMuon; iParam++){
      int ncol = (int)(paramPos->at(iParam)).size();
      vector <double> setParam;
      setParam.clear();
      if(ncol > 0){
        for(int icol = 0; icol < ncol; icol++){
          if(isPhi){
            if((MesPhi->at(iParam)).at(icol)){
              setParam.push_back(fabs(roiPos->at(iMuon)-(paramPos->at(iParam)).at(icol)));
            }
          }
          else{
            if(!(MesPhi->at(iParam)).at(icol)){
              setParam.push_back(fabs(roiPos->at(iMuon)-(paramPos->at(iParam)).at(icol)));
            }
          }
        }
        getParamPos.push_back(setParam);
      }
    }
  }
}
