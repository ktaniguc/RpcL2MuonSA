#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TCollection.h"
#include "THStack.h"
#include "/home/ktaniguc/RootUtils/RootUtils/TLegend_addfunc.h"
#include "/home/ktaniguc/RootUtils/RootUtils/TCanvas_opt.h"


using namespace std;

void RPC_FCBM::DrawFCBM(TString pdf){
  //==================================================================
  //Dear Shiomi and Sumi
  //
  //If you really wanna try the algorithms like below, I recommend you
  //to see "/home/ktaniguc/RpcClusteringHistMaker".
  //It is easier to understand and the construction is 
  //just same with "RpcL2MuonSA".
  //
  //If you copy the codes of "RpcClusteringHistMaker" in 
  //your own directory, go and check the following URL
  //https://gitlab.cern.ch/ktaniguc/RpcClusteringHistMaker
  //
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1080, 700);
  TLegend *legDeltaPt = new TLegend(0.1,0.8,0.3,0.95);
  legDeltaPt->AddEntry(h_DeltaPt_OffSA_CloseBy, "same RoI");
  legDeltaPt->AddEntry(h_DeltaPt_OffSA, "different RoI");
  TLegend *legRpcdistMin = new TLegend(0.1,0.8,0.3,0.95);
  legRpcdistMin->AddEntry(h_DeltaRpcMidMinEta_CloseBy, "Middle");
  legRpcdistMin->AddEntry(h_DeltaRpcOutMinEta_CloseBy, "Outer");
  TLegend *legEtaPhiLay[8];
  for(int i = 0; i < 8; i++){
    legEtaPhiLay[i] = new TLegend(0.8,0.8,0.9,0.9);
    legEtaPhiLay[i]->AddEntry(h_rpcClusEtaSizeLay_CloseBy[i], "#eta strip");
    legEtaPhiLay[i]->AddEntry(h_rpcClusPhiSizeLay_CloseBy[i], "#phi strip");
  }
 
  gStyle->SetOptStat(0);
  // c1->SetGrid();
 // c1->SetRightMargin(0.20);
 // c1->SetLeftMargin(0.23);
  c1->SetTopMargin(0.20);
 // c1->SetBottomMargin(0.10);

  c1 -> Print( pdf + "[", "pdf" );
/*  h_distRoIrpcEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distRoIrpcPhi->Draw();
  c1 -> Print(pdf, "pdf" );*/
  c1->SetLogy();
  h_DeltaPt_OffSA->SetLineColor(kRed);
  h_DeltaPt_OffSA->Draw();
  h_DeltaPt_OffSA_CloseBy->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  h_isSameRpcFitMid_CloseBy->Draw();
  c1 -> Print(pdf, "pdf" );
  
  h_rpcHitSize_CloseBy->Draw();
  h_rpcHitSize->SetLineColor(kRed);
  h_rpcHitSize->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );

  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcHitSizeLay[iLay]->SetLineColor(kRed);
    h_rpcHitSizeLay[iLay]->Draw();
    /*double scalefactor = h_rpcHitSizeLay[iLay]->Integral()/h_rpcHitSizeLay_CloseBy[iLay]->Integral();
    h_rpcHitSizeLay_CloseBy[iLay]->Scale(scalefactor);*/
    h_rpcHitSizeLay_CloseBy[iLay]->Draw("same");
    legDeltaPt->Draw("same");
    c1 -> Print(pdf, "pdf" );
  }

  c1->SetLogy();
  h_distOffEta_CloseBy->Draw();
  h_distOffEta->SetLineColor(kRed);
  h_distOffEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distOffPhi_CloseBy->Draw();
  h_distOffPhi->SetLineColor(kRed);
  h_distOffPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitInnEta_CloseBy->Draw();
  h_distRpcFitInnEta->SetLineColor(kRed);
  h_distRpcFitInnEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitInnPhi_CloseBy->Draw();
  h_distRpcFitInnPhi->SetLineColor(kRed);
  h_distRpcFitInnPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitMidEta_CloseBy->Draw();
  h_distRpcFitMidEta->SetLineColor(kRed);
  h_distRpcFitMidEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitMidPhi_CloseBy->Draw();
  h_distRpcFitMidPhi->SetLineColor(kRed);
  h_distRpcFitMidPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitOutEta_CloseBy->Draw();
  h_distRpcFitOutEta->SetLineColor(kRed);
  h_distRpcFitOutEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcFitOutPhi_CloseBy->Draw();
  h_distRpcFitOutPhi->SetLineColor(kRed);
  h_distRpcFitOutPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_distRpcHitToFitInnEta_CloseBy->Draw();
  h_distRpcHitToFitInnEta->SetLineColor(kRed);
  h_distRpcHitToFitInnEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcHitToFitInnPhi_CloseBy->Draw();
  h_distRpcHitToFitInnPhi->SetLineColor(kRed);
  h_distRpcHitToFitInnPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcHitToFitMidEta_CloseBy->Draw();
  h_distRpcHitToFitMidEta->SetLineColor(kRed);
  h_distRpcHitToFitMidEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcHitToFitMidPhi_CloseBy->Draw();
  h_distRpcHitToFitMidPhi->SetLineColor(kRed);
  h_distRpcHitToFitMidPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcHitToFitOutEta_CloseBy->Draw();
  h_distRpcHitToFitOutEta->SetLineColor(kRed);
  h_distRpcHitToFitOutEta->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_distRpcHitToFitOutPhi_CloseBy->Draw();
  h_distRpcHitToFitOutPhi->SetLineColor(kRed);
  h_distRpcHitToFitOutPhi->Draw("same");
  legDeltaPt->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_DeltaRpcMidMinEta_CloseBy->Draw();
  h_DeltaRpcOutMinEta_CloseBy->SetLineColor(kRed);
  h_DeltaRpcOutMinEta_CloseBy->Draw("same");
  legRpcdistMin->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_DeltaRpcMidMinPhi_CloseBy->Draw();
  h_DeltaRpcOutMinPhi_CloseBy->SetLineColor(kRed);
  h_DeltaRpcOutMinPhi_CloseBy->Draw("same");
  legRpcdistMin->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  h_RoInumDelEta->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_RoInumDelPhi->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  
  //RPC Cluster
  h_rpc2ClusEtaNum_CloseBy->Draw();
  c1 -> Print(pdf, "pdf" );
  h_rpc2ClusPhiNum_CloseBy->Draw();
  c1 -> Print(pdf, "pdf" );
  h_rpcClusEta_CloseBy->Draw();
  c1 -> Print(pdf, "pdf" );
  h_rpcClusPhi_CloseBy->Draw();
  c1 -> Print(pdf, "pdf" );
  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcClusEtaSizeLay_CloseBy[iLay]->SetLineColor(kRed);
    h_rpcClusEtaSizeLay_CloseBy[iLay]->Draw();
    /*double scalefactor = h_rpcHitSizeLay[iLay]->Integral()/h_rpcHitSizeLay_CloseBy[iLay]->Integral();
    h_rpcHitSizeLay_CloseBy[iLay]->Scale(scalefactor);*/
    h_rpcClusPhiSizeLay_CloseBy[iLay]->Draw("same");
    legEtaPhiLay[iLay]->Draw("same");
    c1 -> Print(pdf, "pdf" );
  }

  c1 -> Print( pdf + "]", "pdf" );

  delete c1;
}

void RPC_FCBM::DrawNoCut(TString pdf){
  //==================================================================
  //Set Canvas
  //==================================================================
  //TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1080, 800);
  TFile *fout = new TFile(Form("outroot/DrawNoCut_1layerSearch.root", pdf.Data()), "recreate");
  TCanvas_opt *c1 = new TCanvas_opt();
 // c1->SetGrid();
 // c1->SetRightMargin(0.20);
 // c1->SetLeftMargin(0.23);
  //c1->SetTopMargin(0.10);
 // c1->SetBottomMargin(0.10);
  TLegend *legEtaPhiLay2;
  legEtaPhiLay2 = new TLegend(0.8,0.8,0.9,0.9);
  legEtaPhiLay2->AddEntry(h_nClusLay2Eta, "#eta strip");
  legEtaPhiLay2->AddEntry(h_nClusLay2Phi, "#phi strip");
  
  TLegend *legrpcHitXYLay = new TLegend(0.75, 0.75, 1, 0.95);
  for(int i = 0; i < 8; i++){
    legrpcHitXYLay->AddEntry(h_rpcHitXYLay[i], Form("layer %d",i ));
  }
  TLegend_addfunc *lef_select;
  lef_select = new TLegend_addfunc(6);
  lef_select->AddSelection("Off-p_{T}#geq4GeV");
  lef_select->AddSelection("Off-RoI |#eta|<1.05");
  lef_select->AddSelection("Off-Segment |#eta|<1.05");
  lef_select->AddSelection("Off-|z_{0}|#leq50(mm)");
  lef_select->AddSelection("N_{OffsegMid}#geq1");
  lef_select->AddSelection("N_{OffsegOut}#geq1");

  c1 -> Print( pdf + "[", "pdf" );
  h_distnoCutRoIrpcEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distnoCutRoIrpcPhi->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distnoCutEachRoIrpcEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distnoCutEachRoIrpcPhi->Draw();
  c1 -> Print(pdf, "pdf" );
  h_isCloseBydEdP->Draw();
  c1 -> Print(pdf, "pdf" );
  h_isSameRoI->Draw();
  c1 -> Print(pdf, "pdf" );
  h_isSameRpcFitMid->Draw();
  c1 -> Print(pdf, "pdf" );
  h_OffSegZR->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_OffSegXY->Draw();
  c1 -> Print(pdf, "pdf" );
  h_nOffSegMid->Draw();
  c1 -> Print(pdf, "pdf" );
  h_nOffSegOut->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_OffSegEtaPhi->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  //h_OffSegNEtaPhiLayers->Draw("colz");
  //c1 -> Print(pdf, "pdf" );
  for(int iLay = 0; iLay < 8; iLay++){
    h_OffSegXY_etaLayer[iLay]->Draw();
    c1 -> Print(pdf, "pdf" );
  }
  
  
  
  TLegend *leg_nClusRoad;
  leg_nClusRoad = new TLegend(0.5,0.8,0.7,0.9);
  leg_nClusRoad->AddEntry(h_nClusRoad_withMidInfo, "by middle cluster");
  leg_nClusRoad->AddEntry(h_nClusRoad_withOutInfo, "by outer cluster");
  h_nClusRoad_withOutInfo->SetLineColor(kRed);
  h_nClusRoad_withOutInfo->Draw();
  h_nClusRoad_withMidInfo->Draw("same");
  leg_nClusRoad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusRoadOut_woOutInfo->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusRoadOut_byOuter->SetLineColor(kRed);
  h_nClusRoadOut_byMiddle->Draw();
  h_nClusRoadOut_byOuter->Draw("same");
  leg_nClusRoad->Draw("same");
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );

  h_nClusRoadMid_byOuterTest->SetLineColor(kRed);
  h_nClusRoadMid_byMiddleTest->Draw();
  h_nClusRoadMid_byOuterTest->Draw("same");
  leg_nClusRoad->Draw("same");
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );


  c1->SetLogy();
  h_etamin_MidOffFit_woMidInfo->SetLineColor(kRed);
  h_etamin_MidOffFit_woMidInfo->Draw();
  h_etamin_MidOffFit_wMidInfo->Draw("same");
  leg_nClusRoad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_slopemin_MidOffFit_woMidInfo->SetLineColor(kRed);
  h_slopemin_MidOffFit_woMidInfo->Draw();
  h_slopemin_MidOffFit_wMidInfo->Draw("same");
  leg_nClusRoad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_etamin_OutOffFit->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_slopemin_OutOffFit->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  TLegend *leg_compMidOutOff;
  leg_compMidOutOff = new TLegend(0.5,0.8,0.7,0.9);
  leg_compMidOutOff->SetHeader("compared with:");
  leg_compMidOutOff->AddEntry(h_etamin_MidOffFit_woMidInfo, "OffsegMid");
  leg_compMidOutOff->AddEntry(h_etamin_MidClusOutInfo_OffsegOut, "OffsegOut");
  leg_compMidOutOff->SetTextSize(0.04);
  h_etamin_MidOffFit_woMidInfo->SetLineColor(kRed);
  h_etamin_MidClusOutInfo_OffsegOut->Draw();
  h_etamin_MidOffFit_woMidInfo->Draw("same");


  lef_select->Draw("same");
  leg_compMidOutOff->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_slopemin_MidOffFit_woMidInfo->SetLineColor(kRed);
  h_slopemin_MidClusOutInfo_OffsegOut->Draw();
  h_slopemin_MidOffFit_woMidInfo->Draw("same");
  lef_select->Draw("same");
  leg_compMidOutOff->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_etamin_MidOffFit_strip->Draw();
  c1 -> Print(pdf, "pdf" );
  h_slopemin_MidOffFit_strip->Draw();
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);

  TLegend *leg_ClusFitMid;
  leg_ClusFitMid = new TLegend(0.8,0.8,0.95,0.95);
  leg_ClusFitMid->AddEntry(h_ClusFitMidEta, "clusterRoad ");
  leg_ClusFitMid->AddEntry(h_SARPCFitMidEta, "rpc strip Road");

  h_SARPCFitInnR->Draw();
  c1 -> Print(pdf, "pdf" );
  h_SARPCFitInnEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_SARPCFitMidR->Draw();
  c1 -> Print(pdf, "pdf" );
  h_ClusFitMidEta->SetLineColor(kRed);
  h_ClusFitMidEta->Draw();
  h_SARPCFitMidEta->Draw("same");
  leg_ClusFitMid->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_SARPCFitOutR->Draw();
  c1 -> Print(pdf, "pdf" );
  h_SARPCFitOutEta->Draw();
  h_ClusFitOutEta->SetLineColor(kRed);
  h_ClusFitOutEta->Draw();
  h_SARPCFitOutEta->Draw("same");
  leg_ClusFitMid->Draw("same");
  c1 -> Print(pdf, "pdf" );
  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcHitXYLay[iLay]->Draw();
    c1 -> Print(pdf, "pdf" );
  }
  h_rpcHitXYLay[0]->SetMarkerColor(kBlack);
  h_rpcHitXYLay[0]->SetLineColor(kBlack);
  h_rpcHitXYLay[0]->Draw();
  h_rpcHitXYLay[1]->SetMarkerColor(kGray+2);
  h_rpcHitXYLay[1]->SetLineColor(kGray+2);
  h_rpcHitXYLay[1]->Draw("same");
  h_rpcHitXYLay[2]->SetMarkerColor(kRed);
  h_rpcHitXYLay[2]->SetLineColor(kRed);
  h_rpcHitXYLay[2]->Draw("same");
  h_rpcHitXYLay[3]->SetMarkerColor(kRed-9);
  h_rpcHitXYLay[3]->SetLineColor(kRed-9);
  h_rpcHitXYLay[3]->Draw("same");
  h_rpcHitXYLay[4]->SetMarkerColor(kBlue);
  h_rpcHitXYLay[4]->SetLineColor(kBlue);
  h_rpcHitXYLay[4]->Draw("same");
  h_rpcHitXYLay[5]->SetMarkerColor(kBlue-9);
  h_rpcHitXYLay[5]->SetLineColor(kBlue-9);
  h_rpcHitXYLay[5]->Draw("same");
  h_rpcHitXYLay[6]->SetMarkerColor(kGreen);
  h_rpcHitXYLay[6]->SetLineColor(kGreen);
  h_rpcHitXYLay[6]->Draw("same");
  h_rpcHitXYLay[7]->SetMarkerColor(kGreen-9);
  h_rpcHitXYLay[7]->SetLineColor(kGreen-9);
  h_rpcHitXYLay[7]->Draw("same");
  legrpcHitXYLay->Draw("same");
  c1 -> Print(pdf, "pdf" );

  c1->SetLogy();
  h_distMidOffClusEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusPhi->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusR->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusEtaMin->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusPhiMin->Draw();
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusfitEtaMin->Draw();
  c1 -> Print(pdf, "pdf" );
  TLegend *legDistRoIClusR = new TLegend(0.6,0.6,0.9,0.9);
  legDistRoIClusR->AddEntry(h_distRoIClusR, "no selection");
  legDistRoIClusR->AddEntry(h_distRoIClusR_CloseBy, "2muon in 1RoI");
  legDistRoIClusR->AddEntry(h_distRoIClusR_noCloseBy, "1muon in 1RoI");
  h_distRoIClusR->Draw();
  h_distRoIClusR_CloseBy->SetLineColor(kRed);
  h_distRoIClusR_CloseBy->Draw("same");
  h_distRoIClusR_noCloseBy->SetLineColor(kGreen);
  h_distRoIClusR_noCloseBy->Draw("same");
  legDistRoIClusR->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  /*for(int i = 0; i < 10; i++){
    h_fakeRoadMid[i]->Draw();
    c1 -> Print(pdf, "pdf" );
    h_truthRoadMid[i]->Draw();
    c1 -> Print(pdf, "pdf" );
  }*/
  h_nClusFitMid->Draw();
  c1 -> Print(pdf, "pdf" );
  h_nClusInSet->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_OffsegEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_ClusFitMidEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_ClusFitOutEta->Draw();
  c1 -> Print(pdf, "pdf" );
  TLegend_addfunc *lef_select_outfake;
  lef_select_outfake = new TLegend_addfunc(7);
  lef_select_outfake->AddSelection("Off-p_{T}#geq4GeV");
  lef_select_outfake->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_outfake->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_outfake->AddSelection("Off-|z_{0}|#leq50(mm)");
  lef_select_outfake->AddSelection("N_{OffsegMid}#geq1");
  lef_select_outfake->AddSelection("N_{OffsegOut}#geq1");
  lef_select_outfake->AddSelection("N_{specialLayer} = 0");
  //lef_select_outfake->AddSelection("N_{OutclusInSet}#geq1");
  TLegend_addfunc *lef_option;
  lef_option = new TLegend_addfunc(6, 4);
  lef_option->AddOption(h_fakeRoadMid_opt, "w/ all plane");
  lef_option->AddOption(h_fakeRoadMid_opt23, "(w/o RPC1)&&(w/ RPC2,RPC3)");
  lef_option->AddOption(h_fakeRoadMid_opt13, "(w/o RPC2)&&(w/ RPC1,RPC3)");
  lef_option->AddOption(h_fakeRoadMid_opt12, "(w/o RPC3)&&(w/ RPC1,RPC2)");
  //h_fakeRoadMid->Draw();
  c1->SetLogy(); 
  
  TLegend_addfunc *lef_bug;
  lef_bug = new TLegend_addfunc(7);
  lef_bug->AddSelection("Off-p_{T}#geq4GeV");
  lef_bug->AddSelection("Off-RoI |#eta|<1.05");
  lef_bug->AddSelection("Off-Segment |#eta|<1.05");
  lef_bug->AddSelection("Off-|z_{0}|#leq50(mm)");
  lef_bug->AddSelection("N_{OffsegMid}#geq1");
  lef_bug->AddSelection("N_{OffsegOut}#geq1");
  lef_bug->AddSelection("N_{specialLayer} = 0");
  /*lef_bug->AddSelection("seedlayer > 0"); 
  lef_bug->AddSelection("seedcluster id > 0");*/
  
  h_extdR_bug->SetLineColor(kRed); 
  h_extdR_wobug->Draw();
  h_extdR_bug->Draw("same");
  TLegend_addfunc *lef_bugopt;
  lef_bugopt = new TLegend_addfunc(7, 2);
  lef_bugopt->AddOption(h_extdR_bug, "seedlayer>0 && id>0"); 
  lef_bugopt->AddOption(h_extdR_wobug, "seedlayer=0 || id=0");
  lef_bug->Draw("same");
  lef_bugopt->Draw("same");
  c1 -> Print(pdf, "pdf" );

  c1->SetLogy(0);
  h_RoIEtaPhi_bug->Draw("colZ");
  c1 -> Print(pdf, "pdf" );

  c1->SetLogy();

  h_fakeRoadMid_opt23->SetLineColor(kBlue);
  h_fakeRoadMid_opt13->SetLineColor(kGreen);
  h_fakeRoadMid_opt12->SetLineColor(kRed);
  h_fakeRoadMid_opt12->Draw();
  h_fakeRoadMid_opt23->Draw("same");
  h_fakeRoadMid_opt13->Draw("same");
  h_fakeRoadMid_opt->Draw("same");
  lef_option->Draw("same");
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  //h_truthRoadMid->Draw();
  h_truthRoadMid_opt23->SetLineColor(kBlue);
  h_truthRoadMid_opt13->SetLineColor(kGreen);
  h_truthRoadMid_opt12->SetLineColor(kRed);
  h_truthRoadMid_opt12->Draw();
  h_truthRoadMid_opt23->Draw("same");
  h_truthRoadMid_opt13->Draw("same");
  h_truthRoadMid_opt->Draw("same");
  lef_option->Draw("same");
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_faketruthRoadMid->Draw("colz");
  c1 -> Print(pdf, "pdf" );
 
  //h_fakeRoadOut->Draw();
  //h_fakeRoadOut_opt->SetMarkerColor(kRed);
  h_fakeRoadOut_opt23->SetLineColor(kBlue);
  h_fakeRoadOut_opt13->SetLineColor(kGreen);
  h_fakeRoadOut_opt12->SetLineColor(kRed);
  h_fakeRoadOut_opt12->Draw();
  h_fakeRoadOut_opt23->Draw("same");
  h_fakeRoadOut_opt13->Draw("same");
  h_fakeRoadOut_opt->Draw("same");
  lef_select_outfake->Draw("same");
  lef_option->Draw("same");
  c1 -> Print(pdf, "pdf" );
//  h_truthRoadOut->Draw();
  //h_truthRoadOut_opt->SetMarkerColor(kRed);
  h_truthRoadOut_opt23->SetLineColor(kBlue);
  h_truthRoadOut_opt13->SetLineColor(kGreen);
  h_truthRoadOut_opt12->SetLineColor(kRed);
  h_truthRoadOut_opt12->Draw();
  h_truthRoadOut_opt23->Draw("same");
  h_truthRoadOut_opt13->Draw("same");
  h_truthRoadOut_opt->Draw("same");
  lef_select_outfake->Draw("same");
  lef_option->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_faketruthRoadOut->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_nmaxClustersInSet->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_id_minlayer->SetLineColor(kRed);
  TLegend_addfunc *lef_optID;
  lef_optID = new TLegend_addfunc(6, 2);
  lef_optID->AddOption(h_id_minlayer, "seedLayer > 0");
  lef_optID->AddOption(h_id_minlayer0, "seedLayer = 0");
  h_id_minlayer->Draw();
  h_id_minlayer0->Draw("same");
  lef_select->Draw("same");
  lef_optID->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_seedcluster_setmax->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_minlayerInSetmax->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_iswoRpc1Road->Draw();
  lef_select->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  
  h_fitMidSlopeOffset_cut->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  h_fitMidSlopeOffset_wocut->Draw("colZ");
  c1 -> Print(pdf, "pdf" );

  c1->SetLogz();
  h_distMidOffClus->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  h_distRoIClus->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  h_distMidOffClusMin->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogz(0);
  h_nClusLay2Eta->SetLineColor(kRed);
  h_nClusLay2Eta->Draw();
  h_nClusLay2Phi->Draw("same");
  legEtaPhiLay2->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusLay3Eta->SetLineColor(kRed);
  h_nClusLay3Eta->Draw();
  h_nClusLay3Phi->Draw("same");
  legEtaPhiLay2->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusLay23Eta->SetLineColor(kRed);
  h_nClusLay23Eta->Draw();
  h_nClusLay23Phi->Draw("same");
  legEtaPhiLay2->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusHit23Eta->SetLineColor(kRed);
  h_nClusHit23Eta->Draw();
  h_nClusHit23Phi->Draw("same");
  legEtaPhiLay2->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_nClusLay23EtaPhi->Draw();
  c1 -> Print(pdf, "pdf" );
  h_clusterSizeLay01->Draw();
  c1 -> Print(pdf, "pdf" );
  h_clusterSizeLay23->Draw();
  c1 -> Print(pdf, "pdf" );
  h_clusterSizeLayOut->Draw();
  c1 -> Print(pdf, "pdf" );
  h_clusterSize->Draw();
  c1 -> Print(pdf, "pdf" );
  
  TLegend *legClusSizeLay;
  legClusSizeLay = new TLegend(0.8,0.8,0.9,0.9);
  legClusSizeLay->AddEntry(h_clusterSizeLay01, "RPC1 plane");
  legClusSizeLay->AddEntry(h_clusterSizeLay23, "RPC2 plane");
  legClusSizeLay->AddEntry(h_clusterSizeLayOut, "RPC3 plane");
  h_clusterSizeLayOut->SetLineColor(kGreen);
  h_clusterSizeLayOut->Draw();
  h_clusterSizeLay01->SetLineColor(kRed);
  h_clusterSizeLay01->Draw("same");
  h_clusterSizeLay23->SetLineColor(kBlue);
  h_clusterSizeLay23->Draw("same");
  legClusSizeLay->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_isSameSize_rpcHitEta->Draw();
  c1 -> Print(pdf, "pdf" );
  h_isSameSize_rpcHitPhi->Draw();
  c1 -> Print(pdf, "pdf" );

  TLegend_addfunc *lef_select_Zmumu;
  lef_select_Zmumu = new TLegend_addfunc(4);
  lef_select_Zmumu->AddSelection("2muon in 1RoI");
  lef_select_Zmumu->AddSelection("Off-p_{T}#geq4GeV");
  //lef_select_Zmumu->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_Zmumu->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_Zmumu->AddSelection("Off-|z_{0}|#leq50(mm)");
  //lef_select_Zmumu->AddSelection("N_{OffsegMid}#geq1");
  //lef_select_Zmumu->AddSelection("N_{OffsegOut}#geq1");
//  lef_select_Zmumu->AddSelection("#DeltaR_{OffMu}^{middle}<0.1");
  c1->SetLogy();
  h_Zmin_OffClusRoadMid->SetLineColor(kRed);
  h_Zmin_OffClusRoadMid->Draw();
  h_Zsubmin_OffClusRoadMid->Draw("same");
  lef_select_Zmumu->Draw("same");
  TLegend_addfunc *lopt_Zmin;
  lopt_Zmin = new TLegend_addfunc(7,2);
  lopt_Zmin->AddOption(h_Zmin_OffClusRoadMid, "smallest |#DeltaZ|");
  lopt_Zmin->AddOption(h_Zsubmin_OffClusRoadMid, "second smallest |#DeltaZ|");
  lopt_Zmin->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_Z_plaumin_OffClusRoadMid->SetLineColor(kRed);
  h_Z_plaumin_OffClusRoadMid->Draw();
  h_Z_plausubmin_OffClusRoadMid->Draw("same");
  lef_select_Zmumu->Draw("same");
  TLegend_addfunc *lopt_Z_plaumin;
  lopt_Z_plaumin = new TLegend_addfunc(7,2);
  lopt_Z_plaumin->AddOption(h_Z_plaumin_OffClusRoadMid, "smallest |#DeltaZ|");
  lopt_Z_plaumin->AddOption(h_Z_plausubmin_OffClusRoadMid, "second smallest |#DeltaZ|");
  lopt_Z_plaumin->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  
  h_countClusterRoadMid->Write();
  h_countClusterRoadOut->Write();
  h_countClusterRoadMid->SetLineColor(kRed);
  h_countClusterRoadOut->SetLineColor(kBlue);
  //h_n_clusRoadMid_plau_n->Draw();
  h_countClusterRoadMid->Draw();
  //h_n_clusRoadMid_total->Draw("same");
//  h_countClusterRoadOut->Draw("same");
  //lef_select_Zmumu->Draw("same");
  TLegend_addfunc *lopt_nroad;
  lopt_nroad = new TLegend_addfunc(4,2);
  lopt_nroad->AddOption(h_countClusterRoadMid, "Middle");
  lopt_nroad->AddOption(h_countClusterRoadOut, "Outer");
  //lopt_nroad->AddOption(h_n_clusRoadMid_plau_n, "w/ extended similarRoad cut");
//  lopt_nroad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  hh_num_clusRoadMidvsOut->Draw("colZ");
  c1->Print(pdf, "pdf");
  c1->SetLogz();
  hh_n_CRvsSPCMid->Draw("colZ");
  c1->Print(pdf, "pdf");
  c1->SetLogz(0);
  prof_CRvsSPCMid = hh_n_CRvsSPCMid->ProfileX();
  prof_CRvsSPCMid->SetErrorOption("S");
  prof_CRvsSPCMid->Draw();
  c1->Print(pdf, "pdf");
  
  h_dtheta_total->Draw();
  h_dtheta->SetLineColor(kBlue);
  h_dtheta->Draw("same");
  h_dtheta_sub->SetLineColor(kRed);
  h_dtheta_sub->Draw("same");
  TLegend_addfunc *lopt_dtheta;
  lopt_dtheta = new TLegend_addfunc(4,2);
  lopt_dtheta->AddOption(h_dtheta_total, "sum");
  lopt_dtheta->AddOption(h_dtheta, "most smallest #Delta#theta");
  lopt_dtheta->AddOption(h_dtheta_sub, "the other pair");
  lopt_dtheta->Draw("same");
  c1 -> Print(pdf, "pdf" );


  m_h_deltaExtR_off_tot->Draw();
  m_h_deltaExtR_off_morecand->SetLineColor(kRed);
  m_h_deltaExtR_off_morecand->Draw("same");
  TLegend_addfunc *lef_select_morecand;
  lef_select_morecand = new TLegend_addfunc(3);
  lef_select_morecand->AddSelection("Off-p_{T}#geq4GeV");
  //lef_select_morecand->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_morecand->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_morecand->AddSelection("Off-|z_{0}|#leq50(mm)");
  //lef_select_morecand->AddSelection("N_{OffsegMid}#geq1");
  //lef_select_morecand->AddSelection("N_{OffsegOut}#geq1");
  lef_select_morecand->Draw("same");
  TLegend_addfunc *lopt_eff = new TLegend_addfunc(7,2);
  lopt_eff->AddOption(m_h_deltaExtR_off_tot, "total");
  lopt_eff->AddOption(m_h_deltaExtR_off_morecand, "isMoreCandInRoI = true");
  lopt_eff->Draw("same");
  c1 -> Print(pdf, "pdf" );

  m_h_deltaExtR_off_tot_clus->Draw();
  m_h_deltaExtR_off_2road_clus->SetLineColor(kRed);
  m_h_deltaExtR_off_2road_clus->Draw("same");
  TLegend_addfunc *lef_select_2road_clus;
  lef_select_2road_clus = new TLegend_addfunc(3);
  lef_select_2road_clus->AddSelection("Off-p_{T}#geq4GeV");
  //lef_select_2road_clus->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_2road_clus->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_2road_clus->AddSelection("Off-|z_{0}|#leq50(mm)");
  //lef_select_2road_clus->AddSelection("N_{OffsegMid}#geq1");
  //lef_select_2road_clus->AddSelection("N_{OffsegOut}#geq1");
  lef_select_2road_clus->Draw("same");
  TLegend_addfunc *lopt_eff_clusRoad = new TLegend_addfunc(3,2);
  lopt_eff_clusRoad->AddOption(m_h_deltaExtR_off_tot_clus, "total");
  lopt_eff_clusRoad->AddOption(m_h_deltaExtR_off_2road_clus, "N_{clusterRoadMid}=2");
  lopt_eff_clusRoad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  m_h_deltaExtR_off_mismatch->SetLineColor(kRed);
  m_h_deltaExtR_off_mismatch->Draw();
  lef_select_2road_clus->Draw("same");
  c1 -> Print(pdf, "pdf" );

  TLegend_addfunc *lopt_eff_closeBy = new TLegend_addfunc(6,2);
  lopt_eff_closeBy->AddOption(m_h_deltaExtR_off_tot_clus, "total");
  lopt_eff_closeBy->AddOption(m_h_deltaExtR_off_closeBy, "2muon in 1RoI");
  m_h_deltaExtR_off_closeBy->Draw();
  m_h_deltaExtR_off_tot_clus->Draw("same");
  lef_select_2road_clus->Draw("same");
  lopt_eff_closeBy->Draw("same");
  c1 -> Print(pdf, "pdf" );

  h_countSPMid->SetLineColor(kRed);
  h_countClusterRoadMid->SetLineColor(kBlack);
  TLegend_addfunc *lef_spmid;
  lef_spmid = new TLegend_addfunc(1, 2);
  lef_spmid->AddOption(h_countClusterRoadMid, "N_{clusterRoadMid}"); 
  lef_spmid->AddOption(h_countSPMid, "N_{superpointMid}");
  h_countClusterRoadMid->Draw();
  h_countSPMid->Draw("same");
  lef_spmid->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  h_countSPOut->SetLineColor(kRed);
  h_countClusterRoadOut->SetLineColor(kBlack);
  TLegend_addfunc *lef_spOut;
  lef_spOut = new TLegend_addfunc(1, 2);
  lef_spOut->AddOption(h_countClusterRoadOut, "N_{clusterRoadOut}"); 
  lef_spOut->AddOption(h_countSPOut, "N_{superpointOut}");
  h_countClusterRoadOut->Draw();
  h_countSPOut->Draw("same");
  lef_spOut->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_countSPInn->SetLineColor(kBlue);
  h_countSPMid->SetLineColor(kGreen);
  h_countSPOut->SetLineColor(kRed);
  TLegend_addfunc *lef_sp;
  lef_sp = new TLegend_addfunc(1, 3);
  lef_sp->AddOption(h_countSPInn, "inner"); 
  lef_sp->AddOption(h_countSPMid, "middle");
  lef_sp->AddOption(h_countSPOut, "Outer"); 
  h_countSPMid->Draw();
  h_countSPOut->Draw("same");
  h_countSPInn->Draw("same");
  lef_sp->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  m_h_countPtclus->Draw();
  c1 -> Print(pdf, "pdf" );

  m_h_dSlope_innmid->Draw();
  c1 -> Print(pdf, "pdf" );
  
  h_OffEtaPhi->Draw("colz");
  lef_select_2road_clus->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  c1->SetLogy();
  m_h_SApt->Draw();
//  m_h_SApt_lead->Draw("same");
//  m_h_SApt_sublead->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_Offpt_lead->SetLineColor(kRed);
  h_Offpt_sublead->SetLineColor(kBlue);
  h_Offpt->Draw();
  h_Offpt_lead->Draw("same");
  h_Offpt_sublead->Draw("same");
  c1 -> Print(pdf, "pdf" );
  h_Offpt->SetLineColor(kRed);
  TLegend_addfunc *lef_ptoffclus;
  lef_ptoffclus = new TLegend_addfunc(1, 3);
  lef_ptoffclus->AddOption(h_Offpt, "Offline p_{T}"); 
  lef_ptoffclus->AddOption(m_h_ptclus, "SA p_{T} from cluster");
  m_h_ptclus->Draw();
  h_Offpt->Draw("same");
  lef_ptoffclus->Draw("same");
  c1 -> Print(pdf, "pdf" );
  c1->SetLogy(0);
  h_OffdR->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_deltaExtR_2spmid->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_deltaExtR_2pt->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_deltaExtR_2spout->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_ptclus_vs_off->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  c1 -> Print( pdf + "]", "pdf" );
  fout->Close();
  
  delete c1;
}

void RPC_FCBM::DrawBranchHist(TString pdf){
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1080, 700);
  // c1->SetGrid();
  // c1->SetRightMargin(0.20);
  // c1->SetLeftMargin(0.23);
//  c1->SetTopMargin(0.10);
  // c1->SetBottomMargin(0.10);
  gStyle->SetOptStat(0);

  c1 -> Print( pdf + "[", "pdf" );
  for(int iBranch = 0; iBranch < nBranch; iBranch++)
  {
   // h_Branch[iBranch]->SetTitle((Branch[iBranch]).c_str());
    h_Branch[iBranch]->Draw();
    long long unsigned int HistEntry = h_Branch[iBranch]->Integral();
    TText eventInfo = TText(0.05, 0.02, Form("Entries = %llu", HistEntry));
    eventInfo.SetNDC();
    eventInfo.SetTextSize(0.03);
    eventInfo.Draw("same");
    if(iBranch == 1 || iBranch == 18 || iBranch == 52 || iBranch == 53 || (iBranch >= 40 && iBranch <= 50 ) || iBranch >= 57) //if Branch is about pT contents
    {
      c1->SetLogy();
    }
    c1 -> Print(pdf, "pdf" );
    if(iBranch == 1 || iBranch == 18 || iBranch == 52 || iBranch == 53 || (iBranch >= 40 && iBranch <= 50 ) || iBranch >= 57) //if Branch is about pT contents
    {
      c1->SetLogy(0);
    }
  }
  h_OffEtaPhi->GetXaxis()->SetLabelSize(0.001);
  h_OffEtaPhi->GetYaxis()->SetLabelSize(0.001);
  h_OffEtaPhi->GetXaxis()->SetTitleSize(0.001);
  h_OffEtaPhi->GetYaxis()->SetTitleSize(0.001);
  h_OffEtaPhi->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  c1 -> Print( pdf + "]", "pdf" );
  delete c1;
}

void RPC_FCBM::DrawEfficiency(TString pdf)
{
  TFile *fout = new TFile(Form("outroot/DrawEff_worecovered.root", pdf.Data()), "recreate");
  TCanvas_opt *c1 = new TCanvas_opt();
  //CalcEfficiency(m_h_deltaExtR_off_morecand, m_h_deltaExtR_off_tot, m_h_eff_isMorecand);
  CalcEfficiency(m_h_deltaExtR_2pt_oppcharge, m_h_deltaExtR_off_tot_clus, m_h_eff_clusterRoad);
  CalcEfficiency(m_h_deltaExtR_2pt, m_h_deltaExtR_off_tot_clus, m_h_eff_isMorecand);
//  CalcEfficiency(m_h_deltaExtR_2pt_oppcharge, m_h_deltaExtR_off_tot_clus, m_h_eff_isMorecand);
  c1 -> Print( pdf + "[", "pdf" );
  TLegend_addfunc *lef_select_morecand;
  lef_select_morecand = new TLegend_addfunc(6);
  lef_select_morecand->AddSelection("Off-p_{T}#geq4GeV");
  lef_select_morecand->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_morecand->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_morecand->AddSelection("Off-|z_{0}|#leq50(mm)");
  lef_select_morecand->AddSelection("N_{OffsegMid}#geq1");
  lef_select_morecand->AddSelection("N_{OffsegOut}#geq1");
  //m_h_eff_isMorecand->Draw("P");
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_off_morecand, *m_h_deltaExtR_off_tot)){
    m_peff_morecand = new TEfficiency(*m_h_deltaExtR_off_morecand, *m_h_deltaExtR_off_tot_clus);
    m_peff_morecand->SetTitle("Rate of isMoreCandinRoI=true events !barrel event only!;#DeltaR_{offline}^{mid};rate");
  }
  m_peff_morecand->Draw();
  lef_select_morecand->Draw("same");
  c1 -> Print(pdf, "pdf" );
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2spall, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2spall = new TEfficiency(*m_h_deltaExtR_2spall, *m_h_deltaExtR_off_tot_clus);
    m_peff_2spall->SetTitle("Ratio of N_{spall}=2 events;#DeltaR_{offline}^{all};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2spinn, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2spinn = new TEfficiency(*m_h_deltaExtR_2spinn, *m_h_deltaExtR_off_tot_clus);
    m_peff_2spinn->SetTitle("Ratio of N_{spinn}=2 events;#DeltaR_{offline}^{inn};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2spmid, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2spmid = new TEfficiency(*m_h_deltaExtR_2spmid, *m_h_deltaExtR_off_tot_clus);
    m_peff_2spmid->SetTitle("Ratio of N_{spmid}=2 events;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2spout, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2spout = new TEfficiency(*m_h_deltaExtR_2spout, *m_h_deltaExtR_off_tot_clus);
    m_peff_2spout->SetTitle("Ratio of N_{spout}=2 events;#DeltaR_{offline}^{out};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_3pt, *m_h_deltaExtR_off_tot_clus)){
    m_peff_3pt = new TEfficiency(*m_h_deltaExtR_3pt, *m_h_deltaExtR_off_tot_clus);
    m_peff_3pt->SetTitle(";#DeltaR_{#mu#mu} at MuonSpectrometer;ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2pt, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2pt = new TEfficiency(*m_h_deltaExtR_2pt, *m_h_deltaExtR_off_tot_clus);
    m_peff_2pt->SetTitle(";#DeltaR_{#mu#mu} at MuonSpectrometer;ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_1pt, *m_h_deltaExtR_off_tot_clus)){
    m_peff_1pt = new TEfficiency(*m_h_deltaExtR_1pt, *m_h_deltaExtR_off_tot_clus);
    m_peff_1pt->SetTitle(";#DeltaR_{#mu#mu} at MuonSpectrometer;ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_2pt_oppcharge, *m_h_deltaExtR_off_tot_clus)){
    m_peff_2pt_oppcharge = new TEfficiency(*m_h_deltaExtR_2pt_oppcharge, *m_h_deltaExtR_off_tot_clus);
    m_peff_2pt_oppcharge->SetTitle("Ratio of N_{pT}=2 events;#DeltaR_{offline}^{mid};ratio");
  }
  m_peff_2pt_oppcharge->Write();
  m_peff_2pt->Write();
  m_h_ptres_2off2clus->Write();
  
  THStack *hs_lowpt = new THStack("hs_lowpt", "N_{superpoint} condition when much lower p_{T} was calculated");
  m_h_n_lowptclus->SetBinContent(2, 0);
  m_h_n_lowptclus->SetFillColor(kBlack);
  m_h_n_lowptclus->SetFillStyle(3002);
  hs_lowpt->Add(m_h_n_lowptclus);
  m_h_n_1spstat_inn->SetFillColor(kRed);
  m_h_n_1spstat_inn->SetLineColor(kRed);
  m_h_n_1spstat_inn->SetFillStyle(3002);
  hs_lowpt->Add(m_h_n_1spstat_inn);
  m_h_n_1spstat_mid->SetFillColor(kBlue);
  m_h_n_1spstat_mid->SetLineColor(kBlue);
  m_h_n_1spstat_mid->SetFillStyle(3003);
  hs_lowpt->Add(m_h_n_1spstat_mid);
  m_h_n_1spstat_out->SetFillColor(kGreen);
  m_h_n_1spstat_out->SetLineColor(kGreen);
  m_h_n_1spstat_out->SetFillStyle(3004);
  hs_lowpt->Add(m_h_n_1spstat_out);
  hs_lowpt->Draw();
//  m_h_n_lowptclus->Draw();
//  m_h_n_1spstat_out->Draw("same");
//  m_h_n_1spstat_mid->Draw("same");
//  m_h_n_1spstat_inn->Draw("same");
  c1 -> Print(pdf, "pdf" );

  m_hh_Offpt_2mu1RoI->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  m_hh_Offpt_2mu1RoI_oppcharge->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  if(TEfficiency::CheckConsistency(*m_h_Offpt_2mu1RoI_lead_oppcharge, *m_h_Offpt_2mu1RoI_lead)){
    m_peff_Offpt_lead = new TEfficiency(*m_h_Offpt_2mu1RoI_lead_oppcharge, *m_h_Offpt_2mu1RoI_lead);
    m_peff_Offpt_lead->SetTitle("Ratio of N_{pT}=2 && opposite charge events;Offline p_{T}[GeV];ratio");
  }
  m_peff_Offpt_lead->Draw();
  c1 -> Print(pdf, "pdf" );
  
  THStack *hs_lowpt_patA = new THStack("hs_lowpt_patA", "Why much lower p_{T} was calculated when all station has 2 or more SPs");
  m_h_n_lowptclus_patA_findout->GetYaxis()->SetRangeUser(0, 1200);
  m_h_n_lowptclus_patA_sameIDout->GetYaxis()->SetRangeUser(0, 1200);
  m_h_n_lowptclus_patA->GetYaxis()->SetRangeUser(0, 1200);
  m_h_n_lowptclus_patA_sameIDout->SetFillColor(kBlue);
  hs_lowpt_patA->Add(m_h_n_lowptclus_patA_sameIDout);
  m_h_n_lowptclus_patA_findout->SetFillColor(kGreen);
  m_h_n_lowptclus_patA->SetFillColor(kGray);
  m_h_n_lowptclus_patA->Draw();
  hs_lowpt_patA->Draw("same");
  m_h_n_lowptclus_patA_findout->Draw("same");
  c1 -> Print(pdf, "pdf" );

  m_peff_2spall->Draw();
  c1 -> Print(pdf, "pdf" );
  m_peff_2spinn->Draw();
  c1 -> Print(pdf, "pdf" );
  m_peff_2spmid->Draw();
  c1 -> Print(pdf, "pdf" );
  m_peff_2spout->Draw();
  c1 -> Print(pdf, "pdf" );
  m_peff_2pt_oppcharge->SetLineColor(kRed);
  m_peff_2pt_oppcharge->SetMarkerColor(kRed);
//  lef_pteff->AddOption(m_peff_2pt, "N_{pT}=2"); 
//  lef_pteff->AddOption(m_peff_2pt_oppcharge, "N_{pT}=2&&opposite charge");
  m_peff_3pt->SetLineColor(kViolet+1);
  m_peff_3pt->SetMarkerColor(kViolet+1);
  m_peff_3pt->SetMarkerStyle(21);
  m_peff_3pt->SetMarkerSize(2);
  m_peff_3pt->Draw();
  gPad->Update(); 
  auto graph_2pt = m_peff_3pt->GetPaintedGraph(); 
  graph_2pt->SetMinimum(0);
  graph_2pt->SetMaximum(1.05); 
  gPad->Update(); 
  m_peff_1pt->SetLineColor(kAzure+7);
  m_peff_1pt->SetMarkerColor(kAzure+7);
  m_peff_1pt->SetMarkerStyle(24);
  m_peff_1pt->Draw("same");
//  m_peff_2pt_oppcharge->Draw("same");
//  ATLASLabel(0.2, 0.5, "Work In Progress");
//  myText(0.2, 0.46, 1, "Simulation #sqrt{s}=13TeV");
  ATLASLabel(0.45, 0.9, "Work In Progress");
  myText(0.45, 0.85, 1, "Simulation #sqrt{s}=13TeV");
  TLegend_addfunc *lef_pteff;
  lef_pteff = new TLegend_addfunc(4);
//  lef_pteff = new TLegend_addfunc(0, 0, 4);
  lef_pteff->AddSelection("J/#psi#rightarrow#mu#mu(particle gun)");
  lef_pteff->AddSelection("|#eta_{#mu}| < 1.05");
  lef_pteff->AddEntry(m_peff_3pt, "N_{#mu, mtSA} > 2");
  lef_pteff->AddEntry(m_peff_1pt, "N_{#mu, mtSA} < 2");
//  myText(0.2, 0.42, 1, "J/#psi#rightarrow#mu#mu (particle gun)");
//  myText(0.2, 0.38, 1, "|#eta_{muon}|<1.05 (Barrel)");
  lef_pteff->Draw("same");
  c1 -> Print(pdf, "pdf" );
  m_h_ptres_2off2clus_leading->SetLineColor(kRed);
  m_h_ptres_2off2clus_sub->SetLineColor(kBlue);
  TLegend_addfunc *lef_ptres;
  lef_ptres = new TLegend_addfunc(1, 3);
  lef_ptres->AddOption(m_h_ptres_2off2clus, "total"); 
  lef_ptres->AddOption(m_h_ptres_2off2clus_leading, "leading"); 
  lef_ptres->AddOption(m_h_ptres_2off2clus_sub, "subleading");
  long entries_ptres = m_h_ptres_2off2clus->GetMaximum();
  m_h_ptres_2off2clus->GetYaxis()->SetRangeUser(0, entries_ptres*1.25);
  m_h_ptres_2off2clus->SetLineColor(38);
  m_h_ptres_2off2clus->SetMarkerColor(38);
  m_h_ptres_2off2clus->SetFillStyle(3004);
  m_h_ptres_2off2clus->SetFillColor(38);
  m_h_ptres_2off2clus->Draw();
  ATLASLabel(0.55, 0.9, "Work In Progress");
  myText(0.55, 0.85, 1, "Simulation #sqrt{s}=13TeV");
  myText(0.55, 0.8, 1, "J/#psi#rightarrow#mu#mu (particle gun)");
  myText(0.55, 0.75, 1, "|#eta_{#mu}|<1.05 (Barrel)");
//  m_h_ptres_2off2clus_leading->Draw("same");
//  m_h_ptres_2off2clus_sub->Draw("same");
//  lef_ptres->Draw("same");
  c1 -> Print(pdf, "pdf" );
  m_nclus_vs_npt->SetMarkerSize(1.5);
  m_nclus_vs_npt->SetMarkerColor(kRed);
  m_nclus_vs_npt->Draw("text");
  ATLASLabel(0.55, 0.9, "Work In Progress");
  myText(0.55, 0.85, 1, "Simulation #sqrt{s}=13TeV");
  myText(0.55, 0.8, 1, "J/#psi#rightarrow#mu#mu (particle gun)");
  myText(0.55, 0.75, 1, "|#eta_{#mu}|<1.05 (Barrel)");
  
  c1 -> Print(pdf, "pdf" );
  
  std::cout << "test1muon/1pad" << std::endl; 
  std::cout << "tail total/sa/mpsa = " << m_total << "/" << m_satail << "/" << m_mpsatail << std::endl;
  std::cout << "test2muon/1pad" << std::endl; 
  std::cout << "tail total/sa/mpsa = " << m_total_closeby << "/" << m_satail << "/" << m_mpsatail_closeby << std::endl;
  long maxEntry = m_ptres_default->GetMaximum();
  m_ptres_default->SetLineColor(kAzure+5);
  m_ptres_default->SetMarkerColor(kAzure+5);
  m_ptres_default->GetYaxis()->SetRangeUser(0, maxEntry*1.2);
  m_ptres_default->Draw();
  //m_ptres_default->Fit("gaus");
  m_ptres_mpsa->SetLineStyle(2);
  m_ptres_mpsa->SetLineColor(kRed+1);
  m_ptres_mpsa->SetMarkerColor(kRed+1);
  //m_ptres_mpsa->SetFillColorAlpha(kBlue+5, 0.35);
  m_ptres_mpsa->Draw("same");
  //m_ptres_mpsa->Fit("gaus");
  ATLASLabel(0.55, 0.9, "Work In Progress");
  myText(0.55, 0.85, 1, "Simulation #sqrt{s}=13TeV");
  TLegend_addfunc* leg_resi = new TLegend_addfunc(1,0,4);
  leg_resi->AddSelection("J/#psi#rightarrow#mu#mu (particle gun)");
  leg_resi->AddSelection("|#eta_{muon}|<1.05 (Barrel)");
  leg_resi->AddEntry(m_ptres_default, "MuonSA");
  leg_resi->AddEntry(m_ptres_mpsa, "mtSA");
  leg_resi->Draw("same");
  c1->Print(pdf, "pdf");
  
  m_nclusPerLayer[0]->GetYaxis()->SetRangeUser(0, m_nclusPerLayer[1]->GetMaximum()*1.2);
  m_nclusPerLayer[0]->Draw();
  m_nclusPerLayer[1]->SetLineColor(kRed+1);
  m_nclusPerLayer[1]->SetMarkerColor(kRed+1);
  m_nclusPerLayer[1]->SetLineStyle(2);
  m_nclusPerLayer[1]->Draw("hist same");
  m_nclusPerLayer[2]->SetLineColor(kMagenta+1);
  m_nclusPerLayer[2]->SetMarkerColor(kMagenta+1);
  m_nclusPerLayer[2]->SetLineStyle(3);
  m_nclusPerLayer[2]->Draw("hist same");
  m_nclusPerLayer[3]->SetLineColor(kBlue+1);
  m_nclusPerLayer[3]->SetMarkerColor(kBlue+1);
  m_nclusPerLayer[3]->SetLineStyle(4);
  m_nclusPerLayer[3]->Draw("hist same");
  ATLASLabel(0.45, 0.8, "Work In Progress");
  myText(0.45, 0.75, 1, "Simulation #sqrt{s}=13TeV");
  TLegend_addfunc* leg_hits = new TLegend_addfunc(6);
  leg_hits->AddSelection("J/#psi#rightarrow#mu#mu (particle gun)");
  leg_hits->AddSelection("|#eta_{muon}|<1.05 (Barrel)");
  leg_hits->AddEntry(m_nclusPerLayer[0], "N_{SA muons}=1");
  leg_hits->AddEntry(m_nclusPerLayer[1], "N_{SA muons}=2");
  leg_hits->AddEntry(m_nclusPerLayer[2], "N_{SA muons}=3");
  leg_hits->AddEntry(m_nclusPerLayer[3], "N_{SA muons}=4");
  leg_hits->Draw("same");
  c1 -> Print(pdf, "pdf" );
//  m_h_ptres_2off2clus->Fit("gaus");
  m_h_ptclus_vs_off->Draw("colZ");
  ATLASLabel(0.2, 0.9, "Work In Progress");
  myText(0.2, 0.85, 1, "Simulation #sqrt{s}=13TeV");
  myText(0.2, 0.8, 1, "J/#psi#rightarrow#mu#mu (particle gun)");
  myText(0.2, 0.75, 1, "|#eta_{#mu}|<1.05 (Barrel)");
  c1 -> Print(pdf, "pdf" );
  m_h_ptres_vs_offsegInndR->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  m_h_ptres_vs_offsegMiddR->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  m_h_ptres_vs_offsegOutdR->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  m_h_ptres_vs_offsegMindR->Draw("colZ");
  c1 -> Print(pdf, "pdf" );

  
  TLegend_addfunc *lef_select_clusroad_eff;
  lef_select_clusroad_eff = new TLegend_addfunc(5);
  lef_select_clusroad_eff->AddSelection("Off-p_{T}#geq4GeV");
  //lef_select_clusroad_eff->AddSelection("Off-RoI |#eta|<1.05");
  lef_select_clusroad_eff->AddSelection("Off-Segment |#eta|<1.05");
  lef_select_clusroad_eff->AddSelection("Off-|z_{0}|#leq50(mm)");
  lef_select_clusroad_eff->AddSelection("N_{OffsegMid}#geq1");
  lef_select_clusroad_eff->AddSelection("N_{OffsegOut}#geq1");
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_off_2road_clus, *m_h_deltaExtR_off_tot_clus)){
    m_peff_clusRoad = new TEfficiency(*m_h_deltaExtR_off_2road_clus, *m_h_deltaExtR_off_tot_clus);
    m_peff_clusRoad->SetTitle("Ratio of N_{clusRoad}=2 events;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_off_1road_clus, *m_h_deltaExtR_off_tot_clus)){
    m_peff_1clusRoad = new TEfficiency(*m_h_deltaExtR_off_1road_clus, *m_h_deltaExtR_off_tot_clus);
    m_peff_1clusRoad->SetTitle("Ratio of N_{clusRoad}=2 events;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_off_manyroad_clus, *m_h_deltaExtR_off_tot_clus)){
    m_peff_manyclusRoad = new TEfficiency(*m_h_deltaExtR_off_manyroad_clus, *m_h_deltaExtR_off_tot_clus);
    m_peff_manyclusRoad->SetTitle("Ratio of N_{clusRoad}=2 events;#DeltaR_{offline}^{mid};ratio");
  }
  //m_peff_clusRoad->SetMinimum(0);
  m_peff_clusRoad->SetLineColor(kRed);
  m_peff_1clusRoad->SetLineColor(kBlue);
  m_peff_manyclusRoad->SetLineColor(kGreen);
  m_peff_clusRoad->SetMarkerColor(kRed);
  m_peff_1clusRoad->SetMarkerColor(kBlue);
  m_peff_manyclusRoad->SetMarkerColor(kGreen);
  TLegend_addfunc *lopt_eff_nclus;
  lopt_eff_nclus = new TLegend_addfunc(5, 3);
  lopt_eff_nclus->AddOption(m_peff_clusRoad, "N_{clusRoad} = 2");
  lopt_eff_nclus->AddOption(m_peff_1clusRoad, "N_{clusRoad} = 1");
  lopt_eff_nclus->AddOption(m_peff_manyclusRoad, "N_{clusRoad} > 2");
  //m_h_eff_clusterRoad->Draw("P");
  m_peff_1clusRoad->Draw();
  m_peff_1clusRoad->Write();
  gPad->Update(); 
  auto graph = m_peff_1clusRoad->GetPaintedGraph(); 
  graph->SetMinimum(0);
  graph->SetMaximum(1); 
  gPad->Update(); 
  m_peff_clusRoad->Draw("same");
  m_peff_manyclusRoad->Draw("same");
  lopt_eff_nclus->Draw("same");
  lef_select_clusroad_eff->Draw("same");
  m_peff_clusRoad->Write();
  c1 -> Print(pdf, "pdf" );
//  c1->RedrawAxis();
  m_peff_2spmid->Draw();
  m_peff_clusRoad->Draw("same");
  c1 -> Print(pdf, "pdf" );
  m_peff_2spout->Draw();
  c1 -> Print(pdf, "pdf" );
  
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_canthelp, *m_h_dR_1road_sep)){
    m_peff_canthelp = new TEfficiency(*m_h_dR_1road_canthelp, *m_h_dR_1road_sep);
    m_peff_canthelp->SetTitle("Ratio N_{canthelp!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_1clus, *m_h_dR_1road_canthelp)){
    m_peff_1clus = new TEfficiency(*m_h_dR_1road_1clus, *m_h_dR_1road_canthelp);
    m_peff_1clus->SetTitle("Ratio N_{1clus!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_1phi, *m_h_dR_1road_canthelp)){
    m_peff_1phi = new TEfficiency(*m_h_dR_1road_1phi, *m_h_dR_1road_canthelp);
    m_peff_1phi->SetTitle("Ratio N_{1phi!}/N_{phiRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_phiothers, *m_h_dR_1road_canthelp)){
    m_peff_phiothers = new TEfficiency(*m_h_dR_1road_phiothers, *m_h_dR_1road_canthelp);
    m_peff_phiothers->SetTitle("Ratio N_{phiothers!}/N_{phiRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_outlier, *m_h_dR_1road_sep)){
    m_peff_outlier = new TEfficiency(*m_h_dR_1road_outlier, *m_h_dR_1road_sep);
    m_peff_outlier->SetTitle("Ratio N_{outlier!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_misdire, *m_h_dR_1road_sep)){
    m_peff_misdire = new TEfficiency(*m_h_dR_1road_misdire, *m_h_dR_1road_sep);
    m_peff_misdire->SetTitle("Ratio N_{misdire!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_dR_1road_others, *m_h_dR_1road_sep)){
    m_peff_others = new TEfficiency(*m_h_dR_1road_others, *m_h_dR_1road_sep);
    m_peff_others->SetTitle("Ratio N_{others!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
  }
  TEfficiency* peff_merge;
//  peff_merge = 0;
//  TH1D *hs = new TH1D("hs", "test;#DeltaR;events", 30, 0, 0.3);
//  hs->Add(m_h_dR_1road_outlier);
//  hs->Add(m_h_dR_1road_canthelp);
////  hs->Add(m_h_dR_1road_misdire);
//  hs->Add(m_h_dR_1road_others);
//  if(TEfficiency::CheckConsistency(*hs, *m_h_dR_1road_sep)){
//    peff_merge = new TEfficiency(*hs, *m_h_dR_1road_sep);
//    peff_merge->SetTitle("Ratio N_{canthelp!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
//  }
////  TList* list = new TList();
////  list->Add(m_peff_canthelp);
////  list->Add(m_peff_misdire);
////  list->Add(m_peff_outlier);
////  peff_merge->Merge(list);
////  peff_merge->Add(*m_peff_outlier);
////  peff_merge->Add(*m_peff_misdire);
//  TLegend_addfunc *lopt_eff_1road;
//  lopt_eff_1road = new TLegend_addfunc(5, 5);
//  lopt_eff_1road->AddOption(m_peff_canthelp, "only 1 clus by 1mu");
//  lopt_eff_1road->AddOption(m_peff_outlier, "by finding same clus in lay1");
////  lopt_eff_1road->AddOption(m_peff_misdire, "by ignoring lay2 clus");
//  lopt_eff_1road->AddOption(m_peff_others, "others");
//  lopt_eff_1road->AddOption(peff_merge, "total");
//  m_peff_canthelp->SetLineColor(kViolet);
//  m_peff_outlier->SetLineColor(kRed);
////  m_peff_misdire->SetLineColor(kBlue);
//  m_peff_others->SetLineColor(kGreen);
//  peff_merge->SetLineColor(kBlack);
//  m_peff_canthelp->SetMarkerColor(kViolet);
//  m_peff_outlier->SetMarkerColor(kRed);
//  m_peff_misdire->SetMarkerColor(kBlue);
//  m_peff_others->SetMarkerColor(kGreen);
//  peff_merge->SetMarkerColor(kBlack);
//  m_peff_canthelp->Draw();
//  m_peff_canthelp->Write();
//  m_peff_others->Draw("same");
//  peff_merge->Draw("same");
//  m_peff_canthelp->Write();
//  m_peff_outlier->Draw("same");
//  m_peff_outlier->Write();
//  //m_peff_misdire->Draw("same");
//  lopt_eff_1road->Draw("same");
//  m_peff_misdire->Write();
//  c1 -> Print(pdf, "pdf" );
  
//  TEfficiency* peff_merge_phi;
//  peff_merge_phi = 0;
//  TH1D *hphi = new TH1D("hphi", "test;#DeltaR;events", 30, 0, 0.3);
//  hphi->Add(m_h_dR_1road_1clus);
//  hphi->Add(m_h_dR_1road_1phi);
//  hphi->Add(m_h_dR_1road_phiothers);
//  if(TEfficiency::CheckConsistency(*hphi, *m_h_dR_1road_canthelp)){
//    peff_merge_phi = new TEfficiency(*hphi, *m_h_dR_1road_canthelp);
//    peff_merge_phi->SetTitle("Ratio N_{canthelp!}/N_{clusRoad}=1;#DeltaR_{offline}^{mid};ratio");
//  }
//  m_peff_1clus->SetLineColor(kBlue);
//  m_peff_1clus->SetMarkerColor(kBlue);
//  m_peff_1clus->Draw();
//  gPad->Update(); 
//  auto g_1clus = m_peff_1clus->GetPaintedGraph(); 
//  g_1clus->SetMinimum(0);
//  g_1clus->SetMaximum(1); 
//  gPad->Update(); 
//  m_peff_1phi->SetLineColor(kRed);
//  m_peff_1phi->SetMarkerColor(kRed);
//  m_peff_1phi->Draw("same");
//  m_peff_phiothers->SetLineColor(kGreen);
//  m_peff_phiothers->SetMarkerColor(kGreen);
//  m_peff_phiothers->Draw("same");
//  peff_merge_phi->SetLineColor(kBlack);
//  peff_merge_phi->SetMarkerColor(kBlack);
//  peff_merge_phi->Draw("same");
//  TLegend_addfunc *lopt_eff_1phi;
//  lopt_eff_1phi = new TLegend_addfunc(5, 4);
//  lopt_eff_1phi->AddOption(m_peff_1clus, "1clus per 1lay");
//  lopt_eff_1phi->AddOption(m_peff_1phi, "1 middle phi");
//  lopt_eff_1phi->AddOption(m_peff_phiothers, "others");
//  lopt_eff_1phi->AddOption(peff_merge_phi, "total");
//  lopt_eff_1phi->Draw("same");
//  c1 -> Print(pdf, "pdf" );
  m_hh_RoIEtaPhi_1road->Draw("colZ");
  c1 -> Print(pdf, "pdf" );
  
  if(TEfficiency::CheckConsistency(*m_h_dR_3road_3clus, *m_h_dR_3road)){
    m_peff_3road = new TEfficiency(*m_h_dR_3road_3clus, *m_h_dR_3road);
    m_peff_3road->SetTitle("Ratio N_{3clus in 2layers}/N_{clusRoad}>2;#DeltaR_{offline}^{mid};ratio");
  }
  m_peff_3road->Draw();
  c1 -> Print(pdf, "pdf" );
  
  if(TEfficiency::CheckConsistency(*m_h_2CR_1SPCMid_u3mdt, *m_h_2CR_1SPCMid)){
    m_peff_2CR_1SPCMid = new TEfficiency(*m_h_2CR_1SPCMid_u3mdt, *m_h_2CR_1SPCMid);
    m_peff_2CR_1SPCMid->SetTitle("Ratio N_{< 3 mdts}/N_{clusRoad}=2&&N_{sp}<2;#DeltaR_{offline}^{mid};ratio");
  }
  if(TEfficiency::CheckConsistency(*m_h_2CR_1SPCMid_u3mdt_selected, *m_h_2CR_1SPCMid)){
    m_peff_2CR_1SPCMid_all = new TEfficiency(*m_h_2CR_1SPCMid_u3mdt_selected, *m_h_2CR_1SPCMid);
    m_peff_2CR_1SPCMid_all->SetTitle("Ratio N_{< 3 mdts}/N_{clusRoad}=2&&N_{sp}<2;#DeltaR_{offline}^{Mid_all};ratio");
  }
  m_peff_2CR_1SPCMid_all->SetLineColor(kRed);
  m_peff_2CR_1SPCMid_all->SetMarkerColor(kRed);
  m_peff_2CR_1SPCMid->Draw();
  m_peff_2CR_1SPCMid_all->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  
  if(TEfficiency::CheckConsistency(*m_h_deltaExtR_off_closeBy, *m_h_deltaExtR_off_tot_clus)){
    m_peff_closeBy = new TEfficiency(*m_h_deltaExtR_off_closeBy, *m_h_deltaExtR_off_tot_clus);
    m_peff_closeBy->SetTitle("Ratio of 2muon in 1RoI tower events;#DeltaR_{offline}^{mid};ratio");
  }
  //m_h_eff_clusterRoad->Draw("P");
  //m_peff_closeBy->GetXaxis()->SetRangeUser(0, 0.3):
  m_peff_closeBy->Draw();
  m_peff_clusRoad->SetLineColor(kRed);
  m_peff_clusRoad->SetMarkerColor(kRed);
  m_peff_clusRoad->Draw("same");
  //lef_select_clusroad_eff->Draw("same");
  TLegend_addfunc *lopt_eff_closeBy = new TLegend_addfunc(3,2);
  lopt_eff_closeBy->AddOption(m_peff_closeBy, "2muon in 1RoI evt ratio");
  lopt_eff_closeBy->AddOption(m_peff_clusRoad, "N_{clusRoadMid}=2 evt ratio");
  ATLASLabel(0.6, 0.5, "Work In Progress", 1);
  //lopt_eff_closeBy->Draw("same");
  c1 -> Print(pdf, "pdf" );
  
  m_h_n_spall->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_sector->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_Nclus_npt1->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_countPtclus->Draw();
  c1 -> Print(pdf, "pdf" );
  m_h_offlinePt_res0p7->Draw();
  c1 -> Print(pdf, "pdf" );
  fout->Close();
  c1 -> Print( pdf + "]", "pdf" );
  delete c1;
}

void RPC_FCBM::CalcEfficiency(TH1D* h_num, TH1D* h_den, TH1D* h_set){
  double eff_x, eff_y, eff_xerr, eff_yerr;
  double eff_y_total = 0;
  double closeTotal = 0;
  double road2total = 0;
  int nbins = h_den->GetXaxis()->GetNbins();
  cout << "DrawEfficiency::the efficiency parameter" << endl;
  cout << "eff_y/eff_yerr" << "   " << "numerator/denominator" << endl;
  for(unsigned int ibin = 1; ibin <= nbins; ibin++){
    double denominator = h_den->GetBinContent(ibin);
    double numerator = h_num->GetBinContent(ibin);
    eff_x = h_den->GetBinCenter(ibin);
    eff_y = (denominator != 0) ? numerator/denominator : 0;
    eff_yerr = (denominator != 0) ? sqrt((numerator*(1.-2*eff_y)+(denominator*pow(eff_y, 2)))/pow(denominator, 2)):0;
    cout << eff_y << " / " << eff_yerr << "  " << numerator << "/" << denominator << endl;
    h_set->SetBinContent(ibin, eff_y);
    h_set->SetBinError(ibin, eff_yerr);
    eff_y_total = eff_y_total + eff_y;
    closeTotal = closeTotal + denominator;
    road2total = road2total + numerator;
  }//ibin loop end
  cout << "eff_y total : " << eff_y_total << endl;
  cout << "eff_allregion : " << road2total/closeTotal << endl;
}
