#include "TH1D.h"
#include <iostream>
#include <TH1.h>
#include <TFile.h>
#include "/home/ktaniguc/RootUtils/RootUtils/TCanvas_opt.h"
#include "/home/ktaniguc/RootUtils/RootUtils/TLegend_addfunc.h"
#include "/home/ktaniguc/RootUtils/src/rootlogon.C"
#include "TEfficiency.h"

using namespace std;

void checkFile(){
  rootlogon(); 
  bool checkDelta = false;
  bool checkEfficiency = false;
  bool checkPtres = true;
  bool checkPtNum = false;
//  TFile *_file0 = TFile::Open("rpcCluster_JpsiCloseBy.root");
//  TFile *_file1 = TFile::Open("rpcCluster_JpsiCloseBy_renewAlgo.root");
//  TFile *_file0 = TFile::Open("rpcCluster_JpsiCloseBy_renewAlgo_202007171628.root");
//  TFile *_file1 = TFile::Open("residualWeight_norm30_m0p43_v2/rpcCluster_JpsiCloseBy_renewAlgo.root");
  TFile *_file1 = TFile::Open("separateMDTs/rpcCluster_JpsiCloseBy_renewAlgo.root");
  TFile *_file0 = TFile::Open("test/rpcCluster_JpsiCloseBy_renewAlgo.root");
//  TFile *_file0 = TFile::Open("residualWeight_norm16_m0p4_v3/rpcCluster_JpsiCloseBy_renewAlgo.root");
//  TFile *_file0 = TFile::Open("changerWidth/rpcCluster_JpsiCloseBy_renewAlgo.root");
//  TFile *_file1 = TFile::Open("rpcCluster_JpsiCloseBy_oldL2SA_202007201337.root");
  if(checkPtres){
    TCanvas_opt *c1 = new TCanvas_opt();
    string histname[1];
    //histname[0] = "m_h_SARoadtheta";
    //histname[0] = "m_h_SARoadoffset";
    //histname[0] = "m_h_dtheta_offvsroad";
    //histname[0] = "h_Offpt";
    histname[0] = "m_h_ptres_2off2clus";
    TH1D *h0 = (TH1D*)_file0->Get(histname[0].c_str());
    TH1D *h1 = (TH1D*)_file1->Get(histname[0].c_str());
    //  h0->Sumw2();
    //  h0->Scale(h1->GetEntries()/h0->GetEntries());
    c1->cd();
    //h0->GetYaxis()->SetRangeUser(0, 50);
    //h1->GetYaxis()->SetRangeUser(0, 50);
    h0->SetMarkerStyle(1);
    //c1->SetLogy();
    h1->SetLineColor(kBlue);

    h0->SetLineColor(kGreen+2);
    h0->Draw("hist");
    h1->Draw("same hist");
    TLegend_addfunc *legptres = new TLegend_addfunc(0, 0, 5);
    legptres->AddSelection("J/#psi->#mu#mu");
    legptres->AddSelection("2mu-in-1RoI");
    legptres->AddSelection("muon #eta_{MS} < 1.05");
//    legptres->AddEntry(h0, "Weight A=16 B=-0.4 fixed");
    //legptres->AddEntry(h0, "changerWidth");
//    legptres->AddEntry(h1, "Weight A=30 B=-0.43 fixed");
    legptres->AddEntry(h0, "separateMDTs #Delta#eta matching");
    legptres->AddEntry(h1, "separateMDTs");
    legptres->Draw("same");

    if(checkDelta){
      int nbins = h0->GetNbinsX();
      cout << "nbins = " << nbins << endl;
      double xmin = h0->GetXaxis()->GetXmin();
      double xmax = h0->GetXaxis()->GetXmax();
      cout << "xmin/xmax = " << xmin << "/" << xmax << endl;
      TCanvas_opt *c2 = new TCanvas_opt();
      c2->cd();
      //TH1D *subtr = new TH1D("subtr", "#Delta(events) : (after bugfix - before bugfix);#theta(rad);#Delta(events)", nbins, xmin, xmax);
      //TH1D *subtr = new TH1D("subtr", "#Delta(events) : (after bugfix - before bugfix);offset z(mm);#Delta(events)", nbins, xmin, xmax);
      //TH1D *subtr = new TH1D("subtr", "#Delta(events) : (after bugfix - before bugfix);#Delta#theta(rad);#Delta(events)", nbins, xmin, xmax);
      TH1D *subtr = new TH1D("subtr", "#Delta(events) : (after bugfix - before bugfix);1-p_{T}^{SA}/p_{T}^{offline};#Delta(events)", nbins, xmin, xmax);
      subtr->Add(h0, h1, -1, 1);
      //subtr->GetYaxis()->SetRangeUser(-1, 1);
      subtr->Draw();
    }
  }
  
  if(checkEfficiency){
    TCanvas_opt *c3 = new TCanvas_opt();
    TString num = "m_h_deltaExtR_off_tot_clus;1";
    TString den = "m_h_deltaExtR_2pt;1";
    TH1D *h_num_newL2SA = (TH1D*)_file0->Get(num);
    TH1D *h_den_newL2SA = (TH1D*)_file0->Get(den);
    TH1D *h_num_oldL2SA = (TH1D*)_file1->Get(num);
    TH1D *h_den_oldL2SA = (TH1D*)_file1->Get(den);
    if(h_den_newL2SA != nullptr){
      TEfficiency* m_peff_2pt_newL2SA = new TEfficiency(*h_den_newL2SA, *h_num_newL2SA);
      m_peff_2pt_newL2SA->SetTitle("Ratio of N_{pT}=2 events;#DeltaR_{offline}^{mid};ratio");
      TEfficiency* m_peff_2pt_oldL2SA = new TEfficiency(*h_den_oldL2SA, *h_num_oldL2SA);
      m_peff_2pt_oldL2SA->SetTitle("Ratio of N_{pT}=2 events;#DeltaR_{offline}^{mid};ratio");
      TLegend_addfunc *leg = new TLegend_addfunc("left lower", 0, 5);
      m_peff_2pt_newL2SA->SetMarkerColor(kGreen+2);
      m_peff_2pt_newL2SA->SetLineColor(kGreen+2);
      m_peff_2pt_oldL2SA->SetMarkerColor(kBlue);
      m_peff_2pt_oldL2SA->SetLineColor(kBlue);
      leg->AddSelection("J/#psi->#mu#mu");
      leg->AddSelection("2mu-in-1RoI");
      leg->AddSelection("muon #eta_{MS} < 1.05");
      leg->AddEntry(m_peff_2pt_newL2SA, "Weight A=16 B=-0.4 fixed");
      leg->AddEntry(m_peff_2pt_oldL2SA, "separateMDTs");
      m_peff_2pt_oldL2SA->Draw();
      m_peff_2pt_newL2SA->Draw("same");
      leg->Draw("same");
    }
  }

  if(checkPtNum){
    TCanvas_opt *c4 = new TCanvas_opt();
    TString hist = "m_h_countPtclus;1";
    TH1D *h_hist_newL2SA = (TH1D*)_file0->Get(hist);
    TH1D *h_hist_oldL2SA = (TH1D*)_file1->Get(hist);
    h_hist_oldL2SA->SetLineColor(kBlue);
    h_hist_newL2SA->SetLineColor(kRed);
    TLegend_addfunc *leg = new TLegend_addfunc(5);
    leg->AddSelection("J/#psi->#mu#mu");
    leg->AddSelection("2mu-in-1RoI");
    leg->AddSelection("muon #eta_{MS} < 1.05");
    leg->AddEntry(h_hist_newL2SA, "newL2SA");
    leg->AddEntry(h_hist_oldL2SA, "oldL2SA");
    h_hist_newL2SA->Draw();
    h_hist_oldL2SA->Draw("same");
    leg->Draw("same");
  }
}
  

