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
//  TFile *_file0 = TFile::Open("rpcCluster_JpsiCloseBy.root");
//  TFile *_file1 = TFile::Open("rpcCluster_JpsiCloseBy_renewAlgo.root");
  TFile *_file0 = TFile::Open("rpcCluster_JpsiCloseBy_renewAlgo_202007171628.root");
  TFile *_file1 = TFile::Open("rpcCluster_JpsiCloseBy_oldL2SA_202007201337.root");
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
    h1->Draw("hist e1");

    h0->SetLineColor(kRed);
    h0->Draw("same");
    TLegend_addfunc *legptres = new TLegend_addfunc(0, 0, 5);
    legptres->AddSelection("J/#psi->#mu#mu");
    legptres->AddSelection("2mu-in-1RoI");
    legptres->AddSelection("muon #eta_{MS} < 1.05");
    legptres->AddEntry(h0, "newL2SA");
    legptres->AddEntry(h1, "oldL2SA");
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
      TLegend_addfunc *leg = new TLegend_addfunc("Efficiency", 5);
      m_peff_2pt_newL2SA->SetMarkerColor(kRed);
      m_peff_2pt_newL2SA->SetLineColor(kRed);
      m_peff_2pt_oldL2SA->SetMarkerColor(kBlue);
      m_peff_2pt_oldL2SA->SetLineColor(kBlue);
      leg->AddSelection("J/#psi->#mu#mu");
      leg->AddSelection("2mu-in-1RoI");
      leg->AddSelection("muon #eta_{MS} < 1.05");
      leg->AddEntry(m_peff_2pt_newL2SA, "newL2SA");
      leg->AddEntry(m_peff_2pt_oldL2SA, "oldL2SA");
      m_peff_2pt_newL2SA->Draw();
      m_peff_2pt_oldL2SA->Draw("same");
      leg->Draw("same");
    }
  }

}
  

