#ifndef RPCL2MUONSA_HISTDATA
#define RPCL2MUONSA_HISTDATA

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TFile.h"

namespace RpcL2MuonSA {

class HistData
{
  public :
    HistData(){
      HistInit();
    };
    ~HistData() {};
    void HistInit();
    void HistEnd();

  public : 
    TH2D *m_h_ptOffvsSAclus;
    TH1D *m_h_OffdR;
    TH1D *m_h_deta_SAvsOFF;
    TH1D *m_h_ptres_deta[5];
    TH1D *m_h_NrpcclusfitPhi;
    TH1D *m_h_NclusterPerLayer[8];
    TH1D *m_h_clusdR;
    TH1D *m_h_clusinvMass;
    TH1D *m_h_clusinvMass_rpcclusphi1;
    TH1D *m_h_clusinvMass_rpcclusphiovr1;
    TH1D *m_h_cluscharge;
    TH1D *m_h_n_passedRatio;
    TH1D *m_h_n_passedRatio_exclusive;
    TH1D *m_h_n_passedL1;
    TH1D *m_h_n_passedMoreCand;
    TH1D *m_h_n_passedSA;
    TH1D *m_h_n_passedSAHypo;
    TH1D *m_h_n_passedSAOvRm_forClus;
    TH1D *m_h_n_passedCB;
    TH1D *m_h_n_passedEF;
    TH1D *m_dR_ismore;
    TH1D *m_dR_sapass;
};
}

#endif //RPCL2MUONSA_HISTDATA

#ifdef RPC_FCBM_cxx
#ifndef RPCL2MUONSA_HISTDATA_hhh 
#define RPCL2MUONSA_HISTDATA_hhh 
namespace RpcL2MuonSA {

void HistData::HistInit()
{
  // ===== fill in 2muinBarrel =====
  //dR check
  m_h_OffdR = new TH1D("m_h_OffdR", "Offline extdR;#DeltaR_{#mu#mu} at MuonSpectrometer;events", 50, 0, 0.25);
  m_dR_ismore = new TH1D("m_dR_ismore", "dR_{MS}^{offline} isMoreCand true&L1_MU4 pass;#DeltaR_{MS};events", 40, 0, 0.4);
  m_dR_sapass = new TH1D("m_dR_sapass", "dR_{MS}^{offline} close-by alg success;#DeltaR_{MS};events", 40, 0, 0.4);

  // ===== fill in 2muin1RoI =====
  m_h_ptOffvsSAclus  = new TH2D("m_h_ptOffvsSAclus", "p_{T}^{offline} vs p_{T}^{cluster}:leading & subleading;p_{T}^{SA, cluster}[GeV];p_{T}^{offline}[GeV];events", 200 ,0 ,1000, 250, 0, 500);
  for(int i = 0; i < 8; i++){
    m_h_NclusterPerLayer[i] = new TH1D(Form("m_h_NclusterPerLayer%d", i), Form("the number of clusters at layer %d;N;events", i), 10, -0.5, 9.5);
  }
  for(int i = 0; i < 5; i++){
    m_h_ptres_deta[i] = new TH1D(Form("m_h_ptres_deta%d",i), Form("pt residual:%d;residual;events", i), 100, -1, 1);
    m_h_ptres_deta[i]->Sumw2();
  }
  m_h_deta_SAvsOFF = new TH1D("m_h_deta_SAvsOFF", "|#Delta#eta| : Offline vs SA;|#Delta#eta|;events", 50, 0, 0.1);
  m_h_clusdR = new TH1D("m_h_clusdR", "cluster dR in MS;#DeltaR;events", 80, 0, 0.4);
  m_h_NrpcclusfitPhi = new TH1D("m_h_NrpcclusfitPhi", "N_#rpc cluster fit Phi; N; events", 9, -0.5, 8.5);

  // inv mass check
  m_h_clusinvMass = new TH1D("m_h_clusinvMass", "cluster invMass;Mass (GeV);events", 50, 0, 0.5);
  m_h_clusinvMass_rpcclusphi1 = new TH1D("m_h_clusinvMass_rpcclusphi1", "cluster invMass (N_{rpccluster fit phi} = 1);Mass (GeV);events", 100, 0, 20);
  m_h_clusinvMass_rpcclusphiovr1 = new TH1D("m_h_clusinvMass_rpcclusphiovr1", "cluster invMass (N_{rpccluster fit phi} > 1);Mass (GeV);events", 100, 0, 20);
  //charge 
  m_h_cluscharge = new TH1D("m_h_cluscharge", "is opposite charge;;events", 2, 0, 2);
  m_h_cluscharge->GetXaxis()->SetBinLabel(1, "opposite");
  m_h_cluscharge->GetXaxis()->SetBinLabel(2, "same");
  //trigger passed ratio
  m_h_n_passedRatio = new TH1D("m_h_n_passedRatio", "Trigger pass events HLT_mu26ivm, L1MU20 (2mu-in-1RoI evt only); step; events", 4, 0, 4);
  m_h_n_passedRatio->GetXaxis()->SetBinLabel(1, "L1");
  m_h_n_passedRatio->GetXaxis()->SetBinLabel(2, "isMoreCand");
  m_h_n_passedRatio->GetXaxis()->SetBinLabel(3, "SAHypo");
  m_h_n_passedRatio->GetXaxis()->SetBinLabel(4, "SAOvRm");

  m_h_n_passedRatio_exclusive = new TH1D("m_h_n_passedRatio_exclusive", "Trigger pass events HLT_mu26ivm, L1MU20 (barrel evt only); step; events", 4, 0, 4);
  m_h_n_passedRatio_exclusive->GetXaxis()->SetBinLabel(1, "L1");
  m_h_n_passedRatio_exclusive->GetXaxis()->SetBinLabel(2, "isMoreCand");
  m_h_n_passedRatio_exclusive->GetXaxis()->SetBinLabel(3, "SAHypo");
  m_h_n_passedRatio_exclusive->GetXaxis()->SetBinLabel(4, "SAOvRm");

  //passed events ===> 1bin : 2mu-in-1RoI, 2bin : 2mu-in-barrel(exclusive)
  m_h_n_passedL1 = new TH1D("m_h_n_passedL1", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedMoreCand = new TH1D("m_h_n_passedMoreCand", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedSA = new TH1D("m_h_n_passedSA", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedSAHypo = new TH1D("m_h_n_passedSAHypo", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedSAOvRm_forClus = new TH1D("m_h_n_passedSAOvRm_forClus", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedCB = new TH1D("m_h_n_passedCB", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
  m_h_n_passedEF = new TH1D("m_h_n_passedEF", "Relative Trigger Efficiency (2mu-in-1RoI evt only); step; #epsilon", 3, 0, 3);
}

void HistData::HistEnd()
{
  for(int i = 0; i<5; i++){
    if(m_h_ptres_deta[i]!=0) delete m_h_ptres_deta[i];
  }
  if(m_h_deta_SAvsOFF!=0) delete m_h_deta_SAvsOFF;
  if(m_h_OffdR != 0) delete m_h_OffdR;
  if(m_h_clusdR != 0) delete m_h_clusdR;
  if(m_h_clusinvMass != 0) delete m_h_clusinvMass;
  if(m_h_clusinvMass_rpcclusphi1 != 0) delete m_h_clusinvMass_rpcclusphi1;
  if(m_h_clusinvMass_rpcclusphiovr1 != 0) delete m_h_clusinvMass_rpcclusphiovr1;
  if(m_h_NrpcclusfitPhi != 0 ) delete m_h_NrpcclusfitPhi;
  if(m_h_n_passedRatio != 0) delete m_h_n_passedRatio;
  if(m_h_n_passedRatio_exclusive != 0) delete m_h_n_passedRatio_exclusive;
  if(m_h_n_passedL1 != 0) delete m_h_n_passedL1;
  if(m_h_n_passedMoreCand != 0) delete m_h_n_passedMoreCand;
  if(m_h_n_passedSA != 0) delete m_h_n_passedSA;
  if(m_h_n_passedSAHypo != 0) delete m_h_n_passedSAHypo;
  if(m_h_n_passedSAOvRm_forClus != 0) delete m_h_n_passedSAOvRm_forClus;
  if(m_h_n_passedCB != 0) delete m_h_n_passedCB;
  if(m_h_n_passedEF != 0) delete m_h_n_passedEF;
}
}

#endif //RPCL2MUONSA_HISTDATA
#endif //RPCL2MUONSA_HISTDATA
