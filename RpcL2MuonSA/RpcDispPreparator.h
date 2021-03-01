#ifndef RpcDispPreparator_h
#define RpcDispPreparator_h
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
#include <vector>
#include <iostream>
#include <stdlib.h>

class RpcDispPreparator {
  public :
    RpcDispPreparator(); 
    void SetParameter_rpcFit(//bool isTag, 
                             int NTrigChain, 
                             float* rpcFitInn, 
                             float* rpcFitMid, 
                             float* rpcFitOut,
                             float rpcFitMidPhi,   
                             float rpcFitMidSlope,
                             float rpcFitMidOffset,
                             float rpcFitOutPhi,
                             float rpcFitOutSlope,
                             float rpcFitOutOffset
                             );

    void RpcHitSetter(const Int_t& nRPC, 
                      int& NTrigChain, 
                      vector<vector<double>>*& mesSA_rpcHitStationNumber,
                      vector<vector<double>>*& mesSA_rpcHitR, 
                      vector<vector<double>>*& mesSA_rpcHitEta, 
                      vector<vector<double>>*& mesSA_rpcHitPhi, 
                      vector<vector<uint32_t>>*& mesSA_rpcHitLayer, 
                      vector<vector<uint32_t>>*& mesSA_rpcHitMeasuresPhi, 
                      vector < vector < float > >*& mesSA_rpcHitDistEtaRO, //ktaniguc trial add
                      vector < vector < float > >*& mesSA_rpcHitDistPhiRO,
                      vector < vector < float > >*& mesSA_rpcHitX,
                      vector < vector < float > >*& mesSA_rpcHitY,
                      vector < vector < float > >*& mesSA_rpcHitZ, //ktaniguc trial add end
                      vector<vector<float>>& setRpc1HitEta,
                      vector<vector<float>>& setRpc2HitEta,
                      vector<vector<float>>& setRpc3HitEta,
                      vector<vector<float>>& setRpc1HitPhi,
                      vector<vector<float>>& setRpc2HitPhi,
                      vector<vector<float>>& setRpc3HitPhi,
                      vector < vector < float > >& setRpc1HitDistEtaRO, //ktaniguc trial add
                      vector < vector < float > >& setRpc1HitDistPhiRO,
                      vector < vector < float > >& setRpc2HitDistEtaRO, //ktaniguc trial add
                      vector < vector < float > >& setRpc2HitDistPhiRO,
                      vector < vector < float > >& setRpc3HitDistEtaRO, //ktaniguc trial add
                      vector < vector < float > >& setRpc3HitDistPhiRO,
                      vector<vector<uint32_t>>& setRpc1HitLayer,
                      vector<vector<uint32_t>>& setRpc2HitLayer,
                      vector<vector<uint32_t>>& setRpc3HitLayer,
                      vector<vector<uint32_t>>& setRpc1HitMeasuresPhi,
                      vector<vector<uint32_t>>& setRpc2HitMeasuresPhi,
                      vector<vector<uint32_t>>& setRpc3HitMeasuresPhi );

    void SetEta_rpcFit( int address, 
                        float* rpcFitInn, 
                        float* rpcFitMid, 
                        float* rpcFitOut, 
                        float& rpcFitInnEta, 
                        float& rpcFitMidEta, 
                        float& rpcFitOutEta);
    float SetR_rpcFit(double address);
    
    float calcRpcFitEta(float rpcFitSlope, 
                        float rpcFitOffset, 
                        float rpcFitR);
  
  
    void setGraphRpcStrip(bool isTag, 
                          TGraphErrors& gr_rpcEtaPhi, 
                          vector<float>& rpcHitEta, 
                          vector<float>& rpcHitPhi, 
                          vector<uint32_t>& rpcHitLayer, 
                          vector<uint32_t>& rpcHitMeasuresPhi);
    void setGraphRpcStripLayer(bool isTag, 
                          TGraphErrors& gr_rpcEtaPhi, 
                          vector<float>& rpcHitEta, 
                          vector<float>& rpcHitPhi, 
                          vector<uint32_t>& rpcHitLayer, 
                          vector<uint32_t>& rpcHitMeasuresPhi);
    bool isLowerLayer(uint32_t rpcHitLayer);
};

#endif
