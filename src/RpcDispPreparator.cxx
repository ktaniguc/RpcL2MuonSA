#define RpcDispPreparator_cxx
#include "../RpcL2MuonSA/RPC.h"
#include "../RpcL2MuonSA/RpcDispPreparator.h"
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

RpcDispPreparator::RpcDispPreparator(){ }

void RpcDispPreparator::SetParameter_rpcFit(//bool isTag, 
                                            int NTrigChain, 
                                            float* rpcFitInn, 
                                            float* rpcFitMid, 
                                            float* rpcFitOut,
                                            float rpcFitMidPhi,
                                            float rpcFitMidSlope,
                                            float rpcFitMidOffset,
                                            float rpcFitOutPhi,
                                            float rpcFitOutSlope,
                                            float rpcFitOutOffset)
{
    rpcFitInn[0] = rpcFitMidPhi;       //Inn -> Mid 
    rpcFitInn[1] = rpcFitMidSlope;
    rpcFitInn[2] = rpcFitMidOffset;
    rpcFitMid[0] = rpcFitMidPhi;
    rpcFitMid[1] = rpcFitMidSlope;
    rpcFitMid[2] = rpcFitMidOffset;
    rpcFitOut[0] = rpcFitOutPhi;
    rpcFitOut[1] = rpcFitOutSlope;
    rpcFitOut[2] = rpcFitOutOffset;
}
void RpcDispPreparator::RpcHitSetter(const Int_t& nRPC, 
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
    vector<vector<uint32_t>>& setRpc3HitMeasuresPhi )

{

  vector<float> rpcHitEta_BML, rpcHitEta_BMS, rpcHitEta_BOL, rpcHitEta_BOS;
  vector<float> rpcHitPhi_BML, rpcHitPhi_BMS, rpcHitPhi_BOL, rpcHitPhi_BOS;
  vector<float> rpcHitEta_BME, rpcHitEta_BOE;//Special sector
  vector<float> rpcHitPhi_BME, rpcHitPhi_BOE;//special sector
  vector<float> rpcHitEta_BMF, rpcHitEta_BOF;//special sector
  vector<float> rpcHitPhi_BMF, rpcHitPhi_BOF;//special sector
  vector<float> rpcHitEta_BOG;//special sector 
  vector<float> rpcHitPhi_BOG;//special sector 
  vector<float> rpcHitR_BML, rpcHitR_BMS, rpcHitR_BOL, rpcHitR_BOS;
  vector<uint32_t> rpcHitLayer_BML, rpcHitLayer_BMS, rpcHitLayer_BOL, rpcHitLayer_BOS;
  vector<uint32_t> rpcHitMeasuresPhi_BML, rpcHitMeasuresPhi_BMS, rpcHitMeasuresPhi_BOL, rpcHitMeasuresPhi_BOS;
  vector<float> rpcHitR_BME, rpcHitR_BOE;
  vector<uint32_t> rpcHitLayer_BME, rpcHitLayer_BOE;
  vector<uint32_t> rpcHitMeasuresPhi_BME, rpcHitMeasuresPhi_BOE;
  vector<float> rpcHitR_BMF, rpcHitR_BOF;
  vector<uint32_t> rpcHitLayer_BMF, rpcHitLayer_BOF;
  vector<uint32_t> rpcHitMeasuresPhi_BMF, rpcHitMeasuresPhi_BOF;
  vector<float> rpcHitR_BOG;
  vector<uint32_t> rpcHitLayer_BOG;
  vector<uint32_t> rpcHitMeasuresPhi_BOG;

  vector<float> rpcHitZ_BML, rpcHitZ_BMS, rpcHitZ_BOL, rpcHitZ_BOS;
  vector<float> rpcHitX_BML, rpcHitX_BMS, rpcHitX_BOL, rpcHitX_BOS;
  vector<float> rpcHitY_BML, rpcHitY_BMS, rpcHitY_BOL, rpcHitY_BOS;
  vector<float> rpcHitZ_BME, rpcHitZ_BOE;//Special sector
  vector<float> rpcHitX_BME, rpcHitX_BOE;//special sector
  vector<float> rpcHitY_BME, rpcHitY_BOE;//Special sector
  vector<float> rpcHitZ_BMF, rpcHitZ_BOF;//special sector
  vector<float> rpcHitX_BMF, rpcHitX_BOF;//special sector
  vector<float> rpcHitY_BMF, rpcHitY_BOF;//special sector
  vector<float> rpcHitZ_BOG;//special sector 
  vector<float> rpcHitX_BOG;//special sector 
  vector<float> rpcHitY_BOG;//special sector 
  vector<float> rpcHitDistEtaRO_BML, rpcHitDistEtaRO_BMS, rpcHitDistEtaRO_BOL, rpcHitDistEtaRO_BOS;
  vector<float> rpcHitDistPhiRO_BML, rpcHitDistPhiRO_BMS, rpcHitDistPhiRO_BOL, rpcHitDistPhiRO_BOS;
  vector<float> rpcHitDistEtaRO_BME, rpcHitDistEtaRO_BOE;//Special sector
  vector<float> rpcHitDistPhiRO_BME, rpcHitDistPhiRO_BOE;//special sector
  vector<float> rpcHitDistEtaRO_BMF, rpcHitDistEtaRO_BOF;//special sector
  vector<float> rpcHitDistPhiRO_BMF, rpcHitDistPhiRO_BOF;//special sector
  vector<float> rpcHitDistEtaRO_BOG;//special sector 
  vector<float> rpcHitDistPhiRO_BOG;//special sector 
  for(int iRPC = 0; iRPC < nRPC; iRPC++)
  {
    if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 1)
    {
      rpcHitEta_BML.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BML.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BML.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BML.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BML.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BML.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BML.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BML.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BML.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BML.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 2)
    {
      rpcHitEta_BMS.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BMS.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BMS.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BMS.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BMS.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BMS.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BMS.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BMS.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BMS.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BMS.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 3)
    {
      rpcHitEta_BMF.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BMF.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BMF.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BMF.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BMF.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BMF.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BMF.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BMF.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BMF.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BMF.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 4)
    {
      rpcHitEta_BME.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BME.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BME.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BME.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BME.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BME.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BME.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BME.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BME.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BME.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 5)
    {
      rpcHitEta_BOL.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BOL.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BOL.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BOL.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BOL.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BOL.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BOL.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BOL.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BOL.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BOL.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 6)
    {
      rpcHitEta_BOS.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BOS.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BOS.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BOS.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BOS.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BOS.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BOS.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BOS.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BOS.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BOS.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 7)
    {
      rpcHitEta_BOF.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BOF.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BOF.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BOF.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BOF.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BOF.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BOF.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BOF.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BOF.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BOF.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 8)
    {
      rpcHitEta_BOE.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BOE.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BOE.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BOE.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BOE.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BOE.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BOE.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BOE.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BOE.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BOE.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }
    else if(mesSA_rpcHitStationNumber->at(NTrigChain)[iRPC] == 9)
    {
      rpcHitEta_BOG.push_back(mesSA_rpcHitEta->at(NTrigChain)[iRPC]);
      rpcHitPhi_BOG.push_back(mesSA_rpcHitPhi->at(NTrigChain)[iRPC]);
      rpcHitR_BOG.push_back(mesSA_rpcHitR->at(NTrigChain)[iRPC]);
      rpcHitLayer_BOG.push_back(mesSA_rpcHitLayer->at(NTrigChain)[iRPC]);
      rpcHitMeasuresPhi_BOG.push_back(mesSA_rpcHitMeasuresPhi->at(NTrigChain)[iRPC]);
      rpcHitX_BOG.push_back(mesSA_rpcHitX->at(NTrigChain)[iRPC]);
      rpcHitY_BOG.push_back(mesSA_rpcHitY->at(NTrigChain)[iRPC]);
      rpcHitZ_BOG.push_back(mesSA_rpcHitZ->at(NTrigChain)[iRPC]);
      rpcHitDistEtaRO_BOG.push_back(mesSA_rpcHitDistEtaRO->at(NTrigChain)[iRPC]);
      rpcHitDistPhiRO_BOG.push_back(mesSA_rpcHitDistPhiRO->at(NTrigChain)[iRPC]);
    }else{
      continue;
    }
  }
 
  if(!rpcHitEta_BML.empty()){
    for(int iRpcBML = 0; iRpcBML < rpcHitEta_BML.size(); iRpcBML++)
    {
      if(6600. < rpcHitR_BML.at(iRpcBML) && rpcHitR_BML.at(iRpcBML) < 7000. )
      {
        setRpc1HitEta[0].push_back(rpcHitEta_BML.at(iRpcBML));
        setRpc1HitPhi[0].push_back(rpcHitPhi_BML.at(iRpcBML));
        setRpc1HitDistEtaRO[0].push_back(rpcHitDistEtaRO_BML.at(iRpcBML));
        setRpc1HitDistPhiRO[0].push_back(rpcHitDistPhiRO_BML.at(iRpcBML));
        setRpc1HitLayer[0].push_back(rpcHitLayer_BML.at(iRpcBML));
        setRpc1HitMeasuresPhi[0].push_back(rpcHitMeasuresPhi_BML.at(iRpcBML));
      }else{
        setRpc2HitEta[0].push_back(rpcHitEta_BML.at(iRpcBML));
        setRpc2HitPhi[0].push_back(rpcHitPhi_BML.at(iRpcBML));
        setRpc2HitDistEtaRO[0].push_back(rpcHitDistEtaRO_BML.at(iRpcBML));
        setRpc2HitDistPhiRO[0].push_back(rpcHitDistPhiRO_BML.at(iRpcBML));
        setRpc2HitLayer[0].push_back(rpcHitLayer_BML.at(iRpcBML));
        setRpc2HitMeasuresPhi[0].push_back(rpcHitMeasuresPhi_BML.at(iRpcBML));
      }
    }
  }
  if(!rpcHitEta_BOL.empty()){
    for(int iRpcBOL = 0; iRpcBOL < rpcHitEta_BOL.size(); iRpcBOL++)
    {
      setRpc3HitEta[0].push_back(rpcHitEta_BOL.at(iRpcBOL));
      setRpc3HitPhi[0].push_back(rpcHitPhi_BOL.at(iRpcBOL));
      setRpc3HitDistEtaRO[0].push_back(rpcHitDistEtaRO_BOL.at(iRpcBOL));
      setRpc3HitDistPhiRO[0].push_back(rpcHitDistPhiRO_BOL.at(iRpcBOL));
      setRpc3HitLayer[0].push_back(rpcHitLayer_BOL.at(iRpcBOL));
      setRpc3HitMeasuresPhi[0].push_back(rpcHitMeasuresPhi_BOL.at(iRpcBOL));
    }
  }
  if(!rpcHitEta_BMS.empty()){
    for(int iRpcBMS = 0; iRpcBMS < rpcHitEta_BMS.size(); iRpcBMS++)
    {
      if(7620. < rpcHitR_BMS.at(iRpcBMS) && rpcHitR_BMS.at(iRpcBMS) < 8020. )
      {
        setRpc1HitEta[1].push_back(rpcHitEta_BMS.at(iRpcBMS));
        setRpc1HitPhi[1].push_back(rpcHitPhi_BMS.at(iRpcBMS));
        setRpc1HitDistEtaRO[1].push_back(rpcHitDistEtaRO_BMS.at(iRpcBMS));
        setRpc1HitDistPhiRO[1].push_back(rpcHitDistPhiRO_BMS.at(iRpcBMS));
        setRpc1HitLayer[1].push_back(rpcHitLayer_BMS.at(iRpcBMS));
        setRpc1HitMeasuresPhi[1].push_back(rpcHitMeasuresPhi_BMS.at(iRpcBMS));
      }else{
        setRpc2HitEta[1].push_back(rpcHitEta_BMS.at(iRpcBMS));
        setRpc2HitPhi[1].push_back(rpcHitPhi_BMS.at(iRpcBMS));
        setRpc2HitDistEtaRO[1].push_back(rpcHitDistEtaRO_BMS.at(iRpcBMS));
        setRpc2HitDistPhiRO[1].push_back(rpcHitDistPhiRO_BMS.at(iRpcBMS));
        setRpc2HitLayer[1].push_back(rpcHitLayer_BMS.at(iRpcBMS));
        setRpc2HitMeasuresPhi[1].push_back(rpcHitMeasuresPhi_BMS.at(iRpcBMS));
      }
    }
  }
  if(!rpcHitEta_BOS.empty()){
    for(int iRpcBOS = 0; iRpcBOS < rpcHitEta_BOS.size(); iRpcBOS++)
    {
      setRpc3HitEta[1].push_back(rpcHitEta_BOS.at(iRpcBOS));
      setRpc3HitPhi[1].push_back(rpcHitPhi_BOS.at(iRpcBOS));
      setRpc3HitDistEtaRO[1].push_back(rpcHitDistEtaRO_BOS.at(iRpcBOS));
      setRpc3HitDistPhiRO[1].push_back(rpcHitDistPhiRO_BOS.at(iRpcBOS));
      setRpc3HitLayer[1].push_back(rpcHitLayer_BOS.at(iRpcBOS));
      setRpc3HitMeasuresPhi[1].push_back(rpcHitMeasuresPhi_BOS.at(iRpcBOS));
    }
  }
  if(!rpcHitEta_BMF.empty()){
    for(int iRpcBMF = 0; iRpcBMF < rpcHitEta_BMF.size(); iRpcBMF++)
    {
      if(7600. < rpcHitR_BMF.at(iRpcBMF) && rpcHitR_BMF.at(iRpcBMF) < 8100. )
      {
        setRpc1HitEta[2].push_back(rpcHitEta_BMF.at(iRpcBMF));
        setRpc1HitPhi[2].push_back(rpcHitPhi_BMF.at(iRpcBMF));
        setRpc1HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BMF.at(iRpcBMF));
        setRpc1HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BMF.at(iRpcBMF));
        setRpc1HitLayer[2].push_back(rpcHitLayer_BMF.at(iRpcBMF));
        setRpc1HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BMF.at(iRpcBMF));
      }else{
        setRpc2HitEta[2].push_back(rpcHitEta_BMF.at(iRpcBMF));
        setRpc2HitPhi[2].push_back(rpcHitPhi_BMF.at(iRpcBMF));
        setRpc2HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BMF.at(iRpcBMF));
        setRpc2HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BMF.at(iRpcBMF));
        setRpc2HitLayer[2].push_back(rpcHitLayer_BMF.at(iRpcBMF));
        setRpc2HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BMF.at(iRpcBMF));
      }
    }
  }
  if(!rpcHitEta_BOF.empty()){
    for(int iRpcBOF = 0; iRpcBOF < rpcHitEta_BOF.size(); iRpcBOF++)
    {
      setRpc3HitEta[2].push_back(rpcHitEta_BOF.at(iRpcBOF));
      setRpc3HitPhi[2].push_back(rpcHitPhi_BOF.at(iRpcBOF));
      setRpc3HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BOF.at(iRpcBOF));
      setRpc3HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BOF.at(iRpcBOF));
      setRpc3HitLayer[2].push_back(rpcHitLayer_BOF.at(iRpcBOF));
      setRpc3HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BOF.at(iRpcBOF));
    }
  }
  if(!rpcHitEta_BME.empty()){
    for(int iRpcBME = 0; iRpcBME < rpcHitEta_BME.size(); iRpcBME++)
    {
      if(7600. < rpcHitR_BME.at(iRpcBME) && rpcHitR_BME.at(iRpcBME) < 7900. )
      {
        setRpc1HitEta[2].push_back(rpcHitEta_BME.at(iRpcBME));
        setRpc1HitPhi[2].push_back(rpcHitPhi_BME.at(iRpcBME));
        setRpc1HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BME.at(iRpcBME));
        setRpc1HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BME.at(iRpcBME));
        setRpc1HitLayer[2].push_back(rpcHitLayer_BME.at(iRpcBME));
        setRpc1HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BME.at(iRpcBME));
      }else{
        setRpc2HitEta[2].push_back(rpcHitEta_BME.at(iRpcBME));
        setRpc2HitPhi[2].push_back(rpcHitPhi_BME.at(iRpcBME));
        setRpc2HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BME.at(iRpcBME));
        setRpc2HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BME.at(iRpcBME));
        setRpc2HitLayer[2].push_back(rpcHitLayer_BME.at(iRpcBME));
        setRpc2HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BME.at(iRpcBME));
      }
    }
  }
  if(!rpcHitEta_BOE.empty()){
    for(int iRpcBOE = 0; iRpcBOE < rpcHitEta_BOE.size(); iRpcBOE++)
    {
      setRpc3HitEta[2].push_back(rpcHitEta_BOE.at(iRpcBOE));
      setRpc3HitPhi[2].push_back(rpcHitPhi_BOE.at(iRpcBOE));
      setRpc3HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BOE.at(iRpcBOE));
      setRpc3HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BOE.at(iRpcBOE));
      setRpc3HitLayer[2].push_back(rpcHitLayer_BOE.at(iRpcBOE));
      setRpc3HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BOE.at(iRpcBOE));
    }
  }
  if(!rpcHitEta_BOG.empty()){
    for(int iRpcBOG = 0; iRpcBOG < rpcHitEta_BOG.size(); iRpcBOG++)
    {
      setRpc3HitEta[2].push_back(rpcHitEta_BOG.at(iRpcBOG));
      setRpc3HitPhi[2].push_back(rpcHitPhi_BOG.at(iRpcBOG));
      setRpc3HitDistEtaRO[2].push_back(rpcHitDistEtaRO_BOG.at(iRpcBOG));
      setRpc3HitDistPhiRO[2].push_back(rpcHitDistPhiRO_BOG.at(iRpcBOG));
      setRpc3HitLayer[2].push_back(rpcHitLayer_BOG.at(iRpcBOG));
      setRpc3HitMeasuresPhi[2].push_back(rpcHitMeasuresPhi_BOG.at(iRpcBOG));
    }
  }





 cout << "--Check Tag and Probe muon's RpcHit Strip--" << endl;
 for(auto& itrrpcHitSector : setRpc2HitEta){
   cout << "rpcHit eta at RPC2 ::";
   for(auto& itrRpcHit : itrrpcHitSector){
     cout << itrRpcHit <<", ";
   }
   cout << "\n";
 }
 for(auto& itrrpcHitSector : setRpc2HitPhi){
   cout << "rpcHit Phi at RPC2 ::";
   for(auto& itrRpcHit : itrrpcHitSector){
     cout << itrRpcHit <<", ";
   }
   cout << "\n";
 }
 for(auto& itrrpcHitLayer : setRpc1HitLayer){
   cout << "rpcHit Layer at RPC1 ::";
   for(auto& itrRpcHitLay : itrrpcHitLayer){
     cout << itrRpcHitLay <<", ";
   }
   cout << "\n";
 }
 for(auto& itrrpcHitLayer : setRpc2HitLayer){
   cout << "rpcHit Layer at RPC2 ::";
   for(auto& itrRpcHitLay : itrrpcHitLayer){
     cout << itrRpcHitLay <<", ";
   }
   cout << "\n";
 }
 for(auto& itrrpcHitLayer : setRpc3HitLayer){
   cout << "rpcHit Layer at RPC3 ::";
   for(auto& itrRpcHitLay : itrrpcHitLayer){
     cout << itrRpcHitLay <<", ";
   }
   cout << "\n";
 }

}

float RpcDispPreparator::SetR_rpcFit(double address){

  float rpcFitMidR = 0;

  if (address == 0){
    //Large Sector
    rpcFitMidR = 7478.;

  }
  else if (address == 1){
    //Small Sector
    rpcFitMidR = 8365.;
  }
  else if (address == 2){
    //Large Special Sector (trial)
    rpcFitMidR = 8100.;

  }
  else if (address == 3){
    //Small Special Sector (trial)
    rpcFitMidR = 8400.;
  }
  if(address >= 0){
    return rpcFitMidR;
  }
  else{
    return 0;
  }
}
  
void RpcDispPreparator::SetEta_rpcFit( int address, float* rpcFitInn, float* rpcFitMid, float* rpcFitOut, float& rpcFitInnEta, float& rpcFitMidEta, float& rpcFitOutEta)
{
  float rpcFitInnR = 0;
  float rpcFitMidR = 0;
  float rpcFitOutR = 0;
  if (address == 0){
    //Large Sector
    rpcFitInnR = 6800.;
    rpcFitMidR = 7478.;
    rpcFitOutR = 9832.;

  }
  else if (address == 1){
    //Small Sector
    rpcFitInnR = 7820.;
    rpcFitMidR = 8365.;
    rpcFitOutR = 10229.;
  }
  else if (address == 2){
    //Large Special Sector (trial)
    rpcFitInnR = 7740.;
    rpcFitMidR = 8100.;
    rpcFitOutR = 11000.;

  }
  else if (address == 3){
    //Small Special Sector (trial)
    rpcFitInnR = 7850.;
    rpcFitMidR = 8400.;
    rpcFitOutR = 10400.;
  }

  rpcFitInnEta = calcRpcFitEta(rpcFitMid[1], rpcFitMid[2], rpcFitInnR);
  rpcFitMidEta = calcRpcFitEta(rpcFitMid[1], rpcFitMid[2], rpcFitMidR);
  rpcFitOutEta = calcRpcFitEta(rpcFitOut[1], rpcFitOut[2], rpcFitOutR);

}


float RpcDispPreparator::calcRpcFitEta(float rpcFitSlope, 
                                       float rpcFitOffset, 
                                       float rpcFitR)
{
  float rpcFitEta = 0;
  float theta = atan(rpcFitSlope*rpcFitR/(rpcFitR - rpcFitOffset));
  float abs_tan = fabs(tan(theta/2));
  int sign_tan = tan(theta/2)/abs_tan;
  rpcFitEta = -sign_tan*log(abs_tan);

  
  return rpcFitEta;
}

void RpcDispPreparator::setGraphRpcStrip(bool isTag, 
                                         TGraphErrors& gr_rpcHitEtaPhi, 
                                         vector<float>& rpcHitEta, 
                                         vector<float>& rpcHitPhi, 
                                         vector<uint32_t>& rpcHitLayer, 
                                         vector<uint32_t>& rpcHitMeasuresPhi)
{
  if(isTag){
    gr_rpcHitEtaPhi = TGraphErrors(rpcHitEta.size());
    for(unsigned int iRpcHit = 0; iRpcHit < rpcHitEta.size(); iRpcHit++)
    {
      if(!rpcHitMeasuresPhi.at(iRpcHit))
      {
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.0003, 0.05);


      }else{
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.05, 0.001);
      }
    }
  }
  else{
    gr_rpcHitEtaPhi = TGraphErrors(rpcHitEta.size());
    for(unsigned int iRpcHit = 0; iRpcHit < rpcHitEta.size(); iRpcHit++)
    {
      if(!rpcHitMeasuresPhi.at(iRpcHit))
      {
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.001, 0.05);
      }else{
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.05, 0.003);
      }
    }
  }
}
void RpcDispPreparator::setGraphRpcStripLayer(bool isTag, 
                                         TGraphErrors& gr_rpcHitEtaPhi, 
                                         vector<float>& rpcHitEta, 
                                         vector<float>& rpcHitPhi, 
                                         vector<uint32_t>& rpcHitLayer, 
                                         vector<uint32_t>& rpcHitMeasuresPhi)
{
  if(isTag){
    gr_rpcHitEtaPhi = TGraphErrors(rpcHitEta.size());
    for(unsigned int iRpcHit = 0; iRpcHit < rpcHitEta.size(); iRpcHit++)
    {
      if(!rpcHitMeasuresPhi.at(iRpcHit))
      {
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.0003, 0.05);


      }else{
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.05, 0.001);
      }
    }
  }
  else{
    gr_rpcHitEtaPhi = TGraphErrors(rpcHitEta.size());
    for(unsigned int iRpcHit = 0; iRpcHit < rpcHitEta.size(); iRpcHit++)
    {
      if(!rpcHitMeasuresPhi.at(iRpcHit))
      {
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.001, 0.05);
      }else{
        gr_rpcHitEtaPhi.SetPoint(iRpcHit, rpcHitEta.at(iRpcHit), rpcHitPhi.at(iRpcHit));
        gr_rpcHitEtaPhi.SetPointError(iRpcHit, 0.05, 0.003);
      }
    }
  }
}

bool RpcDispPreparator::isLowerLayer(uint32_t rpcHitLayer){
  if(rpcHitLayer == 0 || rpcHitLayer == 2 || rpcHitLayer == 4 || rpcHitLayer == 6){
    return true;
  }
  if(rpcHitLayer == 1 || rpcHitLayer == 3 || rpcHitLayer == 5 || rpcHitLayer == 7){
    return false;
  }
}





