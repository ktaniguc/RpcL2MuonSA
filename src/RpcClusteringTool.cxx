#define RpcClusteringTool_cxx
#include "../RpcL2MuonSA/RpcClusteringTool.h"
#include <iostream>
#include "TTree.h"
//#include "../RpcL2MuonSA/RPC.h"

using namespace std;

RpcClusteringTool::RpcClusteringTool(){
  //  cout << fChain->GetEntries();
  //  cout << "\nhoge";
}

void RpcClusteringTool::coutSorted(){
  cout << "hoge" << endl;
}


//template <class type> void RpcClusteringRool::sortStrip(vector<type> )

void RpcClusteringTool::sortStrip(TGraphErrors *gr_rpcHitEtaPhi_BL,
    TGraphErrors *gr_rpcHitEtaPhi_BS,
    TGraphErrors *gr_rpcHitEtaPhi_BSP,
    vector<double>& mesSA_rpcHitStationNumber,
    vector<uint32_t>& mesSA_rpcHitMeasuresPhi,
    vector<double>& mesSA_rpcHitEta,
    vector<double>& mesSA_rpcHitPhi,
    vector<uint32_t>& mesSA_rpcHitLayer
    )
{
  vector < vector < float > >    rpc1HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc1HitPhi(9, vector < float > (0));
  vector < vector < float > >    rpc2HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc2HitPhi(9, vector < float > (0));
  vector < vector < float > >    rpc3HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc3HitPhi(9, vector < float > (0));
  vector < vector < uint32_t > > rpcHitLayer(9, vector < uint32_t > (0));
  int nStationNum = mesSA_rpcHitStationNumber.size();

  cout << "RpcClus::Start" << endl;;
  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum) -1;
    if(!mesSA_rpcHitMeasuresPhi.at(iStationNum)){
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
      else{

        rpc3HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
    }
    else{
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
      else{
        rpc3HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
    }
  }

  cout << "RpcClus::beforeinfo" << endl;
  //info for gr_rpcHitEtaPhiLayer

  vector < vector < float > >    rpcHitEta_BL(8, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BL(8, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BL(8, vector < uint32_t > (0));
  vector < vector < float > >    rpcHitEta_BS(8, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BS(8, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BS(8, vector < uint32_t > (0));
  vector < vector < float > >    rpcHitEta_BSP(8, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BSP(8, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BSP(8, vector < uint32_t > (0));

  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum);
    if(StationNum == 1 || StationNum == 5){
      rpcHitEta_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
    }
    else if(StationNum == 2 || StationNum == 6){
      rpcHitEta_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
    }
    else{
      rpcHitEta_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
    }


  }
  cout << "RpcClus::beforelay";
  const int nLayer = 8;
  for(int iLayer = 0; iLayer < nLayer; iLayer++){
    if(!rpcHitEta_BL[iLayer].empty()){
      int nBL = (int)rpcHitEta_BL[iLayer].size();
      gr_rpcHitEtaPhi_BL[iLayer] = TGraphErrors(nBL);
      for(int iBL = 0; iBL < nBL; iBL++){
        if(!rpcHitMeasuresPhi_BL[iLayer].at(iBL)){
          gr_rpcHitEtaPhi_BL[iLayer].SetPoint(iBL, rpcHitEta_BL[iLayer].at(iBL), rpcHitPhi_BL[iLayer].at(iBL));
          gr_rpcHitEtaPhi_BL[iLayer].SetPointError(iBL, 0.001, 0.05);
        }
        else{
          gr_rpcHitEtaPhi_BL[iLayer].SetPoint(iBL, rpcHitEta_BL[iLayer].at(iBL), rpcHitPhi_BL[iLayer].at(iBL));
          gr_rpcHitEtaPhi_BL[iLayer].SetPointError(iBL, 0.05, 0.001);
        }
      }
    }

    cout << "BL end" << endl;
    if(!rpcHitEta_BS[iLayer].empty()){
      int nBS = (int)rpcHitEta_BS[iLayer].size();
      gr_rpcHitEtaPhi_BS[iLayer] = TGraphErrors(nBS);
      for(int iBS = 0; iBS < nBS; iBS++){
        if(!rpcHitMeasuresPhi_BS[iLayer].at(iBS)){
          gr_rpcHitEtaPhi_BS[iLayer].SetPoint(iBS, rpcHitEta_BS[iLayer].at(iBS), rpcHitPhi_BS[iLayer].at(iBS));
          gr_rpcHitEtaPhi_BS[iLayer].SetPointError(iBS, 0.001, 0.05);
        }
        else{
          gr_rpcHitEtaPhi_BS[iLayer].SetPoint(iBS, rpcHitEta_BS[iLayer].at(iBS), rpcHitPhi_BS[iLayer].at(iBS));
          gr_rpcHitEtaPhi_BS[iLayer].SetPointError(iBS, 0.05, 0.001);
        }
      }
    }

    cout << "BS end" << endl;
    
    if(!rpcHitEta_BSP[iLayer].empty()){
      int nBSP = (int)rpcHitEta_BSP[iLayer].size();
      gr_rpcHitEtaPhi_BSP[iLayer] = TGraphErrors(nBSP);
      for(int iBSP = 0; iBSP < nBSP; iBSP++){
        if(!rpcHitMeasuresPhi_BSP[iLayer].at(iBSP)){
          gr_rpcHitEtaPhi_BSP[iLayer].SetPoint(iBSP, rpcHitEta_BSP[iLayer].at(iBSP), rpcHitPhi_BSP[iLayer].at(iBSP));
          gr_rpcHitEtaPhi_BSP[iLayer].SetPointError(iBSP, 0.001, 0.05);
        }
        else{
          gr_rpcHitEtaPhi_BSP[iLayer].SetPoint(iBSP, rpcHitEta_BSP[iLayer].at(iBSP), rpcHitPhi_BSP[iLayer].at(iBSP));
          gr_rpcHitEtaPhi_BSP[iLayer].SetPointError(iBSP, 0.05, 0.001);
        }
      }
    }
    cout << "BSP end" << endl;

  }

  cout << "RpcClus::afterlay";
}

void RpcClusteringTool::sortStrip(TGraphErrors *gr_rpcHitEtaPhi_BL,
                                  TGraphErrors *gr_rpcHitEtaPhi_BS,
                                  TGraphErrors *gr_rpcHitEtaPhi_BSP,
                                  vector<double>& mesSA_rpcHitStationNumber,
                                  vector<uint32_t>& mesSA_rpcHitMeasuresPhi,
                                  vector<float>& mesSA_rpcHitEtaRO,
                                  vector<float>& mesSA_rpcHitPhiRO,
                                  vector<double>& mesSA_rpcHitEta,
                                  vector<double>& mesSA_rpcHitPhi,
                                  vector<uint32_t>& mesSA_rpcHitLayer
                                  )
{
  vector < vector < float > >    rpc1HitEtaRO(9, vector < float > (0));
  vector < vector < float > >    rpc1HitPhiRO(9, vector < float > (0));
  vector < vector < float > >    rpc2HitEtaRO(9, vector < float > (0));
  vector < vector < float > >    rpc2HitPhiRO(9, vector < float > (0));
  vector < vector < float > >    rpc3HitEtaRO(9, vector < float > (0));
  vector < vector < float > >    rpc3HitPhiRO(9, vector < float > (0));
  vector < vector < float > >    rpc1HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc1HitPhi(9, vector < float > (0));
  vector < vector < float > >    rpc2HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc2HitPhi(9, vector < float > (0));
  vector < vector < float > >    rpc3HitEta(9, vector < float > (0));
  vector < vector < float > >    rpc3HitPhi(9, vector < float > (0));
  vector < vector < uint32_t > > rpcHitLayer(9, vector < uint32_t > (0));
  int nStationNum = mesSA_rpcHitStationNumber.size();

  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum) -1;
    if(mesSA_rpcHitMeasuresPhi.at(iStationNum)){
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
      else{

        rpc3HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
    }else{
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
      else{

        rpc3HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
    }
  }
  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum) -1;
    if(!mesSA_rpcHitMeasuresPhi.at(iStationNum)){
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
      else{

        rpc3HitEta[StationNum].push_back(mesSA_rpcHitEta.at(iStationNum));
      }
    }else{
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
      else{

        rpc3HitPhi[StationNum].push_back(mesSA_rpcHitPhi.at(iStationNum));
      }
    }
  }

  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum) -1;
    if(mesSA_rpcHitMeasuresPhi.at(iStationNum)){
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
      else{

        rpc3HitEtaRO[StationNum].push_back(mesSA_rpcHitEtaRO.at(iStationNum));
      }
    }else{
      if(mesSA_rpcHitLayer.at(iStationNum) == 0 || mesSA_rpcHitLayer.at(iStationNum) == 1){
        rpc1HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
      else if(mesSA_rpcHitLayer.at(iStationNum) == 2 || mesSA_rpcHitLayer.at(iStationNum) == 3){

        rpc2HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
      else{

        rpc3HitPhiRO[StationNum].push_back(mesSA_rpcHitPhiRO.at(iStationNum));
      }
    }
  }

//info for gr_rpcHitEtaPhiLayer

  vector < vector < float > >    rpcHitEta_BL(6, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BL(6, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BL(6, vector < uint32_t > (0));
  vector < vector < float > >    rpcHitEta_BS(6, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BS(6, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BS(6, vector < uint32_t > (0));
  vector < vector < float > >    rpcHitEta_BSP(6, vector < float > (0));
  vector < vector < float > >    rpcHitPhi_BSP(6, vector < float > (0));
  vector < vector < uint32_t > >    rpcHitMeasuresPhi_BSP(6, vector < uint32_t > (0));

  for(int iStationNum = 0; iStationNum < nStationNum; iStationNum++){
    if(mesSA_rpcHitStationNumber.at(iStationNum) < 1 || mesSA_rpcHitStationNumber.at(iStationNum) > 9) continue;
    int StationNum = (int)mesSA_rpcHitStationNumber.at(iStationNum);
    if(StationNum == 1 || StationNum == 5){
      rpcHitEta_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BL[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
      }
    else if(StationNum == 2 || StationNum == 6){
      rpcHitEta_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BS[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
    }
    else{
      rpcHitEta_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitEta.at(iStationNum));
      rpcHitPhi_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitPhi.at(iStationNum));
      rpcHitMeasuresPhi_BSP[mesSA_rpcHitLayer.at(iStationNum)].push_back(mesSA_rpcHitMeasuresPhi.at(iStationNum));
      }


  }
  const int nLayer = 6;
  for(int iLayer = 0; iLayer < nLayer; iLayer++){
    cout << "hogeLayer::" << iLayer << endl;
    if(!rpcHitEta_BL[iLayer].empty()){
      int nBL = (int)rpcHitEta_BL[iLayer].size();
      gr_rpcHitEtaPhi_BL[iLayer] = TGraphErrors(nBL);
      for(int iBL = 0; iBL < nBL; iBL++){
        if(!rpcHitMeasuresPhi_BL[iLayer].at(iBL)){
          gr_rpcHitEtaPhi_BL[iLayer].SetPoint(iBL, rpcHitEta_BL[iLayer].at(iBL), rpcHitPhi_BL[iLayer].at(iBL));
          gr_rpcHitEtaPhi_BL[iLayer].SetPointError(iBL, 0.0006, 0.05);
          cout << "hoge::" << rpcHitEta_BL[iLayer].at(iBL) << endl;
        }
        else{
          gr_rpcHitEtaPhi_BL[iLayer].SetPoint(iBL, rpcHitEta_BL[iLayer].at(iBL), rpcHitPhi_BL[iLayer].at(iBL));
          gr_rpcHitEtaPhi_BL[iLayer].SetPointError(iBL, 0.05, 0.0001);
        }
      }
    }
    if(!rpcHitEta_BS[iLayer].empty()){
      int nBS = (int)rpcHitEta_BS[iLayer].size();
      gr_rpcHitEtaPhi_BS[iLayer] = TGraphErrors(nBS);
      for(int iBS = 0; iBS < nBS; iBS++){
        if(!rpcHitMeasuresPhi_BS[iLayer].at(iBS)){
          gr_rpcHitEtaPhi_BS[iLayer].SetPoint(iBS, rpcHitEta_BS[iLayer].at(iBS), rpcHitPhi_BS[iLayer].at(iBS));
          gr_rpcHitEtaPhi_BS[iLayer].SetPointError(iBS, 0.0006, 0.05);
        }
        else{
          gr_rpcHitEtaPhi_BS[iLayer].SetPoint(iBS, rpcHitEta_BS[iLayer].at(iBS), rpcHitPhi_BS[iLayer].at(iBS));
          gr_rpcHitEtaPhi_BS[iLayer].SetPointError(iBS, 0.05, 0.0001);
        }
      }
    }
    if(!rpcHitEta_BSP[iLayer].empty()){
      int nBSP = (int)rpcHitEta_BSP[iLayer].size();
      gr_rpcHitEtaPhi_BSP[iLayer] = TGraphErrors(nBSP);
      for(int iBSP = 0; iBSP < nBSP; iBSP++){
        if(!rpcHitMeasuresPhi_BSP[iLayer].at(iBSP)){
          gr_rpcHitEtaPhi_BSP[iLayer].SetPoint(iBSP, rpcHitEta_BSP[iLayer].at(iBSP), rpcHitPhi_BSP[iLayer].at(iBSP));
          gr_rpcHitEtaPhi_BSP[iLayer].SetPointError(iBSP, 0.0006, 0.05);
        }
        else{
          gr_rpcHitEtaPhi_BSP[iLayer].SetPoint(iBSP, rpcHitEta_BSP[iLayer].at(iBSP), rpcHitPhi_BSP[iLayer].at(iBSP));
          gr_rpcHitEtaPhi_BSP[iLayer].SetPointError(iBSP, 0.05, 0.0001);
        }
      }
    }

  }


 /* for(int i = 0; i < 9; i++){
    cout << "rpc1HitDistTo...RO No." << i+1 << endl;
    if(!rpc1HitEtaRO[i].empty()){
      sort(rpc1HitEtaRO[i].begin(), rpc1HitEtaRO[i].end());
      cout << "EtaRO::";
      for(int t = 0; t < rpc1HitEtaRO[i].size(); t++){
        cout <<  rpc1HitEtaRO[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc1HitPhiRO[i].empty()){
      sort(rpc1HitPhiRO[i].begin(), rpc1HitPhiRO[i].end());
      cout << "PhiRO::";
      for(int t = 0; t < rpc1HitPhiRO[i].size(); t++){
        cout << rpc1HitPhiRO[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  for(int i = 0; i < 9; i++){
    cout << "rpc2HitDistTo...RO No." << i+1 << endl;
    if(!rpc2HitEtaRO[i].empty()){
      sort(rpc2HitEtaRO[i].begin(), rpc2HitEtaRO[i].end());
      cout << "EtaRO::";
      for(int t = 0; t < rpc2HitEtaRO[i].size(); t++){
        cout <<  rpc2HitEtaRO[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc2HitPhiRO[i].empty()){
      sort(rpc2HitPhiRO[i].begin(), rpc2HitPhiRO[i].end());
      cout << "PhiRO::";
      for(int t = 0; t < rpc2HitPhiRO[i].size(); t++){
        cout << rpc2HitPhiRO[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  for(int i = 0; i < 9; i++){
    cout << "rpc3HitDistTo...RO No." << i+1 << endl;
    if(!rpc3HitEtaRO[i].empty()){
      sort(rpc3HitEtaRO[i].begin(), rpc3HitEtaRO[i].end());
      cout << "EtaRO::";
      for(int t = 0; t < rpc3HitEtaRO[i].size(); t++){
        cout <<  rpc3HitEtaRO[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc3HitPhiRO[i].empty()){
      sort(rpc3HitPhiRO[i].begin(), rpc3HitPhiRO[i].end());
      cout << "PhiRO::";
      for(int t = 0; t < rpc3HitPhiRO[i].size(); t++){
        cout << rpc3HitPhiRO[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  
  
  
  for(int i = 0; i < 9; i++){
    cout << "rpc1Hit. No." << i+1 << endl;
    if(!rpc1HitEta[i].empty()){
      sort(rpc1HitEta[i].begin(), rpc1HitEta[i].end());
      cout << "Eta::";
      for(int t = 0; t < rpc1HitEta[i].size(); t++){
        cout <<  rpc1HitEta[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc1HitPhi[i].empty()){
      sort(rpc1HitPhi[i].begin(), rpc1HitPhi[i].end());
      cout << "Phi::";
      for(int t = 0; t < rpc1HitPhi[i].size(); t++){
        cout << rpc1HitPhi[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  for(int i = 0; i < 9; i++){
    cout << "rpc2Hit. No." << i+1 << endl;
    if(!rpc2HitEta[i].empty()){
      sort(rpc2HitEta[i].begin(), rpc2HitEta[i].end());
      cout << "Eta::";
      for(int t = 0; t < rpc2HitEta[i].size(); t++){
        cout <<  rpc2HitEta[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc2HitPhi[i].empty()){
      sort(rpc2HitPhi[i].begin(), rpc2HitPhi[i].end());
      cout << "Phi::";
      for(int t = 0; t < rpc2HitPhi[i].size(); t++){
        cout << rpc2HitPhi[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  for(int i = 0; i < 9; i++){
    cout << "rpc3Hit. No." << i+1 << endl;
    if(!rpc3HitEta[i].empty()){
      sort(rpc3HitEta[i].begin(), rpc3HitEta[i].end());
      cout << "Eta::";
      for(int t = 0; t < rpc3HitEta[i].size(); t++){
        cout <<  rpc3HitEta[i][t] << ", ";
      }
    }
    cout << "\n";
    if(!rpc3HitPhi[i].empty()){
      sort(rpc3HitPhi[i].begin(), rpc3HitPhi[i].end());
      cout << "Phi::";
      for(int t = 0; t < rpc3HitPhi[i].size(); t++){
        cout << rpc3HitPhi[i][t] << ", ";
      }
    }
    cout << "\n";
  }
  calcAdjacentWidth(rpc1HitEta, rpc1HitPhi);
  calcAdjacentWidth(rpc2HitEta, rpc2HitPhi);
  calcAdjacentWidth(rpc3HitEta, rpc3HitPhi);
  calcAdjacentWidth(rpc1HitEtaRO, rpc1HitPhiRO);
  calcAdjacentWidth(rpc2HitEtaRO, rpc2HitPhiRO);
  calcAdjacentWidth(rpc3HitEtaRO, rpc3HitPhiRO);*/
}

void RpcClusteringTool::calcAdjacentWidth(vector < vector < float> >& rpcHitEta, vector < vector < float> >& rpcHitPhi )
{
  const int nStat = 9;
  for(int iStat = 0; iStat < nStat; iStat++){

    cout << "RPC Hit SectorNo..." << iStat << endl;
    if(!rpcHitEta.empty()){
      cout << "rpcHitEta::";
      int nRpcHitEta = rpcHitEta[iStat].size();
      for(int iRpcHitEta = 0; iRpcHitEta < nRpcHitEta -1; iRpcHitEta++){
        float Width = rpcHitEta[iStat].at(iRpcHitEta +1) - rpcHitEta[iStat].at(iRpcHitEta);
        cout << Width << ", ";
      }
    }
    cout << endl;
    if(!rpcHitPhi.empty()){
      cout << "rpcHitPhi::";
      int nRpcHitPhi = rpcHitPhi[iStat].size();
      for(int iRpcHitPhi = 0; iRpcHitPhi < nRpcHitPhi -1; iRpcHitPhi++){
        float Width = rpcHitPhi[iStat].at(iRpcHitPhi +1) - rpcHitPhi[iStat].at(iRpcHitPhi);
        cout << Width << ", ";
      }
    }
    cout << endl;
  }
}









