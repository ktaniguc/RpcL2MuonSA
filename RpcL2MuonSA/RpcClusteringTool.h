#ifndef RpcClusteringTool_h
#define RpcClusteringTool_h
#include "RPC.h"
#include "TChain.h"

class RpcClusteringTool {
  public :
    RpcClusteringTool();
    void coutSorted();
    void sortStrip(TGraphErrors *gr_rpcHitEtaPhi_BL,
                   TGraphErrors *gr_rpcHitEtaPhi_BS,
                   TGraphErrors *gr_rpcHitEtaPhi_BSP,
                   vector<double>& mesSA_rpcHitStationNumber,
                   vector<uint32_t>& mesSA_rpcHitMeasuresPhi,
                   vector<double>& mesSA_rpcHitEta,
                   vector<double>& mesSA_rpcHitPhi,
                   vector<uint32_t>& mesSA_rpcHitLayer
                   );
    void sortStrip(TGraphErrors *gr_rpcHitEtaPhi_BL,
                   TGraphErrors *gr_rpcHitEtaPhi_BS,
                   TGraphErrors *gr_rpcHitEtaPhi_BSP,
                   vector<double>& mesSA_rpcHitStationNumber,
                   vector<uint32_t>& mesSA_rpcHitMeasuresPhi,
                   vector<float>& mesSA_rpcHitEtaRO,
                   vector<float>& mesSA_rpcHitPhiRO,
                   vector<double>& mesSA_rpcHitEta,
                   vector<double>& mesSA_rpcHitPhi,
                   vector<uint32_t>& mesSA_rpcHitLayer
                   );
void calcAdjacentWidth(vector < vector < float> >& rpcHitEta, vector < vector < float> >& rpcHitPhi );
 
};
#endif
