//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 25 13:30:35 2018 by ROOT version 6.12/04
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: /home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root
//////////////////////////////////////////////////////////

#ifndef RPC_FCBM_h
#define RPC_FCBM_h

#include "HistData.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <iostream>
#include <map>
#include "../AtlasStyle/AtlasLabels.h"
#include "../AtlasStyle/AtlasUtils.h"

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class RPC_FCBM {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // fixed size dimensions of array or collections stored in the TTree if any.
    long m_total{0};
    long m_total_closeby{0};
    long m_mpsatail{0};
    long m_satail{0};
    long m_mpsatail_closeby{0};
    const int n4 = 0;
    const int n50 = 0;
    const int ntagproc = 3;
    const int nRpcHitLayer = 7;
    string RPC1 = "RPC1";
    string RPC2 = "RPC2";
    string RPC3 = "RPC3";
    
    
    vector<double> SARPCFitInnR; 
    vector<double> SARPCFitInnEta;
    vector<double> SARPCFitMidR; 
    vector<double> SARPCFitMidEta;
    vector<double> SARPCFitOutR; 
    vector<double> SARPCFitOutEta;
    // declaration of leaf types
    //branch variable
    ULong64_t       EventNumber;
    Int_t           RunNumber;
    Int_t           LumiBlock;
    Float_t         AverageInteractionsPerCrossing;
    string          *mes_name;
    UInt_t          nMuon;
    vector<double>  *OfflinePt;
    vector<double>  *OfflineEta;
    vector<double>  *OfflineExtEta;
    vector<double>  *OfflineExtInnEta;
    vector<double>  *OfflinePhi;
    vector<double>  *OfflineExtPhi;
    vector<double>  *OfflineExtInnPhi;
    vector<double>  *OfflineD0;
    vector<double>  *OfflineZ0;
    vector<double>  *OfflineCharge;
    vector<int>     *OfflineNumSegment;
    vector<vector<double> > *OfflineSegmentX;
    vector<vector<double> > *OfflineSegmentY;
    vector<vector<double> > *OfflineSegmentZ;
    vector<vector<double> > *OfflineSegmentPx;
    vector<vector<double> > *OfflineSegmentPy;
    vector<vector<double> > *OfflineSegmentPz;
    vector<vector<double> > *OfflineSegmentChiSquared;
    vector<vector<double> > *OfflineSegmentNumberDoF;
    vector<vector<double> > *OfflineSegmentSector;
    vector<vector<double> > *OfflineSegmentChamberIndex;
    vector<vector<double> > *OfflineSegmentEtaIndex;
    vector<vector<double> > *OfflineSegmentNPrecisionHits;
    vector<vector<double> > *OfflineSegmentNPhiLayers;
    vector<vector<double> > *OfflineSegmentNTrigEtaLayers;
    vector<int>     *EFTAGPass;
    vector<double>  *EFTAGPt;
    vector<double>  *EFTAGEta;
    vector<double>  *EFTAGPhi;
    vector<int>     *L1nRoI;
    vector<bool>    *L1isMoreCandInRoI;
    vector<int>     *L1Pass;
    vector<bool>     *passedChain;
    vector<double>  *L1Pt;
    vector<double>  *L1Eta;
    vector<double>  *L1Phi;
    vector<int>     *L1RoINumber;
    vector<int>     *L1RoISector;
    vector<int>     *SAHypoPass;
    vector<bool>    *SAOvRmPass;
    vector<bool>    *SAOvRmPass_forClusAlg;
    vector<double>  *SAPt;
    vector<double>  *SAEta;
    vector<double>  *SAPhi;
    vector<double>  *SAEtaMS;
    vector<double>  *SAEtaBE;
    vector<double>  *SAPhiMS;
    vector<double>  *SAPhiBE;
    vector<double>  *SATGCPt;
    vector<double>  *SAPtBarrelRadius;
    vector<double>  *SAPtBarrelSagitta;
    vector<double>  *SAPtEndcapAlpha;
    vector<double>  *SAPtEndcapBeta;
    vector<double>  *SAPtEndcapRadius;
    vector<double>  *SACSCPt;
    vector<double>  *SAsAddress;
    vector<double>  *SASPR_BI;
    vector<double>  *SASPR_BM;
    vector<double>  *SASPR_BO;
    vector<double>  *SASPR_EI;
    vector<double>  *SASPR_EM;
    vector<double>  *SASPR_EO;
    vector<double>  *SASPR_EE;
    vector<double>  *SASPR_CSC;
    vector<double>  *SASPR_BEE;
    vector<double>  *SASPR_BME;
    vector<double>  *SASPZ_BI;
    vector<double>  *SASPZ_BM;
    vector<double>  *SASPZ_BO;
    vector<double>  *SASPZ_EI;
    vector<double>  *SASPZ_EM;
    vector<double>  *SASPZ_EO;
    vector<double>  *SASPZ_EE;
    vector<double>  *SASPZ_CSC;
    vector<double>  *SASPZ_BEE;
    vector<double>  *SASPZ_BME;
    vector<double>  *SASPSlope_BI;
    vector<double>  *SASPSlope_BM;
    vector<double>  *SASPSlope_BO;
    vector<double>  *SASPSlope_EI;
    vector<double>  *SASPSlope_EM;
    vector<double>  *SASPSlope_EO;
    vector<double>  *SASPSlope_EE;
    vector<double>  *SASPSlope_CSC;
    vector<double>  *SASPSlope_BEE;
    vector<double>  *SASPSlope_BME;
    vector<double>  *SASPIntercept_BI;
    vector<double>  *SASPIntercept_BM;
    vector<double>  *SASPIntercept_BO;
    vector<double>  *SASPIntercept_EI;
    vector<double>  *SASPIntercept_EM;
    vector<double>  *SASPIntercept_EO;
    vector<double>  *SASPIntercept_EE;
    vector<double>  *SASPIntercept_CSC;
    vector<double>  *SASPIntercept_BEE;
    vector<double>  *SASPIntercept_BME;
    vector<double>  *SASPChi2_BI;
    vector<double>  *SASPChi2_BM;
    vector<double>  *SASPChi2_BO;
    vector<double>  *SASPChi2_EI;
    vector<double>  *SASPChi2_EM;
    vector<double>  *SASPChi2_EO;
    vector<double>  *SASPChi2_EE;
    vector<double>  *SASPChi2_CSC;
    vector<double>  *SASPChi2_BEE;
    vector<double>  *SASPChi2_BME;
    vector<float>   *SARoIEta;
    vector<float>   *SARoIPhi;
    vector<double>  *SAisRPCFailure;
    vector<double>  *SAisTGCFailure;
    vector<double>  *SABarrelRadius;
    vector<double>  *SABarrelSagitta;
    vector<double>  *SAEtaMap;
    vector<double>  *SAPhiMap;
    vector<unsigned int> *SARoINumber;
    vector<unsigned int> *SARoISector;
    vector<vector<float> > *SARPCHitTime;
    vector<vector<float> > *SARPCHitX;
    vector<vector<float> > *SARPCHitY;
    vector<vector<float> > *SARPCHitZ;
    vector<vector<double> > *SARPCHitR;
    vector<vector<double> > *SARPCHitEta;
    vector<vector<double> > *SARPCHitPhi;
    vector<vector<unsigned int> > *SARPCHitMeasuresPhi;
    vector<vector<unsigned int> > *SARPCHitLayer;
    vector<vector<string> > *SARPCHitStationName;
    vector<vector<double> > *SARPCHitStationNumber;
    
    //RPC Clustering
    vector<vector<vector<int>> > *SAspcSetID; //I'll make an effort not to use 3D vector when I use this value because I hate 3D vector
    vector<vector<double>>  *SAptclus;
    vector<vector<float>>  *SAetaMSclus;
    vector<vector<float>>  *SAchargeclus;
    vector<vector<float>>  *SASPCZ_BI;
    vector<vector<float>>  *SASPCZ_BM;
    vector<vector<float>>  *SASPCZ_BO;
    vector<vector<float>>  *SASPCR_BI;
    vector<vector<float>>  *SASPCR_BM;
    vector<vector<float>>  *SASPCR_BO;
    vector<vector<vector<float>> > *SARPCCluster_zMax;
    vector<vector<vector<float>> > *SARPCCluster_zMin;
    vector<vector<vector<float>> > *SARPCCluster_rMax;
    vector<vector<vector<float>> > *SARPCCluster_rMin;
    vector<vector<float> > *SARPCCluster_gX;
    vector<vector<float> > *SARPCCluster_gY;
    vector<vector<float> > *SARPCCluster_gZ;
    vector<vector<float> > *SARPCCluster_stripWidth;
    vector<vector<float> > *SARPCCluster_loX;
    vector<vector<float> > *SARPCCluster_loY;
    vector<vector<float> > *SARPCCluster_loZ;
    vector<vector<int> > *SARPCCluster_clusterSize;
    vector<vector<int> > *SARPCCluster_clusterLayer;
    vector<vector<int> > *SARPCCluster_clusterMeasPhi;

    vector<vector<float> > *SARPCCluster_fitInnPhi;
    vector<vector<float> > *SARPCCluster_fitMidPhi;
    vector<vector<float> > *SARPCCluster_fitOutPhi;
    vector<vector<float> > *SARPCCluster_fitInnSlope;
    vector<vector<float> > *SARPCCluster_fitMidSlope;
    vector<vector<float> > *SARPCCluster_fitOutSlope;
    vector<vector<float> > *SARPCCluster_fitInnOffset;
    vector<vector<float> > *SARPCCluster_fitMidOffset;
    vector<vector<float> > *SARPCCluster_fitOutOffset;
    vector<vector<bool> > *SARPCCluster_isSuccess;
    vector<vector<bool> > *SARPCCluster_isUsingMidCluster;
    vector<vector<bool> > *SARPCCluster_isPlausibleFitInnMid;
    vector<vector<bool> > *SARPCCluster_isPlausibleFitOut;
    vector<vector<bool> > *SARPCCluster_isPlausiblePhiInnMid;
    vector<vector<bool> > *SARPCCluster_isPlausiblePhiOut;
    vector<vector<vector<int>> > *SARPCCluster_n_foundClusters; //I'll make an effort not to use 3D vector when I use this value because I hate 3D vector
    vector<vector<vector<int>> > *SARPCCluster_id_clustersInSets; //I'll make an effort not to use 3D vector when I use this value because I hate 3D vector
    vector<vector<vector<int>> > *SARPCCluster_id_clustersInPhiSets; //I'll make an effort not to use 3D vector when I use this value because I hate 3D vector

    vector<float>   *SARPCFitInnPhi;
    vector<float>   *SARPCFitInnSlope;
    vector<float>   *SARPCFitInnOffset;
    vector<float>   *SARPCFitMidPhi;
    vector<float>   *SARPCFitMidSlope;
    vector<float>   *SARPCFitMidOffset;
    vector<float>   *SARPCFitOutPhi;
    vector<float>   *SARPCFitOutSlope;
    vector<float>   *SARPCFitOutOffset;
    vector<vector<float> > *SARoadAw;
    vector<vector<float> > *SARoadBw;
    vector<vector<float> > *SAzMin;
    vector<vector<float> > *SAzMax;
    vector<vector<float> > *SArMin;
    vector<vector<float> > *SArMax;
    vector<vector<float> > *SAEtaMin;
    vector<vector<float> > *SAEtaMax;
    vector<vector<int> > *SAMDTHitisOutlier;
    vector<vector<int> > *SAMDTHitChamber;
    vector<vector<float> > *SAMDTHitR;
    vector<vector<float> > *SAMDTHitZ;
    vector<vector<float> > *SAMDTHitPhi;
    vector<vector<float> > *SAMDTHitResidual;
    vector<vector<float> > *SAMDTHitSpace;
    vector<vector<float> > *SAMDTHitSigma;
    vector<vector<float> > *SAMDTHitAllR;
    vector<vector<float> > *SAMDTHitAllZ;
    vector<vector<int> > *SAMDTHitAllisOutlier;
    vector<vector<int> > *SAMDTHitAllclusRoadID;
    vector<int>     *CBHypoPass;
    vector<bool>    *CBOvRmPass;
    vector<double>  *CBPt;
    vector<double>  *CBEta;
    vector<double>  *CBPhi;
    vector<int>     *EFPass;
    vector<double>  *EFPt;
    vector<double>  *EFEta;
    vector<double>  *EFPhi;
    // List of branches
    //branch variable
    TBranch *b_EventNumber;                    
    TBranch *b_RunNumber;                      
    TBranch *b_LumiBlock;                      
    TBranch *b_AverageInteractionsPerCrossing; 
    TBranch *b_mes_name;                       
    TBranch *b_nMuon;                          
    TBranch *b_OfflinePt;                      
    TBranch *b_OfflineEta;                     
    TBranch *b_OfflineExtEta;                  
    TBranch *b_OfflineExtInnEta;               
    TBranch *b_OfflinePhi;                     
    TBranch *b_OfflineExtPhi;                  
    TBranch *b_OfflineExtInnPhi;               
    TBranch *b_OfflineD0;                      
    TBranch *b_OfflineZ0;                      
    TBranch *b_OfflineCharge;                  
    TBranch *b_OfflineNumSegment;              
    TBranch *b_OfflineSegmentX;                
    TBranch *b_OfflineSegmentY;                
    TBranch *b_OfflineSegmentZ;                
    TBranch *b_OfflineSegmentPx;               
    TBranch *b_OfflineSegmentPy;               
    TBranch *b_OfflineSegmentPz;               
    TBranch *b_OfflineSegmentChiSquared;       
    TBranch *b_OfflineSegmentNumberDoF;        
    TBranch *b_OfflineSegmentSector;           
    TBranch *b_OfflineSegmentChamberIndex;     
    TBranch *b_OfflineSegmentEtaIndex;         
    TBranch *b_OfflineSegmentNPrecisionHits;   
    TBranch *b_OfflineSegmentNPhiLayers;       
    TBranch *b_OfflineSegmentNTrigEtaLayers;   
    TBranch *b_EFTAGPass;                      
    TBranch *b_EFTAGPt;                        
    TBranch *b_EFTAGEta;                       
    TBranch *b_EFTAGPhi;                       
    TBranch *b_L1nRoI;                         
    TBranch *b_L1isMoreCandInRoI;                         
    TBranch *b_L1Pass;                         
    TBranch *b_passedChain;                         
    TBranch *b_L1Pt;                           
    TBranch *b_L1Eta;                          
    TBranch *b_L1Phi;                          
    TBranch *b_L1RoINumber;                    
    TBranch *b_L1RoISector;                    
    TBranch *b_SAHypoPass;                     
    TBranch *b_SAOvRmPass;                     
    TBranch *b_SAOvRmPass_forClusAlg;                     
    TBranch *b_SAPt;                           
    TBranch *b_SAEta;                          
    TBranch *b_SAPhi;                          
    TBranch *b_SAEtaMS;                        
    TBranch *b_SAEtaBE;                        
    TBranch *b_SAPhiMS;                        
    TBranch *b_SAPhiBE;                        
    TBranch *b_SATGCPt;                        
    TBranch *b_SAPtBarrelRadius;               
    TBranch *b_SAPtBarrelSagitta;              
    TBranch *b_SAPtEndcapAlpha;                
    TBranch *b_SAPtEndcapBeta;                 
    TBranch *b_SAPtEndcapRadius;               
    TBranch *b_SACSCPt;                        
    TBranch *b_SAsAddress;                     
    TBranch *b_SASPR_BI;                       
    TBranch *b_SASPR_BM;                       
    TBranch *b_SASPR_BO;                       
    TBranch *b_SASPR_EI;                       
    TBranch *b_SASPR_EM;                       
    TBranch *b_SASPR_EO;                       
    TBranch *b_SASPR_EE;                       
    TBranch *b_SASPR_CSC;                      
    TBranch *b_SASPR_BEE;                      
    TBranch *b_SASPR_BME;                      
    TBranch *b_SASPZ_BI;                       
    TBranch *b_SASPZ_BM;                       
    TBranch *b_SASPZ_BO;                       
    TBranch *b_SASPZ_EI;                       
    TBranch *b_SASPZ_EM;                       
    TBranch *b_SASPZ_EO;                       
    TBranch *b_SASPZ_EE;                       
    TBranch *b_SASPZ_CSC;                      
    TBranch *b_SASPZ_BEE;                      
    TBranch *b_SASPZ_BME;                      
    TBranch *b_SASPSlope_BI;                   
    TBranch *b_SASPSlope_BM;                   
    TBranch *b_SASPSlope_BO;                   
    TBranch *b_SASPSlope_EI;                   
    TBranch *b_SASPSlope_EM;                   
    TBranch *b_SASPSlope_EO;                   
    TBranch *b_SASPSlope_EE;                   
    TBranch *b_SASPSlope_CSC;                  
    TBranch *b_SASPSlope_BEE;                  
    TBranch *b_SASPSlope_BME;                  
    TBranch *b_SASPIntercept_BI;               
    TBranch *b_SASPIntercept_BM;               
    TBranch *b_SASPIntercept_BO;               
    TBranch *b_SASPIntercept_EI;               
    TBranch *b_SASPIntercept_EM;               
    TBranch *b_SASPIntercept_EO;               
    TBranch *b_SASPIntercept_EE;               
    TBranch *b_SASPIntercept_CSC;              
    TBranch *b_SASPIntercept_BEE;              
    TBranch *b_SASPIntercept_BME;              
    TBranch *b_SASPChi2_BI;                    
    TBranch *b_SASPChi2_BM;                    
    TBranch *b_SASPChi2_BO;                    
    TBranch *b_SASPChi2_EI;                    
    TBranch *b_SASPChi2_EM;                    
    TBranch *b_SASPChi2_EO;                    
    TBranch *b_SASPChi2_EE;                    
    TBranch *b_SASPChi2_CSC;                   
    TBranch *b_SASPChi2_BEE;                   
    TBranch *b_SASPChi2_BME;                   
    TBranch *b_SARoIEta;                       
    TBranch *b_SARoIPhi;                       
    TBranch *b_SAisRPCFailure;                 
    TBranch *b_SAisTGCFailure;                 
    TBranch *b_SABarrelRadius;                 
    TBranch *b_SABarrelSagitta;                
    TBranch *b_SAEtaMap;                       
    TBranch *b_SAPhiMap;                       
    TBranch *b_SARoINumber;                    
    TBranch *b_SARoISector;                    
    TBranch *b_SARPCHitTime;                      
    TBranch *b_SARPCHitX;                      
    TBranch *b_SARPCHitY;                      
    TBranch *b_SARPCHitZ;                      
    TBranch *b_SARPCHitR;                      
    TBranch *b_SARPCHitEta;                    
    TBranch *b_SARPCHitPhi;                    
    TBranch *b_SARPCHitMeasuresPhi;            
    TBranch *b_SARPCHitLayer;            
    TBranch *b_SARPCHitStationName;            
    TBranch *b_SARPCHitStationNumber;          

    //RPC clustering
    TBranch *b_SAspcSetID;                      
    TBranch *b_SAptclus;                        
    TBranch *b_SAetaMSclus;                        
    TBranch *b_SAchargeclus;                        
    TBranch *b_SASPCZ_BI;                       
    TBranch *b_SASPCZ_BM;                       
    TBranch *b_SASPCZ_BO;                       
    TBranch *b_SASPCR_BI;                       
    TBranch *b_SASPCR_BM;                       
    TBranch *b_SASPCR_BO;                       
    TBranch *b_SARPCCluster_zMax;                      
    TBranch *b_SARPCCluster_zMin;                      
    TBranch *b_SARPCCluster_rMax;                      
    TBranch *b_SARPCCluster_rMin;                      
    TBranch *b_SARPCCluster_gX;                      
    TBranch *b_SARPCCluster_gY;                      
    TBranch *b_SARPCCluster_gZ;                      
    TBranch *b_SARPCCluster_stripWidth;                      
    TBranch *b_SARPCCluster_loX;                      
    TBranch *b_SARPCCluster_loY;                      
    TBranch *b_SARPCCluster_loZ;                      
    TBranch *b_SARPCCluster_clusterSize;                      
    TBranch *b_SARPCCluster_clusterLayer;                      
    TBranch *b_SARPCCluster_clusterMeasPhi;                      
    TBranch *b_SARPCCluster_fitInnPhi;                      
    TBranch *b_SARPCCluster_fitMidPhi;                      
    TBranch *b_SARPCCluster_fitOutPhi;                      
    TBranch *b_SARPCCluster_fitInnSlope;                      
    TBranch *b_SARPCCluster_fitMidSlope;                      
    TBranch *b_SARPCCluster_fitOutSlope;                      
    TBranch *b_SARPCCluster_fitInnOffset;                      
    TBranch *b_SARPCCluster_fitMidOffset;                      
    TBranch *b_SARPCCluster_fitOutOffset;                      
    TBranch *b_SARPCCluster_isSuccess;                      
    TBranch *b_SARPCCluster_isUsingMidCluster;                      
    TBranch *b_SARPCCluster_isPlausibleFitInnMid;                      
    TBranch *b_SARPCCluster_isPlausibleFitOut;                      
    TBranch *b_SARPCCluster_isPlausiblePhiInnMid;                      
    TBranch *b_SARPCCluster_isPlausiblePhiOut;                      
    TBranch *b_SARPCCluster_n_foundClusters;                      
    TBranch *b_SARPCCluster_id_clustersInSets;                      
    TBranch *b_SARPCCluster_id_clustersInPhiSets;                      

    TBranch *b_SARPCFitInnPhi;                 
    TBranch *b_SARPCFitInnSlope;               
    TBranch *b_SARPCFitInnOffset;              
    TBranch *b_SARPCFitMidPhi;                 
    TBranch *b_SARPCFitMidSlope;               
    TBranch *b_SARPCFitMidOffset;              
    TBranch *b_SARPCFitOutPhi;                 
    TBranch *b_SARPCFitOutSlope;               
    TBranch *b_SARPCFitOutOffset;              
    TBranch *b_SARoadAw;                       
    TBranch *b_SARoadBw;                       
    TBranch *b_SAzMin;                         
    TBranch *b_SAzMax;                         
    TBranch *b_SArMin;                         
    TBranch *b_SArMax;                         
    TBranch *b_SAEtaMin;                       
    TBranch *b_SAEtaMax;                       
    TBranch *b_SAMDTHitisOutlier;              
    TBranch *b_SAMDTHitChamber;                
    TBranch *b_SAMDTHitR;                      
    TBranch *b_SAMDTHitZ;                      
    TBranch *b_SAMDTHitPhi;                    
    TBranch *b_SAMDTHitResidual;               
    TBranch *b_SAMDTHitSpace;                  
    TBranch *b_SAMDTHitSigma;                  
    TBranch *b_SAMDTHitAllR;                      
    TBranch *b_SAMDTHitAllZ;                      
    TBranch *b_SAMDTHitAllclusRoadID;                      
    TBranch *b_SAMDTHitAllisOutlier;                      
    TBranch *b_CBHypoPass;                     
    TBranch *b_CBOvRmPass;                     
    TBranch *b_CBPt;                           
    TBranch *b_CBEta;                          
    TBranch *b_CBPhi;                          
    TBranch *b_EFPass;                         
    TBranch *b_EFPt;                           
    TBranch *b_EFEta;                          
    TBranch *b_EFPhi;                          
    // Histgrams

    TH1F *h_distRoIrpcEta;       //!
    TH1F *h_distRoIrpcPhi;       //!
    TH2D *h_SAEtaPhi;
    TH2D *h_SARoIEtaPhi;
    TH2D *h_OffSegZR;
    TH2D *h_OffSegEtaPhi;
    TH2D *h_OffSegNEtaPhiLayers;
    TH2D *h_OffSegXY_etaLayer[8];
    TH2D *h_OffSegXY;
    TH1D *h_nOffSegMid;
    TH1D *h_nOffSegOut;
    TH1D *h_SARpcHitEta;
    TH1D *h_SARpcHitPhi;
    TH1F *h_distnoCutRoIrpcEta;       //!
    TH1F *h_distnoCutRoIrpcPhi;       //!
    TH1D *h_distnoCutEachRoIrpcEta;       //!
    TH1D *h_distnoCutEachRoIrpcPhi;       //!
    TH1D *h_SARPCFitInnR;
    TH1D *h_SARPCFitInnEta;
    TH1D *h_SARPCFitMidR;
    TH1D *h_SARPCFitMidEta;
    TH1D *h_SARPCFitOutR;
    TH1D *h_SARPCFitOutEta;

    //histgrams of branch
    const int nBranch = 61;
    TH1D *h_Branch[61];
    TH2D *h_OffEtaPhi;
    TH2D *h_RoInumDelEta;
    TH2D *h_RoInumDelPhi;

    //histograms of Close by check
    TH1D *h_CloseByCount;
    TH1D *h_DeltaPt_OffSA_CloseBy;
    TH1D *h_DeltaPt_OffSA;
    TH1D *h_isCloseBydEdP;
    TH1D *h_isSameRoI;
    TH1D *h_isSameRpcFitMid;
    TH1D *h_isSameRpcFitMid_CloseBy;
    TH1D *h_distOffEta_CloseBy;
    TH1D *h_distOffPhi_CloseBy;
    TH1D *h_distOffEta;
    TH1D *h_distOffPhi;
    TH1D *h_distRpcFitInnEta_CloseBy;
    TH1D *h_distRpcFitInnPhi_CloseBy;
    TH1D *h_distRpcFitInnEta;
    TH1D *h_distRpcFitInnPhi;
    TH1D *h_distRpcFitMidEta_CloseBy;
    TH1D *h_distRpcFitMidPhi_CloseBy;
    TH1D *h_distRpcFitMidEta;
    TH1D *h_distRpcFitMidPhi;
    TH1D *h_distRpcFitOutEta_CloseBy;
    TH1D *h_distRpcFitOutPhi_CloseBy;
    TH1D *h_distRpcFitOutEta;
    TH1D *h_distRpcFitOutPhi;
    
    TH1D *h_distRpcHitToFitInnPhi;
    TH1D *h_distRpcHitToFitInnEta;
    TH1D *h_distRpcHitToFitMidPhi;
    TH1D *h_distRpcHitToFitMidEta;
    TH1D *h_distRpcHitToFitOutPhi;
    TH1D *h_distRpcHitToFitOutEta;
    TH1D *h_distRpcHitToFitInnPhi_CloseBy;
    TH1D *h_distRpcHitToFitInnEta_CloseBy;
    TH1D *h_distRpcHitToFitMidPhi_CloseBy;
    TH1D *h_distRpcHitToFitMidEta_CloseBy;
    TH1D *h_distRpcHitToFitOutPhi_CloseBy;
    TH1D *h_distRpcHitToFitOutEta_CloseBy;

    TH1D *h_rpcHitSize_CloseBy;
    TH1D *h_rpcHitSize;

    TH1D *h_rpcHitSizeLay_CloseBy[8];
    TH1D *h_rpcHitSizeLay[8];
    TH1D *h_rpcClusEtaSizeLay_CloseBy[8];
    TH1D *h_rpcClusPhiSizeLay_CloseBy[8];
    TH1D *h_DeltaRpcMidMinPhi_CloseBy;
    TH1D *h_DeltaRpcMidMinEta_CloseBy;
    TH1D *h_DeltaRpcOutMinPhi_CloseBy;
    TH1D *h_DeltaRpcOutMinEta_CloseBy;
    TH2D *h_rpcHitXYLay[8];
    TH1D *h_isSameSize_rpcHitEta;
    TH1D *h_isSameSize_rpcHitPhi;
    //RPC Clustrinng
    TH1D *h_rpc2ClusPhiNum_CloseBy;
    TH1D *h_rpc2ClusEtaNum_CloseBy;
    TH1D *h_rpcClusEta_CloseBy;
    TH1D *h_rpcClusPhi_CloseBy;
    TH1D *h_nClusRoad_withMidInfo;
    TH1D *h_nClusRoad_withOutInfo;

    TH2D *h_fitMidSlopeOffset_cut;
    TH2D *h_fitMidSlopeOffset_wocut;
    TH2D *m_nclus_vs_npt;

    TH1D *h_nClusRoadOut_woOutInfo;
    TH1D *h_nClusRoadOut_byMiddle;
    TH1D *h_nClusRoadOut_byOuter;
    TH1D *h_nClusRoadMid_byMiddleTest;
    TH1D *h_nClusRoadMid_byOuterTest;

    TH1D *h_distMidOffClusEta;
    TH1D *h_distMidOffClusPhi;
    TH1D *h_distMidOffClusR;
    TH1D *h_distMidOffClusR_CloseBy;
    TH1D *h_distMidOffClusR_noCloseBy;
    TH1D *h_etamin_MidOffFit_wMidInfo;
    TH1D *h_etamin_MidOffFit_woMidInfo;
    TH1D *h_slopemin_MidOffFit_wMidInfo;
    TH1D *h_slopemin_MidOffFit_woMidInfo;
    
    TH1D *h_etamin_OutOffFit;
    TH1D *h_slopemin_OutOffFit;
    TH1D *h_etamin_MidClusOutInfo_OffsegOut;
    TH1D *h_slopemin_MidClusOutInfo_OffsegOut;
    
    TH1D *h_etamin_MidOffFit_strip;
    TH1D *h_slopemin_MidOffFit_strip;
    //TH1D *h_fakeRoadMid[10];
    //TH1D *h_truthRoadMid[10];
    TH1D *h_nClusInSet;
    TH1D *h_nClusFitMid;
    TH1D *h_OffsegEta;
    TH1D *h_ClusFitMidEta;
    TH1D *h_ClusFitOutEta;
    TH1D *h_fakeRoadMid;
    TH1D *h_truthRoadMid;
    TH2D *h_faketruthRoadMid;
    TH1D *h_fakeRoadOut;
    TH1D *h_truthRoadOut;

    TH1D *h_fakeRoadMid_opt;
    TH1D *h_truthRoadMid_opt;
    TH1D *h_fakeRoadMid_opt23;
    TH1D *h_truthRoadMid_opt23;
    TH1D *h_fakeRoadMid_opt13;
    TH1D *h_truthRoadMid_opt13;
    TH1D *h_fakeRoadMid_opt12;
    TH1D *h_truthRoadMid_opt12;
    TH1D *h_fakeRoadOut_opt;
    TH1D *h_truthRoadOut_opt;
    TH1D *h_fakeRoadOut_opt23;
    TH1D *h_truthRoadOut_opt23;
    TH1D *h_fakeRoadOut_opt13;
    TH1D *h_truthRoadOut_opt13;
    TH1D *h_fakeRoadOut_opt12;
    TH1D *h_truthRoadOut_opt12;

    TH2D *h_faketruthRoadOut;
    TH2D *h_distMidOffClus;
    TH1D *h_distRoIClusR;
    TH1D *h_distRoIClusR_CloseBy;
    TH1D *h_distRoIClusR_noCloseBy;
    TH2D *h_distRoIClus;
    TH1D *h_distMidOffClusEtaMin;
    TH1D *h_distMidOffClusPhiMin;
    TH2D *h_distMidOffClusMin;
    TH1D *h_distMidOffClusfitEtaMin;
    TH1D *h_nClusLay2Phi;
    TH1D *h_nClusLay2Eta;
    TH1D *h_nClusLay3Phi;
    TH1D *h_nClusLay3Eta;
    TH1D *h_nClusLay23Phi;
    TH1D *h_nClusLay23Eta;
    TH1D *h_nClusLay23EtaPhi;
    TH1D *h_nClusHit23Phi;
    TH1D *h_nClusHit23Eta;
    TH1D *h_clusterSizeLay01;
    TH1D *h_clusterSizeLay23;
    TH1D *h_clusterSizeLayOut;
    TH1D *h_clusterSize;

    TH1D *h_nmaxClustersInSet;
    TH1D *h_iswoRpc1Road;
    TH1D *h_id_minlayer;
    TH1D *h_id_minlayer0;
    TH1D *h_seedcluster_setmax;
    TH1D *h_minlayerInSetmax;

    TH1D *h_extdR_bug;
    TH2D *h_RoIEtaPhi_bug;
    TH1D *h_extdR_wobug;
    TH2D *h_RoIEtaPhi_wobug;

    TH1D *h_countClusterRoadMid;
    TH1D *h_countClusterRoadOut;
    TH1D *h_countSPInn;
    TH1D *h_countSPMid;
    TH1D *h_countSPOut;
    TH1D *h_Zmin_OffClusRoadMid;
    TH1D *h_Zsubmin_OffClusRoadMid;
    TH1D *h_Z_plaumin_OffClusRoadMid;
    TH1D *h_Z_plausubmin_OffClusRoadMid;
    TH1D *h_n_clusRoadMid_total;
    TH1D *h_n_clusRoadMid_plau;
    TH1D *h_n_clusRoadMid_plau_n;
    TH1D *h_deltaZR_min_samePlane; 
    
    TH2D *m_hh_RoIEtaPhi_1road;

//  For calculation of efficiency
    TH2D *m_hh_Offpt_2mu1RoI;
    TH2D *m_hh_Offpt_2mu1RoI_oppcharge;
    TH1D *m_h_Offpt_2mu1RoI_lead;
    TH1D *m_h_Offpt_2mu1RoI_lead_oppcharge;
    TEfficiency* m_peff_Offpt_lead;
    TH1D *m_h_deltaExtR_2spall;
    TH1D *m_h_deltaExtR_2spinn;
    TH1D *m_h_deltaExtR_2spmid;
    TH1D *m_h_deltaExtR_2spout;
    TEfficiency* m_peff_2spall;
    TEfficiency* m_peff_2spinn;
    TEfficiency* m_peff_2spmid;
    TEfficiency* m_peff_2spout;
    TH1D *m_h_deltaExtR_3pt;
    TH1D *m_h_deltaExtR_2pt;
    TH1D *m_h_deltaExtR_1pt;
    TH1D *m_h_deltaExtR_2pt_oppcharge;
    TEfficiency* m_peff_3pt;
    TEfficiency* m_peff_2pt;
    TEfficiency* m_peff_1pt;
    TEfficiency* m_peff_2pt_oppcharge;
    TH1D *m_h_deltaExtR_off_tot;
    TH1D *m_h_deltaExtR_off_morecand;
    TH1D *m_h_eff_isMorecand;
    TEfficiency* m_peff_morecand;
    TH1D *m_h_deltaExtR_off_tot_clus;
    TH1D *m_h_deltaExtR_off_2road_clus;
    TH1D *m_h_deltaExtR_off_1road_clus;
    TH1D *m_h_deltaExtR_off_manyroad_clus;
    TH1D *m_h_eff_clusterRoad;
    TEfficiency* m_peff_clusRoad;
    TEfficiency* m_peff_1clusRoad;
    TEfficiency* m_peff_manyclusRoad;
    TH1D *m_h_deltaExtR_off_mismatch;
    TH1D *m_h_deltaExtR_off_closeBy;
    TEfficiency* m_peff_closeBy;
    TH1D *m_h_dR_1road_sep;
    TH1D *m_h_dR_1road_canthelp;
    TH1D *m_h_dR_1road_1clus;
    TH1D *m_h_dR_1road_1phi;
    TH1D *m_h_dR_1road_phiothers;
    TEfficiency* m_peff_1clus;
    TEfficiency* m_peff_1phi;
    TEfficiency* m_peff_phiothers;
    TH1D *m_h_dR_1road_outlier;
    TH1D *m_h_dR_1road_misdire;
    TH1D *m_h_dR_1road_others;
    TH1D *m_h_dR_3road;
    TH1D *m_h_dR_3road_3clus;
    TEfficiency* m_peff_3road;
    TEfficiency* m_peff_canthelp;
    TEfficiency* m_peff_outlier;
    TEfficiency* m_peff_misdire;
    TEfficiency* m_peff_others;

    TH1D* m_nclusPerLayer[4];

    TH2D *hh_num_clusRoadMidvsOut;
    TH2D *hh_n_CRvsSPCMid;
    TProfile *prof_CRvsSPCMid;

    TH1D *m_h_2CR_1SPCMid;
    TH1D *m_h_2CR_1SPCMid_u3mdt;
    TH1D *m_h_2CR_1SPCMid_u3mdt_selected;
    TEfficiency* m_peff_2CR_1SPCMid;
    TEfficiency* m_peff_2CR_1SPCMid_all;

    TH1D *m_h_dSlope_innmid;
    TH1D *m_ptres_default;
    TH1D *m_ptres_mpsa;

    TH1D *m_h_SApt;
    TH1D *h_Offpt;
    TH1D *h_Offpt_lead;
    TH1D *h_Offpt_sublead;
    TH1D *h_OffdR;
    TH1D *m_h_offlinePt_res0p7;
    TH1D *h_dtheta_sub;
    TH1D *h_dtheta_total;
    TH1D *h_dtheta;
    TH1D *m_h_ptclus;
    TH1D *m_h_ptres_2off2clus;
    TH2D *m_h_ptclus_vs_off;
    TH1D *m_h_ptres_2off2clus_sub;
    TH1D *m_h_ptres_2off2clus_leading;
    TH1D *m_h_countPtclus;
    TH1D *m_h_n_lowptclus;
    TH1D *m_h_n_lowptclus_patA;
    TH1D *m_h_n_lowptclus_patA_sameIDout;
    TH1D *m_h_n_lowptclus_patA_findout;
    TH1D *m_h_n_1spstat_inn;
    TH1D *m_h_n_1spstat_mid;
    TH1D *m_h_n_1spstat_out;
    TH1D *m_h_n_spall;
    TH1D *m_h_Nclus_npt1;
    TH1D *m_h_sector;
    TH2D *m_h_ptres_vs_offsegInndR;
    TH2D *m_h_ptres_vs_offsegMiddR;
    TH2D *m_h_ptres_vs_offsegOutdR;
    TH2D *m_h_ptres_vs_offsegMindR;
     
    RPC_FCBM(TChain *tree);
    virtual ~RPC_FCBM();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Loop_FCBM( int Nevents , RpcL2MuonSA::HistData histData);
    virtual void     Display_FCBM(int tap_type, int trig_chain, Long64_t begin_entry, Long64_t limit_entry, TString pdf);
    virtual void     End();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    
    //legend for display
    virtual void setLegend_eventInfo(TText* eventInfo);
    virtual void setLegend_leftSideInfo(TLegend* OffMu1Info, string& plane);
    virtual void setLegend_rightSideInfo(TLegend* OffMu2Info, string& plane);

    //for R-Z event display
    //
    virtual void setFrameRegion(const char* station, float& Zmin, float& Rmin, float& Zmax, float& Rmax);
//    virtual void setFrameRegion(const char* station, double& Zmin, double& Rmin, double& Zmax, double& Rmax);
    virtual void setMdtRegion(unsigned int iMuon, const char* station, float* Z, float* R);
    virtual void setMdtAllRegion(int iClus, const char* station, float* Z, float* R);

    // Initialize
    virtual void Init(TTree *tree);
    virtual void InitHist();
    virtual void InitFCBM();
    virtual void InitNoCut();
    virtual void InitBranchHist();
    virtual void FillFCBM(vector<pair<int, int>>& MUpair);
    virtual void setRPCFitEtaR();
    virtual void setRPCHitRegion(unsigned int nMuon, vector<vector<double>>*& SARPCHitEta, vector<vector<double>>*& SARPCHitPhi, double* rpcHitRegion);
    virtual void FillnoFCBM();
    virtual void FillNoCut();
    virtual void FillClusterInfo();
    virtual void FillEfficiency();
    virtual int  EtaRegion( double eta );
    virtual int  ECWeakRegion( double eta, double phi );
    virtual void FillBranchHist();
    virtual void DrawFCBM(TString pdf);
    virtual void DrawNoCut(TString pdf);
    virtual void DrawBranchHist(TString pdf);
    virtual void CalcEfficiency(TH1D* h_num, TH1D* h_den, TH1D* h_set);
    virtual void DrawEfficiency(TString pdf);
    virtual bool is2muIn1RoI();
    virtual bool is2muInBarrel();
    void RpcClusterSetter(TGraphErrors* grRPC,
        vector<float>& clusgX,
        vector<float>& clusgY,
        vector<float>& clusgZ,
        vector<int>& clusLayer,
        vector<int>& clusMeasPhi );

    double invMass(double m1, double pt1, double eta1, float phi1,
        double m2, double pt2, double eta2, float phi2) const;
    double dR(double eta1, float phi1, double eta2, float phi2) const;
    
    template <class Type> void setFillBranch_v2(TH1D*& h_Branch, vector<vector<Type>>*& BranchName);
    template <class Type> void setFillBranch_v1(TH1D*& h_Branch, vector<Type>*& BranchName);
    template <class Type> void setFillBranch(TH1D*& h_Branch, Type& BranchName);
    virtual bool isCloseByCut(UInt_t nMuon, vector<pair<int, int>>& MUpair);
    virtual bool isCloseBy(int iMuon); 
    virtual void getDistToEachRoI(vector<float>* &roiPos, vector< vector<double>>* &paramPos, vector<vector<double>>& getParamPos,vector<vector<UInt_t>>* &MesPhi, bool isPhi);
    void getCloseBy(auto& param, auto& selectedParam);

    bool isOld{false};
};

#endif

#ifdef RPC_FCBM_cxx
#ifndef RPC_FCBM_hhh
#define RPC_FCBM_hhh
RPC_FCBM::RPC_FCBM(TChain *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");
    }
    f->GetObject("t_tap",tree);

  }
  Init(tree);
}

RPC_FCBM::~RPC_FCBM()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t RPC_FCBM::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t RPC_FCBM::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
    
  }
  return centry;
}

void RPC_FCBM::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  SARPCFitInnR.clear();
  SARPCFitInnEta.clear();
  SARPCFitMidR.clear();
  SARPCFitMidEta.clear();
  SARPCFitOutR.clear();
  SARPCFitOutEta.clear();

  EventNumber = 0;
  RunNumber = 0;
  LumiBlock = 0;
  AverageInteractionsPerCrossing = 0;
  mes_name = 0;
  nMuon = 0;
  OfflinePt = 0;                      
  OfflineEta = 0;                     
  OfflineExtEta = 0;                  
  OfflineExtInnEta = 0;               
  OfflinePhi = 0;                     
  OfflineExtPhi = 0;                  
  OfflineExtInnPhi = 0;               
  OfflineD0 = 0;                      
  OfflineZ0 = 0;                      
  OfflineCharge = 0;                  
  OfflinePhi = 0;                     
  OfflineExtPhi = 0;                  
  OfflineExtInnPhi = 0;               
  OfflineD0 = 0;                      
  OfflineZ0 = 0;                      
  OfflineCharge = 0;                  
  OfflineNumSegment = 0;              
  OfflineSegmentX = 0;                
  OfflineSegmentY = 0;                
  OfflineSegmentZ = 0;                
  OfflineSegmentPx = 0;               
  OfflineSegmentPy = 0;               
  OfflineSegmentPz = 0;               
  OfflineSegmentChiSquared = 0;       
  OfflineSegmentNumberDoF = 0;        
  OfflineSegmentSector = 0;           
  OfflineSegmentChamberIndex = 0;     
  OfflineSegmentEtaIndex = 0;         
  OfflineSegmentNPrecisionHits = 0;   
  OfflineSegmentNPhiLayers = 0;       
  OfflineSegmentNTrigEtaLayers = 0;   
  EFTAGPass = 0;                      
  EFTAGPt = 0;                        
  EFTAGEta = 0;                       
  EFTAGPhi = 0;                       
  L1nRoI = 0;                           
  L1isMoreCandInRoI = 0;                           
  L1Pass = 0;                           
  passedChain = 0;                           
  L1Pt = 0;                             
  L1Eta = 0;                            
  L1Phi = 0;                            
  L1RoINumber = 0;                      
  L1RoISector = 0;                      
  SAHypoPass = 0;                       
  SAOvRmPass = 0;                       
  SAOvRmPass_forClusAlg = 0;                       
  SAPt = 0;                             
  SAEta = 0;                            
  SAPhi = 0;                            
  SAEtaMS = 0;                          
  SAEtaBE = 0;                          
  SAPhiMS = 0;                          
  SAPhiBE = 0;                          
  SATGCPt = 0;                          
  SAPtBarrelRadius = 0;                 
  SAPtBarrelSagitta = 0;                
  SAPtEndcapAlpha = 0;                  
  SAPtEndcapBeta = 0;                   
  SAPtEndcapRadius = 0;                 
  SACSCPt = 0;                          
  SAsAddress = 0;                       
  SASPR_BI = 0;                         
  SASPR_BM = 0;                         
  SASPR_BO = 0;                         
  SASPR_EI = 0;                         
  SASPR_EM = 0;                         
  SASPR_EO = 0;                         
  SASPR_EE = 0;                         
  SASPR_CSC = 0;                        
  SASPR_BEE = 0;                        
  SASPR_BME = 0;                        
  SASPZ_BI = 0;                         
  SASPZ_BM = 0;                         
  SASPZ_BO = 0;                         
  SASPZ_EI = 0;                         
  SASPZ_EM = 0;                         
  SASPZ_EO = 0;                         
  SASPZ_EE = 0;                         
  SASPZ_CSC = 0;                        
  SASPZ_BEE = 0;                        
  SASPZ_BME = 0;                        
  SASPSlope_BI = 0;                     
  SASPSlope_BM = 0;                     
  SASPSlope_BO = 0;                     
  SASPSlope_EI = 0;                     
  SASPSlope_EM = 0;                     
  SASPSlope_EO = 0;                     
  SASPSlope_EE = 0;                     
  SASPSlope_CSC = 0;                    
  SASPSlope_BEE = 0;                    
  SASPSlope_BME = 0;                    
  SASPIntercept_BI = 0;                 
  SASPIntercept_BM = 0;                 
  SASPIntercept_BO = 0;                 
  SASPIntercept_EI = 0;                 
  SASPIntercept_EM = 0;                 
  SASPIntercept_EO = 0;                 
  SASPIntercept_EE = 0;                 
  SASPIntercept_CSC = 0;                
  SASPIntercept_BEE = 0;                
  SASPIntercept_BME = 0;                
  SASPChi2_BI = 0;                      
  SASPChi2_BM = 0;                      
  SASPChi2_BO = 0;                      
  SASPChi2_EI = 0;                      
  SASPChi2_EM = 0;                      
  SASPChi2_EO = 0;                     
  SASPChi2_EE = 0;                     
  SASPChi2_CSC = 0;                    
  SASPChi2_BEE = 0;                    
  SASPChi2_BME = 0;                    
  SARoIEta = 0;                        
  SARoIPhi = 0;                        
  SAisRPCFailure = 0;                 
  SAisTGCFailure = 0;                 
  SABarrelRadius = 0;                 
  SABarrelSagitta = 0;                
  SAEtaMap = 0;                       
  SAPhiMap = 0;                       
  SARoINumber = 0;                    
  SARoISector = 0;                    
  SARPCHitTime = 0;                      
  SARPCHitX = 0;                      
  SARPCHitY = 0;                      
  SARPCHitZ = 0;                      
  SARPCHitR = 0;                      
  SARPCHitEta = 0;                    
  SARPCHitPhi = 0;                    
  SARPCHitMeasuresPhi = 0;            
  SARPCHitLayer = 0;            
  SARPCHitStationName = 0;            
  SARPCHitStationNumber = 0;          
  //RPC clustering
  SAspcSetID = 0;                      
  SAptclus = 0;                          
  SAetaMSclus = 0;                          
  SAchargeclus = 0;                          
  SASPCZ_BI = 0;                         
  SASPCZ_BM = 0;                         
  SASPCZ_BO = 0;                         
  SASPCR_BI = 0;                         
  SASPCR_BM = 0;                         
  SASPCR_BO = 0;                         
  SARPCCluster_zMax = 0;                      
  SARPCCluster_zMin = 0;                      
  SARPCCluster_rMax = 0;                      
  SARPCCluster_rMin = 0;                      
  SARPCCluster_gX = 0;                      
  SARPCCluster_gY = 0;                      
  SARPCCluster_gZ = 0;                      
  SARPCCluster_stripWidth = 0;                      
  SARPCCluster_loX = 0;                      
  SARPCCluster_loY = 0;                      
  SARPCCluster_loZ = 0;                      
  SARPCCluster_clusterSize = 0;                      
  SARPCCluster_clusterLayer = 0;                      
  SARPCCluster_clusterMeasPhi = 0;                      
  SARPCCluster_fitInnPhi = 0;                      
  SARPCCluster_fitMidPhi = 0;                      
  SARPCCluster_fitOutPhi = 0;                      
  SARPCCluster_fitInnSlope = 0;                      
  SARPCCluster_fitMidSlope = 0;                      
  SARPCCluster_fitOutSlope = 0;                      
  SARPCCluster_fitInnOffset = 0;                      
  SARPCCluster_fitMidOffset = 0;                      
  SARPCCluster_fitOutOffset = 0;                      
  SARPCCluster_isSuccess = 0;                      
  SARPCCluster_isUsingMidCluster = 0;                      
  SARPCCluster_isPlausibleFitInnMid = 0;                      
  SARPCCluster_isPlausibleFitOut = 0;                      
  SARPCCluster_isPlausiblePhiInnMid = 0;                      
  SARPCCluster_isPlausiblePhiOut = 0;                      
  SARPCCluster_n_foundClusters = 0;                      
  SARPCCluster_id_clustersInSets = 0;                      
  SARPCCluster_id_clustersInPhiSets = 0;                      

  SARPCFitInnPhi = 0;                 
  SARPCFitInnSlope = 0;               
  SARPCFitInnOffset = 0;              
  SARPCFitMidPhi = 0;                 
  SARPCFitMidSlope = 0;               
  SARPCFitMidOffset = 0;              
  SARPCFitOutPhi = 0;                 
  SARPCFitOutSlope = 0;               
  SARPCFitOutOffset = 0;              
  SARoadAw = 0;                       
  SARoadBw = 0;                       
  SAzMin = 0;                         
  SAzMax = 0;                         
  SArMin = 0;                         
  SArMax = 0;                         
  SAEtaMin = 0;                       
  SAEtaMax = 0;                       
  SAMDTHitisOutlier = 0;              
  SAMDTHitChamber = 0;                
  SAMDTHitR = 0;                      
  SAMDTHitZ = 0;                      
  SAMDTHitPhi = 0;                    
  SAMDTHitResidual = 0;               
  SAMDTHitSpace = 0;                  
  SAMDTHitSigma = 0;                  
  SAMDTHitAllR = 0;                      
  SAMDTHitAllZ = 0;                      
  SAMDTHitAllisOutlier = 0;                      
  SAMDTHitAllclusRoadID = 0;                      
  CBHypoPass = 0;                     
  CBOvRmPass = 0;                     
  CBPt = 0;                           
  CBEta = 0;                          
  CBPhi = 0;                          
  EFPass = 0;                         
  EFPt = 0;                           
  EFEta = 0;                          
  EFPhi = 0;                          
  
  // Set object pointer
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("EventNumber",                       &EventNumber,                       &b_EventNumber                       );
  fChain->SetBranchAddress("RunNumber",                         &RunNumber,                         &b_RunNumber                         );
  fChain->SetBranchAddress("LumiBlock",                         &LumiBlock,                         &b_LumiBlock                         );
  fChain->SetBranchAddress("AverageInteractionsPerCrossing",    &AverageInteractionsPerCrossing,    &b_AverageInteractionsPerCrossing    );
  fChain->SetBranchAddress("mes_name",                          &mes_name,                          &b_mes_name                          );
  fChain->SetBranchAddress("nMuon",                            &nMuon,                            &b_nMuon                            );
  fChain->SetBranchAddress("OfflinePt",                      &OfflinePt,                      &b_OfflinePt                      );
  fChain->SetBranchAddress("OfflineEta",                     &OfflineEta,                     &b_OfflineEta                     );
  fChain->SetBranchAddress("OfflineExtEta",                  &OfflineExtEta,                  &b_OfflineExtEta                  );
  fChain->SetBranchAddress("OfflineExtInnEta",               &OfflineExtInnEta,               &b_OfflineExtInnEta               );
  fChain->SetBranchAddress("OfflinePhi",                     &OfflinePhi,                     &b_OfflinePhi                     );
  fChain->SetBranchAddress("OfflineExtPhi",                  &OfflineExtPhi,                  &b_OfflineExtPhi                  );
  fChain->SetBranchAddress("OfflineExtInnPhi",               &OfflineExtInnPhi,               &b_OfflineExtInnPhi               );
  fChain->SetBranchAddress("OfflineD0",                      &OfflineD0,                      &b_OfflineD0                      );
  fChain->SetBranchAddress("OfflineZ0",                      &OfflineZ0,                      &b_OfflineZ0                      );
  fChain->SetBranchAddress("OfflineCharge",                  &OfflineCharge,                  &b_OfflineCharge                  );
  fChain->SetBranchAddress("OfflineNumSegment",              &OfflineNumSegment,              &b_OfflineNumSegment              );
  fChain->SetBranchAddress("OfflineSegmentX",                &OfflineSegmentX,                &b_OfflineSegmentX                );
  fChain->SetBranchAddress("OfflineSegmentY",                &OfflineSegmentY,                &b_OfflineSegmentY                );
  fChain->SetBranchAddress("OfflineSegmentZ",                &OfflineSegmentZ,                &b_OfflineSegmentZ                );
  fChain->SetBranchAddress("OfflineSegmentPx",               &OfflineSegmentPx,               &b_OfflineSegmentPx               );
  fChain->SetBranchAddress("OfflineSegmentPy",               &OfflineSegmentPy,               &b_OfflineSegmentPy               );
  fChain->SetBranchAddress("OfflineSegmentPz",               &OfflineSegmentPz,               &b_OfflineSegmentPz               );
  fChain->SetBranchAddress("OfflineSegmentChiSquared",       &OfflineSegmentChiSquared,       &b_OfflineSegmentChiSquared       );
  fChain->SetBranchAddress("OfflineSegmentNumberDoF",        &OfflineSegmentNumberDoF,        &b_OfflineSegmentNumberDoF        );
  fChain->SetBranchAddress("OfflineSegmentSector",           &OfflineSegmentSector,           &b_OfflineSegmentSector           );
  fChain->SetBranchAddress("OfflineSegmentChamberIndex",     &OfflineSegmentChamberIndex,     &b_OfflineSegmentChamberIndex     );
  fChain->SetBranchAddress("OfflineSegmentEtaIndex",         &OfflineSegmentEtaIndex,         &b_OfflineSegmentEtaIndex         );
  fChain->SetBranchAddress("OfflineSegmentNPrecisionHits",   &OfflineSegmentNPrecisionHits,   &b_OfflineSegmentNPrecisionHits   );
  fChain->SetBranchAddress("OfflineSegmentNPhiLayers",       &OfflineSegmentNPhiLayers,       &b_OfflineSegmentNPhiLayers       );
  fChain->SetBranchAddress("OfflineSegmentNTrigEtaLayers",   &OfflineSegmentNTrigEtaLayers,   &b_OfflineSegmentNTrigEtaLayers   );
  fChain->SetBranchAddress("EFTAGPass",                      &EFTAGPass,                      &b_EFTAGPass                      );
  fChain->SetBranchAddress("EFTAGPt",                        &EFTAGPt,                        &b_EFTAGPt                        );
  fChain->SetBranchAddress("EFTAGEta",                       &EFTAGEta,                       &b_EFTAGEta                       );
  fChain->SetBranchAddress("EFTAGPhi",                       &EFTAGPhi,                       &b_EFTAGPhi                       );
  fChain->SetBranchAddress("L1nRoI",                         &L1nRoI,                         &b_L1nRoI                         );
  fChain->SetBranchAddress("L1isMoreCandInRoI",                         &L1isMoreCandInRoI,                         &b_L1isMoreCandInRoI                         );
  fChain->SetBranchAddress("L1Pass",                         &L1Pass,                         &b_L1Pass                         );
  fChain->SetBranchAddress("passedChain",                         &passedChain,                         &b_passedChain                         );
  fChain->SetBranchAddress("L1Pt",                           &L1Pt,                           &b_L1Pt                           );
  fChain->SetBranchAddress("L1Eta",                          &L1Eta,                          &b_L1Eta                          );
  fChain->SetBranchAddress("L1Phi",                          &L1Phi,                          &b_L1Phi                          );
  fChain->SetBranchAddress("L1RoINumber",                    &L1RoINumber,                    &b_L1RoINumber                    );
  fChain->SetBranchAddress("L1RoISector",                    &L1RoISector,                    &b_L1RoISector                    );
  fChain->SetBranchAddress("SAHypoPass",                     &SAHypoPass,                     &b_SAHypoPass                     );
  fChain->SetBranchAddress("SAOvRmPass",                     &SAOvRmPass,                     &b_SAOvRmPass                     );
  fChain->SetBranchAddress("SAOvRmPass_forClusAlg",                     &SAOvRmPass_forClusAlg,                     &b_SAOvRmPass_forClusAlg                     );
  fChain->SetBranchAddress("SAPt",                           &SAPt,                           &b_SAPt                           );
  fChain->SetBranchAddress("SAEta",                          &SAEta,                          &b_SAEta                          );
  fChain->SetBranchAddress("SAPhi",                          &SAPhi,                          &b_SAPhi                          );
  fChain->SetBranchAddress("SAEtaMS",                        &SAEtaMS,                        &b_SAEtaMS                        );
  fChain->SetBranchAddress("SAEtaBE",                        &SAEtaBE,                        &b_SAEtaBE                        );
  fChain->SetBranchAddress("SAPhiMS",                        &SAPhiMS,                        &b_SAPhiMS                        );
  fChain->SetBranchAddress("SAPhiBE",                        &SAPhiBE,                        &b_SAPhiBE                        );
  fChain->SetBranchAddress("SATGCPt",                        &SATGCPt,                        &b_SATGCPt                        );
  fChain->SetBranchAddress("SAPtBarrelRadius",               &SAPtBarrelRadius,               &b_SAPtBarrelRadius               );
  fChain->SetBranchAddress("SAPtBarrelSagitta",              &SAPtBarrelSagitta,              &b_SAPtBarrelSagitta              );
  fChain->SetBranchAddress("SAPtEndcapAlpha",                &SAPtEndcapAlpha,                &b_SAPtEndcapAlpha                );
  fChain->SetBranchAddress("SAPtEndcapBeta",                 &SAPtEndcapBeta,                 &b_SAPtEndcapBeta                 );
  fChain->SetBranchAddress("SAPtEndcapRadius",               &SAPtEndcapRadius,               &b_SAPtEndcapRadius               );
  fChain->SetBranchAddress("SACSCPt",                        &SACSCPt,                        &b_SACSCPt                        );
  fChain->SetBranchAddress("SAsAddress",                     &SAsAddress,                     &b_SAsAddress                     );
  fChain->SetBranchAddress("SASPR_BI",                       &SASPR_BI,                       &b_SASPR_BI                       );
  fChain->SetBranchAddress("SASPR_BM",                       &SASPR_BM,                       &b_SASPR_BM                       );
  fChain->SetBranchAddress("SASPR_BO",                       &SASPR_BO,                       &b_SASPR_BO                       );
  fChain->SetBranchAddress("SASPR_EI",                       &SASPR_EI,                       &b_SASPR_EI                       );
  fChain->SetBranchAddress("SASPR_EM",                       &SASPR_EM,                       &b_SASPR_EM                       );
  fChain->SetBranchAddress("SASPR_EO",                       &SASPR_EO,                       &b_SASPR_EO                       );
  fChain->SetBranchAddress("SASPR_EE",                       &SASPR_EE,                       &b_SASPR_EE                       );
  fChain->SetBranchAddress("SASPR_CSC",                      &SASPR_CSC,                      &b_SASPR_CSC                      );
  fChain->SetBranchAddress("SASPR_BEE",                      &SASPR_BEE,                      &b_SASPR_BEE                      );
  fChain->SetBranchAddress("SASPR_BME",                      &SASPR_BME,                      &b_SASPR_BME                      );
  fChain->SetBranchAddress("SASPZ_BI",                       &SASPZ_BI,                       &b_SASPZ_BI                       );
  fChain->SetBranchAddress("SASPZ_BM",                       &SASPZ_BM,                       &b_SASPZ_BM                       );
  fChain->SetBranchAddress("SASPZ_BO",                       &SASPZ_BO,                       &b_SASPZ_BO                       );
  fChain->SetBranchAddress("SASPZ_EI",                       &SASPZ_EI,                       &b_SASPZ_EI                       );
  fChain->SetBranchAddress("SASPZ_EM",                       &SASPZ_EM,                       &b_SASPZ_EM                       );
  fChain->SetBranchAddress("SASPZ_EO",                       &SASPZ_EO,                       &b_SASPZ_EO                       );
  fChain->SetBranchAddress("SASPZ_EE",                       &SASPZ_EE,                       &b_SASPZ_EE                       );
  fChain->SetBranchAddress("SASPZ_CSC",                      &SASPZ_CSC,                      &b_SASPZ_CSC                      );
  fChain->SetBranchAddress("SASPZ_BEE",                      &SASPZ_BEE,                      &b_SASPZ_BEE                      );
  fChain->SetBranchAddress("SASPZ_BME",                      &SASPZ_BME,                      &b_SASPZ_BME                      );
  fChain->SetBranchAddress("SASPSlope_BI",                   &SASPSlope_BI,                   &b_SASPSlope_BI                   );
  fChain->SetBranchAddress("SASPSlope_BM",                   &SASPSlope_BM,                   &b_SASPSlope_BM                   );
  fChain->SetBranchAddress("SASPSlope_BO",                   &SASPSlope_BO,                   &b_SASPSlope_BO                   );
  fChain->SetBranchAddress("SASPSlope_EI",                   &SASPSlope_EI,                   &b_SASPSlope_EI                   );
  fChain->SetBranchAddress("SASPSlope_EM",                   &SASPSlope_EM,                   &b_SASPSlope_EM                   );
  fChain->SetBranchAddress("SASPSlope_EO",                   &SASPSlope_EO,                   &b_SASPSlope_EO                   );
  fChain->SetBranchAddress("SASPSlope_EE",                   &SASPSlope_EE,                   &b_SASPSlope_EE                   );
  fChain->SetBranchAddress("SASPSlope_CSC",                  &SASPSlope_CSC,                  &b_SASPSlope_CSC                  );
  fChain->SetBranchAddress("SASPSlope_BEE",                  &SASPSlope_BEE,                  &b_SASPSlope_BEE                  );
  fChain->SetBranchAddress("SASPSlope_BME",                  &SASPSlope_BME,                  &b_SASPSlope_BME                  );
  fChain->SetBranchAddress("SASPIntercept_BI",               &SASPIntercept_BI,               &b_SASPIntercept_BI               );
  fChain->SetBranchAddress("SASPIntercept_BM",               &SASPIntercept_BM,               &b_SASPIntercept_BM               );
  fChain->SetBranchAddress("SASPIntercept_BO",               &SASPIntercept_BO,               &b_SASPIntercept_BO               );
  fChain->SetBranchAddress("SASPIntercept_EI",               &SASPIntercept_EI,               &b_SASPIntercept_EI               );
  fChain->SetBranchAddress("SASPIntercept_EM",               &SASPIntercept_EM,               &b_SASPIntercept_EM               );
  fChain->SetBranchAddress("SASPIntercept_EO",               &SASPIntercept_EO,               &b_SASPIntercept_EO               );
  fChain->SetBranchAddress("SASPIntercept_EE",               &SASPIntercept_EE,               &b_SASPIntercept_EE               );
  fChain->SetBranchAddress("SASPIntercept_CSC",              &SASPIntercept_CSC,              &b_SASPIntercept_CSC              );
  fChain->SetBranchAddress("SASPIntercept_BEE",              &SASPIntercept_BEE,              &b_SASPIntercept_BEE              );
  fChain->SetBranchAddress("SASPIntercept_BME",              &SASPIntercept_BME,              &b_SASPIntercept_BME              );
  fChain->SetBranchAddress("SASPChi2_BI",                    &SASPChi2_BI,                    &b_SASPChi2_BI                    );
  fChain->SetBranchAddress("SASPChi2_BM",                    &SASPChi2_BM,                    &b_SASPChi2_BM                    );
  fChain->SetBranchAddress("SASPChi2_BO",                    &SASPChi2_BO,                    &b_SASPChi2_BO                    );
  fChain->SetBranchAddress("SASPChi2_EI",                    &SASPChi2_EI,                    &b_SASPChi2_EI                    );
  fChain->SetBranchAddress("SASPChi2_EM",                    &SASPChi2_EM,                    &b_SASPChi2_EM                    );
  fChain->SetBranchAddress("SASPChi2_EO",                    &SASPChi2_EO,                    &b_SASPChi2_EO                    );
  fChain->SetBranchAddress("SASPChi2_EE",                    &SASPChi2_EE,                    &b_SASPChi2_EE                    );
  fChain->SetBranchAddress("SASPChi2_CSC",                   &SASPChi2_CSC,                   &b_SASPChi2_CSC                   );
  fChain->SetBranchAddress("SASPChi2_BEE",                   &SASPChi2_BEE,                   &b_SASPChi2_BEE                   );
  fChain->SetBranchAddress("SASPChi2_BME",                   &SASPChi2_BME,                   &b_SASPChi2_BME                   );
  fChain->SetBranchAddress("SARoIEta",                       &SARoIEta,                       &b_SARoIEta                       );
  fChain->SetBranchAddress("SARoIPhi",                       &SARoIPhi,                       &b_SARoIPhi                       );
  fChain->SetBranchAddress("SAisRPCFailure",                 &SAisRPCFailure,                 &b_SAisRPCFailure                 );
  fChain->SetBranchAddress("SAisTGCFailure",                 &SAisTGCFailure,                 &b_SAisTGCFailure                 );
  fChain->SetBranchAddress("SABarrelRadius",                 &SABarrelRadius,                 &b_SABarrelRadius                 );
  fChain->SetBranchAddress("SABarrelSagitta",                &SABarrelSagitta,                &b_SABarrelSagitta                );
  fChain->SetBranchAddress("SAEtaMap",                       &SAEtaMap,                       &b_SAEtaMap                       );
  fChain->SetBranchAddress("SAPhiMap",                       &SAPhiMap,                       &b_SAPhiMap                       );
  fChain->SetBranchAddress("SARoINumber",                    &SARoINumber,                    &b_SARoINumber                    );
  fChain->SetBranchAddress("SARoISector",                    &SARoISector,                    &b_SARoISector                    );
  fChain->SetBranchAddress("SARPCHitTime",                   &SARPCHitTime,                   &b_SARPCHitTime                   );
  fChain->SetBranchAddress("SARPCHitX",                      &SARPCHitX,                      &b_SARPCHitX                      );
  fChain->SetBranchAddress("SARPCHitY",                      &SARPCHitY,                      &b_SARPCHitY                      );
  fChain->SetBranchAddress("SARPCHitZ",                      &SARPCHitZ,                      &b_SARPCHitZ                      );
  fChain->SetBranchAddress("SARPCHitR",                      &SARPCHitR,                      &b_SARPCHitR                      );
  fChain->SetBranchAddress("SARPCHitEta",                    &SARPCHitEta,                    &b_SARPCHitEta                    );
  fChain->SetBranchAddress("SARPCHitPhi",                    &SARPCHitPhi,                    &b_SARPCHitPhi                    );
  fChain->SetBranchAddress("SARPCHitMeasuresPhi",            &SARPCHitMeasuresPhi,            &b_SARPCHitMeasuresPhi            );
  fChain->SetBranchAddress("SARPCHitLayer",                  &SARPCHitLayer,                  &b_SARPCHitLayer                  );
  fChain->SetBranchAddress("SARPCHitStationName",            &SARPCHitStationName,            &b_SARPCHitStationName            );
  fChain->SetBranchAddress("SARPCHitStationNumber",          &SARPCHitStationNumber,          &b_SARPCHitStationNumber          );
  //RPC clustering
  fChain->SetBranchAddress("SAspcSetID",                &SAspcSetID,                &b_SAspcSetID                );
  fChain->SetBranchAddress("SAptclus",                        &SAptclus,                        &b_SAptclus                        );
  fChain->SetBranchAddress("SAetaMSclus",                        &SAetaMSclus,                        &b_SAetaMSclus                        );
  fChain->SetBranchAddress("SAchargeclus",                        &SAchargeclus,                        &b_SAchargeclus                        );
  fChain->SetBranchAddress("SASPCZ_BI",                       &SASPCZ_BI,                       &b_SASPCZ_BI                       );
  fChain->SetBranchAddress("SASPCZ_BM",                       &SASPCZ_BM,                       &b_SASPCZ_BM                       );
  fChain->SetBranchAddress("SASPCZ_BO",                       &SASPCZ_BO,                       &b_SASPCZ_BO                       );
  fChain->SetBranchAddress("SASPCR_BI",                       &SASPCR_BI,                       &b_SASPCR_BI                       );
  fChain->SetBranchAddress("SASPCR_BM",                       &SASPCR_BM,                       &b_SASPCR_BM                       );
  fChain->SetBranchAddress("SASPCR_BO",                       &SASPCR_BO,                       &b_SASPCR_BO                       );
  fChain->SetBranchAddress("SARPCCluster_zMax",                &SARPCCluster_zMax,                &b_SARPCCluster_zMax                );
  fChain->SetBranchAddress("SARPCCluster_zMin",                &SARPCCluster_zMin,                &b_SARPCCluster_zMin                );
  fChain->SetBranchAddress("SARPCCluster_rMax",                &SARPCCluster_rMax,                &b_SARPCCluster_rMax                );
  fChain->SetBranchAddress("SARPCCluster_rMin",                &SARPCCluster_rMin,                &b_SARPCCluster_rMin                );
  fChain->SetBranchAddress("SARPCCluster_gX",                &SARPCCluster_gX,                &b_SARPCCluster_gX                );
  fChain->SetBranchAddress("SARPCCluster_gY",                &SARPCCluster_gY,                &b_SARPCCluster_gY                );
  fChain->SetBranchAddress("SARPCCluster_gZ",                &SARPCCluster_gZ,                &b_SARPCCluster_gZ                );
  fChain->SetBranchAddress("SARPCCluster_stripWidth",        &SARPCCluster_stripWidth,        &b_SARPCCluster_stripWidth        );
  fChain->SetBranchAddress("SARPCCluster_loX",               &SARPCCluster_loX,               &b_SARPCCluster_loX               );
  fChain->SetBranchAddress("SARPCCluster_loY",               &SARPCCluster_loY,               &b_SARPCCluster_loY               );
  fChain->SetBranchAddress("SARPCCluster_loZ",               &SARPCCluster_loZ,               &b_SARPCCluster_loZ               );
  fChain->SetBranchAddress("SARPCCluster_clusterSize",       &SARPCCluster_clusterSize,       &b_SARPCCluster_clusterSize       );
  fChain->SetBranchAddress("SARPCCluster_clusterLayer",      &SARPCCluster_clusterLayer,      &b_SARPCCluster_clusterLayer      );
  fChain->SetBranchAddress("SARPCCluster_clusterMeasPhi",    &SARPCCluster_clusterMeasPhi,    &b_SARPCCluster_clusterMeasPhi    );
  fChain->SetBranchAddress("SARPCCluster_fitInnPhi",                &SARPCCluster_fitInnPhi,                &b_SARPCCluster_fitInnPhi                );
  fChain->SetBranchAddress("SARPCCluster_fitMidPhi",                &SARPCCluster_fitMidPhi,                &b_SARPCCluster_fitMidPhi                );
  fChain->SetBranchAddress("SARPCCluster_fitOutPhi",                &SARPCCluster_fitOutPhi,                &b_SARPCCluster_fitOutPhi                );
  fChain->SetBranchAddress("SARPCCluster_fitInnSlope",                &SARPCCluster_fitInnSlope,                &b_SARPCCluster_fitInnSlope                );
  fChain->SetBranchAddress("SARPCCluster_fitMidSlope",                &SARPCCluster_fitMidSlope,                &b_SARPCCluster_fitMidSlope                );
  fChain->SetBranchAddress("SARPCCluster_fitOutSlope",                &SARPCCluster_fitOutSlope,                &b_SARPCCluster_fitOutSlope                );
  fChain->SetBranchAddress("SARPCCluster_fitInnOffset",                &SARPCCluster_fitInnOffset,                &b_SARPCCluster_fitInnOffset                );
  fChain->SetBranchAddress("SARPCCluster_fitMidOffset",                &SARPCCluster_fitMidOffset,                &b_SARPCCluster_fitMidOffset                );
  fChain->SetBranchAddress("SARPCCluster_fitOutOffset",                &SARPCCluster_fitOutOffset,                &b_SARPCCluster_fitOutOffset                );
  fChain->SetBranchAddress("SARPCCluster_isSuccess",                &SARPCCluster_isSuccess,                &b_SARPCCluster_isSuccess                );
  fChain->SetBranchAddress("SARPCCluster_isUsingMidCluster",                &SARPCCluster_isUsingMidCluster,                &b_SARPCCluster_isUsingMidCluster                );
  fChain->SetBranchAddress("SARPCCluster_isPlausibleFitInnMid",                &SARPCCluster_isPlausibleFitInnMid,                &b_SARPCCluster_isPlausibleFitInnMid                );
  fChain->SetBranchAddress("SARPCCluster_isPlausibleFitOut",                &SARPCCluster_isPlausibleFitOut,                &b_SARPCCluster_isPlausibleFitOut                );
  fChain->SetBranchAddress("SARPCCluster_isPlausiblePhiInnMid",                &SARPCCluster_isPlausiblePhiInnMid,                &b_SARPCCluster_isPlausiblePhiInnMid                );
  fChain->SetBranchAddress("SARPCCluster_isPlausiblePhiOut",                &SARPCCluster_isPlausiblePhiOut,                &b_SARPCCluster_isPlausiblePhiOut                );
  fChain->SetBranchAddress("SARPCCluster_n_foundClusters",                &SARPCCluster_n_foundClusters,                &b_SARPCCluster_n_foundClusters                );
  fChain->SetBranchAddress("SARPCCluster_id_clustersInSets",                &SARPCCluster_id_clustersInSets,                &b_SARPCCluster_id_clustersInSets                );
  fChain->SetBranchAddress("SARPCCluster_id_clustersInPhiSets",                &SARPCCluster_id_clustersInPhiSets,                &b_SARPCCluster_id_clustersInPhiSets                );

  fChain->SetBranchAddress("SARPCFitInnPhi",                 &SARPCFitInnPhi,                 &b_SARPCFitInnPhi                 );
  fChain->SetBranchAddress("SARPCFitInnSlope",               &SARPCFitInnSlope,               &b_SARPCFitInnSlope               );
  fChain->SetBranchAddress("SARPCFitInnOffset",              &SARPCFitInnOffset,              &b_SARPCFitInnOffset              );
  fChain->SetBranchAddress("SARPCFitMidPhi",                 &SARPCFitMidPhi,                 &b_SARPCFitMidPhi                 );
  fChain->SetBranchAddress("SARPCFitMidSlope",               &SARPCFitMidSlope,               &b_SARPCFitMidSlope               );
  fChain->SetBranchAddress("SARPCFitMidOffset",              &SARPCFitMidOffset,              &b_SARPCFitMidOffset              );
  fChain->SetBranchAddress("SARPCFitOutPhi",                 &SARPCFitOutPhi,                 &b_SARPCFitOutPhi                 );
  fChain->SetBranchAddress("SARPCFitOutSlope",               &SARPCFitOutSlope,               &b_SARPCFitOutSlope               );
  fChain->SetBranchAddress("SARPCFitOutOffset",              &SARPCFitOutOffset,              &b_SARPCFitOutOffset              );
  fChain->SetBranchAddress("SARoadAw",                       &SARoadAw,                       &b_SARoadAw                       );
  fChain->SetBranchAddress("SARoadBw",                       &SARoadBw,                       &b_SARoadBw                       );
  fChain->SetBranchAddress("SAzMin",                         &SAzMin,                         &b_SAzMin                         );
  fChain->SetBranchAddress("SAzMax",                         &SAzMax,                         &b_SAzMax                         );
  fChain->SetBranchAddress("SArMin",                         &SArMin,                         &b_SArMin                         );
  fChain->SetBranchAddress("SArMax",                         &SArMax,                         &b_SArMax                         );
  fChain->SetBranchAddress("SAEtaMin",                       &SAEtaMin,                       &b_SAEtaMin                       );
  fChain->SetBranchAddress("SAEtaMax",                       &SAEtaMax,                       &b_SAEtaMax                       );
  fChain->SetBranchAddress("SAMDTHitisOutlier",              &SAMDTHitisOutlier,              &b_SAMDTHitisOutlier              );
  fChain->SetBranchAddress("SAMDTHitChamber",                &SAMDTHitChamber,                &b_SAMDTHitChamber                );
  fChain->SetBranchAddress("SAMDTHitR",                      &SAMDTHitR,                      &b_SAMDTHitR                      );
  fChain->SetBranchAddress("SAMDTHitZ",                      &SAMDTHitZ,                      &b_SAMDTHitZ                      );
  fChain->SetBranchAddress("SAMDTHitPhi",                    &SAMDTHitPhi,                    &b_SAMDTHitPhi                    );
  fChain->SetBranchAddress("SAMDTHitResidual",               &SAMDTHitResidual,               &b_SAMDTHitResidual               );
  fChain->SetBranchAddress("SAMDTHitSpace",                  &SAMDTHitSpace,                  &b_SAMDTHitSpace                  );
  fChain->SetBranchAddress("SAMDTHitSigma",                  &SAMDTHitSigma,                  &b_SAMDTHitSigma                  );
  fChain->SetBranchAddress("SAMDTHitAllR",                      &SAMDTHitAllR,                      &b_SAMDTHitAllR                      );
  fChain->SetBranchAddress("SAMDTHitAllZ",                      &SAMDTHitAllZ,                      &b_SAMDTHitAllZ                      );
  fChain->SetBranchAddress("SAMDTHitAllclusRoadID",                      &SAMDTHitAllclusRoadID,                      &b_SAMDTHitAllclusRoadID                      );
  fChain->SetBranchAddress("SAMDTHitAllisOutlier",                      &SAMDTHitAllisOutlier,                      &b_SAMDTHitAllisOutlier                      );
  fChain->SetBranchAddress("CBHypoPass",                     &CBHypoPass,                     &b_CBHypoPass                     );
  fChain->SetBranchAddress("CBOvRmPass",                     &CBOvRmPass,                     &b_CBOvRmPass                     );
  fChain->SetBranchAddress("CBPt",                           &CBPt,                           &b_CBPt                           );
  fChain->SetBranchAddress("CBEta",                          &CBEta,                          &b_CBEta                          );
  fChain->SetBranchAddress("EFPass",                     &EFPass,                     &b_EFPass                     );
  InitHist();
  Notify();
}

void RPC_FCBM::InitFCBM(){
  h_distRoIrpcEta  = new TH1F("h_RoIrpcEta", ";|eta_RoI - eta_rpcHit|;events",  150, 0, 0.15);
  h_distRoIrpcPhi  = new TH1F("h_RoIrpcPhi", ";|phi_RoI - phi_rpcHit|;events",  150, 0, 0.15);

  h_CloseByCount = new TH1D("h_CloseByCount", ";;events", 4, 0, 4);
  h_CloseByCount->GetXaxis()->SetBinLabel(1, "#Delta#eta_{RoI-Offline} && #Delta#phi_{RoI-Offline} && same RoI");
  h_CloseByCount->GetXaxis()->SetBinLabel(2, "same RoI");
  h_CloseByCount->GetXaxis()->SetBinLabel(3, "different RoI");
  h_CloseByCount->GetXaxis()->SetBinLabel(4, "single muon");
   
  h_DeltaPt_OffSA_CloseBy = new TH1D("h_CloseByDeltaPt_OffSA", ";p_{t,Offline} - p_{t,SA}[GeV] ;events", 40, -10, 10);
  h_DeltaPt_OffSA = new TH1D("h_noCloseByDeltaPt_OffSA", ";p_{t,Offline} - p_{t,SA}[GeV] ;events", 40, -10, 10);
 
  h_rpcHitSize_CloseBy = new TH1D("h_rpcHitSize_CloseBy", ";number of rpcHitStrip(inclusive);events", 101, -0.5, 100.5); 
  h_rpcHitSize = new TH1D("h_rpcHitSize", ";number of rpcHitStrip(inclusive);events", 101, -0.5, 100.5); 
 
  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcHitSizeLay_CloseBy[iLay] = new TH1D(Form("h_rpcHitSizeLay_CloseBy%d", iLay), Form(";number of rpcHitStrip Layer%d;events", iLay), 51, -0.5, 50.5); 
    h_rpcHitSizeLay[iLay] = new TH1D(Form("h_rpcHitSizeLay%d", iLay), Form(";number of rpcHitStrip Layer%d;events", iLay), 51, -0.5, 50.5); 
  }

  //isSameRpcFitMid
  h_isSameRpcFitMid_CloseBy = new TH1D("h_isSameRpcFitMid_CloseBy", ";isSameRpcFitMid between CloseBy muon;events", 2, 0, 2);
  h_isSameRpcFitMid_CloseBy->GetXaxis()->SetBinLabel(1, "same");
  h_isSameRpcFitMid_CloseBy->GetXaxis()->SetBinLabel(2, "different");
  
  //distance between Offline muon 
  h_distOffEta_CloseBy = new TH1D("h_distOffEta_CloseBy", ";#Delta#eta_{Offline};events", 100, -0.5, 0.5  );
  h_distOffPhi_CloseBy = new TH1D("h_distOffPhi_CloseBy", ";#Delta#phi_{Offline};events", 100, -0.5, 0.5 );
  h_distOffEta = new TH1D("h_distOffEta", ";#Delta#eta_{Offline}, DifferentRoI;events", 100, -0.5, 0.5 );
  h_distOffPhi = new TH1D("h_distOffPhi", ";#Delta#phi_{Offline}, DifferentRoI;events", 100, -0.5, 0.5 );
  
  //distance between rpcFit for each muon
  h_distRpcFitInnEta_CloseBy = new TH1D("h_distRpcFitInnEta_CloseBy", ";#Delta#eta_{rpc1Road};events", 50, -0.5, 0.5 );
  h_distRpcFitInnPhi_CloseBy = new TH1D("h_distRpcFitInnPhi_CloseBy", ";#Delta#phi_{rpc1Road};events", 50, -0.5, 0.5 );
  h_distRpcFitMidEta_CloseBy = new TH1D("h_distRpcFitMidEta_CloseBy", ";#Delta#eta_{rpc2Road};events", 50, -0.5, 0.5 );
  h_distRpcFitMidPhi_CloseBy = new TH1D("h_distRpcFitMidPhi_CloseBy", ";#Delta#phi_{rpc2Road};events", 50, -0.5, 0.5 );
  h_distRpcFitOutEta_CloseBy = new TH1D("h_distRpcFitOutEta_CloseBy", ";#Delta#eta_{rpc3Road};events", 50, -0.5, 0.5 );
  h_distRpcFitOutPhi_CloseBy = new TH1D("h_distRpcFitOutPhi_CloseBy", ";#Delta#phi_{rpc3Road};events", 50, -0.5, 0.5 );
  
  h_distRpcHitToFitInnEta_CloseBy = new TH1D("h_distRpcHitToFitInnEta_CloseBy", ";#eta_{rpc1Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitInnPhi_CloseBy = new TH1D("h_distRpcHitToFitInnPhi_CloseBy", ";#phi_{rpc1Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitMidEta_CloseBy = new TH1D("h_distRpcHitToFitMidEta_CloseBy", ";#eta_{rpc2Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitMidPhi_CloseBy = new TH1D("h_distRpcHitToFitMidPhi_CloseBy", ";#phi_{rpc2Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitOutEta_CloseBy = new TH1D("h_distRpcHitToFitOutEta_CloseBy", ";#eta_{rpc3Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitOutPhi_CloseBy = new TH1D("h_distRpcHitToFitOutPhi_CloseBy", ";#phi_{rpc3Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
 /* 
  h_widedistRpcHitToFitInnEta_CloseBy = new TH1D("h_widedistRpcHitToFitInnEta_CloseBy", ";#eta_{rpc1Road}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitInnPhi_CloseBy = new TH1D("h_widedistRpcHitToFitInnPhi_CloseBy", ";#phi_{rpc1Road}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  h_widedistRpcHitToFitMidEta_CloseBy = new TH1D("h_widedistRpcHitToFitMidEta_CloseBy", ";#eta_{rpc2Road}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitMidPhi_CloseBy = new TH1D("h_widedistRpcHitToFitMidPhi_CloseBy", ";#phi_{rpc2Road}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  h_widedistRpcHitToFitOutEta_CloseBy = new TH1D("h_widedistRpcHitToFitOutEta_CloseBy", ";#eta_{rpc3Road}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitOutPhi_CloseBy = new TH1D("h_widedistRpcHitToFitOutPhi_CloseBy", ";#phi_{rpc3Road}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  */
  h_distRpcFitInnEta = new TH1D("h_distRpcFitInnEta", ";#Delta#eta_{rpc1Road};events", 88, -1.1, 1.1 );
  h_distRpcFitInnPhi = new TH1D("h_distRpcFitInnPhi", ";#Delta#phi_{rpc1Road};events", 128, -3.2, 3.2 );
  h_distRpcFitMidEta = new TH1D("h_distRpcFitMidEta", ";#Delta#eta_{rpc2Road};events", 88, -1.1, 1.1 );
  h_distRpcFitMidPhi = new TH1D("h_distRpcFitMidPhi", ";#Delta#phi_{rpc2Road};events", 128, -3.2, 3.2 );
  h_distRpcFitOutEta = new TH1D("h_distRpcFitOutEta", ";#Delta#eta_{rpc3Road};events", 88, -1.1, 1.1 );
  h_distRpcFitOutPhi = new TH1D("h_distRpcFitOutPhi", ";#Delta#phi_{rpc3Road};events", 128, -3.2, 3.2 );

  h_distRpcHitToFitInnEta = new TH1D("h_distRpcHitToFitInnEta", ";#eta_{rpc1Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitInnPhi = new TH1D("h_distRpcHitToFitInnPhi", ";#phi_{rpc1Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitMidEta = new TH1D("h_distRpcHitToFitMidEta", ";#eta_{rpc2Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitMidPhi = new TH1D("h_distRpcHitToFitMidPhi", ";#phi_{rpc2Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitOutEta = new TH1D("h_distRpcHitToFitOutEta", ";#eta_{rpc3Road}-#eta_{rpcHitStrip};events", 50, -0.1, 0.1 );
  h_distRpcHitToFitOutPhi = new TH1D("h_distRpcHitToFitOutPhi", ";#phi_{rpc3Road}-#phi_{rpcHitStrip};events", 50, -0.1, 0.1 );
  
  h_DeltaRpcMidMinEta_CloseBy = new TH1D("h_DeltaRpcMidabs_minEta_CloseBy", ";abs_min(#eta_{rpcRoad}-#eta_{rpcHitStrip});events", 80, 0, 0.2 );
  h_DeltaRpcMidMinPhi_CloseBy = new TH1D("h_DeltaRpcMidabs_minPhi_CloseBy", ";abs_min(#phi_{rpcRoad}-#phi_{rpcHitStrip});events", 80, 0, 0.2 );
  h_DeltaRpcOutMinEta_CloseBy = new TH1D("h_DeltaRpcOutabs_minEta_CloseBy", ";abs_min(#eta_{rpcRoad}-#eta_{rpcHitStrip});events", 80, 0, 0.2 );
  h_DeltaRpcOutMinPhi_CloseBy = new TH1D("h_DeltaRpcOutabs_minPhi_CloseBy", ";abs_min(#phi_{rpcRoad}-#phi_{rpcHitStrip});events", 80, 0, 0.2 );
  
  /*
  h_widedistRpcHitToFitInnEta = new TH1D("h_widedistRpcHitToFitInnEta", ";#eta_{rpcFitInn}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitInnPhi = new TH1D("h_widedistRpcHitToFitInnPhi", ";#phi_{rpcFitInn}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  h_widedistRpcHitToFitMidEta = new TH1D("h_widedistRpcHitToFitMidEta", ";#eta_{rpcFitMid}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitMidPhi = new TH1D("h_widedistRpcHitToFitMidPhi", ";#phi_{rpcFitMid}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  h_widedistRpcHitToFitOutEta = new TH1D("h_widedistRpcHitToFitOutEta", ";#eta_{rpcFitOut}-#eta_{rpcHitStrip};events", 100, -0.1, 0.1 );
  h_widedistRpcHitToFitOutPhi = new TH1D("h_widedistRpcHitToFitOutPhi", ";#phi_{rpcFitOut}-#phi_{rpcHitStrip};events", 100, -0.1, 3.2 );
  */

  //RPC Clustering
  h_rpc2ClusEtaNum_CloseBy = new TH1D("h_ClusterEtaNum_CloseBy", "nClusters(Eta) at RPC2;N_{Clusters at RPC2};events",21, -0.5, 20.5);
  h_rpc2ClusPhiNum_CloseBy = new TH1D("h_ClusterPhiNum_CloseBy", "nClusters(Phi) at RPC2;N_{Clusters at RPC2};events",21, -0.5, 20.5);
  h_rpcClusEta_CloseBy  = new TH1D("h_rpcClusEta_CloseBy", "SA rpcClusterEta_CloseBy;Eta_CloseBy;Entries",  212, -1.06, 1.06);
  h_rpcClusPhi_CloseBy  = new TH1D("h_rpcClusPhi_CloseBy", "SA rpcClusterPhi_CloseBy;Phi_CloseBy;Entries", 64 , -3.2, 3.2);
  h_distRoIClusR_CloseBy = new TH1D("h_distRoIClusR_CloseBy", "#Delta R(RoI center - cluster center) at RPC2;#Delta R;events", 50, 0, 0.5);
  h_distRoIClusR_noCloseBy = new TH1D("h_distRoIClusR_noCloseBy", "#Delta R(RoI center - cluster center) at RPC2;#Delta R;events", 50, 0, 0.5);
  
  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcClusEtaSizeLay_CloseBy[iLay] = new TH1D(Form("h_rpcClusEtaSizeLay_CloseBy%d", iLay), Form("number of strip clusters at Layer%d (CloseBy);N_{cluster};events", iLay+1), 21, -0.5, 20.5); 
    h_rpcClusPhiSizeLay_CloseBy[iLay] = new TH1D(Form("h_rpcClusPhiSizeLay_CloseBy%d", iLay), Form("number of strip clusters at Layer%d (CloseBy);N_{cluster};events", iLay+1), 21, -0.5, 20.5); 
  }
}

void RPC_FCBM::InitNoCut(){
  h_SAEtaPhi = new TH2D("h_SAEtaPhi", "SA Eta vs Phi;Eta;Phi", 212, -1.06, 1.06, 64, -3.2, 3.2);
  h_SARoIEtaPhi = new TH2D("h_SARoIEtaPhi", "SA RoIEta vs RoIPhi;Eta;Phi", 30, -1.5, 1.5, 62, -3.2, 3.2);
  h_OffSegZR = new TH2D("h_OffSegZR", "OfflineSegment Z vs R;Z[mm];R[mm]", 300, 15000, 15000, 300, -15000, 15000);
  h_OffSegEtaPhi = new TH2D("h_OffSegEtaPhi", "OfflineSegment eta vs phi;eta;phi", 40, -2, 2, 70, -3.5, 3.5);
  for(int iLay = 0; iLay < 8; iLay++){
    h_OffSegXY_etaLayer[iLay] = new TH2D(Form("h_OffSegXY_etaLayer%d", iLay), Form("OfflineSegment position;X[m];Y[m] (NTrigEtaLayers=%d)", iLay), 300, -15, 15, 300, -15, 15); 
  }
  h_OffSegXY = new TH2D("h_OffSegXY", "OfflineSegment position;X[m];Y[m]", 300, -15, 15, 300, -15, 15); 
  h_nOffSegMid = new TH1D("h_nOffSegMid", "the number of OfflineSegments at middle station (for each Offline muon index);N_{OffsegMid};events",8, -0.5, 7.5);
  h_nOffSegOut = new TH1D("h_nOffSegOut", "the number of OfflineSegments at Outer station (for each Offline muon index);N_{OffsegOut};events",8, -0.5, 7.5);
  h_SARpcHitEta  = new TH1D("h_SARpcHitEta", "SA rpcHitEta;Eta;Entries",  212, -1.06, 1.06);
  h_SARpcHitPhi  = new TH1D("h_SARpcHitPhi", "SA rpcHitPhi;Phi;Entries", 64 , -3.2, 3.2);
  h_distnoCutRoIrpcEta  = new TH1F("h_noCutRoIrpcEta", "dist RoI-rpcHitEta noCut;|delta_Eta| ;Entries",  105, 0, 1.05);
  h_distnoCutRoIrpcPhi  = new TH1F("h_noCutRoIrpcPhi", "dist RoI-rpcHitPhi noCut;|delta_Phi|;Entries",  320, 0, 3.2);
  h_distnoCutEachRoIrpcEta = new TH1D("h_EachRoIrpcEta", "dist RoI-rpcHitEta noCut;|delta_Eta| ;Entries",  105, 0, 1.05);
  h_distnoCutEachRoIrpcPhi = new TH1D("h_EachRoIrpcPhi", "dist EachRoI-rpcHitPhi noCut;|delta_Phi|;Entries",  320, 0, 3.2);
  
  
  h_isCloseBydEdP = new TH1D("h_isCloseBydEdP", "#delta#eta#delta#phi;#delta#eta_{SARoI-Offline} < 0.1 && #delta#phi_{SARoI-Offline} < 0.1;events", 2, 0, 2);
  h_isCloseBydEdP->GetXaxis()->SetBinLabel(1, "true");
  h_isCloseBydEdP->GetXaxis()->SetBinLabel(2, "false");
  
  h_isSameRoI = new TH1D("h_isSameRoI", "isSameRoI;isSameRoI between each muon;events", 3, 0, 3);
  h_isSameRoI->GetXaxis()->SetBinLabel(1, "same");
  h_isSameRoI->GetXaxis()->SetBinLabel(2, "different");
  h_isSameRoI->GetXaxis()->SetBinLabel(3, "single muon");
    
  h_SARPCFitInnR = new TH1D("h_SARPCFitInnR", ";RPC1 Fit R[mm];events", 100, 6000, 9000);
  h_SARPCFitInnEta = new TH1D("h_SARPCFitInnEta", ";RPC1 Fit #eta;events", 80, -2, 2);
  h_SARPCFitMidR = new TH1D("h_SARPCFitMidR", ";RPC2 Fit R[mm];events", 100, 6000, 9000);
  h_SARPCFitMidEta = new TH1D("h_SARPCFitMidEta", ";RPC2 Fit #eta;events", 88, -1.1, 1.1);
  h_SARPCFitOutR = new TH1D("h_SARPCFitOutR", ";RPC3 Fit R[mm];events", 100, 6000, 9000);
  h_SARPCFitOutEta = new TH1D("h_SARPCFitOutEta", ";RPC3 Fit #eta;events", 88, -1.1, 1.1);
  h_isSameRpcFitMid = new TH1D("h_isSameRpcFitMid", ";isSameRpcFitMid between each muon;events", 3, 0, 3);
  h_isSameRpcFitMid->GetXaxis()->SetBinLabel(1, "same");
  h_isSameRpcFitMid->GetXaxis()->SetBinLabel(2, "different");
  h_isSameRpcFitMid->GetXaxis()->SetBinLabel(3, "single muon");
  h_extdR_bug = new TH1D("h_extdR_bug", "#DeltaR_{ext} (road algorithm bug event);#DeltaR_{ext}; events", 50, 0, 0.25);
  h_RoIEtaPhi_bug = new TH2D("h_RoIEtaPhi_bug", "#eta_{RoI} vs #phi_{RoI} (road algorithm bug event);#eta_{RoI}; #phi_{RoI}", 22, -1.1, 1.1, 64, -3.2, 3.2);
  h_extdR_wobug = new TH1D("h_extdR_wobug", "#DeltaR_{ext} (correct road algorithm event);#DeltaR_{ext}; events", 50, 0, 0.25);
  h_RoIEtaPhi_wobug = new TH2D("h_RoIEtaPhi_wobug", "#eta_{RoI} vs #phi_{RoI} (correct road algorithm event);#eta_{RoI}; #phi_{RoI}", 22, -1.1, 1.1, 64, -3.2, 3.2);

  h_RoInumDelEta = new TH2D("h_RoInumDelEta", ";#Delta#eta_{offline};RoI state;events", 80, -0.2, 0.2, 2, 0, 2 );
  h_RoInumDelEta->GetYaxis()->SetBinLabel(1, "same");
  h_RoInumDelEta->GetYaxis()->SetBinLabel(2, "different");
  h_RoInumDelPhi = new TH2D("h_RoInumDelPhi", ";#Delta#phi_{offline};RoI state;events", 80, -0.2, 0.2, 2, 0, 2 );
  h_RoInumDelPhi->GetYaxis()->SetBinLabel(1, "same");
  h_RoInumDelPhi->GetYaxis()->SetBinLabel(2, "different");
  for(int iLay = 0; iLay < 8; iLay++){
    h_rpcHitXYLay[iLay] = new TH2D(Form("h_rpcHitXYLay%d", iLay), Form("rpcHit strip position;X[m];Y[m]", iLay), 300, -15, 15, 300, -15, 15); 
    //h_rpcHitXYLay[iLay] = new TH2D(Form("h_rpcHitXYLay%d", iLay), Form("rpcHit strip position(layer%d);X[m];Y[m]", iLay), 300, -15, 15, 300, -15, 15); 
  }

  h_nClusRoad_withMidInfo = new TH1D("h_nClusRoad_withMidInfo", "The number of clusterRoads(middle);N_{clusterRoad(middle)};events", 11, -0.5, 10.5);
  h_nClusRoad_withOutInfo = new TH1D("h_nClusRoad_withOutInfo", "The number of clusterRoads(middle);N_{clusterRoad(middle)};events", 11, -0.5, 10.5);
  
  h_nClusRoadOut_woOutInfo = new TH1D("h_nClusRoadOut_woOutInfo", "The number of clusterRoads(outer) wo outerCluster;N_{clusterRoad(outer)};events", 11, -0.5, 10.5);
  h_nClusRoadOut_byMiddle = new TH1D("h_nClusRoadOut_byMiddle", "The number of clusterRoads(outer) by RPC1 cluster;N_{clusterRoad(outer)};events", 11, -0.5, 10.5);
  h_nClusRoadOut_byOuter = new TH1D("h_nClusRoadOut_byOuter", "The number of clusterRoads(outer);N_{clusterRoad(outer)};events", 11, -0.5, 10.5);
  
  h_nClusRoadMid_byMiddleTest = new TH1D("h_nClusRoadMid_byMiddle", "The number of clusterRoads(Middle) by RPC1 cluster;N_{clusterRoad(Middle)};events", 11, -0.5, 10.5);
  h_nClusRoadMid_byOuterTest = new TH1D("h_nClusRoadMid_byOuter", "The number of clusterRoads(Middle);N_{clusterRoad(Middle)};events", 11, -0.5, 10.5);

  h_distMidOffClusEta = new TH1D("h_distMidOffClusEta", "d=fabs(OffsegMidEta-ClusterEta) at Middle;d;events", 105, 0, 1.05);
  h_distMidOffClusPhi = new TH1D("h_distMidOffClusPhi", "d=fabs(OffsegMidPhi-ClusterPhi) at Middle;d;events", 50, 0, 0.5);
  
  h_distMidOffClusR = new TH1D("h_distMidOffClusR", "#Delta R(OffsegMid - cluster center);#Delta R;events", 150, 0, 1.5);
  h_distMidOffClus = new TH2D("h_distMidOffClus", "distance(OffsegMid-Cluster) at Middle;#Delta#eta;#Delta#phi;events", 80, -0.4, 0.4, 80, -0.4, 0.4);
  h_etamin_MidOffFit_wMidInfo = new TH1D("h_etamin_MidOffFit_wMidInfo", "|#Delta#eta|_{min} : OffsegMid - clusterRoadMid;|#Delta#eta|_{min};events", 100, 0, 0.1);
  h_etamin_MidOffFit_woMidInfo = new TH1D("h_etamin_MidOffFit_woMidInfo", "|#Delta#eta|_{min} : OffsegMid - clusterRoadMid;|#Delta#eta|_{min};events", 100, 0, 0.1);
  h_slopemin_MidOffFit_wMidInfo = new TH1D("h_slopemin_MidOffFit_wMidInfo", "|#Delta(R/Z)|_{min} : OffsegMid - clusterRoadMid;|#Delta(R/Z)|_{min};events", 60, 0, 0.3);
  h_slopemin_MidOffFit_woMidInfo = new TH1D("h_slopemin_MidOffFit_woMidInfo", "|#Delta(R/Z)|_{min} : OffsegMid - clusterRoadMid;|#Delta(R/Z)|_{min};events", 80, 0, 0.4);
  
  h_etamin_OutOffFit = new TH1D("h_etamin_OutOffFit", "|#Delta#eta|_{min} : OffsegOut - clusterRoadOut;|#Delta#eta|_{min};events", 100, 0, 0.1);
  h_slopemin_OutOffFit = new TH1D("h_slopemin_OutOffFit", "|#Delta(R/Z)|_{min} : OffsegOut - clusterRoadOut;|#Delta(R/Z)|_{min};events", 60, 0, 0.3);
  
  h_etamin_MidClusOutInfo_OffsegOut = new TH1D("h_etamin_MidClusOutInfo_OffsegOut", "|#Delta#eta|_{min} : OffsegOut - clusRoadMid(by outer);|#Delta#eta|_{min};events", 100, 0, 0.1);
  h_slopemin_MidClusOutInfo_OffsegOut = new TH1D("h_slopemin_MidClusOutInfo_OffsegOut", "|#Delta(R/Z)|_{min} : OffsegOut - clusRoadMid(by outer);|#Delta(R/Z)|_{min};events", 80, 0, 0.4);
  
  
  h_etamin_MidOffFit_strip = new TH1D("h_etamin_MidOffFit_strip", "|#Delta#eta|_{min} : OffsegMid - stripRoadMid;|#Delta#eta|_{min};events", 100, 0, 0.1);
  h_slopemin_MidOffFit_strip = new TH1D("h_slopemin_MidOffFit_strip", "|#Delta(R/Z)|_{min} : OffsegMid - stripRoadMid;|#Delta(R/Z)|_{min};events", 60, 0, 0.3);
  /*for(int i = 0; i< 10; i++){
    h_fakeRoadMid[i] = new TH1D(Form("h_fakeRoadMid%d", i), Form("N_{fakeRoad} at Middle fake:#Delta#eta>0.01&&#Delta(R/Z)>%f;N_{fakeRoad};events", i*0.1), 13, -0.5, 12.5);

    h_truthRoadMid[i] = new TH1D(Form("h_truthRoadMid%d", i), Form("N_{truthlike-Road} at Middle|truthLike:#Delta#eta<=0.01&&#Delta(R/Z)<=%f ;N_{truthlike-Road};events", i*0.1), 13, -0.5, 12.5);
  }*/
  //EtaMin PhiMin
  h_nClusInSet = new TH1D("h_nClusInSet", "N_{clusters} in a clusterSet;N_{nClusInSet};events", 10, -0.5, 9.5);
  h_nClusFitMid = new TH1D("h_nClusFitMid", "N_{nClusFit} at Middle;N_{nClusFit};events", 13, -0.5, 12.5);
  h_OffsegEta = new TH1D("h_OffsegEta", "Offline Segment #eta;#eta_{offseg};events", 88, -1.1, -1.1);
  h_ClusFitMidEta = new TH1D("h_ClusFitMidEta", "Road(middle) eta;#eta_{Road};events", 88, -1.1, 1.1 );
  h_ClusFitOutEta = new TH1D("h_ClusFitOutEta", "Road(Outer) eta;#eta_{Road};events", 88, -1.1, 1.1 );
  h_fakeRoadMid = new TH1D("h_fakeRoadMid", "N_{fakeRoad} at Middle [fake:#Delta#eta>0.01 && #Delta(R/Z)>0.08];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadMid = new TH1D("h_truthRoadMid", "N_{true-likeRoad} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.08] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_faketruthRoadMid = new TH2D("h_faketruthRoadMid", "N_{fakeRoad} vs N_{true-like Road} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.08];N_{fake};N_{true-like};events", 13, -0.5, 12.5, 13, -0.5, 12.5);

//outer!!!!!!
  h_fakeRoadOut = new TH1D("h_fakeRoadOut", "N_{fakeRoad} at Outer [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadOut = new TH1D("h_truthRoadOut", "N_{true-likeRoad} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_faketruthRoadOut = new TH2D("h_faketruthRoadOut", "N_{fakeRoad} vs N_{true-like Road} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05];N_{fake};N_{true-like};events", 13, -0.5, 12.5, 13, -0.5, 12.5);

  //Middle Outer fake True-like histograms added options
  h_fakeRoadMid_opt = new TH1D("h_fakeRoadMid_opt", "N_{fakeRoad} at Middle [fake:#Delta#eta>0.01 && #Delta(R/Z)>0.08];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadMid_opt = new TH1D("h_truthRoadMid_opt", "N_{true-likeRoad} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.08] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadMid_opt23 = new TH1D("h_fakeRoadMid_opt23", "N_{fakeRoad} at Middle [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadMid_opt23 = new TH1D("h_truthRoadMid_opt23", "N_{true-likeRoad} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadMid_opt13 = new TH1D("h_fakeRoadMid_opt13", "N_{fakeRoad} at Middle [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadMid_opt13 = new TH1D("h_truthRoadMid_opt13", "N_{true-likeRoad} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadMid_opt12 = new TH1D("h_fakeRoadMid_opt12", "N_{fakeRoad} at Middle [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadMid_opt12 = new TH1D("h_truthRoadMid_opt12", "N_{true-likeRoad} at Middle [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadOut_opt = new TH1D("h_fakeRoadOut_opt", "N_{fakeRoad} at Outer [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadOut_opt = new TH1D("h_truthRoadOut_opt", "N_{true-likeRoad} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadOut_opt23 = new TH1D("h_fakeRoadOut_opt23", "N_{fakeRoad} at Outer [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadOut_opt23 = new TH1D("h_truthRoadOut_opt23", "N_{true-likeRoad} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadOut_opt13 = new TH1D("h_fakeRoadOut_opt13", "N_{fakeRoad} at Outer [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadOut_opt13 = new TH1D("h_truthRoadOut_opt13", "N_{true-likeRoad} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);
  h_fakeRoadOut_opt12 = new TH1D("h_fakeRoadOut_opt12", "N_{fakeRoad} at Outer [fake:#Delta#eta>0.01 && #Delta(R/Z) > 0.05];N_{fakeRoad};events", 13, -0.5, 12.5);
  h_truthRoadOut_opt12 = new TH1D("h_truthRoadOut_opt12", "N_{true-likeRoad} at Outer [true-like:#Delta#eta<=0.01 || #Delta(R/Z)<=0.05] ;N_{true-likeRoad};events", 13, -0.5, 12.5);

  h_fitMidSlopeOffset_cut = new TH2D("h_fitMidSlopeOffset_cut", "clusterRoad slope vs offset (middle): cut the RPC1both or RPC2both only; Offset R[m] ; Slope(R/Z)", 40, -8, 8, 50, -5, 5);
  h_fitMidSlopeOffset_wocut = new TH2D("h_fitMidSlopeOffset_wocut", "clusterRoad slope vs offset (middle): wocut the RPC1both or RPC2both only; Offset R[m] ; Slope(R/Z)", 40, -8, 8, 50, -5, 5);

  m_h_countPtclus = new TH1D("m_h_countPtclus", "the number of pT from clusterRoad; N_{pT}; events", 13, -0.5, 12.5);
  h_nmaxClustersInSet = new TH1D("h_nmaxClustersInSet", "N_{max} in cluster sets; N_{max}; events", 8, 1.5, 9.5);
  h_countSPInn = new TH1D("h_countSPInn", "the number of SP from clusterRoad; N_{SP}; events", 13, -0.5, 12.5);
  h_countSPMid = new TH1D("h_countSPMid", "the number of SP from clusterRoad; N_{SP}; events", 13, -0.5, 12.5);
  h_countSPOut = new TH1D("h_countSPOut", "the number of SP from clusterRoad; N_{SP}; events", 13, -0.5, 12.5);
  h_countClusterRoadMid = new TH1D("h_countClusterRoadMid", "the number of clusterRoad; N_{clusterRoadMid}; events", 13, -0.5, 12.5);
  h_countClusterRoadOut = new TH1D("h_countClusterRoadOut", "the number of clusterRoad; N_{clusterRoad}; events", 13, -0.5, 12.5);
  hh_num_clusRoadMidvsOut = new TH2D("h_num_clusRoadMidvsOut", "the number of clusterRoad; N_{middle}; N_{outer}", 6, -0.5, 5.5, 6, -0.5, 5.5);
  hh_n_CRvsSPCMid = new TH2D("h_n_CRvsSPCMid", "N_{superpoint} vs N_{clusterRoad} at middle (2mu-in-1RoI evt); N_{clusterRoad}; N_{superpoint}", 6, -0.5, 5.5, 6, -0.5, 5.5);
  prof_CRvsSPCMid = new TProfile();
  
  m_h_n_lowptclus_patA = new TH1D("m_h_n_lowptclus_patA", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 1, 0, 1);
  m_h_n_lowptclus_patA->GetXaxis()->SetBinLabel(1, "all stats have 2 or more SPs");
  m_h_n_lowptclus_patA_sameIDout = new TH1D("m_h_n_lowptclus_patA_sameIDout", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 1, 0, 1);
  m_h_n_lowptclus_patA_sameIDout->GetXaxis()->SetBinLabel(1, "all stats have 2 or more SPs");
  m_h_n_lowptclus_patA_findout = new TH1D("m_h_n_lowptclus_patA_findout", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 1, 0, 1);
  m_h_n_lowptclus_patA_findout->GetXaxis()->SetBinLabel(1, "all stats have 2 or more SPs");
  
  
  m_h_n_lowptclus = new TH1D("m_h_n_lowptclus", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 3, 0, 3);
  m_h_n_lowptclus->GetXaxis()->SetBinLabel(1, "A");
  m_h_n_lowptclus->GetXaxis()->SetBinLabel(2, "B");
  m_h_n_lowptclus->GetXaxis()->SetBinLabel(3, "C");
  m_h_n_1spstat_inn = new TH1D("m_h_n_1spstat_inn", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 3, 0, 3);
  m_h_n_1spstat_inn->GetXaxis()->SetBinLabel(1, "A");
  m_h_n_1spstat_inn->GetXaxis()->SetBinLabel(2, "B");
  m_h_n_1spstat_inn->GetXaxis()->SetBinLabel(3, "C");
  m_h_n_1spstat_mid = new TH1D("m_h_n_1spstat_mid", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 3, 0, 3);
  m_h_n_1spstat_mid->GetXaxis()->SetBinLabel(1, "A");
  m_h_n_1spstat_mid->GetXaxis()->SetBinLabel(2, "B");
  m_h_n_1spstat_mid->GetXaxis()->SetBinLabel(3, "C");
  m_h_n_1spstat_out = new TH1D("m_h_n_1spstat_out", "N_{superpoint} condition when much lower p_{T} was calculated; condition; events", 3, 0, 3);
  m_h_n_1spstat_out->GetXaxis()->SetBinLabel(1, "A");
  m_h_n_1spstat_out->GetXaxis()->SetBinLabel(2, "B");
  m_h_n_1spstat_out->GetXaxis()->SetBinLabel(3, "C");
  

  m_h_dSlope_innmid = new TH1D("m_h_dSlope_innmid", "superpoint #Delta(Z/R)_{min}: inner-middle;#Delta(Z/R);events", 300, 0, 0.03);
  
  h_iswoRpc1Road = new TH1D("h_iswoRpc1Road", "the number of some patterns' sets (N_{max} set);events", 4, 0, 4);
  h_iswoRpc1Road->GetXaxis()->SetBinLabel(1, "(!RPC1)&&(RPC2)");
  h_iswoRpc1Road->GetXaxis()->SetBinLabel(2, "(RPC1)&&(!RPC2)");
  h_iswoRpc1Road->GetXaxis()->SetBinLabel(3, "(!RPC3)");
  h_iswoRpc1Road->GetXaxis()->SetBinLabel(4, "else");

  h_id_minlayer = new TH1D("h_id_minlayer", "Index of cluster in seed layer; ID; evets", 10, -0.5, 9.5);
  h_id_minlayer0 = new TH1D("h_id_minlayer0", "Index of cluster in seed layer; ID; evets", 10, -0.5, 9.5);
  
  h_seedcluster_setmax = new TH1D("h_seedcluster_setmax", "ID of seed cluster in Set_{max} (seed layer > 0); id; evets", 10, -0.5, 9.5);
  h_minlayerInSetmax = new TH1D("h_minlayerInSetmax", "Seed layer of Set_{max}; layer; evets", 10, -0.5, 9.5);

  h_distMidOffClusEtaMin = new TH1D("h_distMidOffClusEtaMin", "d=min_fabs(OffsegMidEtaMin-ClusterEta) at Middle;d;events", 100, 0, 0.1);
  h_distMidOffClusPhiMin = new TH1D("h_distMidOffClusPhiMin", "d=min_fabs(OffsegMidPhiMin-ClusterPhi) at Middle;d;events", 100, 0, 0.1);
  h_distMidOffClusMin = new TH2D("h_distMidOffClusMin", "Min_distance(OffsegMid-Cluster) at Middle;#Delta#eta_{min};#Delta#phi_{min};events", 80, -0.2, 0.2, 80, -0.2, 0.2);
  h_distMidOffClusfitEtaMin = new TH1D("h_distMidOffClusfitEtaMin", "d=min_fabs(#eta_{Offseg}-#eta_{clusterRoad}) at Middle;d;events", 100, 0, 0.1);
  h_distRoIClusR = new TH1D("h_distRoIClusR", "#Delta R(RoI center - cluster center) at RPC2;#Delta R;events", 50, 0, 0.5);
  h_distRoIClus = new TH2D("h_distRoIClus", "distance(RoI center-Cluster) at RPC2;#Delta#eta;#Delta#phi;events", 80, -0.4, 0.4, 80, -0.4, 0.4);
  h_nClusLay2Eta = new TH1D("h_nClusLay2Eta", "nClusters at Layer2 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay2Phi = new TH1D("h_nClusLay2Phi", "nClusters at Layer2 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay3Eta = new TH1D("h_nClusLay3Eta", "nClusters at Layer3 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay3Phi = new TH1D("h_nClusLay3Phi", "nClusters at Layer3 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay23Eta = new TH1D("h_nClusLay23Eta", "nClusters at Layer2+3 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay23Phi = new TH1D("h_nClusLay23Phi", "nClusters at Layer2+3 (#Delta R < 0.1 clus center-RoI);N_{clusters};events",21, -0.5, 20.5);
  h_nClusLay23EtaPhi = new TH1D("h_nClusLay23EtaPhi", "nClusters at Layer2+3 (#Delta R < 0.1 clus center-RoI);N_{#eta + #phi clusters};events",21, -0.5, 20.5);
  h_nClusHit23Eta = new TH1D("h_nClusHit23Eta", "nHits at Layer2+3 (#Delta R < 0.1 rpcHit center-RoI);N_{hits};events",21, -0.5, 20.5);
  h_nClusHit23Phi = new TH1D("h_nClusHit23Phi", "nHits at Layer2+3 (#Delta R < 0.1 rpcHit center-RoI);N_{hits};events",21, -0.5, 20.5);
  h_clusterSizeLayOut = new TH1D("h_clusterSizeLayOut", "nStrips in cluster;N_{strips};events",21, -0.5, 20.5);
  h_clusterSizeLay01 = new TH1D("h_clusterSizeLay01", "nStrips in cluster || 1;N_{strips};events",21, -0.5, 20.5);
  h_clusterSizeLay23 = new TH1D("h_clusterSizeLay23", "nStrips in cluster at Layer2 || 3;N_{strips};events",21, -0.5, 20.5);
  h_clusterSize = new TH1D("h_clusterSize", "nStrips in cluster at all layer;N_{strips};events",21, -0.5, 20.5);
  m_nclus_vs_npt = new TH2D("m_nclus_vs_npt", ";(N_{RPC #eta cluster}/Layer)_{max};N_{SA muons};events", 15, -0.5, 14.5, 12, -0.5, 11.5);
  m_ptres_default = new TH1D("m_ptres_default", ";p_{T} residual;events", 80, -1, 1);
  m_ptres_mpsa = new TH1D("m_ptres_mpsa", ";p_{T} residual;events", 80, -1, 1);

  h_isSameSize_rpcHitEta = new TH1D("h_isSameSize_rpcHitEta", "rpcHitEta isSameSize for each Off-muon;;events", 2, 0, 2);
  h_isSameSize_rpcHitEta->GetXaxis()->SetBinLabel(1, "same size");
  h_isSameSize_rpcHitEta->GetXaxis()->SetBinLabel(2, "different size");

  h_isSameSize_rpcHitPhi = new TH1D("h_isSameSize_rpcHitPhi", "rpcHitPhi isSameSize for each Off-muon;;events", 2, 0, 2);
  h_isSameSize_rpcHitPhi->GetXaxis()->SetBinLabel(1, "same size");
  h_isSameSize_rpcHitPhi->GetXaxis()->SetBinLabel(2, "different size");

  h_Zmin_OffClusRoadMid = new TH1D("h_Zmin_OffClusRoadMid", "|#DeltaZ| OffsegMid vs clusterRoadMid;|#DeltaZ|[m];events", 100, 0, 2);
  h_Zsubmin_OffClusRoadMid = new TH1D("h_Zsubmin_OffClusRoadMid", "|#DeltaZ| OffsegMid vs clusterRoadMid;|#DeltaZ|[m];events", 100, 0, 2);
  h_Z_plaumin_OffClusRoadMid = new TH1D("h_Z_plaumin_OffClusRoadMid", "|#DeltaZ| OffsegMid vs clusterRoadMid w/ similarRoad cut;|#DeltaZ|[m];events", 100, 0, 2);
  h_Z_plausubmin_OffClusRoadMid = new TH1D("h_Z_plausubmin_OffClusRoadMid", "|#DeltaZ| OffsegMid vs clusterRoadMid w/ similarRoad cut;|#DeltaZ|[m];events", 100, 0, 2);
  h_n_clusRoadMid_total = new TH1D("h_n_clusRoadMid_total", "N_{clusterRoadMid};N_{clusterRoad};events", 13, -0.5, 12.5);
  h_n_clusRoadMid_plau = new TH1D("h_n_clusRoadMid_plau", "N_{clusterRoadMid};N_{clusterRoad};events", 13, -0.5, 12.5);
  h_n_clusRoadMid_plau_n = new TH1D("h_n_clusRoadMid_plau_n", "N_{clusterRoadMid};N_{clusterRoad};events", 13, -0.5, 12.5);
  h_deltaZR_min_samePlane = new TH1D("h_deltaZR_min_samePlane", "|#Delta(Z/R)|_{min} between samePlane's cluster; |#Delta(Z/R)|_{min}x10^{-3}; events", 40, 0, 20);

//Superpoint efficiency

  m_h_deltaExtR_2spall = new TH1D("m_h_deltaExtR_2spall", "#DeltaR_{offline} at middle station (Nspall == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_2spinn = new TH1D("m_h_deltaExtR_2spinn", "#DeltaR_{offline} at middle station (Nspinn == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_2spmid = new TH1D("m_h_deltaExtR_2spmid", "#DeltaR_{offline} at middle station (Nspmid == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_2spout = new TH1D("m_h_deltaExtR_2spout", "#DeltaR_{offline} at middle station (Nspout == 2); #DeltaR;events", 10, 0, 0.1);

  m_h_deltaExtR_3pt = new TH1D("m_h_deltaExtR_3pt", "#DeltaR_{offline} at middle station (Npt > 2); #DeltaR;events", 10, 0, 0.1);
//  m_h_deltaExtR_2pt = new TH1D("m_h_deltaExtR_2pt", "#DeltaR_{offline} at middle station (Npt == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_2pt = new TH1D("m_h_deltaExtR_2pt", "#DeltaR_{offline} at middle station (Npt == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_1pt = new TH1D("m_h_deltaExtR_1pt", "#DeltaR_{offline} at middle station (Npt < 2); #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_2pt_oppcharge = new TH1D("m_h_deltaExtR_2pt_oppcharge", "#DeltaR_{offline} at middle station (Npt == 2); #DeltaR;events", 10, 0, 0.1);
  m_h_ptres_2off2clus = new TH1D("m_h_ptres_2off2clus", "; p_{T} residual; events", 200, -1, 1);
//  m_h_ptres_2off2clus->Sumw2();
  m_h_ptclus_vs_off = new TH2D("m_h_ptclus_vs_off", "; p_{T, offline}[GeV]; p_{T, mtSA}[GeV]; events", 100, 0, 400, 100, 0, 400);
//  m_h_ptres_2off2clus = new TH1D("m_h_ptres_2off2clus", "1/p_{T} residual = (1/p_{T}^{SA}-1/p_{T}^{offline})/(1/p_{T}^{offline}); 1/p_{T} residual; events", 100, -1, 1);
//  m_h_ptres_2off2clus = new TH1D("m_h_ptres_2off2clus", "p_{T} residual = 1-p_{T}^{SAcluster}/p_{T}^{offline} leading&subleading; p_{T} residual; events", 5000, -49, 1);
  m_h_ptres_2off2clus_leading = new TH1D("m_h_ptres_2off2clus_leading", "p_{T} residual = 1-p_{T}^{SAcluster}/p_{T}^{offline} leading; p_{T} residual; events", 100, -1, 1);
  m_h_ptres_2off2clus_sub = new TH1D("m_h_ptres_2off2clus_sub", "p_{T} residual = 1-p_{T}^{SAcluster}/p_{T}^{offline} subleading; p_{T} residual; events", 100, -1, 1);

  m_hh_RoIEtaPhi_1road = new TH2D("m_hh_RoIEtaPhi_1road", "#eta_{RoI} vs #phi_{RoI} (1clusterRoad);#eta_{RoI}; #phi_{RoI}", 44, -1.1, 1.1, 128, -3.2, 3.2);
  m_h_deltaExtR_off_tot = new TH1D("m_h_deltaExtR_off_tot", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_off_morecand = new TH1D("m_h_deltaExtR_off_morecand", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  //m_h_deltaExtR_off_tot = new TH1D("m_h_deltaExtR_off_tot", "#DeltaR_{offline} at middle station; #DeltaR;events", 15, 2, 5);
  //m_h_deltaExtR_off_morecand = new TH1D("m_h_deltaExtR_off_morecand", "#DeltaR_{offline} at middle station; #DeltaR;events", 15, 2, 5);
  m_h_eff_isMorecand = new TH1D("m_h_eff_isMorecand", "R:rate of isMoreCand=true event !barrel event only!;#DeltaR_{offline}^{ext};rate", 10, 0, 0.1);
  //m_h_eff_isMorecand = new TH1D("m_h_eff_isMorecand", "R:rate of isMoreCand=true event !barrel event only!;#DeltaR_{offline}^{ext};rate", 15, 2, 5);
  //m_h_deltaExtR_off_tot_clus = new TH1D("m_h_deltaExtR_off_tot_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 15, 2, 5);
  m_h_deltaExtR_off_tot_clus = new TH1D("m_h_deltaExtR_off_tot_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
//  m_h_deltaExtR_off_tot_clus = new TH1D("m_h_deltaExtR_off_tot_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_off_2road_clus = new TH1D("m_h_deltaExtR_off_2road_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_off_1road_clus = new TH1D("m_h_deltaExtR_off_1road_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_off_manyroad_clus = new TH1D("m_h_deltaExtR_off_manyroad_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_eff_clusterRoad = new TH1D("m_h_eff_clusterRoad", "R:rate of N_{clusterRoadMid} = 2 event;#DeltaR_{offline}^{middle};rate", 10, 0, 0.1);
  //m_h_deltaExtR_off_2road_clus = new TH1D("m_h_deltaExtR_off_2road_clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 15, 2, 5);
  //m_h_eff_clusterRoad = new TH1D("m_h_eff_clusterRoad", "R:rate of N_{clusterRoadMid} = 2 event;#DeltaR_{offline}^{middle};rate", 15, 2, 5);
  m_h_dR_1road_sep = new TH1D("m_h_dR_1road_sep", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_canthelp = new TH1D("m_h_dR_1road_canthelp", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_outlier = new TH1D("m_h_dR_1road_outlier", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_misdire = new TH1D("m_h_dR_1road_misdire", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_others = new TH1D("m_h_dR_1road_others", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_phiothers = new TH1D("m_h_dR_1road_phiothers", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_1clus = new TH1D("m_h_dR_1road_1clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_1road_1phi = new TH1D("m_h_dR_1road_1phi", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_3road = new TH1D("m_h_dR_3road", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_dR_3road_3clus = new TH1D("m_h_dR_3road_3clus", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_2CR_1SPCMid = new TH1D("m_h_2CR_1SPCMid", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_2CR_1SPCMid_u3mdt = new TH1D("m_h_2CR_1SPCMid_u3mdt", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_h_2CR_1SPCMid_u3mdt_selected = new TH1D("m_h_2CR_1SPCMid_u3mdt_selected", "#DeltaR_{offline} at middle station; #DeltaR;events", 10, 0, 0.1);
  m_hh_Offpt_2mu1RoI = new TH2D("m_hh_Offpt_2mu1RoI", "Offline pT : 2mu-in-1RoI; leading p_{T, offline}; subleading p_{T, offline}; events", 100, 0, 500, 50, 0, 250);
  m_hh_Offpt_2mu1RoI_oppcharge = new TH2D("m_hh_Offpt_2mu1RoI_oppcharge", "Offline pT : 2mu-in-1RoI&oppcharge pT calculated; leading p_{T, offline}; subleading p_{T, offline}; events", 100, 0, 500, 50, 0, 250);
  m_h_Offpt_2mu1RoI_lead = new TH1D("m_h_Offpt_2mu1RoI_lead", "Offline pT : 2mu-in-1RoI; leading p_{T, offline}; events", 50, 0, 500);
  m_h_Offpt_2mu1RoI_lead_oppcharge = new TH1D("m_h_Offpt_2mu1RoI_lead_oppcharge", "Offline pT : 2mu-in-1RoI&opposite charge pT calculated; leading p_{T, offline}; events", 50, 0, 500);
  m_h_offlinePt_res0p7 = new TH1D("m_h_offlinePt_res0p7", "Offline pT (pT residual > 0.7); p_{T, offline}[GeV]; events", 250, 0, 500);
  m_peff_Offpt_lead = 0;
  m_peff_2spall = 0;
  m_peff_2spinn = 0;
  m_peff_2spmid = 0;
  m_peff_2spout = 0;
  m_peff_3pt = 0;
  m_peff_2pt = 0;
  m_peff_1pt = 0;
  m_peff_2pt_oppcharge = 0;
  m_peff_morecand = 0;
  m_peff_clusRoad = 0;
  m_peff_1clusRoad = 0;
  m_peff_manyclusRoad = 0;
  m_peff_closeBy = 0;
  m_peff_canthelp = 0;
  m_peff_outlier = 0;
  m_peff_misdire = 0;
  m_peff_others = 0;
  m_peff_1clus = 0;
  m_peff_1phi = 0;
  m_peff_phiothers = 0;
  m_peff_3road = 0;
  m_peff_2CR_1SPCMid = 0;
  m_peff_2CR_1SPCMid_all = 0;
 
  for(int i = 0; i < 4; i++){
    m_nclusPerLayer[i] = new TH1D(Form("m_nclusPerLayer%d", i), ";N_{RPC#eta cluster}/innermost RPC layer;events", 15, -0.5, 14.5);
  }

  m_h_deltaExtR_off_mismatch = new TH1D("m_h_deltaExtR_off_mismatch", "#DeltaR_{offline} at middle station !only mismatching event!; #DeltaR;events", 10, 0, 0.1);
  m_h_deltaExtR_off_closeBy = new TH1D("m_h_deltaExtR_off_closeBy", "#DeltaR_{offline} at middle station !only 2muon in 1RoI event!; #DeltaR;events", 10, 0, 0.1);
  h_OffEtaPhi = new TH2D("h_Offetaphi", ";#Delta#eta_{offline};#Delta#phi_{offline};events", 20, -0.1, 0.1, 20, -0.1, 0.1 );
  h_Offpt  = new TH1D("h_Offpt", ";p_{T}^{offline}[GeV];events", 500 ,0 , 500);
  m_h_SApt  = new TH1D("m_h_SApt", ";p_{T}^{SA}[GeV];events", 100 ,0 , 50);
  h_Offpt_lead  = new TH1D("h_Offpt_lead", ";p_{T}^{offline}[GeV];events", 500 ,0 , 500);
  h_Offpt_sublead  = new TH1D("h_Offpt_sublead", ";p_{T}^{offline}[GeV];events", 500 ,0 , 500);
//  h_Offpt  = new TH1D("h_Offpt", ";p_{T}^{offline}[GeV];events", 500 ,0 , 2000);
  h_OffdR  = new TH1D("h_OffdR", ";#DeltaR_{offline};events", 30 ,0 ,0.3);
  h_dtheta_sub = new TH1D("h_dtheta_sub", "#Delta#theta clusterRoad vs Offline;#Delta#theta(rad);events", 120, -0.3, 0.3);
  h_dtheta = new TH1D("h_dtheta", "#Delta#theta clusterRoad vs Offline;#Delta#theta(rad);events", 120, -0.3, 0.3);
  h_dtheta_total = new TH1D("h_dtheta_total", "#Delta#theta clusterRoad vs Offline;#Delta#theta(rad);events", 120, -0.3, 0.3);
//  m_h_ptclus  = new TH1D("m_h_ptclus", ";p_{T}^{SA, cluster}[GeV];events", 500 ,0 ,50);
  m_h_ptclus  = new TH1D("m_h_ptclus", ";p_{T}^{SA, cluster}[GeV];events", 500 ,0 ,2000);
  m_h_n_spall = new TH1D("m_h_n_spall", "the number of superpoint (N_{pT} != 2); N_{sp}; events", 13, -0.5, 12.5);
  m_h_Nclus_npt1 = new TH1D("m_h_Nclus_npt1", "the number of middle clusterRoad (N_{pT} < 2 & normal sector); N_{clusRoadMid}; events", 13, -0.5, 12.5);
  m_h_sector = new TH1D("m_h_sector", "N_{pT} < 2 sector; condition; events", 2, 0, 2);
  m_h_sector->GetXaxis()->SetBinLabel(1, "normal sec");
  m_h_sector->GetXaxis()->SetBinLabel(2, "special sec");
  m_h_ptres_vs_offsegInndR = new TH2D("m_h_ptres_vs_offsegInndR", "1/p_{T} residual vs offline segment #DeltaR@Inner; 1/p_{T} residual; #DeltaR; events", 100, -1, 4, 60, 0, 0.3);
  m_h_ptres_vs_offsegMiddR = new TH2D("m_h_ptres_vs_offsegMiddR", "1/p_{T} residual vs offline segment #DeltaR@Middle; 1/p_{T} residual; #DeltaR; events", 100, -1, 4, 60, 0, 0.3);
  m_h_ptres_vs_offsegOutdR = new TH2D("m_h_ptres_vs_offsegOutdR", "1/p_{T} residual vs offline segment #DeltaR@Outer; 1/p_{T} residual; #DeltaR; events", 100, -1, 4, 60, 0, 0.3);
  m_h_ptres_vs_offsegMindR = new TH2D("m_h_ptres_vs_offsegMindR", "1/p_{T} residual vs offline segment #DeltaR_{min}; 1/p_{T} residual; #DeltaR; events", 100, -1, 4, 60, 0, 0.3);
}

void RPC_FCBM::InitBranchHist(){
  h_Branch[0]  = new TH1D("h0", ";number of muon;events", 15, -0.5, 14.5);
  h_Branch[1]  = new TH1D("h1", ";p_{t,offline}[GeV];events", 100 ,0 ,100);
  h_Branch[2]  = new TH1D("h2", ";#eta_{offline};events", 120, -1.5, 1.5);
  h_Branch[3]  = new TH1D("h3", ";z_{0,offline}[mm];events", 200, -200, 200);
  h_Branch[4]  = new TH1D("h4", ";#phi_{offline};events", 140, -3.5, 3.5);
  h_Branch[5]  = new TH1D("h5", ";N_{offlineSegment};events", 15, -0.5, 14.5);
  h_Branch[6]  = new TH1D("h6", ";X_{offlineSegment}[mm];events", 150, -15000, 15000);
  h_Branch[7]  = new TH1D("h7", ";Y_{offlineSegment}[mm];events", 150, -15000, 15000);
  h_Branch[8]  = new TH1D("h8", ";Z_{offlineSegment}[mm];events", 200, -20000, 20000);
  h_Branch[9]  = new TH1D("h9" , ";DirectionX of OfflineSegment;events", 60, -1.5, 1.5); 
  h_Branch[10] = new TH1D("h10", ";DirectionY of OfflineSegment;events", 60, -1.5, 1.5);               
  h_Branch[11] = new TH1D("h11", ";DirectionZ of OfflineSegment;events", 60, -1.5, 1.5);               
  h_Branch[12] = new TH1D("h12", ";OfflineSegmentSector;events", 30, -0.5, 29.5);           
  h_Branch[13] = new TH1D("h13", ";p_{t,L1}[GeV];events", 30, 0, 15);   
  h_Branch[14] = new TH1D("h14", ";#eta_{L1};events", 100, -2.5, 2.5);  
  h_Branch[15] = new TH1D("h15", ";#phi_{L1};events", 140, -3.5, 3.5);                           
  h_Branch[16] = new TH1D("h16", ";L1RoINumber;events", 152, -1.5, 150.5);                     
  h_Branch[17] = new TH1D("h17", ";L1RoISector;events", 52, -1.5, 50.5);                     
  h_Branch[18] = new TH1D("h18", ";p_{t,SA}[GeV];events", 500, 0, 500);                            
  h_Branch[19] = new TH1D("h19", ";#eta_{SA};events", 100, -2.5, 2.5);                           
  h_Branch[20] = new TH1D("h20", ";#phi_{SA};events", 140, -3.5, 3.5);                           
  h_Branch[21] = new TH1D("h21", ";#eta_{SAMS};events", 100, -2.5, 2.5);                         
  h_Branch[22] = new TH1D("h22", ";#eta_{SABE};events", 100, -2.5, 2.5);                         
  h_Branch[23] = new TH1D("h23", ";#phi_{SAMS};events", 140, -3.5, 3.5);                         
  h_Branch[24] = new TH1D("h24", ";#phi_{SABE};events", 140, -3.5, 3.5);                         
  h_Branch[25] = new TH1D("h25", ";SABarrelRadius;events", 100, 0, 1); 
  h_Branch[26] = new TH1D("h26", ";SAsAddress;events", 9, -2.5, 6.5);                      
  h_Branch[27] = new TH1D("h27", ";eta_{SARoI};events", 100, -2.5, 2.5);                        
  h_Branch[28] = new TH1D("h28", ";phi_{SARoI};events", 140, -3.5, 3.5);                        
  h_Branch[29] = new TH1D("h29", ";SARoINumber;events", 152, -1.5, 150.5);                    
  h_Branch[30] = new TH1D("h30", ";SARoISector;events", 52, -1.5, 50.5);                    
  h_Branch[31] = new TH1D("h31", ";X_{SARPCHit}[mm];events", 150, -15000, 15000);                      
  h_Branch[32] = new TH1D("h32", ";Y_{SARPCHit}[mm];events", 150, -15000, 15000);                      
  h_Branch[33] = new TH1D("h33", ";Z_{SARPCHit}[mm];events", 150, -15000, 15000);                      
  h_Branch[34] = new TH1D("h34", ";R_{SARPCHit}[mm];events", 75, 5000, 20000);                      
  h_Branch[35] = new TH1D("h35", ";eta_{SARPCHit};events", 100, -2.5, 2.5);                    
  h_Branch[36] = new TH1D("h36", ";phi_{SARPCHit};events", 140, -3.5, 3.5);                    
  h_Branch[37] = new TH1D("h37", ";SARPCHitMeasuresPhi;events", 5, -0.5, 4.5);            
  h_Branch[38] = new TH1D("h38", ";SARPCHitStationName;events", 15, -0.5, 14.5);            
  h_Branch[39] = new TH1D("h39", ";SARPCHitStationNumber;events", 15, -0.5, 14.5);          
  h_Branch[40] = new TH1D("h40", ";phi_{SARPCFitInn};events", 140, -3.5, 3.5);                 
  h_Branch[41] = new TH1D("h41", ";SARPCFitInnSlope;events", 120, -30, 30);               
  h_Branch[42] = new TH1D("h42", ";SARPCFitInnOffset[mm];events", 1000, -10000, 10000);              
  h_Branch[43] = new TH1D("h43", ";phi_{SARPCFitMid};events", 35, -3.5, 3.5);                 
  h_Branch[44] = new TH1D("h44", ";SARPCFitMidSlope;events", 120, -30, 30);               
  h_Branch[45] = new TH1D("h45", ";SARPCFitMidOffset[mm];events", 1000, -10000, 10000);              
  h_Branch[46] = new TH1D("h46", ";phi_{SARPCFitOut};events", 35, -3.5, 3.5);                 
  h_Branch[47] = new TH1D("h47", ";SARPCFitOutSlope;events", 120, -30, 30);               
  h_Branch[48] = new TH1D("h48", ";SARPCFitOutOffset[mm];events", 1000, -10000, 10000);              
  h_Branch[49] = new TH1D("h49", ";SARoadAw;events", 120, -30, 30); 
  h_Branch[50] = new TH1D("h50", ";SARoadBw[mm];events", 100, -100, 100); 
  h_Branch[51] = new TH1D("h51", ";R_{OfflineSegment}[mm];events", 75, 0, 15000); 
  h_Branch[52] = new TH1D("h52", ";p_{t,offline}[GeV];events", 100 ,0 ,50);
  h_Branch[57] = new TH1D("h57", "SARPCHitTime;time(ns);events", 52, -0.5, 25.5);            
  h_Branch[58] = new TH1D("h58", "clusterRoad_slope(middle);slope;events", 200, -100, 100);            
  h_Branch[59] = new TH1D("h59", "clusterRoad_offset(middle);offset;events", 200, -100, 100);            
  h_Branch[60] = new TH1D("h60", "clusterRoad_phi(middle);phi;events", 35, -3.5, 3.5);            
}

void RPC_FCBM::InitHist(){
  InitFCBM();
  InitNoCut();
  InitBranchHist();
}

void RPC_FCBM::End(){
  if(h_distRoIrpcEta != 0) {
    delete h_distRoIrpcEta;
  }
  if(h_distRoIrpcPhi != 0) {
    delete h_distRoIrpcPhi;
  }
  if(h_SAEtaPhi != 0) {
    delete h_SAEtaPhi;
  }
  if(h_SARoIEtaPhi != 0) {
    delete h_SARoIEtaPhi;
  }
  if(h_distnoCutRoIrpcEta != 0) {
    delete h_distnoCutRoIrpcEta;
  }
  if(h_distnoCutRoIrpcPhi != 0) {
    delete h_distnoCutRoIrpcPhi;
  }
  if(h_distnoCutEachRoIrpcEta != 0) {
    delete h_distnoCutEachRoIrpcEta;
  }
  if(h_distnoCutEachRoIrpcPhi != 0) {
    delete h_distnoCutEachRoIrpcPhi;
  }
/*  for(int iBranch = 0; iBranch < nBranch; iBranch++)
  {
    if(h_Branch[iBranch] != 0) {
      delete h_Branch[iBranch];
    }
  }*/
}

Bool_t RPC_FCBM::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void RPC_FCBM::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}


Int_t RPC_FCBM::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


#endif // RPC_FCBM_hhh
#endif // #ifdef RPC_FCBM_cxx

