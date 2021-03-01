#include "../RpcL2MuonSA/RPC_FCBM.h"
#include "TVector3.h"
#include <math.h>
using namespace std;
const double pi = 3.14159265358979323846;

/**
 * this function is used for calculation of efficiency
 * at first, I develop this for the test of "isMoreCandInRoI"
 * second, for clusterRoad 
 * =====strategy of calculation for efficiency=========
 * 1. fill each parameter
 * 2. get each bin content (!bin start from 1!)
 * 3. calculate the efficiency for each bin
 * 4. also calculate the error
 * 5. draw with option "P" using TGraphErrors
 *=====================================================*/
void RPC_FCBM::FillEfficiency()
{
  /*vector<double> deltaExtR_off;
    deltaExtR_off.clear();*/
  bool isDimuon = true;
  for(int iMuon = 0; iMuon < nMuon; iMuon++){
    std::cout << "# muon = " << iMuon << std::endl;
    for(int i=0; i< (int)SAptclus->at(iMuon).size(); i++){
      std::cout << "SAptclus = " << SAptclus->at(iMuon).at(i) << std::endl;
    }
  }
  if(nMuon != 2) isDimuon = false;
  if(isDimuon){
    bool is2mu1RoI = false;
    if(L1nRoI->at(0) == 1 &&
        L1RoINumber->at(0) == L1RoINumber->at(1) &&
        L1RoISector->at(0) == L1RoISector->at(1) )
    {
      is2mu1RoI = true;
      cout << "FillEfficiency::this is 2muon in 1RoI event" << endl;
    }//L1 mismatching cut
    double deltaExtEta = OfflineExtEta->at(0) - OfflineExtEta->at(1);
    double deltaExtPhi = OfflineExtPhi->at(0) - OfflineExtPhi->at(1);
    double deltaExtR_off = sqrt(pow(deltaExtEta, 2) + pow(deltaExtPhi, 2));
    cout << "FillEfficiency::deltaExtR = " << deltaExtR_off << endl;
    m_h_deltaExtR_off_tot->Fill(deltaExtR_off);
    /* clusterRoad analysis start *************************************/
    unsigned int nClus_fit = (SARPCCluster_fitMidSlope->at(1)).size();
    int N_clusRoadInn = SARPCCluster_fitInnSlope->at(0).size();
    int N_clusRoadMid = SARPCCluster_fitMidSlope->at(0).size();
    int N_clusRoadOut = SARPCCluster_fitOutSlope->at(0).size();
    
//    int N_clusRoadInn = 0;
//    int N_clusRoadMid = 0;
//    int N_clusRoadOut = 0;
    
//    for(unsigned int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
//      if(SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit))N_clusRoadInn++;
//      if(SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit))N_clusRoadMid++;
//      if(SARPCCluster_isPlausibleFitOut->at(0).at(iClus_fit))N_clusRoadOut++;
//    }
//    cout << "FillEfficiency::the number of clusterRoadMid = " << N_clusRoadMid << endl;
//    cout << "FillEfficiency::isMoreCand->at(" << 0 << ") = " << L1isMoreCandInRoI->at(0) << endl;
    int N_spinn = 0;
    int N_spmid = 0;
    int N_spout = 0;
    for(int iSPinn=0; iSPinn<SASPCZ_BI->at(0).size(); iSPinn++){
      if( std::fabs(SASPCZ_BI->at(0).at(iSPinn)) < 1e-05 && std::fabs(SASPCR_BI->at(0).at(iSPinn)) < 1e-05) continue;
      N_spinn++;
    }
    for(int iSPmid=0; iSPmid<SASPCZ_BM->at(0).size(); iSPmid++){
      if( std::fabs(SASPCZ_BM->at(0).at(iSPmid)) < 1e-05 && std::fabs(SASPCR_BM->at(0).at(iSPmid)) < 1e-05) continue;
      N_spmid++;
    }
    for(int iSPout=0; iSPout<SASPCZ_BO->at(0).size(); iSPout++){
      if( std::fabs(SASPCZ_BO->at(0).at(iSPout)) < 1e-05 && std::fabs(SASPCR_BO->at(0).at(iSPout)) < 1e-05) continue;
      N_spout++;
    }

    //SA muon6GeV thresholds
    double SAThresholds[4] = {5.92,  5.86,  5.70,  5.64};
    double SAWeakThresholds[2] = {3.91,  2.22};
    int WeakRegion = ECWeakRegion( SAEta->at(0), SAPhi->at(0) );
    int etaRegion = EtaRegion( SAEta->at(0) );
    double threshold = 0;
    //get SA muon eta region
    if( WeakRegion == 0){
      threshold = SAThresholds[etaRegion];
    } else {
      threshold = SAWeakThresholds[WeakRegion];
    }
    std::cout << "SA Hypo threshold = " << threshold << std::endl;
    int N_pt = 0;
    for(int i = 0; i < SAptclus->at(0).size(); i++){
      if(SAptclus->at(0).at(i) == -99999.) continue;
      if(fabs(SAptclus->at(0).at(i)) > threshold) N_pt++;
    }
    int N_overZLIM = 0;
    for(int i=0; i<N_pt; i++){
      if(SAptclus->at(0).at(i) > 1e-5) N_overZLIM++;
    }
    double ptlead = -999999;
    double ptsub = 999999;
    for(int iMuon = 0; iMuon< nMuon; iMuon++){
      if(ptlead < OfflinePt->at(iMuon)) ptlead = OfflinePt->at(iMuon);
      if(ptsub > OfflinePt->at(iMuon)) ptsub = OfflinePt->at(iMuon);
    }
    if(L1RoINumber->at(0) != L1RoINumber->at(1) ||
       L1RoISector->at(0) != L1RoISector->at(1) )
    {
      bool skipEvent = false;
      for(int iMuon=0; iMuon<nMuon; iMuon++){
        if(SAptclus->at(iMuon).size() != 1) skipEvent = true;
      }
      if(!skipEvent){
        for(int iMuon=0; iMuon<nMuon; iMuon++){
          std::cout << "ptcluster/default pt = " << SAptclus->at(iMuon).at(0) << "/" << SAPt->at(iMuon) << std::endl;
          float iptsa = 1/fabs(SAPt->at(iMuon));
          float iptmpsa = 1/fabs(SAptclus->at(iMuon).at(0));
          float iptoff = 1/fabs(OfflinePt->at(iMuon));

          float res_sapt = (iptsa-iptoff)/iptoff;
          float res_mpsapt = (iptmpsa-iptoff)/iptoff;
          m_total++;
          if(res_sapt > 0.1 || res_sapt < -0.1) m_satail++;
          if(res_mpsapt > 0.1 || res_mpsapt < -0.1) m_mpsatail++;
          m_ptres_default->Fill(res_sapt);
          m_ptres_mpsa->Fill(res_mpsapt);
        }
      }
    }

    //for 2muon in 1RoI event (in order to protect the overlap filling)
    if(is2mu1RoI){
      int nclusPerLay[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      int nclusterPerLay = -999;
      for(int i=0; i < SARPCCluster_clusterLayer->at(0).size(); i++){
        if(SARPCCluster_clusterMeasPhi->at(0).at(i)) continue;
        float r = sqrt(SARPCCluster_gX->at(0).at(i)*SARPCCluster_gX->at(0).at(i)+SARPCCluster_gY->at(0).at(i)*SARPCCluster_gY->at(0).at(i));
        float phi = atan(SARPCCluster_gY->at(0).at(i)/SARPCCluster_gX->at(0).at(i));
        cout << "phi = " << phi << endl;
        if (SARPCCluster_gX->at(0).at(i)<0 && SARPCCluster_gY->at(0).at(i)>0) phi = phi + pi;
        if (SARPCCluster_gX->at(0).at(i)<0 && SARPCCluster_gY->at(0).at(i)<0) phi = phi - pi;
        cout << "phi renew = " << phi << endl;
        float l = sqrt(SARPCCluster_gZ->at(0).at(i)*SARPCCluster_gZ->at(0).at(i)+r*r);
        float tan = sqrt( (l-SARPCCluster_gZ->at(0).at(i))/(l+SARPCCluster_gZ->at(0).at(i)) );
        float eta = -log(tan);
        float deta = fabs(SARoIEta->at(0) - eta);
        float dphi = acos(cos( SARoIPhi->at(0) - phi ) );
        std::cout << "deta/dphi = " << deta << "/" << dphi << std::endl;
        //====================================================
        if(deta >= 0.1 && dphi >= 0.1) continue;
        int layer = SARPCCluster_clusterLayer->at(0).at(i);
        nclusPerLay[layer]++;
      }
      for(int i=0; i < 8; i++){
        if(nclusPerLay[i] > nclusterPerLay) nclusterPerLay = nclusPerLay[i];
      }
      m_h_countPtclus->Fill(N_pt);
      m_nclus_vs_npt->Fill(nclusterPerLay, N_pt);
      m_hh_Offpt_2mu1RoI->Fill(ptlead, ptsub);
      m_h_Offpt_2mu1RoI_lead->Fill(ptsub);
      m_h_deltaExtR_off_tot_clus->Fill(deltaExtR_off);
//      m_h_Offpt_2mu1RoI_lead->Fill(ptlead);
      if(N_spinn == 2)m_h_deltaExtR_2spinn->Fill(deltaExtR_off);
      if(N_spmid == 2)m_h_deltaExtR_2spmid->Fill(deltaExtR_off);
      if(N_spout == 2)m_h_deltaExtR_2spout->Fill(deltaExtR_off);
      if(N_spinn == 2 && N_spmid == 2 && N_spout == 2) m_h_deltaExtR_2spall->Fill(deltaExtR_off);
      if(N_pt != 2) m_h_n_spall->Fill(N_spinn+N_spmid+N_spout);
      if(N_pt > 2) m_h_deltaExtR_3pt->Fill(deltaExtR_off);
      if(N_pt < 2){
//        if(N_spmid < 2){
//          if(N_clusRoadMid < 2){m_h_nptless->Fill(0);} 
//          if(N_clusRoadMid == 2 && ){m_h_nptless->Fill(0);} 
//        }
        if(SAsAddress->at(0) > 1) m_h_sector->Fill(1);
        if(SAsAddress->at(0) == 0 || SAsAddress->at(0) == 1){
          m_h_sector->Fill(0);
        }
        if(N_spmid < 2) m_h_Nclus_npt1->Fill(N_clusRoadMid);
        m_h_deltaExtR_1pt->Fill(deltaExtR_off);
      }
      if(N_pt == 1){
        m_nclusPerLayer[0]->Fill(nclusterPerLay);
      }
      if(N_pt == 3) m_nclusPerLayer[2]->Fill(nclusterPerLay);
      if(N_pt == 4) m_nclusPerLayer[3]->Fill(nclusterPerLay);
      if(N_pt == 2){
        m_nclusPerLayer[1]->Fill(nclusterPerLay);
        double ischarge;
        if(isOld){
          ischarge = SAchargeclus->at(0).at(0) * SAchargeclus->at(0).at(1);
        } else {
          if(SAptclus->at(0).at(0) > 0 && SAptclus->at(0).at(1) < 0) ischarge = -1;
          else if(SAptclus->at(0).at(0) < 0 && SAptclus->at(0).at(1) > 0) ischarge = -1;
          else { ischarge = 1; }
        }
        m_h_deltaExtR_2pt->Fill(deltaExtR_off);
        if(ischarge < 0){
          m_h_deltaExtR_2pt_oppcharge->Fill(deltaExtR_off);
          m_hh_Offpt_2mu1RoI_oppcharge->Fill(ptlead, ptsub);
//          m_h_Offpt_2mu1RoI_lead_oppcharge->Fill(ptlead);
          m_h_Offpt_2mu1RoI_lead_oppcharge->Fill(ptsub);
        }
        // compare pt from clusterRoad with OfflinePt
        if(OfflinePt->size() == 2){
          if(std::abs(OfflinePt->at(0) - OfflinePt->at(1)) > 1e-5){
            double ptOff_leading = -9999;
            double ptOff_sub = -9999;
            double ptclus_leading = -9999;
            double ptclus_sub = -9999;
            int leadingID[2] = {0, 0}; //{off, clus}
            for(int i=0; i<2; i++){
              if(std::abs(OfflinePt->at(i)) > ptOff_leading){
                ptOff_leading = std::fabs(OfflinePt->at(i));
                leadingID[0] = i;
              }
              if(std::abs(SAptclus->at(0).at(i)) > ptclus_leading){
                ptclus_leading = std::fabs(SAptclus->at(0).at(i));
                leadingID[1] = i;
              }
            }
            ptOff_sub = (leadingID[0]) ? fabs(OfflinePt->at(0)) : fabs(OfflinePt->at(1));
            ptclus_sub = (leadingID[1]) ? fabs(SAptclus->at(0).at(0)) : fabs(SAptclus->at(0).at(1));
            std::cout << "pt leading(off, clus) = (" << ptOff_leading << "," << ptclus_leading << ")" << std::endl;
            std::cout << "pt sub(off, clus) = (" << ptOff_sub << "," << ptclus_sub << ")" << std::endl;
            double ptres_leading = 1 - ptOff_leading/ptclus_leading;
            double ptres_sub = 1 - ptOff_sub/ptclus_sub;
            m_total_closeby++;
            m_total_closeby++;
            if(ptres_leading > 0.1 || ptres_leading < -0.1) m_mpsatail_closeby++;
            if(ptres_sub > 0.1 || ptres_sub < -0.1) m_mpsatail_closeby++;
//            m_h_ptres_2off2clus->Fill(ptres_leading);
//            m_h_ptres_2off2clus->Fill(ptres_sub);
            m_h_ptres_2off2clus_leading->Fill(ptres_leading);
            m_h_ptres_2off2clus_sub->Fill(ptres_sub);
            m_h_ptclus_vs_off->Fill(ptOff_leading, ptclus_leading);
            m_h_ptclus_vs_off->Fill(ptOff_sub, ptclus_sub);
            if(SASPCZ_BM->at(0).size() == 2){
              for(int i_spmid=0; i_spmid<(int)SASPCZ_BM->at(0).size(); i_spmid++){
                if(fabs(SASPCZ_BM->at(0).at(i_spmid)) > 1e-5 && fabs(SASPCR_BM->at(0).at(i_spmid)) > 1e-5){
                  std::cout << "SuperPoint from cluster R/Z = " << SASPCR_BM->at(0).at(i_spmid) << "/" << SASPCZ_BM->at(0).at(i_spmid) << std::endl;
                  float theta = atan(SASPCR_BM->at(0).at(i_spmid)/fabsf(SASPCZ_BM->at(0).at(i_spmid)));
                  float eta = (tan(theta/2.)!=0.)? -log(tan(theta/2.))*SASPCZ_BM->at(0).at(i_spmid)/fabsf(SASPCZ_BM->at(0).at(i_spmid)): 0.;
                  std::cout << "superpoint (BM) eta = " << eta << std::endl;
                  float deta_min = 999999;
                  float ptres = -99999;
                  int ID = -1;
                  for(int i_mu=0; i_mu < (int)OfflineExtEta->size(); i_mu++){
                    float fdeta = fabs(OfflineExtEta->at(i_mu) - eta);
                    if(fdeta < deta_min){ 
                      deta_min = fdeta;
                      float fiSAptclus = 1/fabs(SAptclus->at(0).at(i_spmid));
                      float fiOfflinePt = 1/fabs(OfflinePt->at(i_mu));
                      ptres = (fiSAptclus - fiOfflinePt)/fiOfflinePt;
                      ID = i_mu;
                    }
                  }
                  std::cout << "deta matching NO/ptres = " << ID << "/" << ptres << std::endl;
                  m_h_ptres_2off2clus->Fill(ptres);
                  if(ptres > 0.7){
                    m_h_offlinePt_res0p7->Fill(fabs(OfflinePt->at(ID)));
                  }
                  std::vector<TVector3> offsegInn, offsegMid, offsegOut;
                  for(int iMuon=0; iMuon < nMuon; iMuon++){
                    int nSeg = OfflineNumSegment->at(iMuon);
                    for(int iSeg = 0; iSeg < nSeg; iSeg++){
                      TVector3 vecOffseg;
                      vecOffseg.SetXYZ((OfflineSegmentX->at(iMuon)).at(iSeg), (OfflineSegmentY->at(iMuon)).at(iSeg), (OfflineSegmentZ->at(iMuon)).at(iSeg));
                      if(fabs(vecOffseg.Eta()) >= 1.05) continue;
                      if(vecOffseg.Perp() <= 6000){
                        bool isDiffseg = true;
                        for(int i_seg=0; i_seg < (int)offsegInn.size(); i_seg++){
                          if(fabs(offsegInn.at(i_seg).Eta() - vecOffseg.Eta()) < 1e-4 && 
                              fabs(offsegInn.at(i_seg).Phi() - vecOffseg.Phi()) < 1e-4) isDiffseg = false;
                        }
                        if(isDiffseg) offsegInn.push_back(vecOffseg);
                      }
                      else if(6500. < vecOffseg.Perp() && vecOffseg.Perp() < 8500.){
                        bool isDiffseg = true;
                        for(int i_seg=0; i_seg < (int)offsegMid.size(); i_seg++){
                          if(fabs(offsegMid.at(i_seg).Eta() - vecOffseg.Eta()) < 1e-4 && 
                              fabs(offsegMid.at(i_seg).Phi() - vecOffseg.Phi()) < 1e-4) isDiffseg = false;
                        }
                        if(isDiffseg) offsegMid.push_back(vecOffseg);
                      }
                      else if(vecOffseg.Perp() >= 8500.){ 
                        bool isDiffseg = true;
                        for(int i_seg=0; i_seg < (int)offsegOut.size(); i_seg++){
                          if(fabs(offsegOut.at(i_seg).Eta() - vecOffseg.Eta()) < 1e-4 && 
                              fabs(offsegOut.at(i_seg).Phi() - vecOffseg.Phi()) < 1e-4) isDiffseg = false;
                        }
                        if(isDiffseg) offsegOut.push_back(vecOffseg);
                      }
                    }
                  }
                  std::cout << "offline segment, N inner/middle/outer = " << offsegInn.size() << "/" << offsegMid.size() << "/" << offsegOut.size() << std::endl;
                  float dR_min = 99999;
                  if(offsegInn.size() == 2){
                    float deta = offsegInn.at(0).Eta() - offsegInn.at(1).Eta();
                    float dphi = offsegInn.at(0).Phi() - offsegInn.at(1).Phi();
                    float dR = sqrt(pow(deta, 2) + pow(dphi,2));
                    m_h_ptres_vs_offsegInndR->Fill(ptres, dR);
                    if(dR_min > dR) dR_min = dR;
                  }
                  if(offsegMid.size() == 2){
                    float deta = offsegMid.at(0).Eta() - offsegMid.at(1).Eta();
                    float dphi = offsegMid.at(0).Phi() - offsegMid.at(1).Phi();
                    float dR = sqrt(pow(deta, 2) + pow(dphi,2));
                    m_h_ptres_vs_offsegMiddR->Fill(ptres, dR);
                    if(dR_min > dR) dR_min = dR;
                  }
                  if(offsegOut.size() == 2){
                    float deta = offsegOut.at(0).Eta() - offsegOut.at(1).Eta();
                    float dphi = offsegOut.at(0).Phi() - offsegOut.at(1).Phi();
                    float dR = sqrt(pow(deta, 2) + pow(dphi,2));
                    m_h_ptres_vs_offsegOutdR->Fill(ptres, dR);
                    if(dR_min > dR) dR_min = dR;
                  }
                  if(dR_min < 99999){
                    m_h_ptres_vs_offsegMindR->Fill(ptres, dR_min);
                  }
                }
              }
            }

            // probe the condition of SuperPoint when ptclus are much less than Offline pT
            if(ptres_sub > 0.7){
              if(N_spinn > 1 && N_spmid > 1 && N_spout > 1){
                m_h_n_lowptclus->Fill("A", 1);
                m_h_n_lowptclus_patA->Fill("all stats have 2 or more SPs", 1);
//                if(SAspcSetID->at(0).at(1).at(0) == SAspcSetID->at(0).at(1).at(1)){
//                  m_h_n_lowptclus_patA_sameIDout->Fill("all stats have 2 or more SPs", 1);
//                
//                  for(int i_spinn=0; i_spinn<(int)SASPCZ_BI->at(0).size(); i_spinn++){
//                    float Z_spinn = SASPCZ_BI->at(0).at(i_spinn);
//                    float R_spinn = SASPCR_BI->at(0).at(i_spinn);
//                    double slope_inn = Z_spinn/R_spinn;
//                    double slopemin_innmid = 9999;
//                    double slopemin_innout = 9999;
//                    for(int i_spmid=0; i_spmid<(int)SASPCZ_BM->at(0).size(); i_spmid++){
//                      float Z_spmid = SASPCZ_BM->at(0).at(i_spmid);
//                      float R_spmid = SASPCR_BM->at(0).at(i_spmid);
//                      double slope_mid = Z_spmid/R_spmid;
//                      if(std::abs(slope_mid - slope_inn) < slopemin_innmid) slopemin_innmid = std::abs(slope_mid - slope_inn);
//                    }
//                    for(int i_spout=0; i_spout<(int)SASPCZ_BO->at(0).size(); i_spout++){
//                      float Z_spout = SASPCZ_BO->at(0).at(i_spout);
//                      float R_spout = SASPCR_BO->at(0).at(i_spout);
//                      double slope_out = Z_spout/R_spout;
//                      if(std::abs(slope_out - slope_inn) < slopemin_innout) slopemin_innout = std::abs(slope_out - slope_inn);
//                    }
//                    if(slopemin_innmid > slopemin_innout) m_h_n_lowptclus_patA_findout->Fill("all stats have 2 or more SPs", 1);
//                  }
//                }
              }
              else if((N_spinn > 1 && N_spmid > 1 && N_spout < 2)||
                 (N_spinn > 1 && N_spmid < 2 && N_spout > 1)||
                 (N_spinn < 2 && N_spmid > 1 && N_spout > 1)){
                m_h_n_lowptclus->Fill("B", 1);
                if(N_spinn < 2) m_h_n_1spstat_inn->Fill("B",1);
                if(N_spmid < 2) m_h_n_1spstat_mid->Fill("B",1);
                if(N_spout < 2) m_h_n_1spstat_out->Fill("B",1);
              }
              else if((N_spinn < 2 && N_spmid < 2 && N_spout > 1)||
                 (N_spinn < 2 && N_spmid > 1 && N_spout < 2)||
                 (N_spinn > 1 && N_spmid < 2 && N_spout < 2)){
                m_h_n_lowptclus->Fill("C", 1);
              }
//              m_h_n_spinn_lowpt->Fill(N_spinn);
//              m_h_n_spmid_lowpt->Fill(N_spmid);
//              m_h_n_spout_lowpt->Fill(N_spout);
            }
          }
        }
      }
      //check the reason why dropping Superpoint despite there are 2 clusterRoads
      if(N_spinn == 1 && N_clusRoadInn == 2){
        m_h_2CR_1SPCMid->Fill(deltaExtR_off);
        int NmdtMin = 999;
        for(int iclus = 0; iclus < 2; iclus++){
          int Nmdt = 0;
          for(int iMdt = 0; iMdt < SAMDTHitAllR->at(0).size(); iMdt++){
            float mdtR = SAMDTHitAllR->at(0).at(iMdt);
//            float minR = SARPCCluster_rMin->at(0).at(1).at(0);
//            float maxR = SARPCCluster_rMax->at(0).at(1).at(0);
            //float minR = SARPCCluster_rMin->at(0).at(0).at(0);
            //float maxR = SARPCCluster_rMax->at(0).at(0).at(0);
            //if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0 && SAMDTHitAllisOutlier->at(0).at(iMdt) == 5) continue;
            //if(minR <= mdtR && mdtR <= maxR 
            //    && SAMDTHitAllclusRoadID->at(0).at(iMdt) == iclus){
            //  cout << "mdt R = " << mdtR << ": mdt rmin/rmax = " << minR << "/" << maxR << endl;
            //  Nmdt++;
            //}
          }
          if(Nmdt < NmdtMin) NmdtMin = Nmdt;
        }
        if(NmdtMin < 3) m_h_2CR_1SPCMid_u3mdt->Fill(deltaExtR_off);
      }
      else if(N_spinn == 0 && N_clusRoadInn == 2){
        m_h_2CR_1SPCMid->Fill(deltaExtR_off);
        int NmdtMax = -9;
        for(int iclus = 0; iclus < 2; iclus++){
          int Nmdt = 0;
          for(int iMdt = 0; iMdt < SAMDTHitAllR->at(0).size(); iMdt++){
            float mdtR = SAMDTHitAllR->at(0).at(iMdt);
            //float minR = SARPCCluster_rMin->at(0).at(0).at(0);
            //float maxR = SARPCCluster_rMax->at(0).at(0).at(0);
//            float minR = SARPCCluster_rMin->at(0).at(1).at(0);
//            float maxR = SARPCCluster_rMax->at(0).at(1).at(0);
            //if(SAMDTHitAllisOutlier->at(0).at(iMdt) == 0 && SAMDTHitAllisOutlier->at(0).at(iMdt) == 5) continue;
            //if(minR <= mdtR && mdtR <= maxR 
            //    && SAMDTHitAllclusRoadID->at(0).at(iMdt) == iclus){
            //  cout << "mdt R = " << mdtR << ": mdt rmin/rmax = " << minR << "/" << maxR << endl;
            //  Nmdt++;
            //}
          }
          if(Nmdt > NmdtMax) NmdtMax = Nmdt;
        }
        if(NmdtMax < 3) m_h_2CR_1SPCMid_u3mdt->Fill(deltaExtR_off);
      }
//      if(N_spmid == 1 && N_clusRoadMid == 2){
      if(N_spinn == 1 && N_clusRoadInn == 2){
        int NmdtMin = 999;
        for(int iclus = 0; iclus < 2; iclus++){
          int Nmdt = 0;
          for(int iMdt = 0; iMdt < SAMDTHitAllR->at(0).size(); iMdt++){
            float mdtR = SAMDTHitAllR->at(0).at(iMdt);
//            float minR = SARPCCluster_rMin->at(0).at(1).at(0);
//            float maxR = SARPCCluster_rMax->at(0).at(1).at(0);
            //float minR = SARPCCluster_rMin->at(0).at(0).at(0);
            //float maxR = SARPCCluster_rMax->at(0).at(0).at(0);
            //if(SAMDTHitAllisOutlier->at(0).at(iMdt) > 2) continue;
            //if(minR <= mdtR && mdtR <= maxR 
            //    && SAMDTHitAllclusRoadID->at(0).at(iMdt) == iclus){
            //  cout << "mdt R = " << mdtR << ": mdt rmin/rmax = " << minR << "/" << maxR << endl;
            //  Nmdt++;
            //}
          }
          if(Nmdt < NmdtMin) NmdtMin = Nmdt;
        }
        if(NmdtMin < 3) m_h_2CR_1SPCMid_u3mdt_selected->Fill(deltaExtR_off);
      }
//      else if(N_spmid == 0 && N_clusRoadMid == 2){
      else if(N_spinn == 0 && N_clusRoadInn == 2){
        int NmdtMax = -9;
//        comment out 2020/07/21
//        for(int iclus = 0; iclus < 2; iclus++){
//          int Nmdt = 0;
//          for(int iMdt = 0; iMdt < SAMDTHitAllR->at(0).size(); iMdt++){
//            float mdtR = SAMDTHitAllR->at(0).at(iMdt);
//            float minR = SARPCCluster_rMin->at(0).at(0).at(0);
//            float maxR = SARPCCluster_rMax->at(0).at(0).at(0);
////            float minR = SARPCCluster_rMin->at(0).at(1).at(0);
////            float maxR = SARPCCluster_rMax->at(0).at(1).at(0);
//            if(SAMDTHitAllisOutlier->at(0).at(iMdt) > 2) continue;
//            if(minR <= mdtR && mdtR <= maxR 
//                && SAMDTHitAllclusRoadID->at(0).at(iMdt) == iclus){
//              cout << "mdt R = " << mdtR << ": mdt rmin/rmax = " << minR << "/" << maxR << endl;
//              Nmdt++;
//            }
//          }
//          if(Nmdt > NmdtMax) NmdtMax = Nmdt;
//        }
        if(NmdtMax < 3) m_h_2CR_1SPCMid_u3mdt_selected->Fill(deltaExtR_off);
      }

      m_h_deltaExtR_off_closeBy->Fill(deltaExtR_off);
      if(N_clusRoadMid == 2 /*|| N_clusRoadOut == 2*/){
        m_h_deltaExtR_off_2road_clus->Fill(deltaExtR_off);
      }
      else if(N_clusRoadMid == 1 /*|| N_clusRoadOut == 1*/){
        m_h_deltaExtR_off_1road_clus->Fill(deltaExtR_off);
        for(int iMuon = 0; iMuon < nMuon; iMuon++){
          //            m_hh_RoIEtaPhi_1road->Fill(SARoIEta->at(0), SARoIPhi->at(0));
          m_hh_RoIEtaPhi_1road->Fill(OfflineExtEta->at(iMuon), OfflineExtPhi->at(iMuon));
        }
      }
      else if(N_clusRoadMid > 2 /*|| N_clusRoadOut == 1*/){
        m_h_deltaExtR_off_manyroad_clus->Fill(deltaExtR_off);
      }
      else{ cout << "what's this event?" << endl; }
      if(N_clusRoadMid == 1 /*&& deltaExtR_off > 0.02*/){
        m_h_dR_1road_sep->Fill(deltaExtR_off);
        int nClus = SARPCCluster_clusterLayer->at(0).size();
        unsigned int nClus_phi = SARPCCluster_fitMidPhi->at(0).size();
        int countPhiSet = 0;
        for(unsigned int iClus_phi = 0; iClus_phi < nClus_phi; iClus_phi++){
          //if(SARPCCluster_isPlausiblePhiInnMid->at(0).at(iClus_phi)){ countPhiSet++; }
        }
        if(countPhiSet == 1){
//          for(int iMuon = 0; iMuon < nMuon; iMuon++){
////            m_hh_RoIEtaPhi_1road->Fill(SARoIEta->at(0), SARoIPhi->at(0));
//            m_hh_RoIEtaPhi_1road->Fill(OfflineExtEta->at(iMuon), OfflineExtPhi->at(iMuon));
//          }
        }

        vector<vector<int>> clus_layer, clus_outlier;
        vector<int> forClear;
        forClear.clear();
        clus_layer.assign(4, forClear);
        clus_outlier.assign(4, forClear);
        for(int iClus = 0; iClus < nClus; iClus++){
          if(SARPCCluster_clusterMeasPhi->at(0).at(iClus)) continue;
          //======check where the test cluster is located======
          float r = sqrt(SARPCCluster_gX->at(0).at(iClus)*SARPCCluster_gX->at(0).at(iClus)+SARPCCluster_gY->at(0).at(iClus)*SARPCCluster_gY->at(0).at(iClus));
          float phi = atan(SARPCCluster_gY->at(0).at(iClus)/SARPCCluster_gX->at(0).at(iClus));
          cout << "phi = " << phi << endl;
          if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)>0) phi = phi + pi;
          if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)<0) phi = phi - pi;
          cout << "phi renew = " << phi << endl;
          float l = sqrt(SARPCCluster_gZ->at(0).at(iClus)*SARPCCluster_gZ->at(0).at(iClus)+r*r);
          float tan = sqrt( (l-SARPCCluster_gZ->at(0).at(iClus))/(l+SARPCCluster_gZ->at(0).at(iClus)) );
          float eta = -log(tan);
          float deta = fabs(SARoIEta->at(0) - eta);
          float dphi = acos(cos( SARoIPhi->at(0) - phi ) );
          std::cout << "deta/dphi = " << deta << "/" << dphi << std::endl;
          //====================================================
          int i_layer = SARPCCluster_clusterLayer->at(0).at(iClus);
          if((i_layer < 4) && deta < 0.1 && dphi < 0.1) clus_layer.at(i_layer).push_back(1);
          if((i_layer < 4) && deta >= 0.1 && dphi >= 0.1) clus_outlier.at(i_layer).push_back(1);
        }//iClus loop
        int countClus = 0;
        int countOutlier = 0;
        for(int iLay = 0; iLay < 4; iLay++){
          if(clus_layer.at(iLay).size() > 1) countClus++;
          if(clus_outlier.at(iLay).size() > 0) countOutlier++;
        }
        if(countClus < 2){
          cout << "this event cant help because of cluster pattern" << endl;
          m_h_dR_1road_canthelp->Fill(deltaExtR_off);
          if(countClus < 1) m_h_dR_1road_1clus->Fill(deltaExtR_off);
          else {
            int countPhiClus = 0;
            for(unsigned int iClus_phi = 0; iClus_phi < nClus_phi; iClus_phi++){
              if(SARPCCluster_isPlausiblePhiInnMid->at(0).at(iClus_phi)){ countPhiClus++; }
            }
            if(countPhiClus < 2) m_h_dR_1road_1phi->Fill(deltaExtR_off);
            else m_h_dR_1road_phiothers->Fill(deltaExtR_off);
          }

        }
        else if(countClus > 1){
//          bool isShared = false;
//          vector<vector<int>> id_clus = SARPCCluster_id_clustersInSets->at(0);
//          for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
//            for(int jClus_fit = 0; jClus_fit < nClus_fit; jClus_fit++){
//              if(jClus_fit == iClus_fit) continue;
//              if(id_clus.at(0).at(iClus_fit) != id_clus.at(0).at(jClus_fit) && id_clus.at(0).at(iClus_fit) > -1 && id_clus.at(0).at(jClus_fit) > -1){
//                if(id_clus.at(1).at(iClus_fit) == id_clus.at(1).at(jClus_fit) && id_clus.at(1).at(iClus_fit) > -1) isShared = true;
//              }
//            }// jClus_fit loop end
//          }//iClus_fit loop end
//          if(isShared) m_h_dR_1road_outlier->Fill(deltaExtR_off);
//          else if(clus_layer.at(2).size() > 1 && (clus_layer.at(0).size() > 1 || clus_layer.at(1).size() > 1)){
//            int layer_base = (clus_layer.at(0).size() < 2) ? 1 : 0;
//            if(!(layer_base == 1 && clus_layer.at(1).size() < 2)){
//              int id_plau, plauRoadID;
//              bool isMisdirection = false;
//              for(int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
//                if(SARPCCluster_isPlausibleFitInnMid->at(0).at(iClus_fit)){
//                  id_plau = id_clus.at(layer_base).at(iClus_fit);
//                  plauRoadID = iClus_fit;
//                }
//              }
//              int id_notplau = -2;
//              for(int iClus_fit = 0; iClus_fit < nClus_fit;  iClus_fit++){
//                if(plauRoadID == iClus_fit) continue;
//                if(id_clus.at(layer_base).at(iClus_fit) != id_plau && id_clus.at(layer_base).at(iClus_fit) > -1) id_notplau = iClus_fit;
//              }
//              cout << "hoge, " << id_notplau << ", plau id = " << id_plau << endl;
//              if(id_notplau > -1){
//                if(id_clus.at(2).at(id_notplau) < 0) isMisdirection = true;
//                cout << "hoge" << endl;
//                int count = 0;
//                for(int ilay = 0; ilay < 4; ilay++){ if(id_clus.at(ilay).at(id_notplau > -1))count++; }
//                if(id_clus.at(2).at(plauRoadID) < 0 && count < 3){ isMisdirection = true; }
//                if(isMisdirection) m_h_dR_1road_others->Fill(deltaExtR_off); //motomoto misdire
//                else if(!isMisdirection) m_h_dR_1road_others->Fill(deltaExtR_off);
//              }
//              else if(id_notplau < 0) m_h_dR_1road_others->Fill(deltaExtR_off);
//            }
//          }//clus_layer if end
//          else{ m_h_dR_1road_others->Fill(deltaExtR_off); }
        }// countClus < 1 if end
      }
      else if(N_clusRoadMid > 2){
        m_h_dR_3road->Fill(deltaExtR_off);
        int nClus = SARPCCluster_clusterLayer->at(0).size();
        vector<vector<int>> clus_layer, clus_outlier;
        vector<int> forClear;
        forClear.clear();
        clus_layer.assign(4, forClear);
        clus_outlier.assign(4, forClear);
        for(int iClus = 0; iClus < nClus; iClus++){
          if(SARPCCluster_clusterMeasPhi->at(0).at(iClus)) continue;
          //======check where the test cluster is located======
          float r = sqrt(SARPCCluster_gX->at(0).at(iClus)*SARPCCluster_gX->at(0).at(iClus)+SARPCCluster_gY->at(0).at(iClus)*SARPCCluster_gY->at(0).at(iClus));
          float phi = atan(SARPCCluster_gY->at(0).at(iClus)/SARPCCluster_gX->at(0).at(iClus));
          cout << "phi = " << phi << endl;
          if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)>0) phi = phi + pi;
          if (SARPCCluster_gX->at(0).at(iClus)<0 && SARPCCluster_gY->at(0).at(iClus)<0) phi = phi - pi;
          cout << "phi renew = " << phi << endl;
          float l = sqrt(SARPCCluster_gZ->at(0).at(iClus)*SARPCCluster_gZ->at(0).at(iClus)+r*r);
          float tan = sqrt( (l-SARPCCluster_gZ->at(0).at(iClus))/(l+SARPCCluster_gZ->at(0).at(iClus)) );
          float eta = -log(tan);
          float deta = fabs(SARoIEta->at(0) - eta);
          float dphi = acos(cos( SARoIPhi->at(0) - phi ) );
          std::cout << "deta/dphi = " << deta << "/" << dphi << std::endl;
          //====================================================
          int i_layer = SARPCCluster_clusterLayer->at(0).at(iClus);
          if((i_layer < 4) && deta < 0.1 && dphi < 0.1) clus_layer.at(i_layer).push_back(1);
          if((i_layer < 4) && deta >= 0.1 && dphi >= 0.1) clus_outlier.at(i_layer).push_back(1);
        }//iClus loop
        int countClus = 0;
        int countOutlier = 0;
        for(int iLay = 0; iLay < 4; iLay++){
          if(clus_layer.at(iLay).size() > 2) countClus++;
          if(clus_outlier.at(iLay).size() > 0) countOutlier++;
        }
        if(countClus > 1) m_h_dR_3road_3clus->Fill(deltaExtR_off);
      }

      if(L1isMoreCandInRoI->at(1)){
        cout << "FillEfficiency::there are more than 2 candidates in 1RoI tower" << endl;
        m_h_deltaExtR_off_morecand->Fill(deltaExtR_off);
      }//isMoreCandInRoI
    }
    //for 1muon in 1RoI event

  }//isDimuon
//  if(isDimuon){
//    for(unsigned int iMuon = 0; iMuon < nMuon; iMuon++){
//      for(unsigned int jMuon = 0; jMuon < nMuon; jMuon++){
//        if(iMuon == jMuon) continue;
//        bool is2mu1RoI = false;
//        if(L1nRoI->at(iMuon) == 1 &&
//            L1RoINumber->at(iMuon) == L1RoINumber->at(jMuon) &&
//            L1RoISector->at(iMuon) == L1RoISector->at(jMuon) )
//        {
//          is2mu1RoI = true;
//          cout << "FillEfficiency::this is 2muon in 1RoI event" << endl;
//        }//L1 mismatching cut
//        double deltaExtEta = OfflineExtEta->at(iMuon) - OfflineExtEta->at(jMuon);
//        double deltaExtPhi = OfflineExtPhi->at(iMuon) - OfflineExtPhi->at(jMuon);
//        double deltaExtR_off = sqrt(pow(deltaExtEta, 2) + pow(deltaExtPhi, 2));
//        cout << "FillEfficiency::deltaExtR = " << deltaExtR_off << endl;
//        m_h_deltaExtR_off_tot->Fill(deltaExtR_off);
//        /* clusterRoad analysis start *************************************/
//        unsigned int nClus_fit = (SARPCCluster_fitMidSlope->at(jMuon)).size();
//        int N_clusRoadMid = 0;
//        int N_clusRoadOut = 0;
//        for(unsigned int iClus_fit = 0; iClus_fit < nClus_fit; iClus_fit++){
//          if(SARPCCluster_isPlausibleFitInnMid->at(jMuon).at(iClus_fit))N_clusRoadMid++;
//          if(SARPCCluster_isPlausibleFitOut->at(jMuon).at(iClus_fit))N_clusRoadOut++;
//        }
//        cout << "FillEfficiency::the number of clusterRoadMid = " << N_clusRoadMid << endl;
//        cout << "FillEfficiency::isMoreCand->at(" << iMuon << ") = " << L1isMoreCandInRoI->at(iMuon) << endl;
//        //for 2muon in 1RoI event (in order to protect the overlap filling)
//        if(is2mu1RoI && jMuon == 0){
//          m_h_deltaExtR_off_tot_clus->Fill(deltaExtR_off);
//          m_h_deltaExtR_off_closeBy->Fill(deltaExtR_off);
//          if(N_clusRoadMid == 2 /*|| N_clusRoadOut == 2*/)m_h_deltaExtR_off_2road_clus->Fill(deltaExtR_off);
//          if(L1isMoreCandInRoI->at(jMuon)){
//            cout << "FillEfficiency::there are more than 2 candidates in 1RoI tower" << endl;
//            m_h_deltaExtR_off_morecand->Fill(deltaExtR_off);
//          }//isMoreCandInRoI
//        }
//        //for 1muon in 1RoI event
//        else if(L1nRoI->at(jMuon) > 1){
//          m_h_deltaExtR_off_tot_clus->Fill(deltaExtR_off);
//          if(N_clusRoadMid == 2/* || N_clusRoadOut == 2*/) m_h_deltaExtR_off_2road_clus->Fill(deltaExtR_off);
//          if(L1isMoreCandInRoI->at(jMuon)){
//            cout << "FillEfficiency::there are more than 2 candidates in 1RoI tower" << endl;
//            m_h_deltaExtR_off_morecand->Fill(deltaExtR_off);
//          }//isMoreCandInRoI
//        }
//
//        /* clusterRoad analysis end * *************************************/
//        /* isMoreCandinRoI analysis start**********************************/
//        /*if(L1isMoreCandInRoI->size() > 0){
//          cout << "L1isMoreCand = " << L1isMoreCandInRoI->at(iMuon) << "/" << L1isMoreCandInRoI->at(jMuon) << endl;
//          if(L1isMoreCandInRoI->at(jMuon)){
//          cout << "FillEfficiency::there are more than 2 candidates in 1RoI tower" << endl;
//          m_h_deltaExtR_off_morecand->Fill(deltaExtR_off);
//          }//isMoreCandInRoI
//          }//L1isMoreCandInRoI size*/
//        /* isMoreCandInRoI analysis end ***********************************/
//
//
//      }//jMuon loop end
//    }//iMuon loop end
//  }//isDimuon
 std::cout << "FillEfficiency end " << std::endl;
}//FillEfficiency

int RPC_FCBM::EtaRegion( double eta )
{
  int region = -1;
  if( fabs(eta) < 1.05 ){
    region = 0;
  } else if( fabs(eta) < 1.5 ){
    region = 1;
  } else if( fabs(eta) < 2.0 ){
    region = 2;
  } else {
    region = 3;
  }
  return region;
}


int RPC_FCBM::ECWeakRegion( double eta, double phi )
{
  int region = 0;
  double absEta = fabs( eta );
  double absPhi = fabs( phi );
  if(  ( 1.3 <= absEta && absEta < 1.45) &&
      ( (0  <= fabs(phi) && fabs(phi) < TMath::Pi()/48. )     ||
        (TMath::Pi()*11./48. <= fabs(phi) && fabs(phi) < TMath::Pi()*13./48. ) ||
        (TMath::Pi()*23./48. <= fabs(phi) && fabs(phi) < TMath::Pi()*25./48. ) ||
        (TMath::Pi()*35./48. <= fabs(phi) && fabs(phi) < TMath::Pi()*37./48. ) ||
        (TMath::Pi()*47./48. <= fabs(phi) && fabs(phi) < TMath::Pi() )
      )
    ) region = 1;
  else if( ( 1.5 <= absEta && absEta < 1.65 ) &&
      ( (TMath::Pi()*3./32.  <= fabs(phi) && fabs(phi) < TMath::Pi()*5./32. ) ||
        (TMath::Pi()*11./32. <= fabs(phi) && fabs(phi) < TMath::Pi()*13./32.) ||
        (TMath::Pi()*19./32. <= fabs(phi) && fabs(phi) < TMath::Pi()*21./32.) ||
        (TMath::Pi()*27./32. <= fabs(phi) && fabs(phi) < TMath::Pi()*29./32.)
      )
      ) region = 2;
  return region;
}
