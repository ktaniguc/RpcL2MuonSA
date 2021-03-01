#include <TLegend.h>
#include "../RpcL2MuonSA/RPC_FCBM.h"

using namespace std;

void RPC_FCBM::setLegend_eventInfo(TText* eventInfo)
{
  *eventInfo = TText(0.05,0.02,Form("EventNumber = %llu, RunNumber = %d, LumiBlock = %d, AverageInteractionsPerCrossing = %5.3f",EventNumber, RunNumber, LumiBlock, AverageInteractionsPerCrossing));
  eventInfo->SetNDC();
  eventInfo->SetTextSize(0.03);
}

void RPC_FCBM::setLegend_leftSideInfo(TLegend* OffMu1Info, string& plane){
  //This function show the information about Offline muon1
  TString label_for_sector;
  if (SAsAddress->at(0) == 0){
    label_for_sector = "Large";
  } else if (SAsAddress->at(0) == 1){
    label_for_sector = "Small";
  } else if (SAsAddress->at(0) == 2){
    label_for_sector = "Large Special";
  } else if (SAsAddress->at(0) == 3){
    label_for_sector = "Small Special";
  } else if (SAsAddress->at(0) == -1){
    label_for_sector = "Endcap";
  }
  OffMu1Info->SetTextSize(0.035);
  OffMu1Info->SetHeader("OfflineMuon1","C");
  if(plane.compare("RPC1") == 0){
    OffMu1Info->AddEntry((TObject*)0, Form("#splitline{RPC1 plane}{%s}", label_for_sector.Data()),"C");
  }
  else if(plane.compare("RPC2") == 0){
    OffMu1Info->AddEntry((TObject*)0, Form("#splitline{RPC2 plane}{%s}", label_for_sector.Data()),"C");
  }
  else if(plane.compare("RPC3") == 0){
    OffMu1Info->AddEntry((TObject*)0, Form("#splitline{RPC3 plane}{%s}", label_for_sector.Data()),"C");
  }
  OffMu1Info->AddEntry((TObject*)0,Form("#splitline{Offline p_{T}}{ : %4.3f [GeV]}", OfflinePt->at(0)),"");
  OffMu1Info->AddEntry((TObject*)0,Form("#splitline{L2SA p_{T}}{ : %4.3f [GeV]}",abs(SAPt->at(0))),"");

}


void RPC_FCBM::setLegend_rightSideInfo(TLegend* OffMu2Info, string& plane){
  //This function show the information about Offline muon2
  TString label_for_sector;
  if (SAsAddress->at(1) == 0){
    label_for_sector = "Large";
  } else if (SAsAddress->at(1) == 1){
    label_for_sector = "Small";
  } else if (SAsAddress->at(1) == 2){
    label_for_sector = "Large Special";
  } else if (SAsAddress->at(1) == 3){
    label_for_sector = "Small Special";
  } else if (SAsAddress->at(1) == -1){
    label_for_sector = "Endcap";
  }
  OffMu2Info->SetTextSize(0.035);
  if(nMuon == 2){
    OffMu2Info->SetHeader("OfflineMuon2","C");
    if(plane.compare("RPC1") == 0){
      OffMu2Info->AddEntry((TObject*)0, Form("#splitline{RPC1 plane}{%s}", label_for_sector.Data()),"C");
    }
    else if(plane.compare("RPC2") == 0){
      OffMu2Info->AddEntry((TObject*)0, Form("#splitline{RPC2 plane}{%s}", label_for_sector.Data()),"C");
    }
    else if(plane.compare("RPC3") == 0){
      OffMu2Info->AddEntry((TObject*)0, Form("#splitline{RPC3 plane}{%s}", label_for_sector.Data()),"C");
    }
    OffMu2Info->AddEntry((TObject*)0,Form("#splitline{Offline p_{T}}{ : %4.3f [GeV]}", OfflinePt->at(1)),"");
    for(int i= 0; i< SAptclus->at(0).size(); i++){
      OffMu2Info->AddEntry((TObject*)0,Form("#splitline{SAclus p_{T}}{ : %4.3f [GeV]}",SAptclus->at(0).at(i)),"");
    }
  }
  else if(nMuon == 1){
    OffMu2Info->AddEntry((TObject*)0, "No content" ,"");
  }
  else if(nMuon >= 3){
    OffMu2Info->AddEntry((TObject*)0, "#splitline{more than 3 muons}{space overflaw} ");
  }

}

