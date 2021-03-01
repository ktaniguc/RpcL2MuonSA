#!/bin/sh

if [ "$1" = "" ]; then
  DATE=$(date '+%Y%m%d%H%M')
else
  DATE="$1"
fi
#PDF_LABEL="rpcCluster_JpsiCloseBy_${DATE}"
#PDF_LABEL="rpcCluster_JpsiCloseBy_oldL2SA"
PDF_LABEL="forMthesis_test"
#PDF_LABEL="plot_forMthesis"
INPUT_NTUPLE="/gpfs/fs7001/ktaniguc/outputfile/OutputCalcEff/RpcCluster/JpsiCollimated_newL2SA_separateMDTs_addetaMS.root"  #new Algorithm     #use for mthesis display

IS_DRAW="false"
IS_EVENTDISPLAY="true"
BEGIN_ENTRY=0
LIMIT_ENTRY=100
TAP_TYPE=3
TRIG_CHAIN=14
IS_CloseByMuon="true"

#IS_DRAW="true"
#IS_EVENTDISPLAY="false"
#BEGIN_ENTRY=0
#LIMIT_ENTRY=-1
##LIMIT_ENTRY=800
#TAP_TYPE=3
#TRIG_CHAIN=0
##IS_CloseByMuon="false"
#IS_CloseByMuon="true"

echo ""
echo "PDF_LABEL: "${PDF_LABEL}
echo "INPUT_NTUPLE: "${INPUT_NTUPLE}
echo ""
mkdir -p ./outroot/${DATE}
mkdir -p ./plot/${DATE}
COMMAND="${setupDir}/build/RpcL2MuonSA.out ${PDF_LABEL} ${INPUT_NTUPLE} ${IS_DRAW} ${IS_EVENTDISPLAY} ${BEGIN_ENTRY} ${LIMIT_ENTRY} ${TAP_TYPE} ${TRIG_CHAIN} ${IS_CloseByMuon} ${DATE}"
LOG="LOGS/log_"${PDF_LABEL}

eval ${COMMAND} 2>&1 | tee ${LOG}

#eof

