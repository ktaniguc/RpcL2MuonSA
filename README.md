To analyze L2MuonSA
- This can be work for CalcEffTool and its specific branch made by Kazuki Yamashita
- the codes in this package are too messy. if you want to know only how to draw EventDisplay, then please check only "RPC_FCBM.cxx"
##setup, compile and run
```sh
source setup.sh
./compile.sh cmake
cd run
# edit run.sh
./run.sh output_file_name #if you added nothing to 1st argument, then output file name will be the time when you run
```
##the meaning of each arguments in run.sh
------------------------------
|arguments name | description |
|:-------------:|:----------:|
|PDF_LABEL|the name of output pdf and root file|
|INPUT_NTUPLE|the input file name|
|IS_DRAW|if true, write the histogram which you draw|
|IS_EVENTDISPLAY|if true, write event display of offine and trigger objects on R-Z plane and middle station of muon spectrometer|
|BEGIN_ENTRY| the entry number started to run|
|LIMIT_ENTRY| the number of events to run from BEGIN_ENTRY|
|TAP_TYPE| the type of tag-and-probe. this explanation in detail is written in README of CalcEffTool. it's not needed to change|
|TRIG_CHAIN| the trigger chain which you want to analyze. this explanation is also in CalcEffTool's README|
|IS_CloseByMuon| if true, switch the mode to read the specific branch of CalcEffTool(made by kayamash)|

