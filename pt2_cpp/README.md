# C++ PT2 Processor

## Setup
Cloning the repository:
```
git clone https://github.com/MatthewDittrich/lst_pt2.git
```
Activate ROOT:
```
cd lst_cpp/
export SCRAM_ARCH=el8_amd64_gcc10
export CMSSW_VERSION=CMSSW_12_5_0
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmssw-el8 --nv
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd - > /dev/null
```
## Input
As the Python version, this processor also uses LST OD Ntuple. Some examples are stored in:
```
/blue/p.chang/aaponteutani/LSTNtuple/
```
In that file **mG** and **PU200** refer to muon gun and 200 pileup, respectively.

## Running the code
If you are working on HiperGator:
```
srun --nodes=1 --ntasks=1 --cpus-per-task=32 --gpus=1 --time=04:00:00 --partition=hpg-b200 --mem=120G --qos=avery --account=avery --pty bash
```
Considering that ROOT is already activated, just run the command:
```
root -l -q 'main_analysis.C'
```
Don't forget to change the location of the input inside **main_analysis.C**.
