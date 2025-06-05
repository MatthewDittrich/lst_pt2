# Python PT2 Processor


## Setup

Clone the repository:

```
git clone https://github.com/MatthewDittrich/lst_pt2.git
```

Make the conda environment:

```
cd lst_pt2/pt2_python/

conda env create -f environment.yml
```

## Input

Input for the processor is the LST OD ntuple. <br>
Example output has been stored in the location below for you to use:

```
/blue/avery/matthew.dittrich/PT2_Studies/LSTNtuple.root
```

## Running the code and producing plots

Assuming you are running on UF's Hipergator:

```
srun -t 600 --qos=avery --account=avery --cpus-per-task=32 --mem=128gb --pty bash -i
```

Activate conda environment:

```
conda activate lst_env
```
