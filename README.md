# qsar
Simple ML model for performing QSAR

### Usage

```code
Usage: 
qsar_model.py train classification --in INFILE_NAME --model MODEL_NAME --cutoff CUTOFF [--flip]
qsar_model.py train regression --in INFILE_NAME --model MODEL_NAME [--log] [--units UNITS]
qsar_model.py predict classification --in INFILE_NAME --model MODEL_NAME --out OUTFILE_NAME 
qsar_model.py predict regression --in INFILE_NAME --model MODEL_NAME --out OUTFILE_NAME [--unlog] [--units UNITS]

--in INFILE_NAME  input SMILES file name, has activity for training
--out OUTFILE_NAME  output csv file
--model MODEL_NAME  model file name
--cutoff CUTOFF classification cutoff, values <= cutoff will be labeled as 1 unless flip is specified
--flip  consider to larger values to be more active, a value of 1 will be assigned to values >= cutoff
--log  convert input activity to log scale
--unlog  convert output data from log scale to uM
--units UNITS units (uM or nM) default is uM
```
### Classification Models
#### Training
```
Training on IC50 data, lower values better 
qsar_model.py train classification --in data/COX-2_train_IC50_uM.smi --model model.model --cutoff 1

Training on pIC50 values, higher values better, note the use of the "--flip" flag 
qsar_model.py train classification --in data/COX-2_train_pIC50.smi --model model.model --cutoff 6 --flip
```
#### Prediction
```
predict classification --in data/COX-2_test.smi --model model.model --out classification_out.csv
```
### Regression Models
#### Training
```
Training on pIC50 data
qsar_model.py train regression --in data/COX-2_train_pIC50.smi --model model.model

Training on IC50 data in uM - the "--log" flag converts data to pIC50
qsar_model.py train regression --in data/COX-2_train_IC50_uM.smi --model model.model --log

Training on IC50 data in uM - the "--log" flag converts data to pIC50 the "--units" data sets unit to nM rather than the default, uM
qsar_model.py train regresssion --in data/COX-2_train_IC50_nM.smi --model model.model --log --units nM
```
#### Prediction
```
Predicting pIC50 values
qsar_model.py predict regression --in data/COX-2_test.smi --out pred.csv --model model.model

Predicting IC50 values in uM
qsar_model.py predict regression --in data/COX-2_test.smi --out pred.csv --model model.model --unlog

Predicting IC50 values in nM
qsar_model.py predict regression --in data/COX-2_test.smi --out pred.csv --model model.model --unlog --units nM
```



```

