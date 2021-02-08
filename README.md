# qsar
Simple ML model for performing QSAR

### Usage

```code
Usage:
qsar_model.py train --in INFILE_NAME --model MODEL_NAME [--log] [--units UNITS]
qsar_model.py predict --in INFILE_NAME --model MODEL_NAME --out OUTFILE_NAME [--unlog] [--units UNITS]

--in INFILE_NAME  input SMILES file name, has activity for training
--out OUTFILE_NAME  output csv file
--model MODEL_NAME  model file name
--log  convert input activity to log scale
--unlog  convert output data from log scale to uM
--units UNITS units (uM or nM) default is uM
```
### Trainining

```
Training on pIC50 data
./qsar_model.py train --in data/COX-2_train_pIC50.smi --model model.model

Training on IC50 data in uM - the "--log" flag converts data to pIC50
./qsar_model.py train --in data/COX-2_train_IC50_uM.smi --model model.model --log

Training on IC50 data in uM - the "--log" flag converts data to pIC50 the "--units" data sets unit to nM rather than the default, uM
./qsar_model.py train --in data/COX-2_train_IC50_nM.smi --model model.model --log --units nM
```

### Predicting
```
Predicting pIC50 values
./qsar_model.py predict --in data/COX-2_test.smi --out pred.csv --model model.model

Predicting IC50 values in uM
./qsar_model.py predict --in data/COX-2_test.smi --out pred.csv --model model.model --unlog

Predicting IC50 values in nM
./qsar_model.py predict --in data/COX-2_test.smi --out pred.csv --model model.model --unlog --units nM
```



```

