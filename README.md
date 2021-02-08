# qsar
Simple ML model for performing QSAR

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

