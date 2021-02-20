#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from xgboost import XGBRegressor, XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, roc_auc_score
from docopt import docopt
import os

cmd_str = """Usage: 
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
"""


def read_data(infile_name, is_pred=False):
    if not os.path.exists(infile_name):
        print(f"Error: {infile_name} not found")
        sys.exit(0)

    cols_names = ["SMILES", "Name", "Activity"]
    if is_pred:
        cols_names = cols_names[0:2]
    df = pd.read_csv(infile_name, names=cols_names, sep=" ")
    print(f"Read {df.shape[0]} records from {infile_name}", file=sys.stderr)
    return df


def calc_fps(df):
    mol_list = []
    for smi in df.SMILES:
        mol = Chem.MolFromSmiles(smi)
        mol_list.append(mol)
    df['mol'] = mol_list
    num_mol_in = df.shape[0]
    df.dropna(subset=['mol'], inplace=True)
    num_mol_out = df.shape[0]
    df['fp'] = [AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in tqdm(df.mol)]
    print(f"Generated fingerprints for {num_mol_out} of {num_mol_in} molecules", file=sys.stderr)


def get_exponent(unit_name):
    exponent_dict = {"uM": -6, "nM": -9}
    exponent = exponent_dict.get(unit_name)
    if exponent is None:
        print(f"{unit_name} is not supported as a unit identifier", file=sys.stderr)
        sys.exit(1)
    return exponent


def train_classification_model(infile_name, model_name, cutoff, flip=False, num_folds=10):
    df = read_data(infile_name)
    if flip:
        print(f"values >= {cutoff} will be assigned a value of 1")
        df['act_class'] = [1 if x >= cutoff else 0 for x in df.Activity]
    else:
        print(f"values <= {cutoff} will be assigned a value of 1")
        df['act_class'] = [1 if x <= cutoff else 0 for x in df.Activity]
    calc_fps(df)
    X = np.asarray(list(df.fp.values))
    y = df.act_class.values
    auc_list = []
    for i in range(0, num_folds):
        X_train, X_test, y_train, y_test = train_test_split(X, y)
        xgb = XGBClassifier()
        xgb.fit(X_train, y_train)
        y_pred = xgb.predict(X_test)
        auc = roc_auc_score(y_test, y_pred)
        auc_list.append(auc)
        print(f"Fold {i + 1:2d} AUC {auc:.2f}", file=sys.stderr)
    print(f"Mean AUC {np.mean(auc_list):.2f}", file=sys.stderr)

    xgb = XGBClassifier()
    xgb.fit(X, y)
    xgb.save_model(model_name)
    print(f"Model saved to {model_name}", file=sys.stderr)


def train_regression_model(infile_name, model_name, num_folds=10, log_scale=False, units="uM"):
    df = read_data(infile_name)
    if log_scale:
        exponent = get_exponent(units)
        df.Activity = -np.log10(df.Activity * pow(10, exponent))
        print(f"Converting data to log scale with exponent {exponent} ", file=sys.stderr)
    calc_fps(df)
    X = np.asarray(list(df.fp.values))
    y = df.Activity.values
    r2_list = []
    for i in range(0, num_folds):
        X_train, X_test, y_train, y_test = train_test_split(X, y)
        xgb = XGBRegressor()
        xgb.fit(X_train, y_train)
        y_pred = xgb.predict(X_test)
        r2 = r2_score(y_test, y_pred)
        r2_list.append(r2)
        print(f"Fold {i + 1:2d} R**2 {r2:.2f}", file=sys.stderr)
    print(f"Mean R**2 {np.mean(r2_list):.2f}",file=sys.stderr)

    xgb = XGBRegressor()
    xgb.fit(X, y)
    xgb.save_model(model_name)
    print(f"Model saved to {model_name}",file=sys.stderr)


def predict_classification(infile_name, model_name, outfile_name):
    df = read_data(infile_name)
    calc_fps(df)
    X = np.asarray(list(df.fp.values))
    xgb = XGBClassifier()
    xgb.load_model(model_name)
    print(f"Loaded model from {model_name}")
    y_pred = xgb.predict(X)
    pred_col = "Pred_Class"
    df[pred_col] = y_pred
    df[["SMILES", "Name", pred_col]].to_csv(outfile_name, index=False, float_format="%0.2f")
    print(f"Predictions written to {outfile_name}")


def predict_regression(infile_name, model_name, outfile_name, un_log=False, units="uM"):
    df = read_data(infile_name)
    calc_fps(df)
    X = np.asarray(list(df.fp.values))
    xgb = XGBRegressor()
    xgb.load_model(model_name)
    print(f"Loaded model from {model_name}")
    y_pred = xgb.predict(X)
    pred_col = "Pred_pIC50"
    if un_log:
        exponent = get_exponent(units)
        pred_col = "Pred_" + units
        print(f"Converting from log scale with exponent {exponent}")
        y_pred = [pow(10, -x) * pow(10, -exponent) for x in y_pred]
    df[pred_col] = y_pred
    df[["SMILES", "Name", pred_col]].to_csv(outfile_name, index=False, float_format="%0.2f")
    print(f"Predictions written to {outfile_name}")


def main(cmd_input_str):
    cmd_input = docopt(cmd_input_str)
    infile_name = cmd_input.get("--in")
    model_name = cmd_input.get("--model")
    # -- Regression params
    units = cmd_input.get("--units") or "uM"
    log_scale = cmd_input.get("--log")
    un_log = cmd_input.get("--unlog")
    # -- Classification params
    cutoff = cmd_input.get("--cutoff")
    if cutoff:
        cutoff = float(cutoff)
    flip = bool(cmd_input.get("--flip"))
    units = units or "uM"
    if cmd_input.get("train"):
        if cmd_input.get("regression"):
            train_regression_model(infile_name, model_name, log_scale=log_scale, units=units)
        else:
            train_classification_model(infile_name, model_name, cutoff=cutoff, flip=flip)
    if cmd_input.get("predict"):
        if not os.path.exists(model_name):
            print(f"Error: Could not open model {model_name}", file=sys.stderr)
            sys.exit(0)
        outfile_name = cmd_input.get("--out")
        if cmd_input.get("regression"):
            predict_regression(infile_name, model_name, outfile_name, un_log=un_log, units=units)
        else:
            predict_classification(infile_name, model_name, outfile_name)


if __name__ == "__main__":
    main(cmd_str)
