import argparse
import chemprop
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from rdkit import Chem
from torch import Tensor
import torch
import numpy as np
from  sirms.files import LoadFragments
from collections import OrderedDict
import os
import optimizer_utils
import re

mol_frag_sep = "###"


def predict_mol(m, model, sclr,i, model_type, frags=None, per_atom_fragments=None, id_field_name=None):
    def pred(mol, model):
        # remove hs for correct work of mpnn
        mol =  Chem.RemoveHs(mol)
        with torch.no_grad():
            pred=model.forward([[mol]])
            if model_type == "reg":
                pred = sclr.inverse_transform(pred)

            elif model_type == "class":
                pred = np.asarray(sigmoid(pred))

        return np.array(pred).squeeze(0)

    mol_dict = OrderedDict()
    if id_field_name is not None:
        nm = m.GetProp(id_field_name)
    elif m.GetProp("_Name") == "":
        nm = 'auto_generated_id_' + str(i + 1)  # 1-based as in sirms.py
    else:
        nm = m.GetProp("_Name")
    res = pred(m, model)
    mol_dict[nm] = res
    if per_atom_fragments:
        counter = 0
        for idx in range(m.GetNumAtoms()):
            if m.GetAtomWithIdx(idx).GetAtomicNum() > 1:
                rw_m = Chem.RWMol(m)
                rw_m.GetAtoms()[idx].SetAtomicNum(0)
                mol_dict[nm + mol_frag_sep + str(idx + 1) + "#" + str(counter)] = pred(rw_m, model)
                counter += 1
    elif frags and nm in frags:
        for k, v in frags[nm].items():
            rw_m = Chem.RWMol(m)
            for idx in sorted(v, reverse=True):  # note we don't check if atom== H (is it ok?)
                rw_m.GetAtoms()[idx - 1].SetAtomicNum(0)
            mol_dict[nm + mol_frag_sep + k] = pred(rw_m, model)
    return mol_dict

def main_params(in_fname, out_fname, model_path, model_type,  multitask, variance_threshold=None,
                frag_fname=None,
                per_atom_fragments=False,
                id_field_name=None,
                save_pred=True, num_frag_id=False):
    # load model
    arguments = [
        '--test_path', '/dev/null',
        '--preds_path', '/dev/null',
        '--checkpoint_dir', model_path
    ]
    args = chemprop.args.PredictArgs().parse_args(arguments)
    model_objects = chemprop.train.load_model(args=args)
    _, __, models, scalers, ___, prop_names = chemprop.train.load_model(args=args)
    sclrs = [i[0] for i in scalers]
    # load sdf and get dict of preds (like sirms dict)
    input_file_extension = in_fname.strip().split(".")[-1].lower()
    if input_file_extension == 'sdf':
        frags = LoadFragments(frag_fname)
        df_lst  = []
        for k,(mod, sclr) in enumerate(zip(models, sclrs)):
            mod.eval()
            mols = None
            mols = OrderedDict()  # key - molname, val- prediction; if frags: key - molname or mol+fragname, val-pred for mol or part b
            for i, m in enumerate(Chem.SDMolSupplier(in_fname, removeHs=False)):

                if m is not None:


                    res=predict_mol(m, mod,sclr, i, model_type, frags=frags,
                                         per_atom_fragments=per_atom_fragments,
                                         id_field_name=id_field_name)

                    mols.update(res)
                #construct df and save to file

            print(mols)

            df = pd.DataFrame.from_dict(mols, orient="index", columns=prop_names)
            df_lst.append(df)
        if variance_threshold is not None:
            df_lst = pd.concat(df_lst).groupby(level=0).agg(['mean', 'var'])
            print(df_lst)
            # Identify all variance columns
            var_cols = [col for col in df_lst.columns if col[1] == 'var']

            # Convert all variance values to boolean based on threshold
            for col in var_cols:
                df_lst[col] = (df_lst[col] <= variance_threshold).astype(int)

            # Flatten multi-level columns and add suffixes
            df_lst.columns = [
                f"{col[0]}" if col[1] == "mean" else f"{col[0]}_bound_box" for col in df_lst.columns
            ]
        else:
            df_lst = pd.concat(df_lst).groupby(level=0).mean()
            print(df_lst)
        df_lst = df_lst.reset_index()
        df_lst = df_lst.rename(columns={'index': 'Compounds'})
        if frags:
            df_lst[['Compounds', 'Fragment']] = df_lst['Compounds'].str.split(mol_frag_sep, expand=True)
            if num_frag_id:  # separate  piece after last # - numerical fragment id
                df_lst[['Fragment', 'Frag_id']] = df_lst.Fragment.str.rsplit("#", n=1, expand=True)
        if multitask:
                out_fname_param = optimizer_utils.retrieve_out_fname_param(out_fname)
                outs_list = []
                for prop_name in prop_names:
                    # todo fix _activity - itss hardcode dependent on chemprop naming
                    out_fname_tmp = os.path.join(os.path.dirname(out_fname),
                                                 re.sub(out_fname_param, re.sub("activity_", "", prop_name),
                                                        os.path.basename(
                                                            out_fname)))  # replace old param name with actually to be written

                    outs_tmp = optimizer_utils.filter_columns_by_keyword(df_lst, prop_name)
                    outs_tmp["Compounds"] = df_lst["Compounds"]
                    outs_tmp = outs_tmp[[outs_tmp.columns[-1]] + list(outs_tmp.columns[:-1])]
                    if frags is not None:
                        outs_tmp["Fragment"] = df_lst["Fragment"]
                        outs_tmp = outs_tmp[
                            [outs_tmp.columns[0]] + [outs_tmp.columns[-1]] + list(outs_tmp.columns[1:-1])]

                        if num_frag_id:  # separate  piece aft
                            outs_tmp["Frag_id"] = df_lst["Frag_id"]
                            outs_tmp = outs_tmp[list(outs_tmp.columns[0:2]) + [outs_tmp.columns[-1]] + list(outs_tmp.columns[2:-1])]

                    outs_tmp["consensus"] = outs_tmp[prop_name]

                    outs_tmp.columns = ["bound_box" if "bound_box" in col else col for col in outs_tmp.columns]
                    if "bound_box" not in  outs_tmp.columns: outs_tmp["bound_box"] = 1 # add fake bb for downstream compatibility

                    outs_list.append(outs_tmp)
                    if save_pred:
                        outs_tmp.to_csv(out_fname_tmp, sep="\t", index=False)
                return outs_list
        else:
            df_lst["consensus"] = df_lst[prop_names[0]]  # there should be only 1 property
            df_lst.columns = ["bound_box" if "bound_box" in col else col for col in df_lst.columns]
            if "bound_box" not in outs_tmp.columns: outs_tmp[
                "bound_box"] = 1  # add fake bb for downstream compatibility

            if save_pred:
                df_lst.to_csv(out_fname, sep="\t", index=False)
            return df_lst

    else:
        print("Input file extension should be SDF Current file has %s. Please check it." %
              input_file_extension.upper())
        return None


def entry_point():
    parser = argparse.ArgumentParser(description='predict molecules')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input file ( sdf with standardized structures')
    parser.add_argument('-o', '--out',  metavar='param_pred.txt', required=True,
                        help='output file with predictions, tsv.')
    parser.add_argument('--model_path',
                        help='Path to model ')
    parser.add_argument('-t', '--model_type', metavar='reg/class', required=True,
                        help='')
    parser.add_argument('-f', '--fragments', metavar='fragments.txt', default=None,
                        help='text file containing list of names of single compounds, fragment names and atom '
                             'indexes of fragment to remove (all values are tab-separated).')
    parser.add_argument('--per_atom_fragments', action='store_true', default=False,
                        help='if set this flag input fragments will be omitted and single atoms will be considered '
                             'as fragments.')
    parser.add_argument('-w', '--id_field_name', metavar='field_name', default=None,
                        help='field name of unique ID for compounds. '
                             'If omitted for sdf molecule titles will be used or auto-generated names')
    parser.add_argument('-m', '--multitask', metavar='True/False', required=True,
                        help='')
    parser.add_argument('--variance_threshold', metavar='', required=False, default=0,
                        help='variance_threshold variance for applicability domain')
    parser.add_argument('--num_frag_id', metavar='', required=False, default=False,
                        help='numerical frag id')
    args = vars(parser.parse_args())

    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "model_path": model_path = v
        if o == "model_type": model_type = v
        if o == "fragments": frag_fname = v
        if o == "per_atom_fragments": per_atom_fragments = v
        if o == "id_field_name": id_field_name = v
        if o == "multitask": multitask = v
        if o == "variance_threshold": variance_threshold = float(v)
        if o == "num_frag_id": num_frag_id = bool(v)

    main_params(in_fname=in_fname, out_fname=out_fname, model_path=model_path, model_type=model_type,  frag_fname=frag_fname,
                per_atom_fragments=per_atom_fragments, id_field_name=id_field_name, multitask=multitask, variance_threshold=variance_threshold, num_frag_id=num_frag_id)


if __name__ == '__main__':
    entry_point()
