""" Create a summary of variables/mip/experiment requirements of recipes

"""

import argparse
import yaml
import os
import datetime as dt
from esmvaltool.preprocessor._derive import get_required


def getfilename(p):
    if not os.path.isfile(p):
        raise Exception("getfilename", "Not a file {0}".format(p))

    return p.split("/")[-1].split('.')[0]


def get_variable_list(infile):
    l = []
    # test if valid path
    if not os.path.isfile(infile):
        raise Exception("getfilename", "Not a file {0}".format(p))
    # load namelist
    with open(infile, 'r') as stream:
        try:
            po = yaml.load(stream)
        except:
            print("Error while parsing {0}. Continue".format(infile))
            raise Exception("getfilename", "Not a file {0}".format(p))
        datasetMip = []
        datasetExp = []
        if "datasets" in po.keys():
            for dataset in po["datasets"]:
                if "mip" in dataset.keys():
                    datasetMip.append(dataset["mip"])
                if "exp" in dataset.keys():
                    datasetExp.append(dataset["exp"])
        datasetMip = list(set(datasetMip))
        datasetExp = list(set(datasetExp))

        # list diagnostic variables (varname, mip, experiment)
        for k, v in po.get('diagnostics').items():

            exp = []
            mip = []

            if datasetMip != []:
                mip += datasetMip
            if datasetExp != []:
                exp += datasetExp

            if 'additional_datasets' in v.keys():
                for l1 in v.get("additional_datasets"):
                    if 'exp' in l1.keys():
                        exp.append(l1['exp'])
                    if 'mip' in l1.keys():
                        mip.append(l1['mip'])
            exp = list(set(exp))
            mip = list(set(mip))

            if 'variables' in v.keys():
                for k1, v1 in v.get("variables").items():
                    if exp == []:
                        exp = [v1.get("exp")]

                    if mip == []:
                        mip = [v1.get("mip")]

                    if v1.get("derive"):
                        for req in get_required(k1, v1.get("field")):
                            for m in mip:
                                for e in exp:
                                    l.append((req[0], m, e))
                    else:
                        for m in mip:
                            for e in exp:
                                l.append((k1, m, e))

    return list(set(l))


def write_csv_file(d, outfile="./out.csv"):
    with open(outfile, "w") as f:
        f.write("This file was created automatically by the ESMValTool")
        f.write("at {0}".format(dt.datetime.now().isoformat()))
        f.write("variable, mip_table, experiment")
        for k, v in d.items():
            pass


def _get_index_of_common_var_mip(t, l):
    idx = 0
    out = list()
    for i in l:
        if i[0:1] == t[0:1]:
            out.append(idx)
        idx += 1
    return out


def _join_experiments(l):
    newlist = list()
    for i in l:
        ci = _get_index_of_common_var_mip(i, l)
        joined_exp = tuple([l[item][2] for item in ci])
        newlist.append((i[0], i[1], joined_exp))
    return list(set(newlist))


def main():
    parser = argparse.ArgumentParser(
        description='Summarizes variable requirements by recipes')
    parser.add_argument(
        'infiles', type=str, nargs='+', help='path to yml recipes')
    parser.add_argument(
        '--uniq',
        action='store_true',
        help='get unique list of variables over all given recipes')
    parser.add_argument(
        '--join',
        action='store_true',
        help='display experiments joined by variable and mip')
    parser.add_argument('--csv', action='store_true', help='output as csv')

    args = parser.parse_args()
    infiles = args.infiles

    if args.uniq:
        out = list()
        for infile in infiles:
            out += get_variable_list(infile)
        out = list(set(out))
        if args.join:
            out = _join_experiments(out)
    else:
        out = dict()
        for infile in infiles:
            if args.join:
                out[getfilename(infile)] = _join_experiments(
                    get_variable_list(infile))
            else:
                out[getfilename(infile)] = get_variable_list(infile)

    print(yaml.dump(out))


if __name__ == '__main__':
    main()
