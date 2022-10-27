import importlib
from argparse import ArgumentParser
import os

def PrintData(Model, logfp):
    keyword_list = list()
    label_list = list()
    DATA_PREFIX = 'Data_'

    for v in sorted(Model.__dict__.keys()):
        if v.startswith(DATA_PREFIX):
            keyword_list.append(v)

    tv = DATA_PREFIX + 'Time'
    if tv in keyword_list:
        i = keyword_list.index(tv)
        keyword_list = [tv] + keyword_list[:i] + keyword_list[i+1:]

    label_list = [v[len(DATA_PREFIX):] for v in keyword_list]

    value_list = [Model.__dict__[k] for k in keyword_list]

    # print header line
    print('\t'.join(label_list), file=logfp)

    # print data points
    for v in zip(*value_list):
        print('\t'.join([str(i) for i in v]), file=logfp)


def RunModel(Model):
    Model.Run()
    return


def LoadModel(modname):
    mod = importlib.import_module(modname)
    if "FModel" in mod.__dict__:
        model = mod.FModel()
    elif "FNetwork" in mod.__dict__:
        model = mod.FNetwork()
    else:
        raise AttributeError()

    return model


def main(args):

    if args.mod_files:
        for fname in args.mod_files:

            modname = os.path.splitext(fname)[0]
            model = LoadModel(modname)
            RunModel(model)

            logfilename = modname + ".out"
            if args.out_dir:
                if not os.path.exists(args.out_dir):
                    os.makedirs(args.out_dir);

                logfilename = os.path.join(args.out_dir, logfilename)

            fp = open(logfilename, "w")
            PrintData(model, fp)
            fp.close()

    return


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-o', '--out-dir', dest='out_dir', type=str,
            help="Write output files to OUT_DIR")
    parser.add_argument('mod_files', nargs='+', help='test filenames')

    args = parser.parse_args()
    main(args)

