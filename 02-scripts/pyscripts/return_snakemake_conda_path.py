#!/usr/bin/env python
########################################################################
# Return the path in which Snakemake would install a conda environment
#
# Alex Huebner, 25/04/23
########################################################################

import argparse
import hashlib
import os


def main():
    ''' Return the path in which Snakemake would install a conda environment
    '''
    md5hash = hashlib.md5()
    if Args['path'].startswith("/"):
        path = Args['path']
    else:
        path = os.getcwd() + "/" + Args['path']
    md5hash.update(path.encode())
    with open(Args['envfile'], "rb") as f:
        md5hash.update(f.read())
    print(md5hash.hexdigest() + "_")


# Argument parser
Parser = argparse.ArgumentParser(description='Return the folder prefix in '
                                 'which Snakemake would install a conda '
                                 'environemnt.')
Parser.add_argument('-e', '--envfile', help='Conda environment YAML file')
Parser.add_argument('-p', '--path', help='path to install conda environment to')
Args = vars(Parser.parse_args())

if __name__ == '__main__':
    main()
