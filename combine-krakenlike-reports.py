#!/usr/bin/env python
"""
NAME: combine-kraken-reports.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

v0.1   2017/01/22    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), sschmeier.com

template version: 1.6 (2016/11/09)
"""

from signal import signal, SIGPIPE, SIG_DFL
import sys
import os
import os.path
import argparse
import csv
import collections
import gzip
import bz2
import zipfile
import time
import re
import numpy as np

try:
    from colorama import init, Fore, Style
except ImportError:
    sys.stderr.write('colorama lib needed.\n')
    sys.exit(1)

# INIT color
# Initialise colours for multi-platform support.
init()

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = 'v0.1'
__date__ = '2017/01/22'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'

# taxids to level
taxid2level = {}

def alert(atype, text, log):
    d = {'success': Fore.GREEN, 'error': Fore.RED, 'warning': Fore.YELLOW, 'info':''}
    textout = '%s [%s] %s\n' % (time.strftime('%Y%m%d-%H:%M:%S'),
                                atype.rjust(7),
                                text)
    log.write('%s%s%s' % (d[atype], textout, Fore.RESET))
    if atype=='error': sys.exit()
def success(text, log=sys.stderr):
    alert('success', text, log)
def error(text, log=sys.stderr):
    alert('error', text, log)
def warning(text, log=sys.stderr):
    alert('warning', text, log)
def info(text, log=sys.stderr):
    alert('info', text, log)    


def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'Read output files from kraken-report and combine into one report for multiple files.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
         'taxfile',
         metavar='TAXFILE',
         help='Two column file with all taxonomy-ids,tax-name in sequence that should be in outputfile.')

    parser.add_argument(
        'list_files',
        metavar='FILE',
        nargs='+',
        help=
        'kraken-report output file.')

    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')

    parser.add_argument('-s',
                        '--scale',
                        action='store_true',
                        dest='scale',
                        default=False,
                        help='Use the mapped reads to scale each value tomwards the mapped reads. [default: "False"]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOAD FILE """
    if filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename, 'rt')
    else:
        filehandle = open(filename, 'rt')
    return filehandle

def parse_krakenreport(filename, dict_taxids):
    fileobj = load_file(filename)    
    reader = csv.reader(fileobj, delimiter = '\t')
    numroot = 0
    for a in reader:
        taxid = int(a[4])
        numreads = int(a[1])
        if taxid == 1:
            numroot = numreads
        level = a[3]
        
        # store the level of the taxid
        taxid2level[taxid] = level        
        try:
            dict_taxids[taxid]=numreads
        except KeyError:
            # do not use this taxid
            continue
    return list(dict_taxids.values()), numroot


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    # parse the taxids that we need to consider
    reader = csv.reader(load_file(args.taxfile), delimiter = '\t')
    taxid2name = {}
    dict_taxids = collections.OrderedDict()
    for a in reader:
        try:
            taxid = int(a[0])
        except ValueError:
            sys.stderr.write('Non-numerical taxid found. Exit.\n')
            sys.exit(1)
        taxid2name[taxid] = a[1]
        dict_taxids[taxid] = 0
    sys.stderr.write('%i taxonomy ids to consider.\n'%len(dict_taxids.keys()))

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')
        
    names = []
    a = []
    roots = []
    for fn in args.list_files:
        #name = os.path.basename(fn)
        name = os.path.dirname(fn).split('/')[-1]
        res, numroot = parse_krakenreport(fn, dict_taxids.copy())
        roots.append(float(numroot))
        names.append(name)
        a.append(res)

    # numpy array from results
    a = np.asarray(a)
    a = a.transpose()
    
    # should we scale?
    if args.scale:
        a = a / roots
    
    outfileobj.write('TAXID\t%s\tLEVEL\tTAXNAME\n' %('\t'.join(names)))
    num = -1    
    for taxid in dict_taxids:
        # tax empty do not consider
        num += 1
        if a[num].sum() == 0:
            continue
        else:
            if args.scale:
                values = '\t'.join(['%.10f'%f for f in a[num]])
            else:
                values = '\t'.join([str(i) for i in a[num]])
            outfileobj.write('%i\t%s\t%s\t%s\n' %(taxid,
                                                  values,
                                                  taxid2level[taxid],
                                                  taxid2name[taxid]))
    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())

