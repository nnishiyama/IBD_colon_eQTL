#!/usr/bin/python
# script to produce covariate matrices with decreasing numbers of RUV factors

import os
import argparse
from argparse import ArgumentParser
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--RNA_covariates', metavar='RNA_COVARIATES_TXT', required=True)
    parser.add_argument('--ruv', metavar='RUV_FACTORS_TXT', required=True)
    parser.add_argument('--out_prefix', metavar='COVARIANCE_OUTPUT_PREFIX', required=True)
    args = parser.parse_args()
    return args

# read in files:
args = parse_args()
# RUV factors: m RUV factors (rows) x n samples (columns)
ruv = pd.read_csv(args.ruv, sep='\t')
# sample covariates: m samples (rows) x n covariates (columns)
meta = pd.read_csv(args.RNA_covariates, sep='\t')

out = str(args.out)
out_dir = out.rsplit('/',1)[0]
# check if output directory exists
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# keep only samples that are in all files
samples = ruv.columns.tolist()

# subset known sample covariates
meta = meta[['sample','sex','batch', 'genotypePC1', 'genotypePC2', 'genotypePC3', 'genotypePC4']]
meta = meta[meta['sample'].isin(samples)]

# transpose data into output format
meta = meta.T
meta.columns = meta.iloc[0]
meta.drop(meta.index[0], inplace=True)
meta.reset_index(inplace=True)

pf = [a + str(b) for a,b in zip((['W']*len(ruv)), list(range(1,len(ruv)+1)))]
ruv.insert(0,'index',pf)

# concatenate matrices into covariate matrix
cov = meta.append(ruv)
cov.reset_index(inplace=True, drop=True)
samples = ['index'] + samples

# loop through to create covariance matrices with decreasing numbers of RUV factors
cov.to_csv('%s.%i.txt' % (out, len(ruv)), sep='\t',index=False, columns=samples)

n = 1
l = len(meta)
for i in range(len(cov)-l):
    cov.drop(cov.tail(n).index, inplace=True)
    p = len(cov)-l
    cov.to_csv('%s_%i.txt' % (out, p), sep='\t',index=False, columns=samples)
