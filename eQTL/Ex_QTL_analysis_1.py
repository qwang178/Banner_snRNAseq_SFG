#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from run_QTL_analysis import run_QTL_analysis
from qtl_utilities import merge_QTL_results
import subprocess
import numpy as np
import pandas as pd
import pytest

def hdf5_results_checking(filename,fun=lambda df: np.mean(df['beta']) ):
    '''For a given hdf5 results file, returns a value derived from the first dataframe
    in the file.'''
    with pd.HDFStore(filename,'r') as h5file:
        feature_id = sorted(h5file.keys())[0]
        df = h5file.select(feature_id)
    return fun(df)

def results_checking(results_checking_dict,error_tolerance=1e-6):
    for f in results_checking_dict.keys():
        check_value = hdf5_results_checking(f)
        print(f,check_value)
        assert(abs(results_checking_dict[f]-check_value)<error_tolerance)


def test_QTL_analysis():
    '''Run a set of test cases'''
    data_path = '/data/breadhea/qwang178/NOMIS/scRNA/Drop8/eQTL/data/'
    covariates_filename = data_path+'expression/Ex/covariates0.tsv'
    geno_prefix = data_path+'genotype/1/1.imputed'
    pheno_filename = data_path+'expression/Ex/e.tsv'
    anno_filename = data_path+'expression/NOMIS.annotation.tsv'
    kinship_filename= data_path+'genotype/1/1.kinship.txt'
    individual2sample_filename = data_path+'NOMIS_gte.txt'
    #feature_filename = data_path+'genotype/11/11.features.txt'
    min_maf = 0.05
    min_hwe_P = 0.001
    min_call_rate = 0.95
    blocksize = 50
    output_dir = data_path+'results/Ex/1/'
    randomSeed = 73
    genetic_range = '1'
    
    ws = 250000
    
    run_QTL_analysis(pheno_filename, anno_filename, geno_prefix, True, output_dir, ws,
                     min_maf, min_hwe_P, min_call_rate, blocksize, cis_mode=True,
                     seed=randomSeed, n_perm=1000, snps_filename=None, 
#                    feature_filename = feature_filename,
                     genetic_range=genetic_range, covariates_filename=covariates_filename, 
                     write_permutations=False, write_feature_top_permutations = True,
                     randomeff_filename=kinship_filename, sample_mapping_filename=individual2sample_filename)

    #results_checking_dict = {output_dir+'qtl_results_1.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)
    
    ##Testing all long names.
    #results_checking_dict = {output_dir+'qtl_results_1.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)
    
    ##Testing all short names + correcting for top SNP per gene (to be added)
    #run again, without specifying chromosome
#    subprocess.call('python run_QTL_analysis2.py '
#                    '-pg {geno_prefix} '
#                    '-af {anno_file} '
#                    '-pf {pheno_file} '
#                    '-od {output_dir} '
#                    '-w {ws} '
#                    '-cf {covariates_file} '
#                    '-rf {kinship_file} '
#                    '-smf {samplemap_file} '
#                    '-c'
#                    '-s {seed}'
#                    '-np {numPerm}'
#                    .format(geno_prefix=geno_prefix,
#                            anno_file=anno_filename,
#                            pheno_file=pheno_filename,
#                            output_dir=output_dir,
#                            ws=ws,
#                            covariates_file=covariates_filename,
#                            kinship_file=kinship_filename,
#                            samplemap_file=individual2sample_filename,
#                            numPerm = 10,
#                            seed = randomSeed
#                            ),
#                    shell=True)

    #results_checking_dict = {output_dir+'qtl_results_all.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)

    

    #run another test case, where we check the subsetting of samples.
    
#    geno_prefix = data_path+'Genotypes/Geuvadis'
#    pheno_filename = data_path+'Expression/Geuvadis_CEU_Expr.txt'
#    output_dir = data_path+'TestOutput/limix_QTL_results/'
    
#    for chromosome in ['1','2']:
#        run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,True,output_dir,ws,min_maf, min_hwe_P,min_call_rate,blocksize,cis_mode=True, seed=randomSeed, n_perm=10,genetic_range=chromosome)

    #results_checking_dict = {output_dir+'qtl_results_1.h5':0.034497, output_dir+'qtl_results_2.h5':0.002150}
    #results_checking(results_checking_dict)


if __name__=='__main__':
    test_QTL_analysis()
