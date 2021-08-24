from ruffus import *
import sys
import os
import sqlite3
import shutil
import cgatcore.experiment as E
from cgatcore import pipeline as P
import re
import glob
import collections


# Pipeline configuration
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

####merged fil####

@follows(mkdir("bedfiles.dir"))
@transform (PARAMS["raw_file.dir"], formatter(), r"bedfiles.dir/{basename[0]}.bed")
def bedile(infile,outfile):
    inputs = P.snip (infile,".ped")
    output = P.snip (outfile,".bed")
    faminfile = output+".b"
    ids = output+".f"
    removefam = output +".bim"
    removeid=output +".fam"
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink --make-bed --file %(inputs)s --out  %(output)s'''
    P.run(statement)


@follows(mkdir("Missing_SNPs.dir"))
@transform (bedile, regex(r"bedfiles.dir/(.+).bed"), r"Missing_SNPs.dir/\1.missingness.bed")

def missingness(infile,outfile):
    inputs = P.snip (infile,".bed")
    module = PARAMS["Plink_module"]
    genotype = PARAMS ["Missing_genotype_threshold"]
    before = "Missing_SNPs.dir/"+"Before_filtering_missigness"
    After = "Missing_SNPs.dir/"+"After_filtering_missigness"
    output = P.snip (outfile,".bed")
    Rmodule=PARAMS["Rmodule"]

    statement = '''module load %(module)s && plink --bfile %(inputs)s 
                                                   --geno %(genotype)s
                                                   --make-bed
                                                   --out %(output)s &&
                                                   plink --bfile %(inputs)s 
                                                   --missing
                                                   --out %(before)s &&
                                                   plink --bfile %(output)s
                                                   --missing
                                                   --out %(After)s &&
                                                   module load %(Rmodule)s &&
                                                   Rscript Missigness_histogram.R &&                                      
                                                   module unload %(Rmodule)s
                                                   '''
    P.run(statement)



@follows (missingness, mkdir ("MAF_filter.dir"))
@transform (missingness, regex(r"Missing_SNPs.dir/(.+).bed"),
            add_inputs(r"bedfiles.dir/*.bed"), 
            r"MAF_filter.dir/\1.MAF.bed")

def allelefreq(infiles,outfile):
    orginal,genotype = infiles
    input1 = P.snip (orginal,".bed")
    input2 = P.snip (genotype,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    Rmodule=PARAMS["Rmodule"]
    before_genotype = "MAF_filter.dir/"+"Before_genotype_filter"
    After_genotye= "MAF_filter.dir/"+"After_genotype_filter"
    MAF_filter = "MAF_filter.dir/"+"After_MAF_filter"
    MAF = PARAMS["MAF_threshold"]
    statement = '''module load %(module)s && plink --bfile %(input2)s  
                   --freq
                   --out %(before_genotype)s &&
                   plink --bfile %(input1)s --maf %(MAF)s --out %(output)s --make-bed &&
                   plink --bfile %(input1)s --freq --out %(After_genotye)s &&
                   plink --bfile %(output)s --freq --out %(MAF_filter)s &&
                   module load %(Rmodule)s && Rscript maf_frequncy_hist.R && module unload %(Rmodule)s'''
                   
    P.run(statement)


@follows (allelefreq,mkdir("HWE_filter.dir"))
@transform (allelefreq, regex(r"MAF_filter.dir/(.+).bed"), r"HWE_filter.dir/\1.bed")

def hwe(infile,outfile):
    inputs = P.snip (infile,".bed")
    module = PARAMS["Plink_module"]
    hwe = PARAMS ["HWE_cutoff"]
    output = P.snip (outfile,".bed")
    statement = '''module load %(module)s && plink --bfile %(inputs)s
                   --hwe %(hwe)s
                   --make-bed
                   --out %(output)s'''
    P.run(statement)


@follows (hwe,mkdir("Samples_missing_GT_rate.dir"))
@transform (hwe, regex("HWE_filter.dir/(.+).MAF.bed"), r"Samples_missing_GT_rate.dir/\1.samples_missing_GT_rate.bed")

def GT (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    GT = PARAMS["missing_GT_rate"]
    Afterfilter = "Samples_missing_GT_rate.dir/" + "After_filter" 
    statement = ''' module load %(module)s && 
                   plink -bfile %(inputs)s --missing --out %(output)s &&
                   plink --bfile %(inputs)s
                   --mind %(GT)s
                   --make-bed
                   --out %(output)s && plink --bfile %(output)s --missing --out %(Afterfilter)s '''
    P.run (statement)





@follows (GT, mkdir("Filter_duplicate_IDs.dir"))
@transform(GT, regex(r"Samples_missing_GT_rate.dir/(.+).bed"), r"Filter_duplicate_IDs.dir/\1.bed")
def removedups (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    plink2 = PARAMS["plink2"]
    Rmodule=PARAMS["Rmodule"]
    statement = '''module load %(plink2)s &&
                   plink2 --bfile %(inputs)s
                         --rm-dup exclude-all
                         --out %(output)s
                         --make-bed
                         '''
    P.run(statement)



@follows (removedups, mkdir("Hetrozygosity.dir"))
@transform(removedups, regex(r"Filter_duplicate_IDs.dir/(.+).bed"), r"Hetrozygosity.dir/\1.hetro.het")
def het (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".het")
    module = PARAMS["Plink_module"]
    Rmodule=PARAMS["Rmodule"]
    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                         --het
                         --out %(output)s && 
                         module load %(Rmodule)s && 
                         Rscript hetro_histo.R'''
    P.run(statement)





@follows (het,mkdir ("Valid_hetrosamples.dir"))
@transform(het,regex(r"Hetrozygosity.dir/(.+).hetro.het"),r"Valid_hetrosamples.dir/\1.Valid_samples")

def keep_valid_het (infile,outfile):
    Rmodule=PARAMS["Rmodule"]
    statement = '''module load %(Rmodule)s && 
                   Rscript Hetro_filter.R %(infile)s %(outfile)s &&
                   module unload %(Rmodule)s'''
    P.run(statement)



@follows (keep_valid_het, mkdir ("Filter_hetro_rate.dir"))
@transform (keep_valid_het, regex (r"Valid_hetrosamples.dir/(.+).Valid_samples"),
                            add_inputs(r"Filter_duplicate_IDs.dir/*.bed"),                    
                            r"Filter_hetro_rate.dir/\1.hetrogygous.bed")

def filterhetro (infiles,outfile):
    bed,samples = infiles
    output = P.snip (outfile,".bed")
    Rmodule=PARAMS["Rmodule"]
    module = PARAMS["Plink_module"]
    inputs = P.snip (samples,".bed")
    statement = '''module load %(module)s &&  plink --bfile %(inputs)s 
                         --keep %(bed)s
                         --out %(output)s 
                         --make-bed'''
    P.run (statement)




@follows (keep_valid_het, mkdir ("unrelated_samples.dir"))
@transform(removedups, regex(r"Filter_duplicate_IDs.dir/(.+).bed"), r"unrelated_samples.dir/king.kin0")

def unrelated (infile,outfile):
    Rmodule=PARAMS["Rmodule"]
    module = PARAMS["Plink_module"]
    King = PARAMS["KING_module"]
    output = P.snip(outfile,".kin0")
    statement = '''module load %(King)s &&  king -b %(infile)s
                         --related
                         --degree 3
                         --prefix %(output)s'''
    P.run (statement)


@follows (removedups, mkdir ("Prunned_file.dir"))
@transform (removedups, regex(r"Filter_duplicate_IDs.dir/(.+).bed"), r"Prunned_file.dir/\1.prunned.bed")

def prunned (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".bed")
    inputs = P.snip(infile,".bed")

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s 
                   --indep-pairwise 200 10 0.1 --out %(output)s
                   --make-bed'''
    P.run (statement)


@follows (prunned, mkdir ("samples_prunnedextract.dir"))
@transform (prunned, regex(r"Prunned_file.dir/(.+).prunned.bed"), r"samples_prunnedextract.dir/\1.prunned.extract.bed")

def extract (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".bed")
    inputs = P.snip(infile,".bed")
    snp_list = inputs + ".prune.in"

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                   --extract %(snp_list)s 
                   --make-bed
                   --out %(output)s'''
    P.run (statement)

@follows (extract, mkdir ("sex_check.dir"))
@transform (hwe, regex("HWE_filter.dir/(.+).MAF.bed"), r"sex_check.dir/\1.sexcheck")

def sexcheck (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".sexcheck")
    inputs = P.snip(infile,".bed")
    statement = ''' module load %(module)s && 
                    plink --bfile %(inputs)s 
                    --check-sex
                    --out %(output)s '''
    
    P.run (statement)



@follows(bedile,allelefreq,hwe,het,keep_valid_het,unrelated,prunned,extract,sexcheck)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))



