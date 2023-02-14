#!/usr/bin/env python
#-*- coding:utf-8 -*-
import os
import sys
import yaml

#usage: snakemake -s ~/SomaticNano.pipeline.v1.smk --configfile ~/SomaticNano.pipeline.v1.yaml -j 20


### Software environment
CondaENV = config["CondaENV"]
CondaBIN = CondaENV + '/bin'

### Snakemake workflow environment
PIPENV_WZK = config["PIPENV_WZK"]
SCRIPT_DIR =  PIPENV_WZK + '/script'
PIPE_DIR = PIPENV_WZK + '/pipeline'
SRC_DIR = PIPENV_WZK + '/src/NanoHub'
RULE_DIR = PIPENV_WZK +  '/rule'

PIPENV = config["PIPENV"]
MY_NGS_DIR = PIPENV + '/NGS'
MY_RULE_DIR = PIPENV + '/rule'
MY_SRC_DIR = PIPENV + '/src/NanoHub'

### General configuation
PROJECT_PATH = config["project_path"]
IN_PATH = PROJECT_PATH + '/RBngs'
Nano_PATH = PROJECT_PATH + '/RBnanopore'
Hifi_PATH = PROJECT_PATH + '/RBpacbio'
THREADS = config["THREADS"]
ThreadFold = config["ThreadFold"]
SAMPLES = config["SAMPLES"]
PAIRS = config["PAIRS"]


### Rules of snakemake  
#### base rule
# include: RULE_DIR + '/databaseFormat.rule.py'
include: RULE_DIR + '/BaseFunction.rule.py'
#### function rule
include: MY_RULE_DIR + "/NanoMapping.rule.v1.py"
include: MY_RULE_DIR + "/NanoCleanForceCall.rule.v1.py"
include: MY_RULE_DIR + "/SV_plus.rule.v1.py"
include: MY_RULE_DIR + "/NanopolishCallMethylation.rule.v1.py"


rule all:
    input:
        ################################## NanoMapping ##############################################
        expand(Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.bam",sample=SAMPLES),
        expand(Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.bam.bai",sample=SAMPLES),
        expand(Nano_PATH + "/mapping/clean_ngmlr/bam/{sample}.bam",sample=SAMPLES),
        expand(Nano_PATH + "/mapping/clean_ngmlr/bam/{sample}.bam.bai",sample=SAMPLES),
        expand(Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/{pair}.filt1plus.5merge_three.flt.depthT.vcf", pair=PAIRS),
        expand(IN_PATH + "/{sample}/converted.{chr}.bam", sample=SAMPLES,chr=CHRS),
        expand(IN_PATH + "/{sample}/converted.{chr}.bam.bai", sample=SAMPLES,chr=CHRS)