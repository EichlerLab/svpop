"""
Handle merging, comparing, and annotating samples from multiple sources.
"""

### Get svpop directory ###

import os

global workflow

PIPELINE_DIR = os.path.dirname(os.path.realpath(workflow.snakefile))


### Library dependencies ###
import sys
sys.path.append(PIPELINE_DIR)
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))  # kanapy: K-mer toolkit
sys.path.append(os.path.join(PIPELINE_DIR, 'dep', 'ply'))  # ply: Python lex/yacc


### Imports ###

import collections
import datetime
import gzip
import os
import re
import shutil
import subprocess
import intervaltree
import itertools
import intervaltree
import pysam
import shutil
import time

import pandas as pd
import numpy as np

import svpoplib

from Bio import SeqIO
import Bio.bgzf

### Fail on warnings ###

#import warnings

#warnings.simplefilter("error")


### Config ###

CONFIG_FILE_NAME = 'config/config.json'
CONFIG_FILE_NAME_INSTALLED = os.path.join(PIPELINE_DIR, 'local', CONFIG_FILE_NAME)

if os.path.isfile(CONFIG_FILE_NAME_INSTALLED):
    configfile: CONFIG_FILE_NAME_INSTALLED

configfile: CONFIG_FILE_NAME

if 'ucsc_ref_name' not in config:
    raise RuntimeWarning('No "ucsc_ref_name" in configuration (Genome name in the UCSC browser): Setting to UCSC_UNKNOWN - will affect downloading track data for annotations')
    config['ucsc_ref_name'] = 'UCSC_UNKNOWN'


UCSC_REF_NAME = config['ucsc_ref_name']


### Setup default configuration items ###

if 'reference_fai' not in config:
    config['reference_fai'] = config['reference'] + '.fai'


### Samples Table ###

SAMPLE_TABLE_FILE_NAME = config.get('variant_table', 'config/samples.tsv')

SAMPLE_TABLE = svpoplib.rules.get_sample_table(SAMPLE_TABLE_FILE_NAME)


### Sample info table ###

SAMPLE_INFO_FILE_NAME = config.get('sample_info_table', 'config/sample_info.tab')

if os.path.isfile(SAMPLE_INFO_FILE_NAME):
    SAMPLE_INFO_TABLE = pd.read_csv(SAMPLE_INFO_FILE_NAME, sep='\t', header=0)
    SAMPLE_INFO_TABLE.set_index(['SAMPLE'], inplace=True, drop=False)
else:
    SAMPLE_INFO_TABLE = pd.DataFrame([], columns=['SAMPLE', 'SEX', 'ETHNICITY'])


### Shell prefix ###

SHELL_PREFIX = 'set -euo pipefail; '

SETENV_SITE = f'{PIPELINE_DIR}/config/setenv.sh'
SETENV_LOCAL = 'config/setenv.sh'

if os.path.isfile(SETENV_SITE):
    SHELL_PREFIX += f'source {SETENV_SITE}; '

if os.path.isfile(SETENV_LOCAL):
    SHELL_PREFIX += f'source {SETENV_LOCAL}; '

shell.prefix(SHELL_PREFIX)


### Global wildcard constraints ###

wildcard_constraints:
    seq_set='[A-Za-z0-9]+'


### Includes ###

## Data ##
include: os.path.join(PIPELINE_DIR, 'rules/data/anno.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/data/homopolymer.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/data/ref.snakefile')


## Variants ##

# Variant tools
include: os.path.join(PIPELINE_DIR, 'rules/variant/intersect.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/intersect_nearest.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/variant_global.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/svset.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/svtypecombined.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bedconversion.snakefile')


# Flag rules
include: os.path.join(PIPELINE_DIR, 'rules/flag/variant.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/flag/intersect.snakefile')

# Variant BED parsers for callers
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/altdup.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/bed.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/pavbed.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/pavbedhap.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/pbsv.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/bed/vcf.snakefile')

# VCF output
include: os.path.join(PIPELINE_DIR, 'rules/variant/vcf.snakefile')

# Annotations
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/altmap.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/homopolymer.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/refseq.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/regions.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/regulation.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/repeats.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/anno/seq_content.snakefile')

# Callerset
include: os.path.join(PIPELINE_DIR, 'rules/variant/callerset.snakefile')

# Sample set
include: os.path.join(PIPELINE_DIR, 'rules/variant/sampleset/bed.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/sampleset/anno.snakefile')
include: os.path.join(PIPELINE_DIR, 'rules/variant/sampleset/tables.snakefile')

# Summary tables
include: os.path.join(PIPELINE_DIR, 'rules/variant/tables.snakefile')

# Tracks
include: os.path.join(PIPELINE_DIR, 'rules/tracks/variant.snakefile')

### Includes ###

## Document and Definitions ##
include: os.path.join(PIPELINE_DIR, 'rules/definitions.snakefile')
