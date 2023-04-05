"""
Handle merging, comparing, and annotating samples from multiple sources.
"""

### Get svpop directory ###

SVPOP_DIR = os.path.dirname(workflow.snakefile)


### Library dependencies ###
import sys
sys.path.append(os.path.join(SVPOP_DIR, 'dep'))  # kanapy: K-mer toolkit
sys.path.append(os.path.join(SVPOP_DIR, 'dep', 'ply'))  # ply: Python lex/yacc


### Imports ###

import collections
import datetime
import gzip
import os
import re
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
CONFIG_FILE_NAME_INSTALLED = os.path.join(SVPOP_DIR, 'local', CONFIG_FILE_NAME)

if os.path.isfile(CONFIG_FILE_NAME_INSTALLED):
    configfile: CONFIG_FILE_NAME_INSTALLED

configfile: CONFIG_FILE_NAME

if 'ucsc_ref_name' not in config:
    config['ucsc_ref_name'] = 'hg38'

UCSC_REF_NAME = config['ucsc_ref_name']


### Setup default configuration items ###

if 'reference_fai' not in config:
    config['reference_fai'] = config['reference'] + '.fai'


### Samples Table ###

SAMPLE_TABLE_FILE_NAME = config.get('variant_table', 'config/samples.tsv')

SAMPLE_TABLE = svpoplib.rules.get_sample_table(SAMPLE_TABLE_FILE_NAME)


### Sample info table ###

SAMPLE_INFO_TABLE_COL_TYPES = {
    'SAMPLE': object,
    'SEX': object,
    'ETHNICITY': object
}

SAMPLE_INFO_FILE_NAME = config.get('sample_info_table', 'config/sample_info.tab')

if os.path.isfile(SAMPLE_INFO_FILE_NAME):
    SAMPLE_INFO_TABLE = pd.read_csv(SAMPLE_INFO_FILE_NAME, sep='\t', header=0, dtype=SAMPLE_INFO_TABLE_COL_TYPES)
    SAMPLE_INFO_TABLE.set_index(['SAMPLE'], inplace=True, drop=False)
else:
    SAMPLE_INFO_TABLE = pd.DataFrame([], columns=['SAMPLE', 'SEX', 'ETHNICITY'])


### Shell prefix ###

SHELL_PREFIX = 'set -euo pipefail; '

SETENV_SITE = f'{SVPOP_DIR}/config/setenv.sh'
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
include: 'rules/data/anno.snakefile'
include: 'rules/data/homopolymer.snakefile'
include: 'rules/data/ref.snakefile'


## Variants ##

# Variant tools
include: 'rules/variant/intersect.snakefile'
include: 'rules/variant/intersect_nearest.snakefile'
include: 'rules/variant/variant_global.snakefile'
include: 'rules/variant/svset.snakefile'
include: 'rules/variant/svtypecombined.snakefile'
include: 'rules/variant/bedconversion.snakefile'

# include: 'rules/variant/intersect_nearest.snakefile'

# Variant BED parsers for callers
include: 'rules/variant/bed/altdup.snakefile'
include: 'rules/variant/bed/bed.snakefile'
include: 'rules/variant/bed/pavbed.snakefile'
include: 'rules/variant/bed/pbsv.snakefile'
#include: 'rules/variant/bed/sniffles.snakefile'  # Moved to VCF parser
include: 'rules/variant/bed/svim.snakefile'
#include: 'rules/variant/bed/svimasm.snakefile'  # Moved to VCF parser
include: 'rules/variant/bed/vcf.snakefile'

# VCF output
include: 'rules/variant/vcf.snakefile'

# Annotations
include: 'rules/variant/anno/altmap.snakefile'
include: 'rules/variant/anno/homopolymer.snakefile'
include: 'rules/variant/anno/refseq.snakefile'
include: 'rules/variant/anno/regions.snakefile'
include: 'rules/variant/anno/regulation.snakefile'
include: 'rules/variant/anno/repeats.snakefile'
include: 'rules/variant/anno/seq_content.snakefile'

# Callerset
include: 'rules/variant/callerset.snakefile'

# Sample set
include: 'rules/variant/sampleset/bed.snakefile'
include: 'rules/variant/sampleset/anno.snakefile'
include: 'rules/variant/sampleset/tables.snakefile'

# Summary tables
include: 'rules/variant/tables.snakefile'

# Tracks
include: 'rules/tracks/variant.snakefile'

### Includes ###

## Document and Definitions ##
include: 'rules/definitions.snakefile'
