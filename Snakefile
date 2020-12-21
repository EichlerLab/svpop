"""
Handle merging, comparing, and annotating samples from multiple sources.
"""

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

import analib

from Bio import SeqIO
import Bio.bgzf

### Fail on warnings ###

#import warnings

#warnings.simplefilter("error")


### Get svpop directory ###

SVPOP_DIR = os.path.dirname(workflow.snakefile)


### Config ###

CONFIG_FILE_NAME = 'config/config.json'

configfile: os.path.join(SVPOP_DIR, CONFIG_FILE_NAME)
configfile: CONFIG_FILE_NAME

if 'ucsc_ref_name' not in config:
    config['ucsc_ref_name'] = 'hg38'

UCSC_REF_NAME = config['ucsc_ref_name']


### Setup default configuration items ###

if 'reference_fai' not in config:
    config['reference_fai'] = config['reference'] + '.fai'


### Samples Table ###

SAMPLE_TABLE_COL_TYPES = {
    'SOURCE': np.object,
    'SAMPLE': np.object,
    'SEQ_SET': np.object,
    'DATA': np.object,
    'VERSION': np.object
}

SAMPLE_TABLE_FILE_NAME = config.get('variant_table', 'config/samples.tsv')

if os.path.isfile(SAMPLE_TABLE_FILE_NAME):
    SAMPLE_TABLE = pd.read_csv(SAMPLE_TABLE_FILE_NAME, sep='\t', header=0, dtype=SAMPLE_TABLE_COL_TYPES)

    if 'SEQ_SET' not in SAMPLE_TABLE.columns:
        SAMPLE_TABLE['SEQ_SET'] = np.nan

    SAMPLE_TABLE.set_index(['SOURCE', 'SAMPLE', 'SEQ_SET'], inplace=True, drop=False)
else:
    SAMPLE_TABLE = pd.DataFrame([], columns=['SOURCE', 'SAMPLE', 'DATA', 'VERSION'])


### Sample info table ###

SAMPLE_INFO_TABLE_COL_TYPES = {
    'SAMPLE': np.object,
    'SEX': np.object,
    'ETHNICITY': np.object
}

SAMPLE_INFO_FILE_NAME = config.get('sample_info_table', 'config/sample_info.tab')

if os.path.isfile(SAMPLE_INFO_FILE_NAME):
    SAMPLE_INFO_TABLE = pd.read_csv(SAMPLE_INFO_FILE_NAME, sep='\t', header=0, dtype=SAMPLE_INFO_TABLE_COL_TYPES)
    SAMPLE_INFO_TABLE.set_index(['SAMPLE'], inplace=True, drop=False)
else:
    SAMPLE_INFO_TABLE = pd.DataFrame([], columns=['SAMPLE', 'SEX', 'ETHNICITY'])


### Shell prefix ###

shell.prefix('. {}/config/setenv.sh; '.format(SVPOP_DIR))


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
include: 'rules/variant/svset.snakefile'
include: 'rules/variant/variant_global.snakefile'

# Variant callers
include: 'rules/variant/caller/bionano/bed.snakefile'
include: 'rules/variant/caller/deepvariant/bed.snakefile'
include: 'rules/variant/caller/extern/bed.snakefile'
include: 'rules/variant/caller/melt/bed.snakefile'
include: 'rules/variant/caller/pbsv/bed.snakefile'
include: 'rules/variant/caller/phasedsv/bed.snakefile'
include: 'rules/variant/caller/std_vcf/bed.snakefile'
include: 'rules/variant/caller/smrtsv/bed.snakefile'
include: 'rules/variant/caller/smrtsv/bed_dup.snakefile'
include: 'rules/variant/caller/sniffles/bed.snakefile'
include: 'rules/variant/caller/svim/bed.snakefile'
include: 'rules/variant/caller/svimASM/bed.snakefile'

include: 'rules/variant/caller/anno/aligndepth.snakefile'
include: 'rules/variant/caller/anno/altmap.snakefile'
include: 'rules/variant/caller/anno/global.snakefile'
include: 'rules/variant/caller/anno/homopolymer.snakefile'
include: 'rules/variant/caller/anno/regions.snakefile'
include: 'rules/variant/caller/anno/refseq.snakefile'
include: 'rules/variant/caller/anno/regulation.snakefile'
include: 'rules/variant/caller/anno/repeats.snakefile'
include: 'rules/variant/caller/anno/seq_content.snakefile'

include: 'rules/variant/caller/global.snakefile'
include: 'rules/variant/caller/tables.snakefile'
include: 'rules/variant/caller/vcf.snakefile'

# Caller set
include: 'rules/variant/callerset/bed.snakefile'
include: 'rules/variant/callerset/anno.snakefile'

# Sample set
include: 'rules/variant/sampleset/bed.snakefile'
include: 'rules/variant/sampleset/anno.snakefile'
include: 'rules/variant/sampleset/tables.snakefile'

# Variant sets (published or variants from other sources)
include: 'rules/variant/varset/varset.snakefile'

#include: 'rules/variant/varset/set/1kgp1.snakefile'
include: 'rules/variant/varset/set/1kgp3.snakefile'
include: 'rules/variant/varset/set/ak1.snakefile'
include: 'rules/variant/varset/set/audano2019.snakefile'
#include: 'rules/variant/varset/set/dbsnp.snakefile'
include: 'rules/variant/varset/set/dbvar.snakefile'
include: 'rules/variant/varset/set/denovodb.snakefile'
#include: 'rules/variant/varset/set/gonl.snakefile'
#include: 'rules/variant/varset/set/hallsv.snakefile'
include: 'rules/variant/varset/set/hgsvc1.snakefile'
include: 'rules/variant/varset/set/huddleston2017.snakefile'
include: 'rules/variant/varset/set/hx1.snakefile'
#include: 'rules/variant/varset/set/kidd2010.snakefile'
#include: 'rules/variant/varset/set/mills2011.snakefile'
include: 'rules/variant/varset/set/sudmant2015a.snakefile'


## Tracks ##

include: 'rules/tracks/variant.snakefile'


## Depth ##

include: 'rules/depth/depth.snakefile'


## Genotyping ##

include: 'rules/gt/genotype.snakefile'


### Includes ###

## Document and Definitions ##
include: 'rules/definitions.snakefile'
