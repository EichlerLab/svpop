"""
Apply filters to sets of variant calls.


A filter spec string specifications is a comma-separated list of filter specifications:

filter_spec_str ::= filter_spec | filter_spec+filter_spec_str


Each specification is a name separated by a colon and it's arguments (colon and arguments are optional)

filter_spec ::= filter_name | filter_name:filter_arg


Filter arguments are a colon separated list of attribute-value pairs delimited by = (value is optional)

filter_arg ::= filter_avp | filter_avp:filter_arg
filter_avp ::= filter_attr | filter_attr=filter_val

For parsing, filter_name may not contain ":", and filter_arg may not contain "=" or ":"


Example: 'notr:ro=any:distance=200:flag,nosd'

```
notr               (filter_name)
  ro = any         (filter_attr = filter_val)
  distance = 200   (filter_attr = filter_val)
  flag             (filter_attr)
nosd               (filter_name)
```
"""

import abc
import operator
import re

import pandas as pd
import numpy as np


######################
### Define Filters ###
######################

### Parent class of filter runners ###

class FilterRunner:
    """
    Parent class of all filters. Defines methods that are available to all filters.
    """

    __metaclass__ = abc.ABCMeta  # Needed for abstract methods (methods that must be defined by subclasses)

    def __init__(self, filter_name, filter_arg, wildcards):
        """
        Create a new filter with name and argument.

        :param filter_name: Name of the filter.
        :param filter_arg: Filter argument.
        :param wildcards: Rule wildcards.
        """

        # Name of the filter and it's argument string
        self.filter_name = filter_name
        self.filter_arg = filter_arg
        self.wildcards = wildcards
        self.filter_arg = filter_arg

        # Overridable by subclasses
        self.default_args = dict()
        self.required_files = dict()

        # Set by init_filter()
        self.args = None
        self.files = None

    def init_filter(self):
        # Initialize argument dictionary
        self.args = self.get_args(self.wildcards, args_to_dict(self.filter_arg))

        # Parse required files
        self.files = dict()

        for key, file_pattern in self.required_files.items():
            try:
                self.files[key] = file_pattern.format(**self.args)

            except KeyError as e:
                raise RuntimeError(
                    'Missing wildcard or argument "{}" while parsing file pattern "{}" ("{}")'.format(
                        e, key, file_pattern
                    )
                )

    @abc.abstractmethod  # Abstract method: must be defined
    def run(self, df):
        """
        Run filter and return a filtered DataFrame.

        :param df: Pandas DataFrame of variant calls.

        :return: Filtered DataFrame.
        """
        pass

    def get_args(self, wildcards, spec_args):

        # Add wildcards to arg_dict
        arg_dict = dict(wildcards)

        # Add default arguments, override wildcards
        if self.default_args is not None:

            for key, val in self.default_args.items():
                if key in arg_dict:
                    raise RuntimeWarning('Default filter argument overrides a wildcard: {}'.format(key))

                arg_dict[key] = val

        # Add arguments from filter specification
        for key, val in spec_args.items():
            if key in wildcards.keys():
                raise RuntimeWarning('Filter argument overrides a wildcard: {}'.format(key))

            arg_dict[key] = val

        # Return the dictionary of arguments
        return arg_dict

    def file_set(self):
        """
        Get a list of files this filter requires.

        :return: A list of required files.
        """

        return set(self.files.values())


### Filter: tr ###

class TRFilterRunner(FilterRunner):
    """
    Filter by annotated tandem repeats (UCSC track).
    """

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(TRFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Default arguments
        self.default_args = {
            'distance': '200',
            'flank': '0',
            'overlap': '50'
        }

        # Required files
        self.required_files = {
            'bed_tr': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/trf/trf_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
        }

        # Init filter
        self.init_filter()

    def run(self, df):

        exclude_set = set(pd.read_csv(self.files['bed_tr'], sep='\t', squeeze=True))

        if self.filter_name == 'notr':
            return df.loc[df['ID'].apply(lambda val: val not in exclude_set)]
        elif self.filter_name == 'intr':
            return df.loc[df['ID'].apply(lambda val: val in exclude_set)]
        else:
            raise RuntimeError('Cannot apply TRFilterRunner to filter with name: {}'.format(self.filter_name))


### Filter: sd ###

class SDFilterRunner(FilterRunner):
    """
    Filter by annotated segmental duplications (UCSC track).
    """

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(SDFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Default arguments
        self.default_args = {
            'distance': '200',
            'flank': '0',
            'overlap': '50'
        }

        # Required files
        self.required_files = {
            'bed_sd': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/sd/sd_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
        }

        # Init filter
        self.init_filter()

    def run(self, df):

        exclude_set = set(pd.read_csv(self.files['bed_sd'], sep='\t', squeeze=True))

        if self.filter_name == 'nosd':
            return df.loc[df['ID'].apply(lambda val: val not in exclude_set)]
        elif self.filter_name == 'insd':
            return df.loc[df['ID'].apply(lambda val: val in exclude_set)]
        else:
            raise RuntimeError('Cannot apply SDFilterRunner to filter with name: {}'.format(self.filter_name))


### Filter: rmsk ###

class RMSKFilterRunner(FilterRunner):
    """
    Filter by annotated repeat-mask annotated regions (UCSC track).
    """

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(RMSKFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Default arguments
        self.default_args = {
            'distance': '200',
            'flank': '0',
            'overlap': '50'
        }

        # Required files
        self.required_files = {
            'bed_rmsk': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-any-all_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
        }

        # Init filter
        self.init_filter()

    def run(self, df):

        exclude_set = set(pd.read_csv(self.files['bed_rmsk'], sep='\t', squeeze=True))

        if self.filter_name == 'normsk':
            return df.loc[df['ID'].apply(lambda val: val not in exclude_set)]

        elif self.filter_name == 'inrmsk':
            return df.loc[df['ID'].apply(lambda val: val in exclude_set)]

        else:
            raise RuntimeError('Cannot apply SDFilterRunner to filter with name: {}'.format(self.filter_name))


### Filter: trsd ###

class TRSDFilterRunner(FilterRunner):
    """
    Filter by annotated segmental duplications (UCSC track).
    """

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(TRSDFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Default arguments
        self.default_args = {
            'distance': '200',
            'flank': '0',
            'overlap': '50'
        }

        # Required files
        self.required_files = {
            'bed_sd': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/sd/sd_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz',
            'bed_tr': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/trf/trf_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
        }

        # Init filter
        self.init_filter()

    def run(self, df):

        exclude_set = set(
            pd.read_csv(self.files['bed_sd'], sep='\t', squeeze=True)
        ).union(
            pd.read_csv(self.files['bed_tr'], sep='\t', squeeze=True)
        )

        if self.filter_name == 'notrsd':
            return df.loc[df['ID'].apply(lambda val: val not in exclude_set)]
        elif self.filter_name == 'intrsd':
            return df.loc[df['ID'].apply(lambda val: val in exclude_set)]
        else:
            raise RuntimeError('Cannot apply TRSDFilterRunner to filter with name: {}'.format(self.filter_name))


### Filter: trsdrmsk ###

class TRSDRMSKFilterRunner(FilterRunner):
    """
    Filter by annotated tandem repeats, segmental duplications, or RepeatMasker (UCSC track).
    """

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(TRSDRMSKFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Default arguments
        self.default_args = {
            'distance': '200',
            'flank': '0',
            'overlap': '50'
        }

        # Required files
        self.required_files = {
            'bed_sd': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/sd/sd_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz',
            'bed_tr': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/trf/trf_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz',
            'bed_rmsk': 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-any-all_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
        }

        # Init filter
        self.init_filter()

    def run(self, df):

        exclude_set = set(
            pd.read_csv(self.files['bed_sd'], sep='\t', squeeze=True)
        ).union(
            pd.read_csv(self.files['bed_tr'], sep='\t', squeeze=True)
        ).union(
            pd.read_csv(self.files['bed_rmsk'], sep='\t', squeeze=True)
        )

        if self.filter_name == 'notrsdrmsk':
            return df.loc[df['ID'].apply(lambda val: val not in exclude_set)]
        elif self.filter_name == 'intrsdrmsk':
            return df.loc[df['ID'].apply(lambda val: val in exclude_set)]
        else:
            raise RuntimeError('Cannot apply TRSDRMSKFilterRunner to filter with name: {}'.format(self.filter_name))


### Filter: svlen ###

class SVLenFilterRunner(FilterRunner):

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(SVLenFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Init filter
        self.init_filter()

    def run(self, df):

        # Check for SVLEN
        if 'SVLEN' not in df.columns:
            raise RuntimeError('SV length filter cannot be applied to a DataFrame with no SVLEN column')

        # Get min and max
        min = None
        max = None

        if 'min' in self.args:
            min = np.int32(self.args['min'])

        if 'max' in self.args:
            max = np.int32(self.args['max'])

        if 'range' in self.args:

            if 'min' in self.args or 'max' in self.args:
                raise RuntimeError(
                    'SV length "range" cannot be specified with "min" or "max": {}'.format(self.filter_arg)
                )

            tok = re.split('\s*-\s*', self.args['range'].strip(), 1)

            if len(tok) != 2:
                raise RuntimeError(
                    'SV length filter with range argument requires a range with a dash (e.g. n-m): {}'.format(
                        self.filter_arg
                    )
                )

            min = np.int32(tok[0])
            max = np.int32(tok[1])

        if min is None and max is None:
            raise RuntimeError(
                'SV length filter was not given min, max, or range: {}'.format(self.filter_arg)
            )

        # Filter
        if min is not None:
            df = df.loc[df['SVLEN'] >= min]

        if max is not None:
            df = df.loc[df['SVLEN'] <= max]

        # Return
        return df


### Filter: field ###

# Filter by a field value in the DataFrame (field, value, and operation are set as attributes to the filter)

class FieldFilterRunner(FilterRunner):

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(FieldFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Init filter
        self.init_filter()

    def run(self, df):

        # Check field
        if 'name' not in self.args:
            raise RuntimeError('Field filter missing required argument: name (field name)')

        if 'value' not in self.args:
            raise RuntimeError('Field filter missing required argument: value (field value)')

        if self.args['name'] not in df.columns:
            raise RuntimeError('DataFrame has no column "{}": {}'.format(self.args['name'], self.filter_arg))

        # Get operator
        if 'op' not in self.args:
            op_str = 'eq'
        else:
            op_str = self.args['op']

        if op_str == 'eq':
            op = operator.eq
        elif op_str == 'ne':
            op = operator.ne
        else:
            raise RuntimeError('Unrecognized attribute: {} ({})'.format(op_str, self.filter_arg))

        # Do filter
        filter_val = self.args['value']

        df = df.loc[df[self.args['name']].apply(lambda val: op(val, filter_val))]

        # Return
        return df

### Filter: autosome ###

# Filter out X and Y

class AutosomeFilterRunner(FilterRunner):

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(AutosomeFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Init filter
        self.init_filter()

    def run(self, df):

        df = df.loc[df['#CHROM'].apply(lambda val: re.match('^chr[XY](_.*)?.*', val) is None)]

        # Return
        return df

### Filter: all ###

# Null filter: Return an unaltered DataFrame. Used to support filter operations such as "notrf_vs_all", which
# variant intersections may do ("all" will filter the second sample with the null filter)

class NullFilterRunner(FilterRunner):

    def __init__(self, filter_name, filter_arg, wildcards):

        # Init parent
        super(NullFilterRunner, self).__init__(filter_name, filter_arg, wildcards)

        # Init filter
        self.init_filter()

    def run(self, df):
        return df


### List of built-in filters ###

filter_runner_dict = {
    'notr': TRFilterRunner,
    'intr': TRFilterRunner,
    'nosd': SDFilterRunner,
    'insd': SDFilterRunner,
    'inrmsk': RMSKFilterRunner,
    'normsk': RMSKFilterRunner,
    'svlen': SVLenFilterRunner,
    'notrsd': TRSDFilterRunner,
    'intrsd': TRSDFilterRunner,
    'notrsdrmsk': TRSDRMSKFilterRunner,
    'intrsdrmsk': TRSDRMSKFilterRunner,
    'field': FieldFilterRunner,
    'autosome': AutosomeFilterRunner,
    'all': NullFilterRunner
}


#########################
### Utility Functions ###
#########################

def args_to_dict(filter_arg):
    """
    Turn a colon-separated list of arguments into a dictionary of elements. Each colon-separated item is an
    attribute-value pair (separated by '=') with a default value of `None` for missing values.

    Example: 'ro=any:distance=200:flag' becomes:
    ```
    {
       'ro': 'any',
       'distance': '200',
       'flag': None
    }
    ```

    Note: 'flag' was added to illustrate the `None` default. It's probably not going to be in real a specification.

    :param filter_arg: Specification string.

    :return: A dictionary of elements.
    """

    spec_args = dict()

    if filter_arg is None:
        filter_arg = ''

    filter_arg = filter_arg.strip()

    # Parse arguments as a colon separated list where each element may have an assignment
    for filter_avp in re.split('\s*:\s*', filter_arg):

        # Split on =
        tok = re.split('\s*=\s*', filter_avp, 1)

        filter_attr = tok[0]
        filter_val = tok[1] if len(tok) > 1 else None

        # Assign to arg_dict
        if filter_attr in spec_args:
            raise RuntimeError('Duplicate attribute in specification: {}: spec={}'.format(filter_attr, filter_arg))

        spec_args[filter_attr] = filter_val

    # Return
    return spec_args


def filter_config_def(filter_config_name, config):
    """
    Check configuration for a filter definition alias.

    :param filter_config_name:
    :param config:
    :return:
    """

    # Check for config section
    if 'filter_def' not in config:
        return None

    if filter_config_name not in config['filter_def']:
        return None

    return config['filter_def'][filter_config_name]


def get_filter_spec_list(filter_spec_str):

    if filter_spec_str is None:
        raise RuntimeError('Filter spec string is None: None')

    filter_spec_str = filter_spec_str.strip()

    if not filter_spec_str:
        raise RuntimeError('Filter spec string is empty')

    filter_spec_list = [
        re.split('\s*:\s*', spec.strip(), 1) for spec in filter_spec_str.split('+') if spec.strip()
    ]

    # Do not allow empty filter specs
    for spec_element in filter_spec_list:
        if not spec_element[0]:
            raise RuntimeError(
                'Filter spec contains at least one element with no filter name (e.g. "++" or "+:"): {}'.format(filter_spec_str)
            )

    # Add a "None" element to spec elements that did not have an argument
    # For example, ('notr') becomes ('notr', None)

    filter_spec_list = [(spec_element[0], None) if len(spec_element) == 1 else tuple(spec_element) for spec_element in filter_spec_list ]

    # Return list
    return filter_spec_list


def get_filter_list(filter_spec_str, wildcards, config):
    """
    Get a list of instantiated filters to be applied.

    :param filter_spec_str: String specifying filters.
    :param wildcards: Rule wildcards.
    :param config: Snakemake config.

    :return: A list of filters.
    """

    # Get filter alias config
    filter_config = config.get('filter_def', {})

    # Create filter list
    filter_list = list()

    for filter_name, filter_arg in get_filter_spec_list(filter_spec_str):

        # Process alias
        if filter_name in filter_config:

            # Alias cannot have arguments
            if filter_arg is not None:
                raise RuntimeError('Specified filter alias with non-empty arguments: {}'.format(filter_name))

            # Add each filter in the alias to the list
            for filter_name_sub, filter_arg_sub in get_filter_spec_list(filter_config[filter_name]):

                # Check for a known filter
                if filter_name_sub not in filter_runner_dict:
                    raise RuntimeError('Filter alias {} contains a definition for an unknown filter: {}'.format(
                        filter_name, filter_name_sub
                    ))

                filter_list.append(
                    filter_runner_dict[filter_name_sub](filter_name_sub, filter_arg_sub, wildcards)
                )

        else:

            # Check for a known filter
            if filter_name not in filter_runner_dict:
                raise RuntimeError('Unknown filter: {}'.format(filter_name))

            # Add to the filter list
            filter_list.append(
                filter_runner_dict[filter_name](filter_name, filter_arg, wildcards)
            )

    # Return list of filters
    return filter_list


####################
### Apply Filter ###
####################

def apply_svset_filter(df, filter_spec_str, wildcards, config):
    """
    Apply filter(s) to a dataframe.

    :param df: DataFrame to filter.
    :param filter_spec_str: String specifying filters and parameters.
    :param wildcards: Rule wildcards.
    :param config: Snakemake config.

    :return: Filtered DataFrame.
    """

    # Check input
    if df is None:
        raise RuntimeError('Cannot filter dataframe: None')

    # Get filter list
    filter_list = get_filter_list(filter_spec_str, wildcards, config)

    for filter_runner in filter_list:
        df = filter_runner.run(df)

    # Return filtered dataframe
    return df


def get_filter_input_files(filter_spec_str, wildcards, config):
    """
    Get a list of input files. Useful for Snakemake rules.

    :param filter_spec_str: Filter specification string (likely `wildcards.svset`)
    :param wildcards: Rule wildcards.
    :param config: Snakemake config.

    :return: A list of input files.
    """

    # Get filter list
    filter_list = get_filter_list(filter_spec_str, wildcards, config)

    file_set = set()

    for filter_runner in filter_list:
        file_set = file_set.union(filter_runner.file_set())

    return sorted(file_set)
