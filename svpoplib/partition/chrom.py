"""
Balance sequences into a set number of partitions.

Implementations solves a partition problem by Largest Differencing Method (LDM) using the Karmarkar–Karp algorithm.
"""

import pandas as pd

import svpoplib.ds.max_heap


class SeqFai:
    """
    Contains a sequence and it's length. Comparable and summable.
    """

    def __init__(self, chrom, len):
        self.chrom = chrom
        self.chrom_len = len

    def __lt__(self, other):
        return self.chrom_len < other.chrom_len

    def __cmp__(self, other):
        return cmp(self.chrom_len, other.chrom_len)

    def __add__(self, other):
        if issubclass(other, SeqFai):
            return self.chrom_len + other.chrom_len
        else:
            return self.chrom_len + other

    def __repr__(self):
        return f'{self.chrom}:{self.chrom_len}'


class Partition:
    """
    One partition containing k bins (k is the number of portitions). Each bin is a list of sequence SeqFai objects
    split into each partition. A collection of partitions is collapsed to one iteratively using merge(), and the final
    Partition's `val_list` object contains a list of SeqFai objects for each partition.
    """

    def __init__(self, val_list):

        if len(val_list) == 0:
            raise ValueError('Cannot create partition from an empty list')

        self.val_list = sorted(val_list, key=lambda vals: sum([val.chrom_len for val in vals]), reverse=True)  # Lists of values in this partition
        self.weight = sum([val.chrom_len for val in val_list[0]]) - sum([val.chrom_len for val in val_list[-1]])  # Largest difference among weights in this partition

    def merge(self, other):
        """
        Merge this partition with another and return a new partition object.

        :param other: Other partition.

        :return: New partition object with `self` and `other` merged optimally.
        """
        return Partition([
                val_l + val_r for val_l, val_r in zip(
                    self.val_list, other.val_list[::-1]
                )
            ])

    def __lt__(self, other):
        return self.weight < other.weight

    def __repr__(self):
        repr_str = '[{'

        repr_str += '}, {'.join(
            [
                ', '.join([str(val) for val in vals]) for vals in self.val_list
            ]
        )

        repr_str += '}]'

        return repr_str


def partition(chrom_series, partitions):
    """
    Split records evently into a set number of partitions.

    :param chrom_series: Pandas Series object with chromosome or contig names as the index and lengths as the values.
    :param partitions: Number of partitions to split sequence records into.

    :return: A list (`partitions` elements long) where each element is a sorted tuple of chromosoeme/contig names
        assigned to one partition.
    """

    # Check arguments
    if not isinstance(partitions, int):
        raise ValueError(f'partition(): Argument partitions is not type int: {type(partitions)}')

    if partitions < 2:
        raise ValueError(f'partition(): Argument partitions must be 2 or greater: {partitions}')

    if chrom_series is None:
        raise ValueError(f'partition(): Argument chrom_series is None')

    if not isinstance(chrom_series, pd.Series):
        raise ValueError(f'partition(): Argument chrom_series is not type pd.Series: {type(chrom_series)}')

    if chrom_series.shape[0] == 0:
        return [()] * partitions

    # Construct a max-heap of initial Partition objects (one sequence per partition)
    heap = list()

    for chrom, chrom_len in chrom_series.items():
        val_list = list()

        val_list.append([SeqFai(chrom, chrom_len)])

        for i in range(partitions - 1):
            val_list.append(list())

        heap.append(Partition(val_list))

    svpoplib.ds.max_heap.heapify(heap)

    # Iteratively collapse partitions until there is one left
    while len(heap) > 1:
        svpoplib.ds.max_heap.insert(
            heap,
            svpoplib.ds.max_heap.pop(heap).merge(
                svpoplib.ds.max_heap.pop(heap)
            )
        )

    return [
        tuple(sorted([val.chrom for val in vals])) for vals in heap[0].val_list
    ]
