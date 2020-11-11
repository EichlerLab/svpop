import pandas as pd
import numpy as np

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3

mpl.rcParams['pdf.fonttype'] = 42


def make_venn(
        set_a, set_b,
        out_file,
        name_a='A', name_b='B',
        title='',
        len_stat=True,
        cb_palette=("#e41a1c", "#377eb8", "#4daf4a")
):
    """
    Make a Venn diagram of variants in two variant BED files.

    :param a_file: List of elements in set A.
    :param b_file: List of elements in set B.
    :param out_file: Output file or list of output files.
    :param name_a: Name of set A.
    :param name_b: Name of set B.
    :param title: Figure title.
    :param len_stat: Include variant size statistics extracted from the sequence ID.
    :param cb_palette: Color palette as a list (Set A, Set B, Set AB).
    """

    names = (name_a, name_b)

    fontsize = 12

    label_font_size = fontsize
    sublabel_font_size = fontsize

    font = {'size': fontsize}
    mpl.rc('font', **font)

    sets = [set_a, set_b]

    # Make plot
    f = plt.figure(figsize=(10, 8))

    # Make venn
    v = venn2(sets, names)

    named_sets = {
        "10": sets[0] - sets[1],
        "01": sets[1] - sets[0],
        "11": sets[0].intersection(sets[1])
    }

    for name, named_set in named_sets.items():

        if len_stat:
            set_lengths = sorted([int(i.split('-')[3]) for i in named_set], reverse=True)

            if len(set_lengths) > 0:
                new_label_text = "{:,d}\n({:,d} bp,\n{:,d} bp)".format(
                    len(named_set),
                    int(np.ceil(np.mean(set_lengths))),
                    int(np.ceil(np.median(set_lengths))),
                )

            else:
                new_label_text = '{:,d}'.format(len(named_set))


        else:
            new_label_text = '{:,d}'.format(len(named_set))

        this_label = v.get_label_by_id(name)\

        if this_label is not None:
            this_label.set_text(new_label_text)

    # Set font and style
    for l in v.subset_labels:
        if l is not None:
            l.set_fontsize(sublabel_font_size)

    for l in v.set_labels:
        if l is not None:
            l.set_fontsize(label_font_size)

    for p in v.patches:
        if p is not None:
            p.set_color(cb_palette[v.patches.index(p)])
            p.set_edgecolor("black")

    # Set title
    plt.title(title, fontsize=14)

    # Write
    if not (isinstance(out_file, list) or isinstance(out_file, tuple)):
        out_file = [out_file]

    for out_file_name in out_file:
        plt.savefig(out_file_name, bbox_inches='tight')

    # Close plots
    plt.close('all')


def make_venn_3way(
        set_a, set_b, set_c,
        out_file,
        name_a='A', name_b='B', name_c='C',
        title='',
        len_stat=True,
        cb_palette=("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")
):
    """
    Make a Venn diagram of variants in two variant BED files.

    :param set_a: Elements in set A.
    :param set_b: Elements in set B.
    :param set_c: Elements in set C.
    :param out_file: Output file or list of output files.
    :param name_a: Name of set A.
    :param name_b: Name of set B.
    :param name_c: Name of set C.
    :param title: Figure title.
    :param cb_palette: Color palette as a list (Set A, Set B, Set AB, ...).
    """

    names = (name_a, name_b, name_c)

    fontsize = 12

    label_font_size = fontsize
    sublabel_font_size = fontsize

    font = {'size': fontsize}
    mpl.rc('font', **font)

    sets = [set_a, set_b, set_c]

    # Make plot
    f = plt.figure(figsize=(10, 8))

    # Make plot
    f = plt.figure(figsize=(10, 8))

    # Make venn
    v = venn3(sets, names)

    named_sets = {
        '001': (sets[2] - sets[1]) - sets[0],
        '010': (sets[1] - sets[0]) - sets[2],
        '011': (sets[1] & sets[2]) - sets[0],
        '100': (sets[0] - sets[1]) - sets[2],
        '101': (sets[0] & sets[2]) - sets[1],
        '110': (sets[0] & sets[1]) - sets[2],
        '111': sets[0] & sets[1] & sets[2]
    }

    for name, named_set in named_sets.items():
        if len_stat:
            set_lengths = sorted([int(i.split('-')[3]) for i in named_set], reverse=True)

            if len(set_lengths) > 0:
                new_label_text = "{:,d}\n({:,d} bp,\n{:,d} bp)".format(
                    len(named_set),
                    int(np.ceil(np.mean(set_lengths))),
                    int(np.ceil(np.median(set_lengths))),
                )

            else:
                new_label_text = '{:,d}'.format(len(named_set))

        else:
            new_label_text = '{:,d}'.format(len(named_set))

        this_label = v.get_label_by_id(name)\

        if this_label is not None:
            this_label.set_text(new_label_text)


    # Set font and style
    for l in v.subset_labels:
        if l is not None:
            l.set_fontsize(sublabel_font_size)

    for l in v.set_labels:
        if l is not None:
            l.set_fontsize(label_font_size)

    for p in v.patches:
        if p is not None:
            p.set_color(cb_palette[v.patches.index(p)])
            p.set_edgecolor("black")

    # Set title
    plt.title(title, fontsize=14)

    # Write
    if not (isinstance(out_file, list) or isinstance(out_file, tuple)):
        out_file = [out_file]

    for out_file_name in out_file:
        plt.savefig(out_file_name, bbox_inches='tight')

    # Close plots
    plt.close('all')
