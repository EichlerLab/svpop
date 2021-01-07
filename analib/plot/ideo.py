"""
Generate ideogram plots.
"""

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import analib

IDEO_HIST_PARAMS = {
    'x_adjust': 0.1,  # Adjust x in by this proportion to make space for y-axis labels
    'y_pad': 0.3,     # Add this proportion to the total height and distribute between plots as padding
    'band_color': {   # Color for bands keyed by "gieStain" is chromosome BED (from UCSC).
        'gneg': '0.875',
        'gpos100': '0.00',
        'gpos75': '0.16',
        'gpos50': '0.32',
        'gpos25': '0.48',
        'acen': 'darkred',
        'gvar': '0.00',
        'stalk': '0.64'
    },
    'chroms': [  # Chromosome order. This order minimizes dead space in multi-faceted plot with two columns (hg38)
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chrY',
        'chrX', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'
    ],
    'anno_fig_height': 0.25,  # Proportion of original figure to expand below for annotation space
    'band_prop': 0.70,        # Proportion of annotation space for chromosome bands
    'gap_prop': 0.20,         # Proportion of annotation space for gap lines through chromosome bands
    'sd_color': 'goldenrod',       # SD overlay color
    'sd_alpha_range': (0.2, 0.8),  # The range of alpha values by SD identity are scaled within this range
    'tr_color': 'darkcyan',  # Tandem repeat overlay color
    'tr_alpha': 0.8,         # Tandem repeat overlay alpha
    'tr_min': 5000           # Tandem repeat minimum length
}

class IdeoHistogram:
    """
    Object returned from ideo_hist. Contains the figure and a dictionary of axes with tuple keys with the same
    indexes as a matrix of chromosome names (matrix_chr_name). This object makes it possible to write the figure as
    well as modify it (e.g. add additional annotations).
    """

    def __init__(self, fig, ax_dict, matrix_chr_name):
        self.fig = fig
        self.ax_dict = ax_dict
        self.matrix_chr_name = matrix_chr_name

    def close(self):
        """
        Close figure.
        """
        return plt.close(self.fig)

    def savefig(self, *args, **kwargs):
        return self.fig.savefig(*args, **kwargs)


def ideo_hist(
        df, fai_file_name,
        df_band, df_gap, df_sd=None, df_tr=None,
        width=10, height=12, dpi=300,
        label_col='SVTYPE',
        label_order=['INS', 'DEL', 'INV', 'SNV'],
        label_color={'INS': 'blue', 'DEL': 'red', 'INV': 'green', 'SNV': 'black'},
        plot_params=None,
        ideo_patch_z_start=100,
        ylim_dict=None,
        cb_func=None,
        verbose=False
    ):
    """
    Create an ideogram with a panel for each chromosome.

    :param df: DataFrame of values to plot. Must have "#CHROM", "POS", "END", and `label_col'. May be `None` if all
        plotting is handled by a callback function.
    :param fai_file_name: FAI file name.
    :param df_band: DataFrame of chromosome bands.
    :param df_gap: DataFrame of assembly gaps.
    :param df_sd: DataFrame of SDs with identity.
    :param df_tr: DataFrame of tandem repeat locations.
    :param width: Figure width in inches.
    :param height: Figure height in inches.
    :param dpi: Figure DPI.
    :param label_col: Labels to plot against. Stacked bar colors are separated on this column in `df`.
    :param label_order: Order of labels. Labels not in this list are ignored.
    :param label_color: Label color (color of stacked bars). May be a list of the same length as `label_order` or a
        dict of values with `label_order` items as keys. Must be valid matplotlib colors.
    :param plot_params: Plot parameters as a dictionary. `IDEO_HIST_PARAMS` is the default set of parameters.
    :param ideo_patch_z_start: Relative order of patches (chromosome ideo pieces). Set to 100 by default should make
        it appear behind other plot elements if they overlap.
    :param ylim_dict: Optional y-limits as a dict keyed on chromosomes.
    :param cb_func: If not `None`, call this function for each chromosome. Function signature is
        `cb_func(df, chrom, ax, fig)` where `df` is the DataFrame given to this function (not subset for the
        chromosome), `chrom` is the chromosome name, 'ax` is the matplotlib axes object for this chromosome, and `fig`
        is the matplotlib figure object.
    :param verbose: If `True`, print verbose information while plotting.

    :return: A `IdeoHistogram` object with the figure, a dict of axes, and a matrix of chromosome names.
    """

    # Init parameters
    if plot_params is None:
        plot_params = IDEO_HIST_PARAMS

    plot_params = plot_params.copy()

    chroms = plot_params['chroms']

    if ylim_dict is None:
        ylim_dict = dict()

    if issubclass(label_color.__class__, dict):
        label_color = [label_color[val] for val in label_order]

    ### Read ###

    # Get variant midpoints
    if df is not None:
        df = df[['#CHROM', 'POS', 'END', label_col]].copy()

        df['POS'] = (df['POS'] + df['END']) // 2

    # Read FAI
    df_fai = analib.ref.get_df_fai(fai_file_name)

    if df_sd is not None:
        df_sd = df_sd.copy()

        # Scale identity to span SD_ALPHA_RANGE
        sd_min = np.min(df_sd['MATCH'])

        df_sd['ALPHA'] = (df_sd['MATCH'] - sd_min) / (1 - sd_min) * (plot_params['sd_alpha_range'][1] - plot_params['sd_alpha_range'][0]) + plot_params['sd_alpha_range'][0]

    # Read TR
    if df_tr is not None:
        df_tr = df_tr.loc[(df_tr['END'] - df_tr['POS']) >= plot_params['tr_min']].copy()

    ### Assign chrom matrix ###
    # Note: assumes len(chroms) is an even number (update to handle more general cases)

    # Get names
    matrix_chr_name = np.array(
        [
            chroms[:len(chroms) // 2],
            chroms[-1:len(chroms) // 2 - 1:-1]
        ]
    ).T

    # Get chromosome lengths
    matrix_chr_len = np.array(
        list(map(lambda val: df_fai[val], matrix_chr_name))
    )

    # Find max horizontal pair
    max_pair_bp = np.max(
        np.apply_along_axis(np.sum, 1, matrix_chr_len)
    ) * (
        (1 + plot_params['x_adjust'] * 2)
    )

    # Get height per subplot
    subplot_height = 1 / matrix_chr_name.shape[0]  # Proportion for each row
    subplot_height *= (1 - plot_params['y_pad'])  # Shrink by total pad height

    pad_height = plot_params['y_pad'] / (matrix_chr_name.shape[0])  # Pad between each row

    ### Plot ###

    fig = plt.figure(figsize=(width, height), dpi=dpi)

    ax_dict = dict()

    for i in range(matrix_chr_name.shape[0]):
        for j in range(matrix_chr_name.shape[1]):

            chrom = matrix_chr_name[i, j]

            # Make axes
            ax_len = matrix_chr_len[i, j] / max_pair_bp

            ax = fig.add_axes((
                plot_params['x_adjust'] if j == 0 else 1 - ax_len,  # Left align if column 0, right align if column 1
                1 - (subplot_height * (i + 1) + pad_height * i),
                ax_len,
                subplot_height
            ))

            ax_dict[(i, j)] = ax

            # Histogram
            if df is not None:
                if verbose:
                    print('SVTYPE hist: {}'.format(chrom))

                ax.hist(
                    [
                        df.loc[
                            (df['#CHROM'] == chrom) & (df[label_col] == label),
                            'POS'
                        ] for label in label_order
                    ],
                    histtype='bar',
                    stacked=True,
                    bins=np.int32(matrix_chr_len[i, j] // 1e6),
                    color=label_color,
                    label=label_order
                )

            # Callback plot function
            if cb_func is not None:
                if verbose:
                    print('Callback: {}'.format(chrom))

                cb_func(df, chrom, ax, fig)

            if chrom in ylim_dict:
                ax.set_ylim(ylim_dict[chrom])

            # Scale x to Mbp
            ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: str(int(x // 1e6))))

            # Adjust x-axis (0 to chromosome max)
            ax.set_xlim((-(df_fai[chrom] * 0.01), df_fai[chrom] * 1.01))

            # Remove spines
            for loc in ('top', 'right'):
                ax.spines[loc].set_visible(False)

            # Title
            ax.text(
                0.5, 0.95,
                matrix_chr_name[i, j],
                horizontalalignment='center',
                verticalalignment='top',
                transform=ax.transAxes,
                bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.5}
            )

            # Aestetics
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)

    # Add bands, gaps, SDs, and TRs
    if verbose:
        print('Adding ideo bands')

    patch_z = ideo_patch_z_start

    for i in range(matrix_chr_name.shape[0]):
        for j in range(matrix_chr_name.shape[1]):

            ax = ax_dict[(i, j)]
            chrom = matrix_chr_name[i, j]

            y_min, y_max = ax.get_ylim()

            # Create annotation space below histogram
            anno_height = (y_max - y_min) * plot_params['anno_fig_height']

            ax.set_ylim(
                np.min([y_min, -anno_height]),
                np.max([y_max, anno_height])
            )

            # Set positions for band annotation
            band_y_pos = - anno_height * ((1 - plot_params['band_prop']) / 2)
            band_y_height = anno_height * plot_params['band_prop']

            # Set positions for gap annotations
            gap_y_pos = - anno_height * ((1 - plot_params['gap_prop']) / 2)
            gap_y_height = anno_height * plot_params['gap_prop']

            # Add Bands
            for index, row in df_band.loc[df_band['#chrom'] == chrom].iterrows():

                if row['gieStain'] == 'acen':
                    continue  # Skip centromeres (plot over all other annotations)

                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['start'], band_y_pos),
                        row['end'] - row['start'],
                        -band_y_height,
                        facecolor=plot_params['band_color'][row['gieStain']],
                        zorder=patch_z
                    )
                )

                patch_z += 1

            # Add Gaps
            for index, row in df_gap.loc[df_gap['#CHROM'] == chrom].iterrows():
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['START'], gap_y_pos),
                        row['END'] - row['START'],
                        -gap_y_height,
                        facecolor='black',
                        zorder=patch_z
                    )
                )

                patch_z += 1

            # Add TR
            if df_tr is not None:
                for index, row in df_tr.loc[df_tr['#CHROM'] == chrom].iterrows():
                    ax.add_patch(
                        mpl.patches.Rectangle(
                            (row['POS'], band_y_pos),
                            row['END'] - row['POS'],
                            -band_y_height,
                            facecolor=plot_params['tr_color'],
                            alpha=plot_params['tr_alpha'],
                            zorder=patch_z
                        )
                    )

                    patch_z += 1

            # Add SD
            if df_sd is not None:
                for index, row in df_sd.loc[df_sd['#CHROM'] == chrom].iterrows():
                    ax.add_patch(
                        mpl.patches.Rectangle(
                            (row['POS'], band_y_pos),
                            row['END'] - row['POS'],
                            -band_y_height,
                            facecolor=plot_params['sd_color'],
                            alpha=row['ALPHA'],
                            zorder=patch_z
                        )
                    )

                    patch_z += 1

            # Add centromere bands
            for index, row in df_band.loc[(df_band['#chrom'] == chrom) & (df_band['gieStain'] == 'acen')].iterrows():
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['start'], band_y_pos),
                        row['end'] - row['start'],
                        -band_y_height,
                        facecolor=plot_params['band_color'][row['gieStain']],
                        zorder=patch_z
                    )
                )

                patch_z += 1

    # Return
    return IdeoHistogram(fig, ax_dict, matrix_chr_name)
