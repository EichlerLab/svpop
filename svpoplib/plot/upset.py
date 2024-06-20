"""
Create an upset plot from a dataframe.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

import collections
import numpy as np
import os
import pandas as pd


def get_upset_fig(
    df,
    list_column='CALLERSET_LIST',
    caller_order=None,
    caller_order_subset=True,
    color_column=None, color=None, color_label=None, color_order=None,
    title=None,
    legend=True,
    legend_title=None,
    caller_label_dict=None,
    n_bar=None,
    bubble_size=70, bubble_label_size=10,
    width=7, height=7, dpi=300):
    """
    Create an upset plot figure.

    :param df: DataFrame to create the figure from.
    :param list_column: Column of comma-separated lists of callers.
    :param caller_order: Order of callers or None for a default order.
    :param caller_order_subset: If True and caller_order is set, subset `df` to include only records with these callers
        and do not count callers not in caller_order in the upset plot (treat other callers as if they are not present).
    :param color_column: Column of values indicating the colors of stacked bars.
    :param color: Dictionary mapping values in `color_column` to colors.
    :param color_label: Dictionary mapping values in `color_column` to legend labels for each color..
    :param color_order: Order of the legend labels (values of `color_label`).
    :param title: Title of the figure or None for no title.
    :param legend: If True, show the legend.
    :param legend_title: Title of the legend or None for no title.
    :param caller_label_dict: Dictionary mapping caller names in the comma-separated name list `list_column` to
        upset bubble labels for each caller.
    :param n_bar: Maximum number of bars.
    :param bubble_size: Size of bubbles.
    :param bubble_label_size: Size of caller labels to the left of bubbles.
    :param width: Width of the figure.
    :param height: Height of the figure.
    :param dpi: Figure resolution in DPI.

    :return: A figure object.
    """

    # Check arguments
    if not isinstance(df, pd.DataFrame):
        raise ValueError(f'df must be a DataFrame: {type(df)}')

    if list_column is None or not isinstance(list_column, str):
        raise ValueError(f'list_column must be a string: {type(list_column)}')

    if list_column not in df.columns:
        raise ValueError(f'list_column must be a column in df: {list_column}')

    # Subset by caller_order
    if caller_order is not None and caller_order_subset:

        if not isinstance(caller_order, (list, tuple)):
            raise ValueError(f'caller_order must be a list or tuple: {type(caller_order)}')

        caller_set = set(caller_order)

        df[list_column] = df[list_column].apply(lambda vals: ','.join([val for val in vals.split(',') if val in caller_set]))

        df = df.loc[df[list_column].str.len() > 0]

    # Count caller sets
    set_counts = df.groupby(list_column)[list_column].count()

    set_counts = set_counts.sort_values(ascending=False)

    if n_bar is not None:
        set_counts = set_counts[:n_bar]

    # Get counts for each stacked bar
    if color_column is not None:

        if color_order is None:
            color_order = []
        else:
            color_order = color_order.copy()

        # Normalize, ensure all values in color_column are represented once
        color_order_seen = set()

        color_order = [
            val for val in color_order if not (val in color_order_seen or color_order_seen.add(val))
        ]

        unordered_color = sorted([val for val in set(df[color_column]) if val not in color_order_seen])

        color_order = color_order + unordered_color

        # Get counts per color
        set_counts_color = {
            color_val: df.loc[
                df[color_column] == color_val
                ].groupby(
                list_column
            )[list_column].count(
            ).reindex(
                set_counts.index
            ).fillna(
                0
            ).astype(int).values
            for color_val in color_order
        }

        # set_counts_color = {
        #     color_val: data_table.loc[
        #         data_table[color_column] == color_val
        #     ].groupby(list_column)[list_column].count()[
        #         set_counts.index
        #     ].fillna(0).astype(int).values
        #         for color_val in color_order
        # }


        # Set all colors and labels
        if color_label is None:
            color_label = dict()
        else:
            color_label = color_label.copy()

        for val in set_counts_color.keys():
            if val not in color_label:
                color_label[val] = val

        if color is None:
            color = dict()
        else:
            color = color.copy()

        for val in set_counts_color.keys():
            if val not in color:
                color[val] = 'black'

    else:
        color_order = ['ALL']
        color_label = {'ALL': 'All'}
        color = {'ALL': 'Black'}

        set_counts_color = {'ALL': set_counts.values}

    bar_labels = np.array(set_counts.index)


    # Count total variants per caller (among all caller groups)

    df_caller_count = pd.Series(
        collections.Counter(
            [caller for caller_set in df[list_column] for caller in caller_set.split(',')]
        )
    ).sort_values(ascending=False)
    
    # Set order of variant calls
    if caller_order is None:
        df_caller_order = list(df_caller_count.index)
    else:
        df_caller_order = list(reversed([val for val in caller_order if val in df_caller_count.index] + [val for val in df_caller_count.index if val not in caller_order]))
    
    # Set caller labels
    if caller_label_dict != None:
        df_caller_label = [caller_label_dict.get(val, val) for val in df_caller_order]
    else:
        df_caller_label = df_caller_order


    # Make upset bubble matrix points

    set_items = [
        set(caller_set.split(',')) for caller_set in set_counts.index
    ]

    bubble_x = np.asarray(
        list(range(bar_labels.shape[0])) * len(df_caller_order)
    )

    bubble_y = np.asarray(
        [
            val for val_list in [[n] * bar_labels.shape[0]
                for n in range(len(df_caller_order))]
                for val in val_list
        ]
    )

    bubble_col = np.asarray(
        [
            'black' if df_caller_order[y] in set_items[x] else 'lightgray'
                for x, y in zip(bubble_x, bubble_y)
        ]
    )


    ### Make plot ###

    fig = plt.figure(figsize=(width, height), dpi=dpi)

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], figure=fig)

    if title is not None:
        fig.suptitle(title, fontsize=16)


    ## Variant count per caller combination bar plot (ax1) ##

    ax1 = plt.subplot(gs[0])

    # Make bars
    bar_bottom = np.zeros(bar_labels.shape[0]).astype(int)  # Bottom of next stacked bar

    for val in color_order:
        ax1.bar(
            bar_labels, set_counts_color[val],
            bottom=bar_bottom,
            color=color[val], label=color_label[val],

        )

        bar_bottom += set_counts_color[val]

    # Float counts over bars
    count_shift = 0.02 * np.max(set_counts.values)

    for index in range(bar_labels.shape[0]):
        cat_label = bar_labels[index]
        ax1.text(
            index,
            set_counts[cat_label] + count_shift,
            '{:,}'.format(set_counts[cat_label]),
            horizontalalignment='left',
            rotation=45,
            rotation_mode='anchor'
        )

    # Format axis
    ax1.get_xaxis().set_visible(False)

    ax1.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax1.set_ylabel('Variants (n)', fontsize=12, labelpad=20)

    for loc in ('top', 'right'):
        ax1.spines[loc].set_visible(False)

    # Legend
    if legend:
        ax1.legend(prop={'size': 18}, title=legend_title, title_fontsize=16)

    # Shift limits (make y-space for bar counts)
    y_lim = ax1.get_ylim()

    ax1.set_ylim(y_lim[0], y_lim[1] * 1.04)


    ## Caller combination dots (ax2) ##

    ax2 = plt.subplot(gs[1], sharex=ax1)

    ax2.scatter(
        'x', 'y', c='c', s=bubble_size,
        data={'x': bubble_x, 'y': bubble_y, 'c': bubble_col}
    )

    # Format axis & remove spines (lines around bubbles)
    ax2.get_xaxis().set_visible(False)

    for loc in ('top', 'bottom', 'left', 'right'):
        ax2.spines[loc].set_visible(False)

    # Set labels
    ax2.set_yticks(np.arange(len(df_caller_label)))
    ax2.set_yticklabels(np.array(df_caller_label), fontsize=bubble_label_size)

    ax2.tick_params('y', length=0)

    # Return figure
    return fig
