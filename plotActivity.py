#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 22:44:36 2020

@author: mishugeb
"""

import json
import math
import random

import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.font_manager
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['figure.constrained_layout.use'] = True

colors = [
    "tab:green",
    "tab:purple",
    "tab:pink",
    "tab:orange",
    "tab:olive",
    "tab:brown",
    "tab:red",
    "tab:cyan",
    "deeppink",
    "orangered",
    "blueviolet",
    "chocolate",
    "darkgreen",
    "dodgerblue",
    "mediumvioletred",
    "salmon",
    "magenta",
    "sandybrown",
    "forestgreen",
    "royalblue",
    "orchid",
    "indigo",
    "darkseagreen",
    "blue",
    "palevioletred",
    "darkslateblue",
    "olivedrab",
    "cyan",
    "hotpink",
    "rebeccapurple",
    "lime",
]

color_code = pd.DataFrame(
    np.array(
        [
            ["SBS1", "#ebf5bc"],
            ["SBS5", "#63d69e"],
            ["SBS2", "#f8b6b3"],
            ["SBS13", "#f17fb2"],
            ["SBS3", "#c4abc4"],
            ["SBS4", "#bcf2f5"],
            ["SBS7a", "#b5d7f5"],
            ["SBS7b", "#9ecef7"],
            ["SBS7c", "#84bdf0"],
            ["SBS7d", "#6cb2f0"],
            ["SBS8", "#dfc4f5"],
            ["SBS9", "#acf2d0"],
            ["SBS10a", "#f2aeae"],
            ["SBS10b", "#f08080"],
            ["SBS17a", "#d9f7b0"],
            ["SBS17b", "#8cc63f"],
            ["SBS40", "#c4c4f5"],
            ["SBS6", "#faf1dc"],
            ["SBS14", "#faecca"],
            ["SBS15", "#fcebc2"],
            ["SBS20", "#fae4af"],
            ["SBS21", "#fae1a5"],
            ["SBS26", "#fcde97"],
            ["SBS44", "#fad682"],
            ["SBS89", "darkslateblue"],
            ["SBS88", "olivedrab"],
            ["SBS87", "cyan"],
            ["SBS32", "hotpink"],
            ["SBS30", "rebeccapurple"],
            ["SBS25", "lime"],
            ["SBS19", "forestgreen"],


        ],
        dtype=object,
    ),
    columns=["signature", "color"],
)

artifacts_sigs = pd.DataFrame(
    np.array(
        [
            ["SBS27", "#C8C8C8"],
            ["SBS43", "#C0C0C0"],
            ["SBS45", "#BEBEBE"],
            ["SBS46", "#B8B8B8"],
            ["SBS47", "#B0B0B0"],
            ["SBS48", "#A9A9A9"],
            ["SBS49", "#A8A8A8"],
            ["SBS50", "#A0A0A0"],
            ["SBS51", "#989898"],
            ["SBS52", "#909090"],
            ["SBS53", "#888888"],
            ["SBS54", "#808080"],
            ["SBS55", "#787878"],
            ["SBS56", "#707070"],
            ["SBS57", "#696969"],
            ["SBS58", "#686868"],
            ["SBS59", "#606060"],
            ["SBS60", "#585858"],
        ],
        dtype=object,
    ),
    columns=["signature", "color"],
)


########## If read from file ############
# color_code = pd.read_table('color_dic.txt')
# artifacts_sigs = pd.read_table('artifacts_sigs.txt')
# with open('additional_colors.txt', 'r') as f: colors = f.read().splitlines()
##########################################
# function to plot Sample Activities
def plotActivity(
    inputDF, output_file="Activity_in_samples.pdf", log=False):
#    inputDF = inputDF.loc[:, (inputDF != 0).any(axis=0)]
    all_sig = list(inputDF.columns.values)
    s1 = color_code[color_code["signature"].isin(all_sig)]["signature"].tolist()
    s3 = artifacts_sigs[artifacts_sigs["signature"].isin(all_sig)]["signature"].tolist()
    a = [x for x in all_sig if x not in s1]
    s2 = [x for x in a if x not in s3]
    signature_list = s1 + s2 + s3
    if len(s2) <= len(colors):
        # print("Haha! we have defined all colors")
        color_list = (
            color_code[color_code["signature"].isin(all_sig)]["color"].tolist()
            + colors[: len(s2)]
            + artifacts_sigs[artifacts_sigs["signature"].isin(all_sig)][
                "color"
            ].tolist()
        )
    else:
        # print("Generating random colors...")
        rand_colors_list = []
        for i in range(0, len(s2) - len(colors)):
            rand_color = "#%06x" % random.randint(0, 0xFFFFFF)
            rand_colors_list += [rand_color]
        color_list = (
            color_code[color_code["signature"].isin(all_sig)]["color"].tolist()
            + colors[: len(s2)]
            + rand_colors_list
            + artifacts_sigs[artifacts_sigs["signature"].isin(all_sig)][
                "color"
            ].tolist()
        )

    # Start plotting
    nsamples = len(inputDF)
    names = inputDF.index.values.tolist()
    x_pos1 = list(map(lambda x: x + 0.2, list(range(0, len(names)))))
    x_pos = list(range(0, len(names)))
    barWidth = 1
 
    plt.rc("axes", edgecolor="gray")
    bottom_bars = []
    plot_list = []
    plt.xlim([-0.5, len(names) - 0.5])
#    if inputDF.values.max() > 1:
#        plt.ylim([0,6000])
#    else: plt.ylim([0,1])
    sig_activity_list = []
    for i in range(0, len(signature_list)):
        sig_activity_list.append(inputDF[signature_list[i]].tolist())

    for i in range(0, len(signature_list)):
        if i == 0:
            plot_list.append(
                plt.bar(
                    x_pos,
                    sig_activity_list[i],
                    color=color_list[i],
                    edgecolor="white",
                    width=barWidth,
                )
            )
        elif i == 1:
            bottom_bars = sig_activity_list[0]
            plot_list.append(
                plt.bar(
                    x_pos,
                    sig_activity_list[i],
                    bottom=bottom_bars,
                    color=color_list[i],
                    edgecolor="white",
                    width=barWidth,
                )
            )
        else:
            bottom_bars = np.add(bottom_bars, sig_activity_list[i - 1]).tolist()
            plot_list.append(
                plt.bar(
                    x_pos,
                    sig_activity_list[i],
                    bottom=bottom_bars,
                    color=color_list[i],
                    edgecolor="white",
                    width=barWidth,
                )
            )


    pos_df = pd.DataFrame()
    pos_df['pos'] = x_pos
    pos_df['name'] = names
    pos_df_grp = pos_df.groupby('name',sort=False)['pos'].min().reset_index()

    plt.xticks(pos_df_grp['pos'], pos_df_grp['name'], rotation=90, ha="right", rotation_mode="anchor")
    plt.legend(
        plot_list,
        signature_list,
        bbox_to_anchor=(1.05, 0.5),
        loc="center left",
        borderaxespad=0.0,
        frameon=False
    )
    plt.ylim([0,8000])
    plt.ylabel("Number of mutations in each signature")
    plt.savefig(output_file)
    plt.show()
    plt.close()
