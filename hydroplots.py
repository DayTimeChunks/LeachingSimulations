import pylab
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def stackdata9(time, l1, l2, l3, l4, l5, l6 , l7, l8, l9):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6 , l7, l8, l9))
    return data

def stackdata6(time, l1, l2, l3, l4, l5, l6 ):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6))
    return data

def add_margin(ax, x=0.05, y=0.05):
    # This will, by default, add 5% to the x and y margins. You
    # can customise this using the x and y arguments when you call it.

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xmargin = (xlim[1]-xlim[0])*x
    ymargin = (ylim[1]-ylim[0])*y

    ax.set_xlim(xlim[0]-xmargin, xlim[1]+xmargin)
    ax.set_ylim(ylim[0]-ymargin, ylim[1]+ymargin)


def hydroplot(time,
              y1, y2, y3,
              y4, y5, y6,
              y7, y8, y9):
    sns.set(style="darkgrid")
    # sns.set(style="whitegrid")
    intesities = ['P(135 mm/h)', 'P(55 mm/h)', 'P(30 mm/h)',
                  ]
    color_sequence = ['#d62728', '#2ca02c', '#1f77b4']
    fig, ax1 = plt.subplots()
    ax1.plot(time, y1, color_sequence[0], label=intesities[0])
    ax1.plot(time, y2, color_sequence[1], label=intesities[1])
    ax1.plot(time, y3, color_sequence[2], label=intesities[2])
    ax1.plot([6.0], [28.35], color_sequence[0], marker='o')
    ax1.plot([6.0], [28.35], color_sequence[1], marker='o')
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Volume (mL)')

    plt.legend(loc='upper left')
    plt.show()


def hydroplot2(data):
    """
    :param data: see: leach_hydrology.py. Need to ensure data
    data[:, 0] = time
    data[:, 1] = cum infiltration (135 mm/h)
    data[:, 2] = cum infiltration (55 mm/h)
    data[:, 3] = cum infiltration (30 mm/h)
    data[:, 4] = cum leached (135 mm/h)
    data[:, 5] = cum leached (55 mm/h)
    data[:, 6] = cum leached (30 mm/h)

    :return:
    Plot with cumu. inf. and cum. leached volumes [mL]
    """
    sns.set(style="whitegrid")
    intesities = ['Inf. (135 mm/h)', 'Inf. (55 mm/h)', 'Inf. (30 mm/h)',
                  'Leach (135 mm/h)', 'Leach (55 mm/h)', 'Leach (30 mm/h)']
    color_sequence = ['#d62728', '#2ca02c', '#1f77b4']

    filled_markers = ['d', 'd', 'd',
                      'o', 'o', 'o',
                      # '.', '.', '.',
                      ]

    fig, ax1 = plt.subplots()
    c = 0
    for i in range(1, len(data[2])):
        if c < 3:
            ax1.plot(data[:, 0], data[:, i], color_sequence[c], marker=filled_markers[i-1], label=intesities[i-1])
            c += 1
        else:
            c = 0
            ax1.plot(data[:, 0], data[:, i], color_sequence[c], marker=filled_markers[i-1], label=intesities[i-1])
            c += 1

    # ax1.plot([6.0], [28.35], color_sequence[0], marker='s')

    """ Lab results """
    # Minutes
    six = np.array([6, 6, 6, 6, 6, 6, 6, 6])
    twelve = np.array([12, 12, 12, 12, 12, 12, 12, 12])
    thirty = np.array([30, 30, 30, 30, 30, 30, 30, 30])

    leach_high_6min = np.array([16.253, 12.952, 13.609, 13.610, 17.536, 14.29, 17.676, 17.705])  # all at 6 min
    leach_med_12min = np.array([10.089, 5.902, 13.787, 11.112, 13.981, 10.602, 11.858, 11.294])  # all at 12 min
    leach_med_30min = np.array([49.197, 40.402, 48.185, 46.402, 45.772, 47.201, 48.164, 47.032])  # all at 30min
    leach_low_30min = np.array([20.037, 17.508, 22.595, 19.082, 22.376, 20.085, 21.285, 20.871])  # all at 30min

    ax1.plot(six, leach_high_6min, color='k', marker='s', linestyle='None')
    ax1.plot(twelve, leach_med_12min, color='k', marker='^', linestyle='None')
    ax1.plot(thirty, leach_med_30min, color='k', marker='^', linestyle='None')
    ax1.plot(thirty, leach_low_30min, color='k', marker='8', linestyle='None')

    # plt.axis((0, 30.1, 0, 160))

    add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Volume (mL)')

    plt.legend(loc='upper left')
    plt.show()


def copperplot(data):
    """
    :param data:
     index 0 = time
      index 1 to 3 = Sterile, Kd_min
       index 4 to 6 = Untreated, Kd_max
    :return: Plot showing cumulative mass discharge
    """
    sns.set(style="whitegrid")

    treatment = ['Sterile', 'Untreated']
    intesities = ['(135 mm/h)', '(55 mm/h)', '(30 mm/h)']

    treat_intes = ['Sterile (135 mm/h)', 'Sterile (55 mm/h)', 'Sterile (30 mm/h)',
                   'Untreated (135 mm/h)', 'Untreated (55 mm/h)', 'Untreated (30 mm/h)']

    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',
                      '#d62728', '#2ca02c', '#1f77b4']

    filled_markers = ['d', '^', '^', 'o','d', '^', '^', 'o']

    fig, ax1 = plt.subplots()

    #  Plot data, which is in numpy array format, with:
    #  index 0 = time
    #  index 1:end =
    c = 0
    for i in range(1, len(data[0])):
        if i < 4:
            ax1.plot(data[:, 0], data[:, i], color_sequence[c], label=treat_intes[i - 1])
            c += 1
        else:
            ax1.plot(data[:, 0], data[:, i], color_sequence[c], label=treat_intes[i - 1], linestyle='dashed')
            c += 1

    """ Lab results """
    cu_time = np.array([6, 12, 30, 30])  # Minutes
    cu_sol_sterile = np.array([11.29, 11.63, 306.80, 21.08])
    cu_sol_untreat = np.array([0.5, 1.405, 37.0, 1])

    ax1.plot(cu_time[0], cu_sol_sterile[0], color_sequence[0], marker='s', linestyle='None', label=treat_intes[0])
    ax1.plot(cu_time[1], cu_sol_sterile[1], color_sequence[1], marker='o', linestyle='None', label=treat_intes[1])
    ax1.plot(cu_time[2], cu_sol_sterile[2], color_sequence[1], marker='o', linestyle='None')
    ax1.plot(cu_time[3], cu_sol_sterile[3], color_sequence[2], marker='h', linestyle='None', label=treat_intes[2])

    ax1.plot(cu_time[0], cu_sol_untreat[0], color_sequence[0], marker='^', linestyle='None', label=treat_intes[3])
    ax1.plot(cu_time[1], cu_sol_untreat[1], color_sequence[1], marker='v', linestyle='None', label=treat_intes[4])
    ax1.plot(cu_time[2], cu_sol_untreat[2], color_sequence[1], marker='v', linestyle='None')
    ax1.plot(cu_time[3], cu_sol_untreat[3], color_sequence[2], marker='*', linestyle='None', label=treat_intes[5])

    # plt.axis((0, 30, 0, 400))
    # Update the limits using set_xlim and set_ylim
    add_margin(ax1, x=0.01, y=0.01)  # Call this after plt.subbplot

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Cu (mu.g)')

    plt.legend(loc='upper left')
    plt.show()
