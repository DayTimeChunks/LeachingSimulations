import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def stackdata9(time, l1, l2, l3, l4, l5, l6, l7, l8, l9):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6, l7, l8, l9))
    return data


def stackdata6(time, l1, l2, l3, l4, l5, l6):
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


def hydroplot(data, leach_high_6min, leach_med_12min, leach_med_30min, leach_low_30min):
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
                  'Leach Vol. (135 mm/h)', 'Leach Vol. (55 mm/h)', 'Leach Vol. (30 mm/h)']
    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',
                      '#d62728', '#2ca02c', '#1f77b4']

    filled_markers = ['d', 'd', 'd',
                      'o', 'o', 'o',
                      # '.', '.', '.',
                      ]
    line_styles = ['dashed', 'dashed', 'dashed']

    fig, ax1 = plt.subplots()
    c = 0
    for i in range(1, len(data[2])):
        if c < 3:
            # ax1.plot(data[:, 0], data[:, i], color_sequence[c], label=intesities[i-1])
            c += 1
        else:
            # c = 0
            ax1.plot(data[:, 0], data[:, i], color_sequence[c], linestyle='dashed', label=intesities[i-1])
            c += 1
            # marker=filled_markers[i-1]

    # ax1.plot([6.0], [28.35], color_sequence[0], marker='s')

    """ Lab results """
    # Minutes
    six = np.array([6, 6, 6, 6])
    twelve = np.array([12, 12, 12, 12])
    thirty = np.array([30, 30, 30, 30])

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


def pestiplot(data, obs_sol_sterile, obs_sol_untreat, title):
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

    treat_intens = ['Sim. Ster. (135 mm/h)', 'Sim. Ster. (55 mm/h)', 'Sim. Ster. (30 mm/h)',
                   'Sim. Unt. (135 mm/h)', 'Sim. Unt. (55 mm/h)', 'Sim. Unt. (30 mm/h)']
    obs_intens = ['Obs. Ster. (135 mm/h)', 'Obs. Ster. (55 mm/h)', 'Obs. Ster. (30 mm/h)',
                   'Obs. Unt. (135 mm/h)', 'Obs. Unt. (55 mm/h)', 'Obs. Unt. (30 mm/h)']

    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',
                      '#d62728', '#2ca02c', '#1f77b4']

    filled_markers = ['d', '^', '^', 'o', 'd', '^', '^', 'o']

    fig, ax1 = plt.subplots()

    #  Plot data, which is in numpy array format, with:
    #  index 0 = time
    #  index 1:end =
    c = 0
    for i in range(1, len(data[0][0])):
        if i < 4:
            ax1.plot(data[0][:, 0], data[0][:, i], color_sequence[c], label=treat_intens[i - 1])
            c += 1
        else:
            ax1.plot(data[0][:, 0], data[0][:, i], color_sequence[c], label=treat_intens[i - 1], linestyle='dashed')
            c += 1

    """ Lab results """
    time = np.array([6, 12, 30, 30])  # Minutes

    ax1.plot(time[0], obs_sol_sterile[0], color_sequence[0], marker='o', linestyle='None', label=obs_intens[0])
    ax1.plot(time[1], obs_sol_sterile[1], color_sequence[1], marker='o', linestyle='None', label=obs_intens[1])
    ax1.plot(time[2], obs_sol_sterile[2], color_sequence[1], marker='o', linestyle='None')
    ax1.plot(time[3], obs_sol_sterile[3], color_sequence[2], marker='o', linestyle='None', label=obs_intens[2])

    ax1.plot(time[0], obs_sol_untreat[0], color_sequence[0], marker='v', linestyle='None', label=obs_intens[3])
    ax1.plot(time[1], obs_sol_untreat[1], color_sequence[1], marker='v', linestyle='None', label=obs_intens[4])
    ax1.plot(time[2], obs_sol_untreat[2], color_sequence[1], marker='v', linestyle='None')
    ax1.plot(time[3], obs_sol_untreat[3], color_sequence[2], marker='v', linestyle='None', label=obs_intens[5])

    # plt.axis((0, 30, 0, 400))
    # Update the limits using set_xlim and set_ylim
    add_margin(ax1, x=0.01, y=0.01)  # Call this after plt.subbplot

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('mu.g')
    plt.title(title)

    plt.legend(loc='upper left')
    plt.show()


def pestiplot_inst(data, obs_sol_sterile, obs_sol_untreat, title):
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

    treat_intes = ['Sim. Ster. (135 mm/h)', 'Sim. Ster. (55 mm/h)', 'Sim. Ster. (30 mm/h)',
                   'Sim. Unt. (135 mm/h)', 'Sim. Unt. (55 mm/h)', 'Sim. Unt. (30 mm/h)']

    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',
                      '#d62728', '#2ca02c', '#1f77b4']

    filled_markers = ['d', '^', '^', 'o', 'd', '^', '^', 'o']

    fig, ax1 = plt.subplots()

    #  Plot data, which is in numpy array format, with:
    #  index 0 = time
    #  index 1:end =
    c = 0
    for i in range(1, len(data[1][0])):
        if i < 4:
            ax1.plot(data[1][:, 0], data[1][:, i], color_sequence[c], label=treat_intes[i - 1])
            c += 1
        else:
            ax1.plot(data[1][:, 0], data[1][:, i], color_sequence[c], label=treat_intes[i - 1], linestyle='dashed')
            c += 1

    """ Lab results """
    time = np.array([6, 12, 30, 30])  # Minutes
    """
    ax1.plot(time[0], obs_sol_sterile[0], color_sequence[0], marker='s', linestyle='None', label='Cum. Ster. (135mm/h)')
    ax1.plot(time[1], obs_sol_sterile[1], color_sequence[1], marker='o', linestyle='None', label='Cum. Ster. (55mm/h)')
    ax1.plot(time[2], obs_sol_sterile[2], color_sequence[1], marker='o', linestyle='None')
    ax1.plot(time[3], obs_sol_sterile[3], color_sequence[2], marker='h', linestyle='None', label='Cum. Ster. (30mm/h)')

    ax1.plot(time[0], obs_sol_untreat[0], color_sequence[0], marker='^', linestyle='None', label='Cum. Unt. (135mm/h)')
    ax1.plot(time[1], obs_sol_untreat[1], color_sequence[1], marker='v', linestyle='None', label='Cum. Unt. (55mm/h)')
    ax1.plot(time[2], obs_sol_untreat[2], color_sequence[1], marker='v', linestyle='None')
    ax1.plot(time[3], obs_sol_untreat[3], color_sequence[2], marker='*', linestyle='None', label='Cum. Unt. (30mm/h)')
    """
    # plt.axis((0, 30, 0, 400))
    # Update the limits using set_xlim and set_ylim
    add_margin(ax1, x=0.01, y=0.01)  # Call this after plt.subbplot

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('mu.g/dt')
    plt.title(title)

    plt.legend(loc='upper left')
    plt.show()
