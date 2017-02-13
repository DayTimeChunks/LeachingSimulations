
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import mpld3
from mpld3 import plugins
sns.set_context(rc={'lines.markeredgewidth': 0.1})


def add_margin(ax, x=0.05, y=0.05):
    # This will, by default, add 5% to the x and y margins. You
    # can customise this using the x and y arguments when you call it.

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xmargin = (xlim[1]-xlim[0])*x
    ymargin = (ylim[1]-ylim[0])*y

    ax.set_xlim(xlim[0]-xmargin, xlim[1]+xmargin)
    ax.set_ylim(ylim[0]-ymargin, ylim[1]+ymargin)


def hydroplot(data,
              y1, y2, y3,
              obs_high_6min,
              obs_med_12min, obs_med_30min,
              obs_low_30min,
              title):
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

    y_var = [y1, y2, y3]

    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',  # red, green, blue
                      '#d62728', '#2ca02c', '#1f77b4']  # red, green, blue

    # Convert from mm3 to cm3 all volume hydro data.
    fig, ax1 = plt.subplots()
    ax1.plot(data[:, 0], data[:, 1]/10**3, color_sequence[0], linestyle='dashed', label=y_var[1 - 1])
    ax1.plot(data[:, 0], data[:, 2]/10**3, color_sequence[1], linestyle='dashed', label=y_var[2 - 1])
    ax1.plot(data[:, 0], data[:, 3]/10**3, color_sequence[2], linestyle='dashed', label=y_var[3 - 1])

    """ Lab results """

    soil_modal = ['Sterile', 'Untreated',
                  'Ster. Aged', 'Untr. Aged']
    # Minutes
    # [sterile, untreat, sterile_aged, untreat_aged]
    six1 = np.array([5.8])
    six = np.array([6.1])

    twelve1 = np.array([11.8])
    twelve = np.array([12.1])

    thirty1 = np.array([29.8])
    thirty = np.array([30])

    # [sterile, untreat, sterile_aged, untreat_aged]
    ax1.plot(six, obs_high_6min[0], color='#d62728', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(six, obs_high_6min[1], color='#d62728', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(six1, obs_high_6min[2], color='#d62728', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(six1, obs_high_6min[3], color='#d62728', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(twelve, obs_med_12min[0], color='#2ca02c', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(twelve, obs_med_12min[1], color='#2ca02c', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(twelve1, obs_med_12min[2], color='#2ca02c', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(twelve1, obs_med_12min[3], color='#2ca02c', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(thirty, obs_med_30min[0], color='#2ca02c', marker='v', linestyle='None')
    ax1.plot(thirty1, obs_med_30min[1], color='#2ca02c', marker='o', linestyle='None')
    ax1.plot(thirty, obs_med_30min[2], color='#2ca02c', marker='^', linestyle='None')
    ax1.plot(thirty1, obs_med_30min[3], color='#2ca02c', marker='s', linestyle='None')

    ax1.plot(thirty, obs_low_30min[0], color='#1f77b4', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty1, obs_low_30min[1], color='#1f77b4', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(thirty, obs_low_30min[2], color='#1f77b4', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(thirty1, obs_low_30min[3], color='#1f77b4', marker='s', linestyle='None', label=soil_modal[3])

    add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Volume (mL)')

    plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # fig.subplots_adjust(right=0.7)
    plt.title(title)
    plt.show()


def hydroplot3(
        data,
        y1, y2, y3, y4,
        obs_high_6min, obs_med_12min, obs_med_30min, obs_low_30min,
        title):
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

    y_var = [y1, y2, y3, y4]

    color_sequence = ['#d62728', 'darkviolet', '#2ca02c', '#1f77b4',  # red, green, blue
                      '#d62728', 'darkviolet', '#2ca02c', '#1f77b4']  # red, green, blue

    # Convert from mm3 to cm3 all volume hydro data.
    fig, ax1 = plt.subplots()
    ax1.plot(data[:, 0], data[:, 1]/10**3, color_sequence[0], linestyle='solid', label=str(y_var[0]+ " (SF)"))
    ax1.plot(data[:, 0], data[:, 2]/10**3, color_sequence[1], linestyle='solid', label=y_var[1]+ " (SF)")
    ax1.plot(data[:, 0], data[:, 3]/10**3, color_sequence[2], linestyle='solid', label=y_var[2]+ " (SF)")
    ax1.plot(data[:, 0], data[:, 4]/10 ** 3, color_sequence[3], linestyle='solid', label=y_var[3]+ " (SF)")

    ax1.plot(data[:, 0], data[:, 5] / 10 ** 3, color_sequence[0], linestyle='dashed', label=y_var[0] + " (SA)")
    ax1.plot(data[:, 0], data[:, 6] / 10 ** 3, color_sequence[1], linestyle='dashed', label=y_var[1] + " (SA)")
    ax1.plot(data[:, 0], data[:, 7] / 10 ** 3, color_sequence[2], linestyle='dashed', label=y_var[2] + " (SA)")
    ax1.plot(data[:, 0], data[:, 8] / 10 ** 3, color_sequence[3], linestyle='dashed', label=y_var[3] + " (SA)")

    ax1.plot(data[:, 0], data[:, 9] / 10 ** 3, color_sequence[0], linestyle='dashdot', label=y_var[0] + " (LF)")
    ax1.plot(data[:, 0], data[:, 10] / 10 ** 3, color_sequence[1], linestyle='dashdot', label=y_var[1] + " (LF)")
    ax1.plot(data[:, 0], data[:, 11] / 10 ** 3, color_sequence[2], linestyle='dashdot', label=y_var[2] + " (LF)")
    ax1.plot(data[:, 0], data[:, 12] / 10 ** 3, color_sequence[3], linestyle='dashdot', label=y_var[3] + " (LF)")

    ax1.plot(data[:, 0], data[:, 13] / 10 ** 3, color_sequence[0], linestyle='dotted', label=y_var[0] + " (LA)")
    ax1.plot(data[:, 0], data[:, 14] / 10 ** 3, color_sequence[1], linestyle='dotted', label=y_var[1] + " (LA)")
    ax1.plot(data[:, 0], data[:, 15] / 10 ** 3, color_sequence[2], linestyle='dotted', label=y_var[2] + " (LA)")
    ax1.plot(data[:, 0], data[:, 16] / 10 ** 3, color_sequence[3], linestyle='dotted', label=y_var[3] + " (LA)")

    """ Lab results """
    soil_modal = ['Sterile Fresh', 'Live Fresh',
                  'Ster. Aged', 'Live Aged']
    # Minutes
    # [sterile, untreat, sterile_aged, untreat_aged]
    six1 = np.array([5.8])
    six = np.array([6.1])

    twelve1 = np.array([11.8])
    twelve = np.array([12.1])

    thirty1 = np.array([29.8])
    thirty = np.array([30])

    ax1.plot(six, obs_high_6min[0], color='#d62728', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(six, obs_high_6min[1], color='#d62728', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(twelve, obs_med_12min[0], color='darkviolet', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(twelve, obs_med_12min[1], color='darkviolet', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(thirty, obs_med_30min[0], color='#2ca02c', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty1, obs_med_30min[1], color='#2ca02c', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(thirty, obs_low_30min[0], color='#1f77b4', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty1, obs_low_30min[1], color='#1f77b4', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(six1, obs_high_6min[2], color='#d62728', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(six1, obs_high_6min[3], color='#d62728', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(twelve1, obs_med_12min[2], color='darkviolet', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(twelve1, obs_med_12min[3], color='darkviolet', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(thirty, obs_med_30min[2], color='#2ca02c', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(thirty1, obs_med_30min[3], color='#2ca02c', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(thirty, obs_low_30min[2], color='#1f77b4', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(thirty1, obs_low_30min[3], color='#1f77b4', marker='s', linestyle='None', label=soil_modal[3])

    add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Volume (mL)')

    # plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    ax1.legend(bbox_to_anchor=(1.65, 0.9), fancybox=True, framealpha=0.7, ncol=2)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # fig.subplots_adjust(right=0.7)

    plt.show()


def hydroplot2(
        data,
        y1, y2, y3, y4,
        obs_high_6min,
        obs_med_12min, obs_med_30min,
        obs_low_30min,
        title,
        AGED = False):
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

    y_var = [y1, y2, y3, y4]

    color_sequence = ['#d62728', 'darkviolet', '#2ca02c', '#1f77b4',  # red, green, blue
                      '#d62728', 'darkviolet', '#2ca02c', '#1f77b4']  # red, green, blue

    # Convert from mm3 to cm3 all volume hydro data.
    fig, ax1 = plt.subplots()
    ax1.plot(data[:, 0], data[:, 1]/10**3, color_sequence[0], linestyle='dashed', label=y_var[1 - 1])
    ax1.plot(data[:, 0], data[:, 2]/10**3, color_sequence[1], linestyle='dashed', label=y_var[2 - 1])
    ax1.plot(data[:, 0], data[:, 3]/10**3, color_sequence[2], linestyle='dashed', label=y_var[3 - 1])
    ax1.plot(data[:, 0], data[:, 4]/10 ** 3, color_sequence[3], linestyle='dashed', label=y_var[4 - 1])

    """ Lab results """

    soil_modal = ['Sterile Fresh', 'Live Fresh',
                  'Ster. Aged', 'Live Aged']
    # Minutes
    # [sterile, untreat, sterile_aged, untreat_aged]
    six1 = np.array([5.8])
    six = np.array([6.1])

    twelve1 = np.array([11.8])
    twelve = np.array([12.1])

    thirty1 = np.array([29.8])
    thirty = np.array([30])

    if not AGED:
        ax1.plot(six, obs_high_6min[0], color='#d62728', marker='v', linestyle='None', label=soil_modal[0])
        ax1.plot(six, obs_high_6min[1], color='#d62728', marker='o', linestyle='None', label=soil_modal[1])

        ax1.plot(twelve, obs_med_12min[0], color='darkviolet', marker='v', linestyle='None', label=soil_modal[0])
        ax1.plot(twelve, obs_med_12min[1], color='darkviolet', marker='o', linestyle='None', label=soil_modal[1])

        ax1.plot(thirty, obs_med_30min[0], color='#2ca02c', marker='v', linestyle='None', label=soil_modal[0])
        ax1.plot(thirty1, obs_med_30min[1], color='#2ca02c', marker='o', linestyle='None', label=soil_modal[1])

        ax1.plot(thirty, obs_low_30min[0], color='#1f77b4', marker='v', linestyle='None', label=soil_modal[0])
        ax1.plot(thirty1, obs_low_30min[1], color='#1f77b4', marker='o', linestyle='None', label=soil_modal[1])
    else:
        ax1.plot(six1, obs_high_6min[2], color='#d62728', marker='^', linestyle='None', label=soil_modal[2])
        ax1.plot(six1, obs_high_6min[3], color='#d62728', marker='s', linestyle='None', label=soil_modal[3])

        ax1.plot(twelve1, obs_med_12min[2], color='darkviolet', marker='^', linestyle='None', label=soil_modal[2])
        ax1.plot(twelve1, obs_med_12min[3], color='darkviolet', marker='s', linestyle='None', label=soil_modal[3])

        ax1.plot(thirty, obs_med_30min[2], color='#2ca02c', marker='^', linestyle='None', label=soil_modal[2])
        ax1.plot(thirty1, obs_med_30min[3], color='#2ca02c', marker='s', linestyle='None', label=soil_modal[3])

        ax1.plot(thirty, obs_low_30min[2], color='#1f77b4', marker='^', linestyle='None', label=soil_modal[2])
        ax1.plot(thirty1, obs_low_30min[3], color='#1f77b4', marker='s', linestyle='None', label=soil_modal[3])

    add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Volume (mL)')

    plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # fig.subplots_adjust(right=0.7)
    if AGED:
        plt.title(title + " (Aged)")
    else:
        plt.title(title + " (Fresh)")
    plt.show()


def pestiplot_all(data,
                  obs_high_6min,
                  obs_med_12min, obs_med_30min,
                  obs_low_30min,
                  title,
                  pest_name):
    """
    :param data: length 7, where:
     index 0 = time
      index 1 to 3 = Sterile, Kd_min
       index 4 to 6 = Untreated, Kd_max
    :return: Plot showing mass output (2 curves per soil)
    """
    sns.set(style="whitegrid")

    treat_intens = ['Sim. Ster. (135 mm/h)', 'Sim. Ster. (55 mm/h)', 'Sim. Ster. (30 mm/h)',
                    'Sim. Unt. (135 mm/h)', 'Sim. Unt. (55 mm/h)', 'Sim. Unt. (30 mm/h)']
    obs_intens = ['Obs. Ster. (135 mm/h)', 'Obs. Ster. (55 mm/h)', 'Obs. Ster. (30 mm/h)',
                  'Obs. Unt. (135 mm/h)', 'Obs. Unt. (55 mm/h)', 'Obs. Unt. (30 mm/h)']

    color_sequence = ['#d62728', '#2ca02c', '#1f77b4',
                      '#d62728', '#2ca02c', '#1f77b4']

    fig, ax1 = plt.subplots()

    #  Plot data, which is in numpy array format, with:
    c = 0
    for i in range(1, len(data[0])):
        if i < 4:
            ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[c], label=treat_intens[i - 1])
            c += 1
        else:
            ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[c], label=treat_intens[i - 1], linestyle='dashed')
            c += 1

    """ Lab results """
    soil_modal = ['Sterile', 'Untreated',
                  'Ster. Aged', 'Untr. Aged']

    six = np.array([6])
    twelve = np.array([12])
    thirty = np.array([30])

    six1 = np.array([5.8])
    twelve1 = np.array([11.8])
    thirty1 = np.array([29.8])

    ax1.plot(six, obs_high_6min[0], color='#d62728', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(six, obs_high_6min[1], color='#d62728', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(six1, obs_high_6min[2], color='#d62728', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(six1, obs_high_6min[3], color='#d62728', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(twelve, obs_med_12min[0], color='#2ca02c', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(twelve, obs_med_12min[1], color='#2ca02c', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(twelve1, obs_med_12min[2], color='#2ca02c', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(twelve1, obs_med_12min[3], color='#2ca02c', marker='s', linestyle='None', label=soil_modal[3])

    ax1.plot(thirty, obs_med_30min[0], color='#2ca02c', marker='v', linestyle='None')
    ax1.plot(thirty1, obs_med_30min[1], color='#2ca02c', marker='o', linestyle='None')
    ax1.plot(thirty, obs_med_30min[2], color='#2ca02c', marker='^', linestyle='None')
    ax1.plot(thirty1, obs_med_30min[3], color='#2ca02c', marker='s', linestyle='None')

    ax1.plot(thirty, obs_low_30min[0], color='#1f77b4', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty1, obs_low_30min[1], color='#1f77b4', marker='o', linestyle='None', label=soil_modal[1])
    ax1.plot(thirty, obs_low_30min[2], color='#1f77b4', marker='^', linestyle='None', label=soil_modal[2])
    ax1.plot(thirty1, obs_low_30min[3], color='#1f77b4', marker='s', linestyle='None', label=soil_modal[3])

    obs_high_6min,
    obs_med_12min, obs_med_30min,
    obs_low_30min,

    # plt.axis((0, 30, 0, 400))
    # Update the limits using set_xlim and set_ylim
    add_margin(ax1, x=0.01, y=0.01)  # Call this after plt.subbplot

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel(pest_name + ' [' + r'$\mu$' + 'g]')
    plt.title(title)

    plt.legend(loc='upper left', fancybox=True, framealpha=0.8)
    plt.show()


def pestiplot_condition(
        data,
        pest_dict,
        pest_name,
        soil_type,
        cycle,
        LEACH = True,
        STERILE = True
):
    """

    :param data:
        np.array with [time,
                        fresh_high, aged_high,
                        fresh_med12, aged_med12..
    :param pest_dict:
    :param pest_name:
    :param soil_type:
    :param LEACH:
    :param STERILE:
    :return:
    """

    sns.set(style="whitegrid")

    treat_intens = ['Fresh (135 mm/h) - A', 'Aged (135 mm/h) - A',
                    'Fresh (55 mm/h) - B', 'Aged - A (55 mm/h) - B',
                    'Fresh (55 mm/h) - C', 'Aged - B (55 mm/h) - C',
                    'Fresh (30 mm/h) - D', 'Aged (30 mm/h) - D'
                    ]

    obs_intens = ['Obs. Ster. (135 mm/h)', 'Obs. Ster. (55 mm/h)', 'Obs. Ster. (30 mm/h)',
                  'Obs. Unt. (135 mm/h)', 'Obs. Unt. (55 mm/h)', 'Obs. Unt. (30 mm/h)']

    color_sequence = ['#d62728', '#d62728',
                      'darkviolet', 'darkviolet',
                      '#2ca02c', '#2ca02c',
                      '#1f77b4', '#1f77b4']

    fig, ax1 = plt.subplots()

    #  Plot data, which is in numpy array format, with:

    if cycle == '1st pulse':
        pattern = 0
        for i in range(1, len(data[0])):
            pattern += 1
            if i % 2 != 0 and (pattern == 3):
                ax1.plot(data[:][:, 0], data[:][:, i]+1, color_sequence[i - 1], label=treat_intens[i - 1])
            elif i % 2 != 0:
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1])
            elif i % 2 == 0 and (pattern == 4):
                ax1.plot(data[:][:, 0], data[:][:, i]+1, color_sequence[i - 1], label=treat_intens[i - 1], linestyle='dashed')
            else:
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1], linestyle='dashed')
    else:
        for i in range(1, len(data[0])):
            if i % 2 != 0:
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1])
            else:
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1],
                         linestyle='dashed')



    """ Lab results """
    soil_modal = ['Fresh', 'Aged']

    six = np.array([6])
    twelve = np.array([12])

    six1 = np.array([5.9])
    twelve1 = np.array([11.9])

    thirty = np.array([30])
    thirty2 = np.array([29.5])
    thirty1 = np.array([29.6])
    thirty3 = np.array([29.9])

    if LEACH:
        dict_index = 1
    else:
        dict_index = 2

    x_values = [six, six1, twelve, twelve1, thirty, thirty1, thirty2, thirty3]
    y_values = []
    y_errors = []
    for index, (name, (mass_i, leach_obs, pond_obs, initial_mass_error, leach_error, pond_error)) in enumerate(sorted(pest_dict.items())):
        if LEACH:
            y_values.append(float(leach_obs))
            # if float(leach_obs) - float(leach_error) < 0:
            #    initial_mass_error = float(leach_error)
            y_errors.append(float(leach_error))
        else:
            y_values.append(float(pond_obs))
            y_errors.append(float(pond_error))

    for x, y, err, color in zip(x_values, y_values, y_errors, color_sequence):
        ax1.errorbar(x, y, err, lw=2, capsize=5, capthick=2, color=color)

        # ax1.errorbar(x_values, y_values, yerr=y_errors, linestyle='None')

    ax1.plot(six, float(pest_dict['a_high_0d'][dict_index]),
             color='#d62728', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(six1, float(pest_dict['b_high_1d'][dict_index]),
             color='#d62728', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(twelve, float(pest_dict['c_med12_0d'][dict_index]),
             color='darkviolet', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(twelve1, float(pest_dict['d_med12_1d'][dict_index]),
             color='darkviolet', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(thirty, float(pest_dict['e_med30_0d'][dict_index]),
             color='#2ca02c', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty1, float(pest_dict['f_med30_1d'][dict_index]),
             color='#2ca02c', marker='o', linestyle='None', label=soil_modal[1])

    ax1.plot(thirty2, float(pest_dict['g_low_0d'][dict_index]),
             color='#1f77b4', marker='v', linestyle='None', label=soil_modal[0])
    ax1.plot(thirty3, float(pest_dict['h_low_1d'][dict_index]),
             color='#1f77b4', marker='o', linestyle='None', label=soil_modal[1])

    # plt.axis((0, 30, 0, 400))
    # Update the limits using set_xlim and set_ylim
    add_margin(ax1, x=0.01, y=0.01)  # Call this after plt.subbplot

    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel(pest_name + ' [' + r'$\mu$' + 'g]')

    if LEACH and STERILE:
        plt.title('Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)')
    elif LEACH and not STERILE:
        plt.title('Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)')
    elif not LEACH and STERILE:
        plt.title('Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)')
    elif not LEACH and not STERILE:
        plt.title('Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)')
    else:
        print("Title error")

    # fig.plugins = [plugins.PointLabelTooltip(x_values, y_values)]
    plt.legend(loc='upper left', fancybox=True, framealpha=0.8)
    plt.show()


def stackdata36(time,
                a, b, c, d, e, f, g, h, i, j,
                k, l, m, n, o, p, q, r, s, t,
                u, v, w, x, y, z, aa, ab, ac, ad,
                af, ag, ah, ai, aj, ak):
    data = np.column_stack((
        time,
        a, b, c, d, e, f, g, h, i, j, k, l, m, n,
        o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad,
        af, ag, ah, ai, aj, ak))
    return data


def stackdata28(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o, p,
                q, r, s, t,
                u, v, w, x,
                y, z, aa, ab):
    data = np.column_stack((time,
                            a, b, c, d, e, f, g, h, i, j, k, l, m, n,
                            o, p, q, r, s, t, u, v, w, x, y, z, aa, ab))
    return data

def stackdata21(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o,
                p, q, r,
                s, t, u):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u))
    return data


def stackdata18(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o,
                p, q, r):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r))
    return data


def stackdata16(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o, p):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p))
    return data


def stackdata15(time,
                a, b, c, d,
                e, f, g, h,
                i, j, k, l,
                m, n, o):
    data = np.column_stack((time, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o))
    return data


def stackdata9(time, l1, l2, l3, l4, l5, l6, l7, l8, l9):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6, l7, l8, l9))
    return data


def stackdata8(time, l1, l2, l3, l4, l5, l6, l7, l8):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6, l7, l8))
    return data

def stackdata6(time, l1, l2, l3, l4, l5, l6):
    data = np.column_stack((time, l1, l2, l3, l4, l5, l6))
    return data


def stackdata4(time, a, b, c, d):
    data = np.column_stack((time, a, b, c, d))
    return data

def stackdata3(time, a, b, c):
    data = np.column_stack((time, a, b, c))
    return data

#######################
# OLD PLOTS, not in use
########################
def pestiplot(data, obs_sol_sterile, obs_sol_untreat, title):
    """
    :param data: length 7, where:
     index 0 = time
      index 1 to 3 = Sterile, Kd_min
       index 4 to 6 = Untreated, Kd_max
    :return: Plot showing mass output (2 curves per soil)
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
    c = 0
    for i in range(1, len(data[0])):
        if i < 4:
            ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[c], label=treat_intens[i - 1])
            c += 1
        else:
            ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[c], label=treat_intens[i - 1], linestyle='dashed')
            c += 1

    """ Lab results """
    soil_modal = ['Sterile', 'Untreated',
                  'Ster. Aged', 'Untr. Aged']

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
