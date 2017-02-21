import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from stackdata import *
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
    plt.savefig('../Figures/' + str(title) + '.png', dpi=600)
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

    plt.title(title)

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.3,  0.9), fancybox=True, framealpha=0.7, ncol=2)
    # bbox_to_anchor=(0.5, -0.1)
    # ax1.legend(bbox_to_anchor=(1.65, 0.9), fancybox=True, framealpha=0.7, ncol=2)

    # plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # fig.subplots_adjust(right=0.7)
    plt.savefig('../Figures/' + str(title) + '.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
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
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1])
            elif i % 2 != 0:
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1])
            elif i % 2 == 0 and (pattern == 4):
                ax1.plot(data[:][:, 0], data[:][:, i], color_sequence[i - 1], label=treat_intens[i - 1], linestyle='dashed')
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

    # fig.plugins = [plugins.PointLabelTooltip(x_values, y_values)]
    plt.legend(loc='upper left', fancybox=True, framealpha=0.8)

    if LEACH and STERILE:
        plt.title('Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)')
        plt.savefig('../Figures/' + str(pest_name) + '/' + 'Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)' + '.png', dpi=600)
    elif LEACH and not STERILE:
        plt.title('Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)')
        plt.savefig('../Figures/' + str(pest_name) + '/' + 'Leached mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)' + '.png', dpi=600)
    elif not LEACH and STERILE:
        plt.title('Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)')
        plt.savefig('../Figures/' + str(pest_name) + '/' + 'Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Sterile)' + '.png', dpi=600)
    elif not LEACH and not STERILE:
        plt.title('Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)')
        plt.savefig('../Figures/' + str(pest_name) + '/' + 'Ponded mass, ' + pest_name + ' ' + cycle + ' ' + soil_type + ' (Living)' + '.png', dpi=600 )
    else:
        print("Title error")

    plt.show()


def extract_and_plothydro(
        water_data,
        series_name1, series_name2, series_name3, series_name4,
        title,
        soil,
        isFirstCycle,
        isPercolation
):
    if soil == 'Alteck':
        if isFirstCycle:
            # Observed percolation
            # fresh, fresh, aged, aged
            # all at 6 min, high inetnesity
            leach_high_6min = np.array([16.253, 12.958, 17.536, 14.29])
            # all at 12 min, med intensity
            leach_med_12min = np.array([10.089, 5.902, 13.981, 10.602])
            # all at 30min, med intensity
            leach_med_30min = np.array([49.197, 40.402, 45.772, 47.201])
            # all at 30min, low intensity
            leach_low_30min = np.array([20.037, 17.508, 22.376, 20.085])


            # Time
            cum_time_30min = water_data[:, 0]

            # Cummulative infiltration
            cum_inf_135mmh = water_data[:, 4]
            cum_inf_55mmh = water_data[:, 5]
            cum_inf_30mmh = water_data[:, 6]

            # Cummulative leaching
            cum_leach_135mmh = water_data[:, 7]
            cum_leach_55mmh = water_data[:, 8]
            cum_leach_30mmh = water_data[:, 9]

            # Ponding
            roff_135mmh = water_data[:, 10]
            roff_55mmh = water_data[:, 11]
            roff_30mmh = water_data[:, 12]

            # Cummulative ponding
            cum_roff_135mmh = water_data[:, 13]
            cum_roff_55mmh = water_data[:, 14]
            cum_roff_30mmh = water_data[:, 15]

            infil_135mmh = water_data[:, 16]
            infil_55mmh = water_data[:, 17]
            infil_30mmh = water_data[:, 18]

            percol_data1 = stackdata3(cum_time_30min,
                                      cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

            runoff_data1 = stackdata3(cum_time_30min,
                                      cum_roff_135mmh, cum_roff_55mmh, cum_roff_30mmh)

            infil_data1 = stackdata3(cum_time_30min,
                                     infil_135mmh, infil_55mmh, infil_30mmh)

            time_size_135mmh = water_data[:, 19]
            time_size_55mmhA = water_data[:, 20]
            time_size_55mmhB = water_data[:, 20]
            time_size_30mmh = water_data[:, 21]

            time_sizes1 = [time_size_135mmh, time_size_135mmh,
                           time_size_55mmhA, time_size_55mmhA,
                           time_size_55mmhB, time_size_55mmhB,
                           time_size_30mmh, time_size_30mmh]

            return hydroplot(percol_data1,
                  series_name1, series_name2, series_name3,
                  leach_high_6min,
                  leach_med_12min, leach_med_30min,
                  leach_low_30min,
                  title)

        else:
            if isPercolation:
                # Observed percolation
                # Order if array is:
                #  [sterile, untreat, sterile_aged, untreat_aged]

                # At 6 min, high inetnesity
                leach_high_6min = np.array([14.192, 8.245, 2.410, 5.469])
                # At 12 min, med intensity
                leach_med_12min = np.array([18.672, 19.0, 0.830, 11.407])
                # At 30min, med intensity
                leach_med_30min = np.array([12.697, 2.473, 3.52, 20.291])
                # At 30min, low intensity
                leach_low_30min = np.array([29.656, 9.375, 0.409, 3.385])

                # Time axis
                cum_time_30min = water_data[0][:, 0]

                # Cumulative leachate
                cum_leach_135mmh_SF = water_data[0][:, 2]
                cum_leach_135mmh_SA = water_data[0][:, 4]
                cum_leach_135mmh_LF = water_data[0][:, 6]
                cum_leach_135mmh_LA = water_data[0][:, 8]

                cum_leach_55mmhA_SF = water_data[0][:, 10]
                cum_leach_55mmhA_SA = water_data[0][:, 12]
                cum_leach_55mmhA_LF = water_data[0][:, 14]
                cum_leach_55mmhA_LA = water_data[0][:, 16]

                cum_leach_55mmhB_SF = water_data[0][:, 18]
                cum_leach_55mmhB_SA = water_data[0][:, 20]
                cum_leach_55mmhB_LF = water_data[0][:, 22]
                cum_leach_55mmhB_LA = water_data[0][:, 24]

                cum_leach_30mmh_SF = water_data[0][:, 26]
                cum_leach_30mmh_SA = water_data[0][:, 28]
                cum_leach_30mmh_LF = water_data[0][:, 30]
                cum_leach_30mmh_LA = water_data[0][:, 32]

                # Group each compartment for graphing
                percol_data2 = stackdata16(
                    cum_time_30min,
                    cum_leach_135mmh_SF, cum_leach_55mmhA_SF, cum_leach_55mmhB_SF, cum_leach_30mmh_SF,
                    cum_leach_135mmh_SA, cum_leach_55mmhA_SA, cum_leach_55mmhB_SA, cum_leach_30mmh_SA,
                    cum_leach_135mmh_LF, cum_leach_55mmhA_LF, cum_leach_55mmhB_LF, cum_leach_30mmh_LF,
                    cum_leach_135mmh_LA, cum_leach_55mmhA_LA, cum_leach_55mmhB_LA, cum_leach_30mmh_LA)

                time_size_135mmh = water_data[0][:, 33]
                time_size_55mmhA = water_data[0][:, 34]
                time_size_55mmhB = water_data[0][:, 35]
                time_size_30mmh = water_data[0][:, 36]

                time_sizes2 = [time_size_135mmh, time_size_135mmh,
                               time_size_55mmhA, time_size_55mmhA,
                               time_size_55mmhB, time_size_55mmhB,
                               time_size_30mmh, time_size_30mmh]

                return hydroplot3(
                    percol_data2,
                    series_name1, series_name2, series_name3, series_name4,
                    leach_high_6min, leach_med_12min, leach_med_30min, leach_low_30min,
                    title
                )
            else:
                # Observed ponding
                # [sterile, untreat, sterile_aged, untreat_aged]
                # all at 6 min, high inetnesity
                roff_high_6min = np.array([10.824, 20.935, 24.75, 19.041])

                # all at 12 min, med intensity
                roff_med_12min = np.array([0, 3.907, 19.436, 7.313])

                # all at 30min, med intensity
                roff_med_30min = np.array([43.764, 28.911, 51.964, 33.478])

                # all at 30min, low intensity
                roff_low_30min = np.array([0, 22.618, 28.598, 27.314])

                # Time axis
                cum_time_30min = water_data[0][:, 0]

                # Ponding cumulative
                cum_roff_135mmh_SF = water_data[0][:, 1]
                cum_roff_135mmh_SA = water_data[0][:, 3]
                cum_roff_135mmh_LF = water_data[0][:, 5]
                cum_roff_135mmh_LA = water_data[0][:, 7]

                cum_roff_55mmhA_SF = water_data[0][:, 9]
                cum_roff_55mmhA_SA = water_data[0][:, 11]
                cum_roff_55mmhA_LF = water_data[0][:, 13]
                cum_roff_55mmhA_LA = water_data[0][:, 15]

                cum_roff_55mmhB_SF = water_data[0][:, 17]
                cum_roff_55mmhB_SA = water_data[0][:, 19]
                cum_roff_55mmhB_LF = water_data[0][:, 21]
                cum_roff_55mmhB_LA = water_data[0][:, 23]

                cum_roff_30mmh_SF = water_data[0][:, 25]
                cum_roff_30mmh_SA = water_data[0][:, 27]
                cum_roff_30mmh_LF = water_data[0][:, 29]
                cum_roff_30mmh_LA = water_data[0][:, 31]

                runoff_data2 = stackdata16(
                    cum_time_30min,
                    cum_roff_135mmh_SF, cum_roff_55mmhA_SF, cum_roff_55mmhB_SF, cum_roff_30mmh_SF,
                    cum_roff_135mmh_SA, cum_roff_55mmhA_SA, cum_roff_55mmhB_SA, cum_roff_30mmh_SA,
                    cum_roff_135mmh_LF, cum_roff_55mmhA_LF, cum_roff_55mmhB_LF, cum_roff_30mmh_LF,
                    cum_roff_135mmh_LA, cum_roff_55mmhA_LA, cum_roff_55mmhB_LA, cum_roff_30mmh_LA)

                time_size_135mmh = water_data[0][:, 33]
                time_size_55mmhA = water_data[0][:, 34]
                time_size_55mmhB = water_data[0][:, 35]
                time_size_30mmh = water_data[0][:, 36]

                time_sizes2 = [time_size_135mmh, time_size_135mmh,
                               time_size_55mmhA, time_size_55mmhA,
                               time_size_55mmhB, time_size_55mmhB,
                               time_size_30mmh, time_size_30mmh]

                return hydroplot3(
                    runoff_data2,
                    series_name1, series_name2, series_name3, series_name4,
                    roff_high_6min, roff_med_12min, roff_med_30min, roff_low_30min,
                    title)

    elif soil == 'Rouff':
        if isFirstCycle:
            leach_high_6min = np.array([13.609, 13.610, 17.676, 17.705])  # all at 6 min
            leach_med_12min = np.array([13.787, 11.112, 11.858, 11.294])  # all at 12 min
            leach_med_30min = np.array([48.185, 46.402, 48.164, 47.032])  # all at 30min
            leach_low_30min = np.array([22.595, 19.082, 21.285, 20.871])  # all at 30min

            # Time
            cum_time_30min = water_data[:, 0]

            # Cummulative infiltration
            cum_inf_135mmh = water_data[:, 4]
            cum_inf_55mmh = water_data[:, 5]
            cum_inf_30mmh = water_data[:, 6]

            # Cummulative leaching
            cum_leach_135mmh = water_data[:, 7]
            cum_leach_55mmh = water_data[:, 8]
            cum_leach_30mmh = water_data[:, 9]

            # Ponding
            roff_135mmh = water_data[:, 10]
            roff_55mmh = water_data[:, 11]
            roff_30mmh = water_data[:, 12]

            # Cummulative ponding
            cum_roff_135mmh = water_data[:, 13]
            cum_roff_55mmh = water_data[:, 14]
            cum_roff_30mmh = water_data[:, 15]

            infil_135mmh = water_data[:, 16]
            infil_55mmh = water_data[:, 17]
            infil_30mmh = water_data[:, 18]

            percol_data1 = stackdata3(cum_time_30min,
                                      cum_leach_135mmh, cum_leach_55mmh, cum_leach_30mmh)

            runoff_data1 = stackdata3(cum_time_30min,
                                      cum_roff_135mmh, cum_roff_55mmh, cum_roff_30mmh)

            infil_data1 = stackdata3(cum_time_30min,
                                     infil_135mmh, infil_55mmh, infil_30mmh)

            time_size_135mmh = water_data[:, 19]
            time_size_55mmhA = water_data[:, 20]
            time_size_55mmhB = water_data[:, 20]
            time_size_30mmh = water_data[:, 21]

            time_sizes1 = [time_size_135mmh, time_size_135mmh,
                           time_size_55mmhA, time_size_55mmhA,
                           time_size_55mmhB, time_size_55mmhB,
                           time_size_30mmh, time_size_30mmh]

            return hydroplot(percol_data1,
                             series_name1, series_name2, series_name3,
                             leach_high_6min,
                             leach_med_12min, leach_med_30min,
                             leach_low_30min,
                             title)

        else:
            if isPercolation:
                leach_high_6min = np.array([13.309, 0., 7.394, 6.549])
                leach_med_12min = np.array([0.958, 3.669, 16.06, 12.988])
                leach_med_30min = np.array([0.941, 18.601, 51.834, 29.232])
                leach_low_30min = np.array([10.157, 26.737, 27.533, 6.197])

                # Time axis
                cum_time_30min = water_data[0][:, 0]

                # Cumulative leachate
                cum_leach_135mmh_SF = water_data[0][:, 2]
                cum_leach_135mmh_SA = water_data[0][:, 4]
                cum_leach_135mmh_LF = water_data[0][:, 6]
                cum_leach_135mmh_LA = water_data[0][:, 8]

                cum_leach_55mmhA_SF = water_data[0][:, 10]
                cum_leach_55mmhA_SA = water_data[0][:, 12]
                cum_leach_55mmhA_LF = water_data[0][:, 14]
                cum_leach_55mmhA_LA = water_data[0][:, 16]

                cum_leach_55mmhB_SF = water_data[0][:, 18]
                cum_leach_55mmhB_SA = water_data[0][:, 20]
                cum_leach_55mmhB_LF = water_data[0][:, 22]
                cum_leach_55mmhB_LA = water_data[0][:, 24]

                cum_leach_30mmh_SF = water_data[0][:, 26]
                cum_leach_30mmh_SA = water_data[0][:, 28]
                cum_leach_30mmh_LF = water_data[0][:, 30]
                cum_leach_30mmh_LA = water_data[0][:, 32]

                # Group each compartment for graphing
                percol_data2 = stackdata16(
                    cum_time_30min,
                    cum_leach_135mmh_SF, cum_leach_55mmhA_SF, cum_leach_55mmhB_SF, cum_leach_30mmh_SF,
                    cum_leach_135mmh_SA, cum_leach_55mmhA_SA, cum_leach_55mmhB_SA, cum_leach_30mmh_SA,
                    cum_leach_135mmh_LF, cum_leach_55mmhA_LF, cum_leach_55mmhB_LF, cum_leach_30mmh_LF,
                    cum_leach_135mmh_LA, cum_leach_55mmhA_LA, cum_leach_55mmhB_LA, cum_leach_30mmh_LA)

                time_size_135mmh = water_data[0][:, 33]
                time_size_55mmhA = water_data[0][:, 34]
                time_size_55mmhB = water_data[0][:, 35]
                time_size_30mmh = water_data[0][:, 36]

                time_sizes2 = [time_size_135mmh, time_size_135mmh,
                               time_size_55mmhA, time_size_55mmhA,
                               time_size_55mmhB, time_size_55mmhB,
                               time_size_30mmh, time_size_30mmh]

                return hydroplot3(
                    percol_data2,
                    series_name1, series_name2, series_name3, series_name4,
                    leach_high_6min, leach_med_12min, leach_med_30min, leach_low_30min,
                    title
                )

            else:
                roff_high_6min = np.array([8.991, 26.633, 15.720, 19.350])
                roff_med_12min = np.array([21.193, 17.731, 0.756, 8.025])
                roff_med_30min = np.array([54.633, 39.350, 0., 23.688])
                roff_low_30min = np.array([13.973, 3.717, 0., 22.827])

                # Time axis
                cum_time_30min = water_data[0][:, 0]

                # Ponding cumulative
                cum_roff_135mmh_SF = water_data[0][:, 1]
                cum_roff_135mmh_SA = water_data[0][:, 3]
                cum_roff_135mmh_LF = water_data[0][:, 5]
                cum_roff_135mmh_LA = water_data[0][:, 7]

                cum_roff_55mmhA_SF = water_data[0][:, 9]
                cum_roff_55mmhA_SA = water_data[0][:, 11]
                cum_roff_55mmhA_LF = water_data[0][:, 13]
                cum_roff_55mmhA_LA = water_data[0][:, 15]

                cum_roff_55mmhB_SF = water_data[0][:, 17]
                cum_roff_55mmhB_SA = water_data[0][:, 19]
                cum_roff_55mmhB_LF = water_data[0][:, 21]
                cum_roff_55mmhB_LA = water_data[0][:, 23]

                cum_roff_30mmh_SF = water_data[0][:, 25]
                cum_roff_30mmh_SA = water_data[0][:, 27]
                cum_roff_30mmh_LF = water_data[0][:, 29]
                cum_roff_30mmh_LA = water_data[0][:, 31]

                runoff_data2 = stackdata16(
                    cum_time_30min,
                    cum_roff_135mmh_SF, cum_roff_55mmhA_SF, cum_roff_55mmhB_SF, cum_roff_30mmh_SF,
                    cum_roff_135mmh_SA, cum_roff_55mmhA_SA, cum_roff_55mmhB_SA, cum_roff_30mmh_SA,
                    cum_roff_135mmh_LF, cum_roff_55mmhA_LF, cum_roff_55mmhB_LF, cum_roff_30mmh_LF,
                    cum_roff_135mmh_LA, cum_roff_55mmhA_LA, cum_roff_55mmhB_LA, cum_roff_30mmh_LA)

                time_size_135mmh = water_data[0][:, 33]
                time_size_55mmhA = water_data[0][:, 34]
                time_size_55mmhB = water_data[0][:, 35]
                time_size_30mmh = water_data[0][:, 36]

                time_sizes2 = [time_size_135mmh, time_size_135mmh,
                               time_size_55mmhA, time_size_55mmhA,
                               time_size_55mmhB, time_size_55mmhB,
                               time_size_30mmh, time_size_30mmh]

                return hydroplot3(
                    runoff_data2,
                    series_name1, series_name2, series_name3, series_name4,
                    roff_high_6min, roff_med_12min, roff_med_30min, roff_low_30min,
                    title)







