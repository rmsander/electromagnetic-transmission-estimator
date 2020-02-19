#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 14:56:19 2018

@author: rachelmorgan
"""

import numpy as np
import matplotlib.pyplot as plt


def radiation_pattern(phaseDifference, distance, power_tx, gain_tx,
                      frequencyOfOperationInHertz=2. * 10 ** 9,
                      numberOfAntennas=2):
    '''
    Input: phase difference between monopoles.
    Output: r and theta variables calculated for
    '''

    # Speed of light and distance
    c = 3. * 10 ** (8)  # speed of light in meters per second
    D = .075  # in meters

    # Electromagnetic wavelength
    wavelength = c / frequencyOfOperationInHertz

    # Plotting parameter
    numberOfPointsToConsiderForPlotting = 360

    gains = []
    predicted_insertion_losses = []
    thetas = np.arange(0, 2 * np.pi,
                       2 * np.pi / numberOfPointsToConsiderForPlotting)

    # create list of gains corresponding to each theta value:
    for theta in thetas:
        sumOfComplexExponentials = 0.
        # sum phases from each antenna to find gain at each angle:
        for i in range(numberOfAntennas):
            total_phase = -1j * i * ((2 * np.pi / wavelength) * D * np.sin(
                theta) - np.pi * phaseDifference / 180.)
            sumOfComplexExponentials += np.exp(total_phase)
        gain_rx = np.abs(sumOfComplexExponentials) ** 2
        gains.append(gain_rx)
        predicted_prx = predict_power_rx(power_tx, gain_tx, gain_rx, distance)
        predicted_IL = -10 * np.log10(predicted_prx)
        if predicted_IL > 60.:
            predicted_IL = 60.
        predicted_insertion_losses.append(predicted_IL)

    return gains, predicted_insertion_losses, thetas


def plot_radiation_patterns(phis):
    '''
    Computes and plots the radiation patterns for three phase offsets in one
    polar plot.
    '''

    plt.figure()
    ax = plt.subplot(111, projection='polar')
    labels = []
    for phi in phis:
        r, IL, theta = radiation_pattern(phi, 0.2, 1., 1.)  # note: extra
        # arguments not needed for plot

        labels.append('Phase offset: {}'.format(str(phi)))
        ax.plot(theta, r)

    # ax.set_rmax(max(sumOfComplexExponentials))
    # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
    ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
    ax.grid(True)
    ax.set_title("Theoretical Radiation Patterns", va='bottom')
    plt.legend(labels)
    plt.show();


def plot_results(phis, angle_predictions):
    '''
    Computes and plots the radiation patterns for three phase offsets in one
    polar plot.
    '''

    plt.figure(figsize=(5, 5))
    ax = plt.subplot(111, projection='polar')
    labels = []
    for phi in phis:
        r, IL, theta = radiation_pattern(phi, 0.2, 1., 1.)  # note: extra
        # arguments not needed for plot

        labels.append('Phase: {}'.format(str(phi)))

        # only plot on right half:
        theta_plot = np.append(theta[0:91], theta[270:-1])

        r_plot = np.append(r[0:91], r[270:359])

        ax.plot(theta_plot, r_plot)

    # add lines for predicted angle:
    radial_line = np.linspace(0, 4, 10)

    for angle in angle_predictions:
        angle_line = [angle * np.pi / 180. for x in radial_line]
        ax.plot(angle_line, radial_line, 'r')

    # ax.set_rmax(max(sumOfComplexExponentials))
    # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
    ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
    ax.grid(True)
    ax.set_title("Predicted Transmitter Angle: {}".format(
        np.round(angle_predictions[0], 2)), va='bottom')
    labels.append('Predicted Tx Location')
    plt.legend(labels)
    plt.show();


def friis_solve(power_rx, power_tx, gain_tx, r, frequency=2. * 10 ** 9):
    # solve for receiver gain from received power using Friis Transmisison
    # Equation
    wavelength = (3. * 10 ** 8) / frequency

    # solve Friis:
    gain_rx = power_rx / ((power_tx / (
                4 * np.pi * r ** 2)) * gain_tx * wavelength ** 2 / (4 * np.pi))
    return gain_rx


def predict_angles(gain_rx, gains, thetas, error_cutoff):
    # returns sorted list of candidate angles based on input power recieved
    # and gain pattern

    # find angles with gains close to theoretical and return list of (index,
    # angle, error) tuple for each:
    candidate_angles = np.array(
        [(i, x, np.abs(gains[i] - gain_rx)) for i, x in enumerate(thetas) if
         np.abs(gains[i] - gain_rx) < error_cutoff],
        dtype=[('ind', int), ('theta', float), ('err', float)])

    # sort angles by total errors:
    sorted_angles = np.sort(candidate_angles, order=['err'])
    return sorted_angles


def find_best_angle(sorted_anglesA, sorted_anglesB, sorted_anglesC):
    # returns angle that shows up earliest in both lists

    # unpack input tuples so it's easier to compare and add:
    angles_in_a = [angle[1] for angle in sorted_anglesA]
    angles_in_b = [angle[1] for angle in sorted_anglesB]
    angles_in_c = [angle[1] for angle in sorted_anglesC]
    errors_a = [angle[2] for angle in sorted_anglesA]
    errors_b = [angle[2] for angle in sorted_anglesB]
    errors_c = [angle[2] for angle in sorted_anglesC]

    # compare lists of angles and sum up errors for angles appearing in all
    # three lists:
    common_angles = np.array([(angle, errors_a[i] +
                               errors_b[angles_in_b.index(angle)] +
                               errors_c[angles_in_c.index(angle)])
                              for i, angle in enumerate(angles_in_a) if
                              angle in angles_in_b and angle in angles_in_c],
                              dtype=[('angle', float), ('error', float)])

    # sort angles by total error:
    sorted_candidates = np.sort(common_angles, order=['error'])

    #    print(sorted_candidates)

    return sorted_candidates


def find_best_angle_four(sorted_anglesA, sorted_anglesB, sorted_anglesC,
                         sorted_anglesD):
    # returns angle that shows up earliest in both lists

    # unpack input tuples so it's easier to compare and add:
    angles_in_a = [angle[1] for angle in sorted_anglesA]
    angles_in_b = [angle[1] for angle in sorted_anglesB]
    angles_in_c = [angle[1] for angle in sorted_anglesC]
    angles_in_d = [angle[1] for angle in sorted_anglesD]
    errors_a = [angle[2] for angle in sorted_anglesA]
    errors_b = [angle[2] for angle in sorted_anglesB]
    errors_c = [angle[2] for angle in sorted_anglesC]
    errors_d = [angle[2] for angle in sorted_anglesD]

    # compare lists of angles and sum up errors for angles appearing in all
    # three lists:
    common_angles = np.array([(angle, errors_a[i] + errors_b[
        angles_in_b.index(angle)] + errors_c[angles_in_c.index(angle)] +
                               errors_d[angles_in_d.index(angle)]) for i, angle
                              in enumerate(angles_in_a) if
                              angle in angles_in_b and angle in angles_in_c
                              and angle in angles_in_d],
                             dtype=[('angle', float), ('error', float)])

    # sort angles by total error:
    sorted_candidates = np.sort(common_angles, order=['error'])
    return sorted_candidates


def predict_power_rx(power_tx, gain_tx, gain_rx, r, frequency=2. * 10 ** 9):
    # calculate expected received power from Friis Transmission Equation:

    wavelength = (3. * 10 ** 8) / frequency

    power_rx = gain_rx * ((power_tx / (
                4 * np.pi * r ** 2)) * gain_tx * wavelength ** 2 / (4 * np.pi))

    return power_rx


# simulate system to test:
# set angle of transmitter for predicting IL values (int between 0 and 90,
# 270 and 360), corresponds to index in list of angles:

# plot radiation pattern:


# note: this section needs to be moved to below p_rx calculation for simulation

# calc expected received power for each phase setting:


'''
Input Parameters
'''

# input phase settings:
phi1 = 0
phi2 = 90
phi3 = 180
# phi4 = 225

# set angles for plot of radiation pattern for phase offsets:
phis = [phi1, phi2, phi3]
# phis = [ 0, 100]
plot_radiation_patterns(phis)

# set parameters for transmission
power_tx = 1.  # note: this is assumed 1 so that p_rx variable is now ratio
# of p_rx/p_tx
distance = 0.2
error_cutoff = 300.  # this is the cutoff to find angles with similar gains
# to theoretical

'''
Calibrate System
'''

# calculate theoretical insertion loss for just one monopole (for calibration):
monopole_gain = 3. / 2

# input calibration measurement value:
cal_IL_meas = 19.5

p_rx_cal = 10 ** (-cal_IL_meas / 10)

# solve for g_tx from one monopole calibration measurement assuming g_rx is 1.5
g_tx = friis_solve(p_rx_cal, power_tx, monopole_gain, distance)
print('calculated TX gain from calibration: ', np.round(g_tx, 2))

'''
calculate theoretical radiation pattern and insertion losses
'''

gains1, predicted_IL1raw, thetas1raw = radiation_pattern(phi1, distance,
                                                         power_tx, g_tx)
gains2, predicted_IL2raw, thetas2raw = radiation_pattern(phi2, distance,
                                                         power_tx, g_tx)
gains3, predicted_IL3raw, thetas3raw = radiation_pattern(phi3, distance,
                                                         power_tx, g_tx)
# gains4, predicted_IL4raw, thetas4raw = radiation_pattern(phi4, distance,
# power_tx, g_tx)

# get rid of angles that are in left half plane:
thetas1 = np.append(thetas1raw[0:91], thetas1raw[270:-1])
thetas2 = np.append(thetas2raw[0:91], thetas2raw[270:-1])
thetas3 = np.append(thetas3raw[0:91], thetas3raw[270:-1])
# thetas4 = np.append(thetas4raw[0:91], thetas4raw[270:-1])

predicted_IL1 = np.append(predicted_IL1raw[0:91], predicted_IL1raw[270:-1])
predicted_IL2 = np.append(predicted_IL2raw[0:91], predicted_IL2raw[270:-1])
predicted_IL3 = np.append(predicted_IL3raw[0:91], predicted_IL3raw[270:-1])
# predicted_IL4 = np.append(predicted_IL4raw[0:91], predicted_IL4raw[270:-1])


'''
predict insertion losses for given angle
'''

# tx_angle_index = 300
# power_rx_1 = predict_power_rx(power_tx, g_tx, gains1[tx_angle_index],
# distance)
# power_rx_2 = predict_power_rx(power_tx, g_tx, gains2[tx_angle_index],
# distance)
# power_rx_3 = predict_power_rx(power_tx, g_tx, gains3[tx_angle_index],
# distance)
##power_rx_4 = predict_power_rx(power_tx, g_tx, gains4[tx_angle_index],
# distance)
##convert to IL to compare directly to collected data
# theoretical_IL1 = -10*np.log10(power_rx_1)
# theoretical_IL2 = -10*np.log10(power_rx_2)
# theoretical_IL3 = -10*np.log10(power_rx_3)
# theoretical_IL4 = -10*np.log10(power_rx_4)

# print('theoretical IL 1:', theoretical_IL1)
# print('theoretical IL 2:', theoretical_IL2)
# print('theoretical IL 3:', theoretical_IL3)
# print('theoretical IL 4:', theoretical_IL4)


'''
predict location from measured IL values:
'''

# 0 degree phase offset measurement:
IL1 = 22.3
# 90 degree phase offset measurement:
IL2 = 20.4
# 180 degree phase offset measurement:
IL3 = 22.4
# IL4 = 24


# predict angles by comparing measured IL to theoretical calculations
predicted_angles1 = predict_angles(IL1, predicted_IL1, thetas1, error_cutoff)
predicted_angles2 = predict_angles(IL2, predicted_IL2, thetas2, error_cutoff)
predicted_angles3 = predict_angles(IL3, predicted_IL3, thetas3, error_cutoff)
# predicted_angles4 = predict_angles(IL4, predicted_IL4, thetas4, error_cutoff)


# find best angle:
list_candidates = find_best_angle(predicted_angles1, predicted_angles2,
                                  predicted_angles3)  # , predicted_angles4)

predicted_angles = [np.degrees(angle[0]) for angle in list_candidates]

# make list of top ten predictions:
output_prediction = predicted_angles[0:10]

print('predicted angles: ', output_prediction)

'''
Plot Results
'''
plot_results(phis, output_prediction)
# print(list_candidates)
