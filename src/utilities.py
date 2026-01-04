import numpy as np
from scipy.optimize import fsolve
import datetime
import json
from nupack import *
from matplotlib import pyplot as plt
import re
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.optimize import leastsq
from scipy.integrate import odeint
from matplotlib.patches import Patch
from scipy.stats import linregress



#Function to convert time in HH:MM:SS format to seconds
def time_to_seconds(time):
    """
    Convvert single HH:MM:SS time string to seconds
    """
    pc = re.split(':',time)
    return int(pc[0])*3600 + int(pc[1])*60 + int(pc[2])

def timelist_to_seconds(time_series):
    """
    Convert a list of time strings in HH:MM:SS format to seconds
    """
    # Split the time string into components (hours, minutes, seconds)
    split_string = lambda timestring: re.split(':', timestring) 
    time = [float(split_string(value)[0])*3600 + float(split_string(value)[1])*60 + float(split_string(value)[2]) for value in time_series]
    return time

def lookup_condition(column_header, conds):
    """
    Given a column header string and a dictionary of conditions,
    return the condition key that contains the column header in its value list.
    Assumes that at most one condition contains the column header.
    """
    count = 0
    out = ""
    for e in conds.values():
        if column_header in e:
            count += 1
            out = list(conds.keys())[list(conds.values()).index(e)]
    assert(count <= 1)
    return out

def get_endpoint(level, x, data, conditions):
    """
    Calculate the average of the last 10 data points for each column in the data
    that matches the specified level in the conditions dictionary.
    Returns a DataFrame with the averages, row average, and row standard deviation.
    """
    averages = []
    for i, c in enumerate(data.columns):
        col = c
        if lookup_condition(col, conditions) == level:
            y = [float(i) for i in data[col].values]
            avg = np.average(y[-10:])
            averages.append(avg)
    df = pd.DataFrame(averages).T
    df.columns = ['Average '+str(i+1) for i in range(len(averages))]
    df.insert(0, "Level", level)  # Insert the level name as the first column
    df['Row Average'] = df.iloc[:, 1:].mean(axis=1)  # Compute the mean excluding the "Level" column
    df['Row Std'] = df.iloc[:, 1:].std(axis=1)  # Compute the std excluding the "Level" column
    return df

def get_initial(level, x, data, conditions):
    """
    Calculate the average of the first 10 data points for each column in the data
    that matches the specified level in the conditions dictionary.
    Returns a DataFrame with the averages, row average, and row standard deviation.
    """
    averages = []

    for i, c in enumerate(data.columns):
        col = c
        if lookup_condition(col, conditions) == level:
            y = [float(i) for i in data[col].values]
            avg = np.average(y[0:10])
            averages.append(avg)

    df = pd.DataFrame(averages).T
    df.columns = ['Average '+str(i+1) for i in range(len(averages))]
    df.insert(0, "Level", f"{level} init")  # Insert the level name as the first column
    df['Row Average'] = df.iloc[:, 1:].mean(axis=1)  # Compute the mean excluding the "Level" column
    df['Row Std'] = df.iloc[:, 1:].std(axis=1)  # Compute the std excluding the "Level" column

    return df

def plot_calibrated_averages(level, x, ax, data, conditions, fit, color_dict, format_dict=None):
    """
    Plot the calibrated averages and standard deviations for a given level, and best linear fit point estimate.
    """
    data_list = []
    for i, c in enumerate(data.columns):
        col = c 
        if i>= 2 and lookup_condition(col, conditions) == level:
            ynew = [(float(x) - fit[0])/fit[1] for x in data[col].values]
            data_list.append(ynew)
            
    data_array = np.array(data_list)
    y_mean = np.mean(data_array, axis=0)
    y_std = np.std(data_array, axis=0)
    lab = 'Av' + level
    color = color_dict[level] if level in color_dict else 'blue' # set default color to blue

    if format_dict is None:
        format_dict = {'point_size': 1, 'point_opacity': 0.6, 'fill_opacity': 0.5}

    point_size = format_dict.get('point_size', 1)
    point_opacity = format_dict.get('point_opacity', 0.6)
    fill_opacity = format_dict.get('fill_opacity', 0.5)

    ax.scatter(x, y_mean, label=lab, color=color, alpha=point_opacity, s=[point_size]*len(y_mean))
    ax.fill_between(x, y_mean - y_std, y_mean + y_std, color=color, alpha=fill_opacity)
    
def completion_at_time(data,
                      baseline_data,
                      triggered_data,
                      time_list_seconds,
                      time_in_seconds,
                      conditions,
                      baseline_conditions,
                      triggered_conditions,
                      experimental_level,
                      base_line_level,
                      triggered_level,
                      verbose = False):
  """
  Calculate the completion level at a specific time point based on baseline and triggered data.
  """

  def progegate_uncertainty_averaged_stddev(std_dev_list):
    return np.sqrt(sum(s**2 for s in std_dev_list)/(len(std_dev_list)**2))

  def average_last_n_with_std(n, average_series, std_series):
    av = average_series[-n:].mean()
    std = progegate_uncertainty_averaged_stddev(std_series[-n:].values)
    return av, std
  def average_first_n_with_std(n, average_series, std_series):
    av = average_series[:n].mean()
    std = progegate_uncertainty_averaged_stddev(std_series[:n].values)
    return av, std
  def calculate_completion_level(V, T, B):
    return (float(V) - float(B))/(float(T) - float(B))
  def series_av_and_std(dat, cond, level_name):
    target_columns = cond[level_name]
    dat[target_columns] = dat[target_columns].apply(pd.to_numeric)
    average_by_point = dat[target_columns].mean(axis=1)
    std_by_point = dat[target_columns].std(axis=1)
    return average_by_point, std_by_point

  #Get the x-values in seconds
  x = time_list_seconds
  #Fine the index of the first datapoint greater than time in seconds
  desired_time_index = min([i for i, time in enumerate(x) if time >= time_in_seconds])
  #Fist calculate baseline average and uncertainty at LAST AVAILABLE 5 DATAPOINTS
  base_line_average, base_line_std = average_first_n_with_std(5, *series_av_and_std(baseline_data, baseline_conditions, base_line_level))
  #Next calculate the triggering average and uncerainty at LAST 5 AVAILABLE DATAPOINTS
  triggered_average, triggered_std = average_last_n_with_std(5, *series_av_and_std(triggered_data, triggered_conditions, triggered_level))
  completion_levels = []
  for l_name in conditions[experimental_level]:
    completion_levels += [calculate_completion_level(data[l_name][desired_time_index], triggered_average, base_line_average)]

  if verbose:
    print(f"{experimental_level} : {np.average(completion_levels)} +/- {np.std(completion_levels)}")
  return np.average(completion_levels), np.std(completion_levels)



def average_and_std(curves):
    """
    Calculate the pointwise mean and standard deviation of a list of curves.
    """
    
    curves_array = np.array(curves)
    mean_curve = np.mean(curves_array, axis=0)
    std_curve = np.std(curves_array, axis=0)
    return mean_curve.tolist(), std_curve.tolist()

def get_average_curve(startTime, offset, level, x, data, conditions, fit):
    """
    Calculate the average curve and standard deviation from data columns based on a condition label, and a dictionary assigning conditions lavels to replicate column lavels.
    """
    y_curves = []

    for i, c in enumerate(data.columns):
        if lookup_condition(c, conditions) == level:
            ynew = [(float(data[c].values[z]) - fit[0]) / fit[1] for z in range(len(data[c].values)) if float(x[z]) >= startTime]
            y_curves.append(ynew)

    y_avg, y_std = average_and_std(y_curves)
    x_vals = [i + offset - startTime for i in x if i >= startTime]

    return x_vals, y_avg, y_std 

def get_triplicates_jackknife(startTime, offset, level, x, data, conditions, fit):
    """
    Leave-one-out version of get_average_curve.
    Returns 3 average curves and their standard deviations, each leaving out one replicate.
    """
    y_all = []
    for i, c in enumerate(data.columns):
        if lookup_condition(c, conditions) == level:
            ynew = [(float(data[c].values[z]) - fit[0]) / fit[1] for z in range(len(data[c].values)) if float(x[z]) >= startTime]
            y_all.append(ynew)

    if len(y_all) != 3:
        raise ValueError(f"Expected 3 replicates, found {len(y_all)} for level {level}")

    x_vals = [i + offset - startTime for i in x if i >= startTime]

    y_curves = []
    y_stds = []

    for i in range(3):
        subset = [y for j, y in enumerate(y_all) if j != i]
        avg, std = average_and_std(subset)
        y_curves.append(avg)
        y_stds.append(std)

    return x_vals, y_curves, y_stds