import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq

def model_one_step(kf, tArray, fixed_params):
    """
    One-step model of a strand displacement reaction.
    """
    t = tArray
    # Unpack the input parameters
    y0, scale = fixed_params['y0'], fixed_params['scale']

    y0 = [i*1e-9 for i in y0]
    scale = 1e-9*scale

    def translator_model(y, t, kf, scale):
        """
        Differential equations for the one-step strand displacement reaction model.
        y: list of current concentrations [Y1, P, W]
        t: time
        kf: forward rate constant
        scale: initial concentration of reactant R1
        Returns the derivatives [dY1dt, dPdt, dWdt]
        
        implementation of the model:
        R1 + Y1 --kf--> P + W
        """
        Y1, P, W = y

        R1 = float(scale - P)
        dY1dt = -float(kf*(R1)*(Y1))
        dPdt = -dY1dt
        dWdt = dPdt

        return (dY1dt, dPdt, dWdt)

    #Solve transalot model based on initial conditions and fit
    y1 = [i[-1] for i in odeint(translator_model, y0, t, args=(kf, scale), rtol = 1e-12, atol = 1e-12)]
    if (kf < 0):
      return [np.inf]*len(y1)

    return y1

def residuals(k, model_func, t, y_list, fixed_params_list):
    """
    Loop through each set of y values and fixed parameters and calculate the residuals between the data and the model.
    """
    total_residuals = []
    n=0
    for y, fixed_params in zip(y_list, fixed_params_list):
        model_y = model_func(k, t[n], fixed_params)
        total_residuals.append(y - model_y)
        n += 1
        
    # Concatenate all the residuals into one array
    return np.concatenate(total_residuals)

def fit_one_step(tArray, yList, initial_k, fixed_params_list):
    """
    Fit the one-step strand displacement model to the data.
    tArray: list of time arrays for each dataset
    yList: list of data arrays for each dataset
    initial_k: initial guess for the rate constant
    fixed_params_list: list of dictionaries of fixed parameters for each dataset
    Returns the result of the least squares fitting.
    """
    # Fit the model to the data
    result  = leastsq(residuals, initial_k, args=(model_one_step, tArray, yList, fixed_params_list))

    return result