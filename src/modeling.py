import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq

def model_one_step(kf, tArray, fixed_params):
    """
    One-step model of a strand displacement reaction.
    Using nM scale internally for better fitting stability.
    But solver reports kf in 1/(M*s)
    """
    t = tArray
    # Unpack the input parameters
    y0, scale = fixed_params['y0'], fixed_params['scale']
    kf = float(kf)*1e-9  # in 1/(nM*s)

    def onestep_model(y, t, kf, scale):
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
        P, Y1, W = y

        R1 = float(scale - P)
        dY1dt = -float(kf*(R1)*(Y1))
        dPdt = -dY1dt
        dWdt = dPdt

        return (dPdt, dY1dt, dWdt)

    #Solve transalot model based on initial conditions and fit
    y1 = [i[0] for i in odeint(onestep_model, y0, t, args=(kf, scale), rtol = 1e-12, atol = 1e-12)]
    if (kf < 0):
      return [np.inf]*len(y1)

    return y1

def model_two_step(kf, tArray, fixed_params):
    """
    Cascade model of a strand displacement reaction.
    Implements:
    (F0t + F1 <kr-kf> F1t + W1)
    (F1t + F2 <kr2-kf2> Y1 + W2)
    Y1 + R1 -krep> P
    """

    t = tArray
    # Unpack the input parameters
    y0, krep, scale, dG = fixed_params['y0'], fixed_params['krep'], fixed_params['scale'], fixed_params['dG']

    y0 = [i*1e-9 for i in y0]
    scale = 1e-9*scale

    krep = krep
    
    kf = float(kf)
    kr = float(kf / np.exp(-dG/(0.592)))  # Calculate kr from kf and dG
    

    def two_step_model(y, t, krep, kf, kr, scale):
        P, Y1, F1t, F2, W2 = y
        # CRN:
        # (F1t + F2 <kr-kf> Y1 + W2)
        # Y1 + R1 -kr> P

        #Rate of reporting
        R1 = scale - P
        dPdt = krep*(R1)*(Y1)
        dY1dt = kf*F1t*F2 - kr*Y1*W2 - dPdt

        #Rate of change of input and output

        dF1tdt = -kf*F1t*F2 + kr*Y1*W2 
        dF2dt  = -kf*F1t*F2 + kr*Y1*W2

        dW2dt  =  kf*F1t*F2 - kr*W2*Y1

        return (dPdt, dY1dt, dF1tdt, dF2dt, dW2dt)

    #Solve transalot model based on initial conditions and fit

    y1 = [i[0] for i in odeint(two_step_model, y0, t, args=(krep, kf, kr, scale), rtol = 1e-12, atol = 1e-12)]

    if kf < 0:
      y1 = [np.inf for _ in range(len(y1))]

    return y1


def model_three_step(kf,  tArray, fixed_params):
    """
    Cascade model of a strand displacement reaction.
    Implements:
    (X1 + F0 <kr-kf> F0t + W0)
    (F0t + F1 <kr1-kf1> F1t + W1)
    (F1t + F2 <kr2-kf2> Y1 + W2)
    Y1 + R1 -krep> P
    """

    t = tArray
    # Unpack the input parameters
    y0, kf2, kr2, krep, scale, dG = fixed_params['y0'],fixed_params['kf2'], fixed_params['kr2'], fixed_params['krep'],  fixed_params['scale'], fixed_params['dG']

    y0 = [i*1e-9 for i in y0]
    scale = 1e-9*scale

    krep = krep

    kf = float(kf)
    kr = float(kf / np.exp(-dG/(0.592)))  # Calculate kr from kf and dG
    

    def three_step_model(y, t, krep, kf2 , kr2, kf, kr, scale):
        P, Y1, F1t, F0t, F2, F1, W1, W2 = y
        # CRN:
        # (F0t + F1 <kr1-kf1> F1t + W1)
        # (F1t + F2 <kr2-kf2> Y1 + W2)
        # Y1 + R1 -krep> P


        #Rate of reporting
        R1 = scale - P
        dPdt = krep*(R1)*(Y1)

        dY1dt = kf2*F1t*F2 - kr2*Y1*W2 - dPdt

        #Rate of change of gate concentrations
        dF1tdt =  - kf2*F1t*F2 + kr2*Y1*W2 + kf*F0t*F1 - kr*F1t*W1
        dF0tdt = - kf*F0t*F1 + kr*F1t*W1
        dF1dt = -kf*F0t*F1 + kr*F1t*W1  
        dF2dt = -kf2*F1t*F2 + kr2*Y1*W2 

        #Rate of change of waste concentrations
        dW1dt = + kf*F0t*F1 - kr*W1*F1t
        dW2dt = + kf2*F1t*F2 - kr2*W2*Y1

        return (dPdt, dY1dt, dF1tdt, dF0tdt, dF2dt, dF1dt, dW1dt, dW2dt)

    #Solve transalot model based on initial conditions and fit

    y1 = [i[0] for i in odeint(three_step_model, y0, t, args=(krep, kf2, kr2, kf, kr, scale), rtol = 1e-12, atol = 1e-12)]

    if kf < 0:
      y1 = [np.inf for _ in range(len(y1))]

    return y1

def model_cascade(kf,  tArray, fixed_params):
    """
    Cascade model of a strand displacement reaction.
    Implements:
    (X1 + F0 <kr0-kf0> F0t + W0)
    (F0t + F1 <kr1-kf1> F1t + W1)
    (F1t + F2 <kr2-kf2> Y1 + W2)
    Y1 + R1 -krep> P
    """

    t = tArray
    # Unpack the input parameters
    y0, kf1, kf2, kr1, kr2, krep, scale, dG = fixed_params['y0'],fixed_params['kf1'], fixed_params['kf2'], fixed_params['kr1'],fixed_params['kr2'], fixed_params['krep'],  fixed_params['scale'], fixed_params['dG']

    y0 = [i*1e-9 for i in y0]
    scale = 1e-9*scale

    krep = krep
    
    
    kf = float(kf)
    kr = float(kf / np.exp(-dG/(0.592)))  # Calculate kr from kf and dG

    def cascade_model(y, t, krep, kf1 ,kf2 ,kr1,kr2, kf, kr, scale):
        P, Y1,F1t, F0t, X1, F2, F1, F0, W0, W1, W2 = y
        # CRN:
        # X1 <krev-kf> Y1
        # (X1 + F0 <krev-kf> F0t + W0)
        # (F0t + F1 <krev-kf> F1t + W1)
        # (F1t + F2 <krev-kf> Y1 + W2)
        # Y1 + R1 -kr> P


        #Rate of reporting
        R1 = scale - P
        dPdt = krep*(R1)*(Y1)

        #Rate of change of input and output
        dX1dt = kr*F0t*W0 - kf*X1*F0
        dY1dt = kf2*F1t*F2 - kr2*Y1*W2 - dPdt

        #Rate of change of gate concentrations
        dF1tdt = kf1*F0t*F1 - kr1*F1t*W1 - kf2*F1t*F2 + kr2*Y1*W2
        dF0tdt = kf*X1*F0 - kr*F0t*W0 - kf1*F0t*F1 + kr1*F1t*W1
        dF0dt = -kf*X1*F0 + kr*F0t*W0
        dF1dt = -kf1*F0t*F1 + kr1*F1t*W1
        dF2dt = -kf2*F1t*F2 + kr2*Y1*W2

        #Rate of change of waste concentrations
        dW0dt =  kf*X1*F0 - kr*W0*F0t
        dW1dt = kf1*F0t*F1 - kr1*W1*F1t
        dW2dt = kf2*F1t*F2 - kr2*W2*Y1

        return (dPdt, dY1dt, dF1tdt, dF0tdt, dX1dt, dF2dt, dF1dt, dF0dt, dW0dt, dW1dt, dW2dt)

    #Solve transalot model based on initial conditions and fit

    y1 = [i[0] for i in odeint(cascade_model, y0, t, args=(krep, kf1, kf2, kr1,kr2,kf,kr, scale), rtol = 1e-12, atol = 1e-12)]

    if kf < 0:
      y1 = [np.inf for _ in range(len(y1))]

    return y1


def model_occlusion(kb_off, tArray, fixed_params):
    """
    Occlusion model of a strand displacement reaction.
    Implements:
    (Y1 + R1 -krep> P + W)
    (Y1 + B <-kb_off kb_on->  S)
    Using nM scale internally for better fitting stability.
    But solver reports kb_off in 1/s
    """
    t = tArray
    # Unpack the input parameters
    y0, kf, scale = fixed_params['y0'], fixed_params['kf'], fixed_params['scale']
    y0 = y0[0:1] + [scale] + y0[1:]
    kb_off = float(kb_off) #Unimolecular off-rate in 1/s
    kf = float(kf)*1e-9  # in 1/(nM*s)
    def occlusion_model(y, t, kb_off, kf):
        kb_off = float(kb_off)

        Y1, R1, P, W, B, S = y
        # CRN:
        # Y1 + R1 -krep> P + W
        # Y1 + B <-kb_off kb_on->  S
        #Rate of reporting
        kb_on = 1e7*1e-9 # in 1/(nM*s)
        dR1dt = -float(kf*(R1)*(Y1)) -float(kb_on*(B)*(R1)) + float(kb_off*S)
        dY1dt = -float(kf*(R1)*(Y1)) 
        dPdt = float(kf*(R1)*(Y1))
        dWdt = dPdt
        dBdt = -float(kb_on*(B)*(R1)) + float(kb_off*S)
        dSdt = float(kb_on*(B)*(R1)) - float(kb_off*S)

        return (dY1dt, dR1dt, dPdt, dWdt, dBdt, dSdt)

    #Solve transalot model based on initial conditions and fit
    y1 = [i[-3] for i in odeint(occlusion_model, y0, t, args=(kb_off, kf), rtol = 1e-12, atol = 1e-12)]
    if (kb_off < 0):
      return [np.inf]*len(y1)

    return y1



def residuals(k, model_func, t, y_list, fixed_params_list):
    """
    Loop through each set of y values and fixed parameters and calculate the residuals between the data and the model.
    """
    total_residuals = []
    n=0
    for y, fixed_params in zip(y_list, fixed_params_list):
        #this logic is required so we can give t0 at time 0
        tfull = np.arange(0, t[n][-1]+(t[n][-1]-t[n][-2]),t[n][-1]-t[n][-2])
        mindex = min([i for i in range(len(tfull)) if tfull[i]>= t[n][0]])
        model_y = model_func(k, tfull, fixed_params)[mindex:]
        total_residuals.append(np.array(y) - np.array(model_y))
        n += 1
        
    # Concatenate all the residuals into one array
    return np.concatenate(total_residuals)

def fit_model(model, tArray, yList, initial_k, fixed_params_list):
    """
    Fit the one-step strand displacement model to the data.
    tArray: list of time arrays for each dataset
    yList: list of data arrays for each dataset
    initial_k: initial guess for the rate constant
    fixed_params_list: list of dictionaries of fixed parameters for each dataset
    Returns the result of the least squares fitting.
    """
    # Fit the model to the data
    result  = leastsq(residuals, initial_k, args=(model, tArray, yList, fixed_params_list))

    return result

