# config.py serves to store model parameters and reference statistics

# 1. Reference statistics from Habib's data analysis
# Alberta Historical Wildfire Data 2006-2025
# Calculated from ALberta wildfire assessment records.
# 
# These values are used to standardlize wind and humidity.
# Z_W = (W - WIND_MEAN) / WIND_STD
# Z_H = (H - HUMIDITY_MEAN) / HUMIDITY_STD
# 
# Positive anomaly: Above reference mean
# Negative anomaly: Below reference mean
#
# Larger STD means observations are more spread out.
# These are measured statistics, not tuning parameters.


WIND_MEAN = 8.920230            # km/h
WIND_STD = 8.415991             # km/h

HUMIDITY_MEAN = 45.203347       # percent %
HUMIDITY_STD = 18.674584        # percentage points 


# 2. Ornstein-Uhlenbeck (OU) Parameters
# dZ_t = -kappa * Z_t dt + sigma * dB_t
#
# These values are illustrative, not estimated from data
#
# KAPPA:
#   Controls the speed of mean reversion
#   Higher kappa  -> faster return toward 0
#   Lower kappa   -> slower return toward zero
#
# SIGMA: 
#   Controls the strength of random fluctuations
#   Higher sigma  -> greater variability
#   Lower sigma   -> smaller random fluctuations

KAPPA_W = 0.8
SIGMA_W = 1.0

KAPPA_H = 0.5
SIGMA_H = 0.7

# 3. Simulation settings
# T:
#   Total simulation duration
#
# DT:
#   Time step between observations of the simulation
#   Smaller DT  -> more time points and more computation
# 
# SEED:
#   Controls reproducibility of random number generation
#   Changing SEED changes the random realization
#
# Time units are illustrative until a consistent physical time scale has been selected

T = 10.0
DT = 0.01
SEED = 42

# 4. Wildfire Intensity Parameters
# All values below are illustrative assumptions
# They are not calibrated from observed fire intensity
# 
# Fuel category: C2 (Boreal Spruce)
#
# ALPHA_C:
#   Baseline intensity growth rate
#
# BETA_W:
#   Effect of standardlized wind anomaly
#
# BETA_H:
#   Effect of standardlized humidity anomaly
#
# GAMMA:
#   Wind-humidity interaction coefficient
#   GAMMA = 0 disables the interaction
#
# SIGMA_I:
#   Strength of random intensity fluctuations
#
# I0:
#   Initial wildfire intensity in kW/m
#   This is assumed starting value

FUEL_TYPE = "C2"
ALPHA_C = 0.05
BETA_W = 0.12
BETA_H = 0.10
GAMMA = -0.04
SIGMA_I = 0.15
I0 = 500.0


# Monte Carlo simulation settings
N_PATHS = 10000


