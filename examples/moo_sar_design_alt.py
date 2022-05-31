""" This script performs SAR instrument design based on a multi-objective optimization algorithm (NSGA-II), enabled by the 
    `pymoo 0.4.2` package. A baseline instrument is defined, and some of the instrument parameters are kept as optimization variables.
    Multiple objectives are defined in terms of the data-metrics produced by the InstruPy package.
    The result is a multi-dimensional Pareto-curve illustrating the trade-offs between different instrument designs. Any point on the
    Pareto curve is optimal and a suitable point can be selected by the user for implementation.

    The Pareto Curve (instrument variables and corresponding data-metrics) in a csv file in the folder ``moo_results``. Various plots
    such as the design-space, objective-space, etc are saved (as svg files) and displayed. If you want faster runtime (with potential 
    reduction in optimality), reduce the number of generations and population size in the NSGA-II configuration.

    Baseline instrument: Synthetic Aperture Radar with following fixed parameters:
                            * Operational altitude is 350km
                            * Dual-polarization (SMAP pulse configuration with default pulse-separation)
                            * L-band
                            * Operating points to be considered at incidence angles of 30 deg, 45 deg and 60 deg
                            * 1kW transmit power
                            * Full swath configuration (processed swath = illuminated swath)
                            * Losses, Noise: Antenna efficiency of 60%, 2dB radar loss, 2dB system noise figure

    Variables: (1) Daz [m] (antenna height in meters) range: 0.5m to 15m
               (2) Delv [m] (antenna width in meters) range: 0.5m to 5m
               (3) Chirp BW [MHz] (chirp bandwidth in Mega Hertz) range: 0.1MHz to 40MHz
               (4) Pulse width [us] (transmit pulse width in micro-seconds) range: 1us to 1000us
    
    (Minimization) Objectives: 
                (1) Antenna area: Antenna height * Antenna width
                
                (2) -1* Swath Width
                
                (3) NESZ [dB]
                
                (4) Speckle reduction [dB] for 1km2 pixel:
                    Number of Looks (N) = processed pixel area (1km2) / unprocessed pixel area (along-track resolution * cross-track resolution)
                    The speckle reduction in decibels is given by 10 log_10(1/N) (See 5.8, Pg 187, David Gardner Long , Fawwaz T. Ulaby - Microwave Radar and Radiometric Remote Sensing)
    
    Constraints: Valid operating point (not all designs are valid for SAR since a suitable operational PRF may not be available.)
                 The operations at incidence angles of 35 deg, 45 deg and 55 deg are considered.
"""
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
import pickle

from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_problem, get_sampling, get_crossover, get_mutation, get_termination
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.util.misc import stack
from pymoo.model.problem import Problem
from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel

c = 299792458 # m/s
inc = [35, 45, 55] # [deg]

out_dir = os.path.dirname(os.path.realpath(__file__)) + '/moo_results/'
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.makedirs(out_dir)

class MyProblem(Problem):
    """     # Variables: (1) Daz [m], (2) Delv [m], (3), Chirp BW [MHz] (4) Pulse width [us]
            # Objectives: (1) Antenna area (2) -1* Swath Width (3) NESZ [dB] (4) # Speckle reduction [dB] for 1km2 pixel
            # Constraints: Valid operating point
    """
    def __init__(self):
        super().__init__(n_var=5,
                         n_obj=3,
                         n_constr=1,
                         xl=np.array([0.1, 0.1, 1, 1, 300]),
                         xu=np.array([15, 15, 40, 1000, 1000]),
                         elementwise_evaluation=True)

    @staticmethod
    def define_test_sar(daz_m, delv_m, chirpbw_Mhz, pulse_w_us):

        test_sar = SyntheticApertureRadarModel.from_json(   '{"@type": "Synthetic Aperture Radar",'
                                                            '"orientation": {'
                                                            '    "referenceFrame": "SC_BODY_FIXED",'
                                                            '    "convention": "SIDE_LOOK",'
                                                            '    "sideLookAngle": 45' # not of significance in this problem
                                                            '},'
                                                            '"pulseWidth":'+ str(pulse_w_us*1e-6) +','
                                                            '"antenna":{"shape": "RECTANGULAR", "height":'+ str(daz_m) +',"width":'+ str(delv_m) +',' 
                                                            '            "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},'
                                                            '"operatingFrequency": 1.2757e9,' 
                                                            '"peakTransmitPower": 1000,' 
                                                            '"chirpBandwidth":'+ str(chirpbw_Mhz*1e6) +','     
                                                            '"minimumPRF": 1,' 
                                                            '"maximumPRF": 20000,' 
                                                            '"radarLoss": 2,' 
                                                            '"systemNoiseFigure": 2,'
                                                            '"swathConfig": {'
                                                            '   "@type": "fixed",'
                                                            '   "fixedSwathSize": 25'
                                                            '},'
                                                            '"polarization": {'
                                                            '   "@type": "dual",'
                                                            '   "pulseConfig": {'
                                                            '        "@type": "SMAP"''}' 
                                                            '}'                                                   
                                                            '}')
        return test_sar            


    def _evaluate(self, x, out, *args, **kwargs):
        
        daz_m = x[0]
        delv_m = x[1]
        chirpbw_Mhz = x[2]
        pulse_w_us = x[3]
        alt_km = x[4]
        h = alt_km*1e3
        Re = 6378.137e3 # [m]
        orb_speed = np.sqrt(3.986004418e14/(Re + h)) # [m/s]

        test_sar = MyProblem.define_test_sar(daz_m, delv_m, chirpbw_Mhz, pulse_w_us)
      
        # Calculate the metrics at the 3 incidence angles. A valid operation point must be found for all of them.
        obsv_metrics_1 = test_sar.calc_data_metrics(alt_km = h*1e-3, sc_speed_kmps = orb_speed*1e-3, sc_gnd_speed_kmps = orb_speed*1e-3*(Re/(Re+h)), 
                                                        inc_angle_deg = inc[0], 
                                                        instru_look_angle_from_target_inc_angle = True)
        
        obsv_metrics_2 = test_sar.calc_data_metrics(alt_km = h*1e-3, sc_speed_kmps = orb_speed*1e-3, sc_gnd_speed_kmps = orb_speed*1e-3*(Re/(Re+h)), 
                                                        inc_angle_deg = inc[1], 
                                                        instru_look_angle_from_target_inc_angle = True)

        obsv_metrics_3 = test_sar.calc_data_metrics(alt_km = h*1e-3, sc_speed_kmps = orb_speed*1e-3, sc_gnd_speed_kmps = orb_speed*1e-3*(Re/(Re+h)), 
                                                        inc_angle_deg = inc[2], 
                                                        instru_look_angle_from_target_inc_angle = True)

        # initialize to an invalid observation point
        g1 = 1
        f1 = 1e9
        f2 = 1e9
        f3 = 1e9
            
        if(obsv_metrics_1["NESZ [dB]"] and obsv_metrics_2["NESZ [dB]"] and obsv_metrics_3["NESZ [dB]"] ):

            if not np.isnan(obsv_metrics_1["NESZ [dB]"]) and not np.isnan(obsv_metrics_2["NESZ [dB]"]) and not np.isnan(obsv_metrics_3["NESZ [dB]"]):
            
                # Include below condition in case of fixed-swath configuration. If this condition is not true it means that the illuminated swath size < the fixed swath size
                if(obsv_metrics_1["swath-width [km]"] == 25 and obsv_metrics_2["swath-width [km]"] == 25 and obsv_metrics_3["swath-width [km]"] == 25):
                    
                    # valid observation point
                    g1 = -1

                    f1 = obsv_metrics_3["NESZ [dB]"]
                    f2 = obsv_metrics_3["ground pixel along-track resolution [m]"]
                    N_looks = 1e6 / (obsv_metrics_3["ground pixel along-track resolution [m]"] * obsv_metrics_3["ground pixel cross-track resolution [m]"])
                    f3 = -N_looks
                    #print(f1,f2,f3,f4)
            



        out["F"] = np.column_stack([f1, f2, f3])
        out["G"] = np.column_stack([g1])

problem = MyProblem()

algorithm = NSGA2(
    pop_size=48,
    sampling=get_sampling("real_random"),
    crossover=get_crossover("real_sbx", prob=0.9, eta=20),
    mutation=get_mutation("real_pm", eta=20),
    eliminate_duplicates=True
)

termination = get_termination("n_eval", 4800)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               pf=problem.pareto_front(use_cache=False),
               save_history=True,
               verbose=False)

#pickle.dump(res, open(out_dir+ "data.p", "wb" ) ) # uncomment to save as pickle file

# Save Pareto curve as csv file
import pandas as pd 
df = pd.DataFrame({"Daz [m]" : res.X[:,0], "Delv [m]" : res.X[:,1], "Chirp BW [MHz]" : res.X[:,2], "Pulse width [us]" : res.X[:,3], "Altitude [m]" : res.X[:,4], "NESZ [dB]" : res.F[:,0], "along-track-res" : res.F[:,1], "N_looks" : res.F[:,2]})
df.to_csv(out_dir+"pareto_front_030122.csv")

# Plot Design Space
# x1: Daz [m], x2: Delv [m], x3: Chirp BW [MHz], x4: Pulse width [us]
plot = Scatter(title = "Design Space", axis_labels="x", tight_layout=True)
plot.add(res.X, s=30, facecolors='none', edgecolors='c')
plot.do()
plt.savefig(out_dir+"design_space.svg")
#plot.apply(lambda ax: ax.set_xlim(0.5, 0.5, 0.1, 1))
#plot.apply(lambda ax: ax.set_ylim(20, 20, 80, 1000))
plot.show()

# Plot Objective Space
# f1: Antenna area [m2], f2: -Swath width [km], f3: Sigma NESZ [dB], f4: Speckle reduction [dB]
plot = Scatter(title = "Objective Space", tight_layout=True)
plot.add(res.F)
plot.do()
plt.savefig(out_dir+"obj_space.svg")
plot.show()

# In[ ]:
'''
plt.plot(n_evals, hv, '-o')
plt.title("Convergence")
plt.xlabel("Function Evaluations")
plt.ylabel("Hypervolume")
plt.savefig("hypervolume.svg")
plt.show()
'''


# In[ ]:
from pymoo.visualization.pcp import PCP
plot = PCP().add(res.F)
plot.add(res.F[0], linewidth=5, color="green")
plot.add(res.F[8], linewidth=5, color="cyan")
plot.add(res.F[163], linewidth=5, color="red")
plot.add(res.F[190], linewidth=5, color="black")
plot.add(res.F[519], linewidth=5, color="yellow")
plot.show()
plt.savefig(out_dir+"pcp.svg")


# In[ ]:
from pymoo.factory import get_performance_indicator
A = res.F[10]
hv = get_performance_indicator("hv", ref_point=np.array([15*15, 0, 0, 0]))
print("hv", hv.calc(A))


'''
# Debug statements.

# In[ ]:


res.F[10]


# In[ ]:


res.X


# In[ ]:


res.F[1]


# In[ ]:


from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

np.set_printoptions(threshold=np.inf)
res.F[res.F[:,2] < -30]


# In[ ]:


np.where(res.F[:,2] < -30)


# In[ ]:


res.F[84]


# In[ ]:


res.X[84]


# In[ ]:


np.size(res.F)


# In[ ]:

'''


