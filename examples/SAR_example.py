""" The following examples illustrate the SAR models with different possible set of configurations.
Please refer to
* html docs -> Instrument Models -> Synthetic Aperture Radar Model 
* html docs -> API Reference -> Synthetic Aperture Radar Model

"""
import json
import numpy
import sys, os


from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel

h = 500e3 # [m]
Re = 6378.137e3 # [m]
c = 299792458 # m/s

orb_speed = numpy.sqrt(3.986004418e14/(Re + h)) # [m/s] orbital speed
print('orb_speed is {} m/s'.format(orb_speed))

gnd_spd = orb_speed*(Re/(Re+h))
print('Ground speed is {} m/s'.format(gnd_spd))

# L- band SAR, Dual polarization (SMAP pulse configuration), Fixed swath size (25km)
sar1_dict = {  "@type": "Synthetic Aperture Radar",
                "name": "L-Band SAR #1",
                "orientation": {
                    "referenceFrame": "SC_BODY_FIXED",
                    "convention": "SIDE_LOOK",
                    "sideLookAngle": 45
                },
                "pulseWidth": 14.16e-6,
                "antenna":{"shape": "RECTANGULAR", "height": 14.38, "width": 1.48, 
                    "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                "operatingFrequency": 1275.7e6, 
                "peakTransmitPower": 1000, 
                "chirpBandwidth": 0.86e6,      
                "minimumPRF": 1, 
                "maximumPRF": 20000, 
                "radarLoss": 2, 
                "systemNoiseFigure": 2,
                "swathConfig": {
                    "@type": "fixed",
                    "fixedSwathSize": 25
                },
                "polarization": {
                    "@type": "dual",
                    "pulseConfig": {
                        "@type": "SMAP"   
                } 
                }                                                   
            }
sar1 = SyntheticApertureRadarModel.from_dict(sar1_dict)

""" Some points.

    *   By setting the `instru_look_angle_from_target_inc_angle` to `True`, the specified instrument look-angle is ignored. 
        Instead the instrument look angle is calculated internally from the input target incidence angle.

    *   Valid operating point has to be found over the range of angles which the instrument will slew. Below we test to see
        if operating points are available incidence of 30 deg, 45 deg and 60 deg. Sometimes the antenna sizes and other instrument parameters 
        may make imaging at certain angles not possible.

    *   Number of looks = processed pixel area/ unprocessed pixel area.
        The speckle reduction in decibels is given by 10 log_10(1/N) (See 5.8, Pg 187, David Gardner Long , Fawwaz T. Ulaby - Microwave Radar and Radiometric Remote Sensing)

"""
print("L- band SAR, Dual polarization (SMAP pulse configuration), Fixed swath size (25km)")
##
inc_deg = 35
obsv_metrics = sar1.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"])  
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 45
obsv_metrics = sar1.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 55
obsv_metrics = sar1.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

####################################################################################################################################################
# L-band SAR, Dual polarization (AIRSAR pulse configuration), Full swath
sar2_dict = {  "@type": "Synthetic Aperture Radar",
                "name": "L-Band SAR #2",
                "orientation": {
                    "referenceFrame": "SC_BODY_FIXED",
                    "convention": "SIDE_LOOK",
                    "sideLookAngle": 45
                },
                "pulseWidth": 1.39e-6,
                "antenna":{"shape": "RECTANGULAR", "height": 7.90, "width": 10.83, 
                    "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                "operatingFrequency": 1275.7e6, 
                "peakTransmitPower": 1000, 
                "chirpBandwidth": 0.70e6,      
                "minimumPRF": 1, 
                "maximumPRF": 20000, 
                "radarLoss": 2, 
                "systemNoiseFigure": 2,
                "swathConfig": {
                    "@type": "full"
                },
                "polarization": {
                    "@type": "dual",
                    "pulseConfig": {
                        "@type": "AIRSAR"   
                } 
                }                                                   
            }
sar2 = SyntheticApertureRadarModel.from_dict(sar2_dict)

print("L- band SAR, Dual polarization (AIRSAR pulse configuration), Full swath")
##
inc_deg = 35
obsv_metrics = sar2.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"])  
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 45
obsv_metrics = sar2.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 55
obsv_metrics = sar2.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

####################################################################################################################################################
# P-band SAR, Single  polarization, Fixed swath size (50km)
sar3_dict = {   "@type": "Synthetic Aperture Radar",
                "name": "P-Band SAR",
                "orientation": {
                    "referenceFrame": "SC_BODY_FIXED",
                    "convention": "SIDE_LOOK",
                    "sideLookAngle": 45
                },
                "pulseWidth": 22.14e-6,
                "antenna":{"shape": "RECTANGULAR", "height": 9.8, "width": 9.04, 
                    "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                "operatingFrequency": 435e6, 
                "peakTransmitPower": 1000, 
                "chirpBandwidth": 6e6,      
                "minimumPRF": 1, 
                "maximumPRF": 20000, 
                "radarLoss": 2, 
                "systemNoiseFigure": 2,
                "swathConfig": {
                    "@type": "fixed",
                    "fixedSwathSize": 50
                },
                "polarization": {
                    "@type": "single"
                }                                                   
            }
sar3 = SyntheticApertureRadarModel.from_dict(sar3_dict)

print("P-band SAR, Single  polarization, Fixed swath size (50km)")
##
inc_deg = 35
obsv_metrics = sar3.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"])  
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 45
obsv_metrics = sar3.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))

##
inc_deg = 55
obsv_metrics = sar3.calc_data_metrics(alt_km=h*1e-3, sc_speed_kmps=orb_speed*1e-3, sc_gnd_speed_kmps=gnd_spd*1e-3, inc_angle_deg=inc_deg, 
                                              instru_look_angle_from_target_inc_angle=True)
print(obsv_metrics)
N_looks = 1e6 / (obsv_metrics["ground pixel along-track resolution [m]"] * obsv_metrics["ground pixel cross-track resolution [m]"]) 
print('Speckle noise improvement is {} decibels'.format(10*numpy.log10(1/numpy.sqrt(N_looks))))
