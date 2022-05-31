""" The following examples illustrate the SAR models with different possible set of configurations.
Please refer to
* html docs -> Instrument Models -> Synthetic Aperture Radar Model 
* html docs -> API Reference -> Synthetic Aperture Radar Model

"""
import json
import numpy
import sys, os
from flask import Flask, request


from instrupy.radiometer_model import RadiometerModel

Re = 6378.137e3 # [m]
c = 299792458 # m/s 

def evaluate_rad(daz_m, delv_m, alt):
    rad_dict = {"@type": "Radiometer",
                "name": "Test radiometer",
                "orientation": {
                   "referenceFrame": "NADIR_POINTING",
                   "convention": "REF_FRAME_ALIGNED"
                },
                "antenna": {"shape": "RECTANGULAR", "height": str(daz_m), "width": str(delv_m), "radiationEfficiency": 0.6,
                           "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM", "phyTemp": 50},
                "operatingFrequency": 1275.7e6,
                "system": {
                    "@type": "TOTAL_POWER",
                    "integrationTime": 0.1,
                    "bandwidth": 40e6,
                    "predetectionGain": 83,
                    "predetectionInpNoiseTemp": 200,
                    "predetectionGainVariation": 2,
                    "integratorVoltageGain": 1
                },
                "scan":{
                    "@type": "FIXED"
                },
                "targetBrightnessTemp": 290
                }
    test_rad = RadiometerModel.from_dict(rad_dict)
    h = float(alt)*1e3
    orb_speed = numpy.sqrt(3.986004418e14/(Re + h)) # [m/s] orbital speed
    gnd_spd = orb_speed*(Re/(Re+h))
    obsv_metrics = test_rad.calc_data_metrics(alt_km = h*1e-3, gnd_spd = orb_speed*(Re/(Re+h)))
    return obsv_metrics

app = Flask(__name__)

@app.route('/', methods=['GET','POST'])
def index():
	if request.method == 'POST':
		print(request.values)
		h = request.form.get("height")
		w = request.form.get("width")
		alt = request.form.get("altitude")
		metrics = evaluate_rad(h,w,alt)
		print(metrics)
		return metrics
	else:
		return "hello world"

if __name__ == '__main__':
	app.run(port='5000')
