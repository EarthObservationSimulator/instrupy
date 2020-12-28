""" 
.. module:: basic_sensor_model

:synopsis: *Module to handle sensor type with the most basic attributes.*

"""
import json
import numpy
import copy
import pandas, csv
from instrupy.util import Entity, Orientation, FieldOfView, MathUtilityFunctions, Constants, SensorGeometry, Maneuverability, SyntheticDataConfiguration
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
import scipy.interpolate
import metpy.interpolate
import astropy.time
import logging
logger = logging.getLogger(__name__)

class BasicSensorModel(Entity):
    """A basic sensor class. 
      
        :ivar name: Full name of the instrument.
        :vartype name: str
        
        :ivar acronym: Acronym, initialism, or abbreviation.
        :vartype acronym: str

        :ivar mass: Total mass (kg) of this entity.
        :vartype mass: float

        :ivar volume: Total volume (m3) of this entity.
        :vartype volume: float
        
        :ivar power: Nominal operating power (W).
        :vartype power: float

        :ivar orientation: Orientation of the instrument with respect to Nadir-frame.
        :vartype orientation: :class:`instrupy.util.Orientation`

        :ivar fieldOfView: Field of view specification of instrument. 
        :vartype fieldOfView: :class:`instrupy.util.FieldOfView`   

        :ivar maneuver: Maneuvarability specification of the instrument
        :vartype maneuver: :class:`instrupy.util.Maneuverability`  
       
        :ivar dataRate: Rate of data recorded (Mbps) during nominal operations.
        :vartype dataRate: float  

        :ivar bitsPerPixel: Bits encoded per pixel of image
        :vartype bitsPerPixel: int

        :ivar fieldOfRegard: Field of regard calculated taking into account FOV and manueverability of the payload.
        :vartype fieldOfRegard: :class:`instrupy.util.FieldOfView`  
   
    """

    def __init__(self, name=None, acronym=None, mass=None,
            volume=None, power=None,  orientation=None,
            fieldOfView=None, maneuver=None, dataRate=None, 
            syntheticDataConfig=None, bitsPerPixel=None, 
            numberDetectorRows=None, numberDetectorCols=None, _id=None):
        """Initialization

        """
        self.name = str(name) if name is not None else None
        self.acronym = str(acronym) if acronym is not None else self.name
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None else None
        self.fieldOfView = copy.deepcopy(fieldOfView) if fieldOfView is not None else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None else None
        self.dataRate = float(dataRate) if dataRate is not None else None
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None        
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.fieldOfView)
        self.syntheticDataConfig = copy.deepcopy(syntheticDataConfig) if syntheticDataConfig is not None else None
        self.numberDetectorRows = int(numberDetectorRows) if numberDetectorRows is not None else 4
        self.numberDetectorCols = int(numberDetectorCols) if numberDetectorCols is not None else 4

        super(BasicSensorModel,self).__init__(_id, "Basic Sensor")

    @staticmethod
    def from_dict(d):
        """Parses an instrument from a normalized JSON dictionary."""
        default_fov = dict({'sensorGeometry': 'CONICAL', 'fullConeAngle':25}) # default fov is a 25 deg conical
        default_orien = dict({"convention": "NADIR"}) #  default orientation = Nadir pointing
        default_manuv = dict({"@type": "FIXED"}) #  default maneuverability = Nadir pointing
        default_synobsconf = dict({"sourceFilePaths": None, "environVar": None, "interplMethod":None})
        return BasicSensorModel(
                name = d.get("name", None),
                acronym = d.get("acronym", None),
                mass = d.get("mass", None),
                volume = d.get("volume", None),
                power = d.get("power", None),
                orientation = Orientation.from_json(d.get("orientation", default_orien)),
                fieldOfView =  FieldOfView.from_json(d.get("fieldOfView", default_fov)),
                maneuver =  Maneuverability.from_json(d.get("maneuverability", default_manuv)),
                dataRate = d.get("dataRate", None),
                bitsPerPixel = d.get("bitsPerPixel", None),
                syntheticDataConfig = SyntheticDataConfiguration.from_json(d.get("syntheticDataConfig", default_synobsconf)),
                numberDetectorRows = d.get("numberDetectorRows", None),
                numberDetectorCols = d.get("numberDetectorCols", None),
                _id = d.get("@id", None)
                )

    def to_dict(self):
        return dict({
                "@type": "Basic Sensor",
                "name":self.name,
                "acronym":self.acronym,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "fieldOfView":self.fieldOfView.to_dict(),
                "orientation":self.orientation.to_dict(),
                "manuverability":self.maneuver.to_dict(),
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "syntheticDataConfig": self.syntheticDataConfig.to_dict(),
                "numberDetectorRows": self.numberDetectorRows,
                "numberDetectorCols": self.numberDetectorCols,
                "@id": self._id
                })

    def calc_typ_data_metrics(self, SpacecraftOrbitState, TargetCoords):
        ''' Calculate typical observation data metrics.

            :param SpacecraftOrbitState: Spacecraft state at the time of observation. This is approximately taken to be the middle (or as close as possible to the middle) of the access interval.

                               Dictionary keys are: 
                               
                               * :code:`Time[JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
                               * :code:`x[km]` (:class:`float`), :code:`y[km]` (:class:`float`), :code:`z[km]` (:class:`float`), Cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial-plane frame at the time of observation.
                               * :code:`vx[km/s]` (:class:`float`), :code:`vy[km/s]` (:class:`float`), :code:`vz[km/s]` (:class:`float`), Velocity of spacecraft in Earth Centered Inertial frame with equatorial-plane frame at the time of observation.
            :paramtype SpacecraftOrbitState: dict

            
            :param TargetCoords: Location of the observation.

                               Dictionary keys are: 
                                
                               * :code:`Lat [deg]` (:class:`float`), :code:`Lon [deg]` (:class:`float`), indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
            :paramtype TargetCoords: dict

            :returns: Typical calculated observation data metrics.
                      
                      Dictionary keys are: 
                    
                      * :code:`Coverage [T/F]` (:class:`bool`) indicating if observation was possible during the access event.
                      * :code:`Incidence angle [deg]` (:class:`float`) Incidence angle in degrees at target point calculated assuming spherical Earth.
                      * :code:`Look angle [deg]` (:class:`float`) Look angle in degrees at target point calculated assuming spherical Earth.
                      * :code:`Observation Range [km]` (:class:`float`) Distance in kilometers from satellite to ground-point during the observation.
                      * :code:`Solar Zenith [deg]` (:class:`float`) Solar Zenith angle in degrees during observation.

            :rtype: dict                      
        '''
        # Observation time in Julian Day UT1
        tObs_JDUT1 = SpacecraftOrbitState["Time[JDUT1]"]

        # Calculate Target position in ECI-frame
        TargetPosition_km = MathUtilityFunctions.geo2eci([TargetCoords["Lat [deg]"], TargetCoords["Lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position in Cartesian coordinates ECI-frame
        SpacecraftPosition_km = numpy.array([SpacecraftOrbitState["x[km]"], SpacecraftOrbitState["y[km]"], SpacecraftOrbitState["z[km]"]])  
        SpacecraftPosition_vel_kmps = numpy.array([SpacecraftOrbitState["vx[km/s]"], SpacecraftOrbitState["vy[km/s]"], SpacecraftOrbitState["vz[km/s]"]])  

        alt_km = numpy.linalg.norm(SpacecraftPosition_km) - Constants.radiusOfEarthInKM

        #  Calculate actual range between spacecraft and POI (Target)
        range_vector_km = TargetPosition_km - SpacecraftPosition_km

        range_km = numpy.linalg.norm(range_vector_km)

        # Calculate look angle
        look_angle = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(SpacecraftPosition_km)))
        look_angle_deg = numpy.rad2deg(look_angle)
        
        # Look angle to corresponding incidence angle conversion for spherical Earth
        incidence_angle = numpy.arcsin(numpy.sin(look_angle)*(Constants.radiusOfEarthInKM +alt_km)/Constants.radiusOfEarthInKM)
        incidence_angle_deg =  numpy.rad2deg(incidence_angle)

        # Solar zenith angle
        [solar_zenith, solar_distance] = MathUtilityFunctions.compute_sun_zenith(tObs_JDUT1, TargetPosition_km)
        if solar_zenith is not None:
            solar_zenith_deg =  numpy.rad2deg(solar_zenith)
        else:
            solar_zenith_deg = numpy.nan

        # assign sign to look-angle. positive sign => look is in opposite direction of the orbit-plane (i.e. the negative of the orbit plane normal vector) and vice-versa
        neg_orbit_normal = numpy.cross(SpacecraftPosition_vel_kmps, -1*SpacecraftPosition_km)
        sgn = numpy.sign(numpy.dot(range_vector_km, neg_orbit_normal))
        if(sgn==0):
            sgn = 1

        # For a basic-sensor, coverage is true for all access-events
        isCovered = True 
    
        obsv_metrics = {}
        obsv_metrics["Observation Range [km]"] = range_km
        obsv_metrics["Look angle [deg]"] = sgn*look_angle_deg
        obsv_metrics["Incidence angle [deg]"] = incidence_angle_deg
        obsv_metrics["Solar Zenith [deg]"] = solar_zenith_deg
        obsv_metrics["Coverage [T/F]"] = isCovered

        return obsv_metrics
    
    def synthesize_observation(self, time_JDUT1, pixel_center_pos):        

        source_fps = self.syntheticDataConfig.sourceFilePaths
        envvar = self.syntheticDataConfig.environVar
        interpl_method = self.syntheticDataConfig.interplMethod
        
        # get the source data file closest to the input time
        dataset_i = 0
        forecast_time_hrs = []
        for _src_fp in source_fps:
            dataset  = Dataset(_src_fp, "r", format="NETCDF4")
            _initial_time = dataset[envvar].initial_time
            _forecast_time = dataset[envvar].forecast_time
            _forecast_time_units = dataset[envvar].forecast_time_units
            if(_forecast_time_units != 'hours'):
                raise Exception('Forecast time units must be of hours')            
            forecast_time_hrs.append(_forecast_time)
            dataset_i = dataset_i + 1
        
        # extract the time in UTC from the string (example string is: 12/11/2020 (12:00))
        t = astropy.time.Time.strptime(_initial_time, '%m/%d/%Y (%H:%M)')
        initial_time_JDUT1 = astropy.time.Time(t, scale='utc').jd # TODO UTC / UT1?

        required_forecast_time_hrs = (time_JDUT1 - initial_time_JDUT1)*24
        # find dataset closest to the required forecast time
        idx = (np.abs(forecast_time_hrs - required_forecast_time_hrs)).argmin()

        dataset  = Dataset(source_fps[idx], "r", format="NETCDF4")
        lons = dataset.variables['lon_0'][:]
        lats = dataset.variables['lat_0'][:]
        var = dataset.variables[envvar][:]

        #print(lons)
        #print(lats)
        #print(temperature)

        # interpolate the variable data to the pixel center positions       
        if interpl_method == 'scipy.interpolate.linear':
            f = scipy.interpolate.interp2d( lons, lats, var, kind='cubic')
            interpl_data = []
            for _pix_p in pixel_center_pos:
                lon = _pix_p['lon[deg]']
                lat = _pix_p['lat[deg]']
                interpl_data.append(f(lon, lat)[0]) # [0] is needed to convert from the single element np.array to float
            logger.info("Inaccuracy around longitude=0 deg")
            

        elif interpl_method == 'metpy.interpolate.linear':
            (X,Y) = np.meshgrid(lons, lats)
            X = X.reshape((np.prod(X.shape),))
            Y = Y.reshape((np.prod(Y.shape),))
            coords = np.dstack((X,Y))[0]
            temperature = temperature.flatten()
            
            interpl_data = metpy.interpolate.interpolate_to_points(coords, temperature, pixel_center_pos, interp_type='linear', 
                                                    minimum_neighbors=3, gamma=0.25, kappa_star=5.052, 
                                                    search_radius=None, rbf_func='linear', rbf_smooth=0)


        #print(result)

        '''
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()

        plt.contourf(lons, lats, temperature,60, transform=ccrs.PlateCarree(), cmap=get_cmap("jet"))
        plt.colorbar(ax=ax, shrink=.98)
        plt.show()
        '''
        dataset.close()

        # return the interpoalted variable data
        return [pixel_center_pos, interpl_data, envvar]


