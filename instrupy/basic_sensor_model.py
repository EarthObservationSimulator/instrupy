""" 
.. module:: basic_sensor_model

:synopsis: *Module to handle sensor model with the most basic attributes.*

"""
import copy
from instrupy.util import Entity, Orientation, SphericalGeometry, Maneuver, ViewGeometry, MathUtilityFunctions, GeoUtilityFunctions, Constants, SyntheticDataConfiguration
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
import astropy.time
import uuid
from collections import namedtuple
import logging
logger = logging.getLogger(__name__)

class BasicSensorModel(Entity):
    """A basic sensor model class. 
      
        :ivar name: Full name of the instrument.
        :vartype name: str

        :ivar mass: Total mass (kg) of this instrument.
        :vartype mass: float

        :ivar volume: Total volume (m3) of this instrument.
        :vartype volume: float
        
        :ivar power: Nominal operating power (W) of this instrument.
        :vartype power: float

        :ivar orientation: Orientation of the instrument.
        :vartype orientation: :class:`instrupy.util.Orientation`

        :ivar fieldOfView: Field of view of instrument specification (SphericalGeometry and Orientation).
        :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

        :ivar sceneFieldOfView: Scene field of view specification (SphericalGeometry and Orientation).
        :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

        :ivar maneuver: Maneuver specification of the instrument. 
        :vartype maneuver: :class:`instrupy.util.Maneuver`  

        :ivar fieldOfRegard: Field of regard of the instrument taking into account the sceneFOV and manueverability of the satellite-sensor system. 
                             Note that this shall be a list in order to accommodate non-intersecting view-geometries.
                             TODO: Modify behavior to have FOR = sceneFOV when no maneuver is specified (hence fixed pointing). Currently FOR = None if manuever is not specified.
        :vartype fieldOfRegard: list, :class:`instrupy.util.ViewGeometry`  
       
        :ivar pointingOption: List of ``Orientation`` objects which specify the orientations into which the instrument-axis can be maneuvered. 
                               The orientations must be specified in the NADIR_POINTING frame.
        :vartype pointingOption: list, :class:`orbitpy.util.Orientation`

        :ivar dataRate: Rate of data recorded (Mega bits per sec) during nominal operations.
        :vartype dataRate: float  

        :ivar syntheticDataConfig: Synthetic data configuration for the instrument.
        :vartype syntheticDataConfig: :class:`instrupy.util.SyntheticDataConfiguration`  

        :ivar bitsPerPixel: Number of bits encoded per pixel of image.
        :vartype bitsPerPixel: int        

        :ivar numberDetectorRows: Number of detector rows (along the Y-axis of the SENOR_BODY_FIXED frame). If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the along-track direction.
        :vartype numberDetectorRows: int

        :ivar numberDetectorCols: Number of detector columns (along the X-axis of the SENOR_BODY_FIXED frame). If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the cross-track direction.
        :vartype numberDetectorCols: int

        :ivar _id: Unique instrument identifier.
        :vartype _id: str or int
   
    """

    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None,
                 fieldOfViewGeometry=None, sceneFieldOfViewGeometry=None, maneuver=None, pointingOption=None, dataRate=None, syntheticDataConfig=None, bitsPerPixel=None, 
                 numberDetectorRows=None, numberDetectorCols=None, _id=None):
        """ Initialization. All except the below two parameters have identical description as that of the corresponding class instance variables.

        :param fieldOfViewGeometry: Instrument field-of-view spherical geometry.
        :paramtype fieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :param sceneFieldOfViewGeometry: Scene field-of-view spherical geometry.
        :paramtype sceneFieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :return: None
        :rtype: None

        """
        self.name = str(name) if name is not None else None
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None and isinstance(orientation, Orientation) else None
        self.fieldOfView = ViewGeometry(orien=self.orientation, sph_geom=fieldOfViewGeometry) if self.orientation is not None and fieldOfViewGeometry is not None and isinstance(fieldOfViewGeometry, SphericalGeometry) else None
        self.sceneFieldOfView = ViewGeometry(orien=self.orientation, sph_geom=sceneFieldOfViewGeometry) if self.orientation is not None and sceneFieldOfViewGeometry is not None and isinstance(sceneFieldOfViewGeometry, SphericalGeometry) else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None and isinstance(maneuver, Maneuver) else None
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.sceneFieldOfView.sph_geom) if (self.maneuver is not None and self.sceneFieldOfView is not None) else None
        self.pointingOption = None
        if isinstance(pointingOption, list):
            if all(isinstance(x, Orientation) for x in pointingOption):
                self.pointingOption = pointingOption 
        elif isinstance(pointingOption, Orientation):
            self.pointingOption = [pointingOption] # make into single-element list
        self.dataRate = float(dataRate) if dataRate is not None else None
        self.syntheticDataConfig = copy.deepcopy(syntheticDataConfig) if syntheticDataConfig is not None and isinstance(syntheticDataConfig, SyntheticDataConfiguration) else None
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None                
        self.numberDetectorRows = int(numberDetectorRows) if numberDetectorRows is not None else None
        self.numberDetectorCols = int(numberDetectorCols) if numberDetectorCols is not None else None
        self._id = _id if _id is not None else None

        super(BasicSensorModel,self).__init__(_id, "Basic Sensor")

    @staticmethod
    def from_dict(d):
        """ Parses a ``BasicSensorModel`` object from a normalized JSON dictionary. 
            Refer to :ref:`basic_sensor_model_desc` for description of the accepted key/value pairs.
        
        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            orientation, Referenced and aligned to the SC_BODY_FIXED frame.
            fieldOfViewGeometry, CIRCULAR shape with 25 deg diameter.
            sceneFieldOfViewGeometry, fieldOfViewGeometry
            numberDetectorRows, 4
            numberDetectorRows, 4
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: ``BasicSensorModel`` object initialized with the input specifications.
        :rtype: :class:`instrupy.basic_sensor_model.BasicSensorModel`
            
        """
        instru_fov_geom = d.get("fieldOfViewGeometry", {'shape': 'CIRCULAR', 'diameter':25})  # default instrument FOV geometry is a 25 deg diameter circular shape
        scene_fov_geom = d.get("sceneFieldOfViewGeometry", instru_fov_geom)  # default sceneFOV geometry is the instrument FOV geometry
        default_orien = dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}) #  default orientation = referenced and aligned to the SC_BODY_FIXED frame.
        
        # parse the pointing options as a list of Orientation objects.
        pnt_opt_dict = d.get("pointingOption", None)
        _pointing_option = None
        if pnt_opt_dict:
            # translate to a list of Orientation objects
            if isinstance(pnt_opt_dict, list):
                _pointing_option = [Orientation.from_dict(x) for x in pnt_opt_dict]
            else:
                _pointing_option = [Orientation.from_dict(pnt_opt_dict)]
        
        return BasicSensorModel(
                name = d.get("name", None),
                mass = d.get("mass", None),
                volume = d.get("volume", None),
                power = d.get("power", None),
                orientation = Orientation.from_json(d.get("orientation", default_orien)),
                fieldOfViewGeometry =  SphericalGeometry.from_json(instru_fov_geom),
                sceneFieldOfViewGeometry = SphericalGeometry.from_json(scene_fov_geom),
                maneuver =  Maneuver.from_json(d.get("maneuver", None)),
                dataRate = d.get("dataRate", None),
                bitsPerPixel = d.get("bitsPerPixel", None),
                syntheticDataConfig = SyntheticDataConfiguration.from_json(d.get("syntheticDataConfig", None)),
                numberDetectorRows = d.get("numberDetectorRows", 4),
                numberDetectorCols = d.get("numberDetectorCols", 4),
                pointingOption = _pointing_option,
                _id = d.get("@id", str(uuid.uuid4()))
                )

    def to_dict(self):
        """ Translate the ``BasicSensorModel`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: ``BasicSensorModel`` specifications as python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfView.sph_geom.to_dict() if self.fieldOfView is not None and isinstance(self.fieldOfView, ViewGeometry) else None
        sceneFieldOfViewGeometry_dict = self.sceneFieldOfView.sph_geom.to_dict() if self.sceneFieldOfView is not None and isinstance(self.sceneFieldOfView, ViewGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        syntheticDataConfig_dict = self.syntheticDataConfig.to_dict() if self.syntheticDataConfig is not None and isinstance(self.syntheticDataConfig, SyntheticDataConfiguration) else None
        pointing_opt_dict = [Orientation.to_dict(x) for x in self.pointingOption] if self.pointingOption is not None else None
        return dict({
                "@type": "Basic Sensor",
                "name":self.name,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "fieldOfViewGeometry":fieldOfViewGeometry_dict,
                "sceneFieldOfViewGeometry": sceneFieldOfViewGeometry_dict,
                "orientation":orientation_dict,
                "maneuver":maneuver_dict,
                "pointingOption": pointing_opt_dict,
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "syntheticDataConfig": syntheticDataConfig_dict,
                "numberDetectorRows": self.numberDetectorRows,
                "numberDetectorCols": self.numberDetectorCols,
                "@id": self._id
                })

    def __repr__(self):
        return "BasicSensorModel.from_dict({})".format(self.to_dict())

    def calc_data_metrics(self, sc_orbit_state, target_coords):
        """ Calculate typical observation data metrics. Refer to :ref:`basic_sensor_model_desc` for description.            

        :param sc_orbit_state: Spacecraft state at the time of observation.

                            Dictionary keys are: 
                            
                            * :code:`time [JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
                            * :code:`x [km]` (:class:`float`), :code:`y [km]` (:class:`float`), :code:`z [km]` (:class:`float`), Cartesian coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
                            * :code:`vx [km/s]` (:class:`float`), :code:`vy [km/s]` (:class:`float`), :code:`vz [km/s]` (:class:`float`), Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        
        :paramtype sc_orbit_state: dict

        
        :param target_coords: Location of the observation. Also sometimes the Point-Of-Interest (POI).

                            Dictionary keys are: 
                            
                            * :code:`lat [deg]` (:class:`float`), :code:`lon [deg]` (:class:`float`) indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
        
        :paramtype target_coords: dict

        :returns: Calculated observation data metrics.
                    
                    Dictionary keys are: 
                
                    * :code:`incidence angle [deg]` (:class:`float`) Incidence angle in degrees at the target point calculated assuming spherical Earth.
                    * :code:`look angle [deg]` (:class:`float`) Look angle in degrees at the target point calculated assuming spherical Earth. Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
                    * :code:`observation range [km]` (:class:`float`) Distance in kilometers from satellite to ground-point during the observation.
                    * :code:`solar zenith [deg]` (:class:`float`) Solar Zenith angle in degrees during observation.
                
                    .. todo:: Include AT, CT footprint size calculations.

        :rtype: dict                      

        """
        # Observation time in Julian Day UT1
        tObs_JDUT1 = sc_orbit_state["time [JDUT1]"]

        # Calculate Target cartesian position in EARTH_CENTERED_INERTIAL frame
        target_pos = GeoUtilityFunctions.geo2eci([target_coords["lat [deg]"], target_coords["lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position in Cartesian coordinates in the EARTH_CENTERED_INERTIAL frame
        sc_pos = np.array([sc_orbit_state["x [km]"], sc_orbit_state["y [km]"], sc_orbit_state["z [km]"]])  
        sc_vel = np.array([sc_orbit_state["vx [km/s]"], sc_orbit_state["vy [km/s]"], sc_orbit_state["vz [km/s]"]])  

        alt_km = np.linalg.norm(sc_pos) - Constants.radiusOfEarthInKM # altitude

        #  Calculate the range vector between spacecraft and POI (Target)
        range_vector_km = target_pos - sc_pos

        range_km = np.linalg.norm(range_vector_km)

        # Calculate look angle
        # below snippet required to prevent invalid 'x' due to numerical errors
        x = np.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(sc_pos))
        if abs(x-1) < 1e-7:
            x=1
        look_angle = np.arccos(x)
        look_angle_deg = np.rad2deg(look_angle)
        
        # Look angle to corresponding incidence angle conversion for spherical Earth
        incidence_angle = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
        incidence_angle_deg =  np.rad2deg(incidence_angle)

        # Solar zenith angle
        [solar_zenith, solar_distance] = GeoUtilityFunctions.compute_sun_zenith(tObs_JDUT1, target_pos)
        if solar_zenith is not None:
            solar_zenith_deg =  np.rad2deg(solar_zenith)
        else:
            solar_zenith_deg = np.nan

        # assign sign to look-angle. 
        # positive sign => Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
        orbit_normal = np.cross(sc_pos, sc_vel)
        sgn = np.sign(np.dot(range_vector_km, orbit_normal))
        if(sgn==0):
            sgn = 1
    
        obsv_metrics = {}
        obsv_metrics["observation range [km]"] = round(range_km,1)
        obsv_metrics["look angle [deg]"] = round(sgn*look_angle_deg, 2)
        obsv_metrics["incidence angle [deg]"] = round(incidence_angle_deg, 2)
        obsv_metrics["solar zenith [deg]"] = round(solar_zenith_deg, 2)

        return obsv_metrics

    def synthesize_observation(self, time_JDUT1, pixel_center_pos):
        """ Compute the value of the geophysical variable (associated with the instrument) at the 
            input pixel (center) positions. 
            
            The interpolation (in space) of the source data is done using the interpolation method
            indicated in the synthetic data configuration of the instrument. The data file closest 
            to the input time of observation is chosen as the source data file.

        :param time_JDUT1: Observation time in Julian Date UT1.
        :paramtype time_JDUT1: float

        :param pixel_center_pos: Array of 'N' pixel center positions (longitudes, latitudes) of the pixels in the observation.
        :paramtype pixel_center_pos: array_like, shape (N, 2), float
        
        :returns: The interpolated data along with the corresponding pixel-center positions. The geophysical variable is also returned. Contents of the
                  namestuple are: 
                  
                  1. Pixel center positions (lon, lat) 
                  2. Interpolated data at the Pixel center positions
                  3. Geophysical variable
        :rtype: namedtuple, (array_like, shape (N, 2), float), (array_like, shape (N, 1), float), (str)

        """        
        source_fps = self.syntheticDataConfig.sourceFilePaths
        geophy_var = self.syntheticDataConfig.geophysicalVar
        interpl = self.syntheticDataConfig.get_interpolator()
        
        # get the source data file closest to the input time
        # collect all the datasets at the different forecast times
        forecast_time_hrs = [] # list of (relative wrt initial_time) forecast times
        for _src_fp in source_fps:
            dataset  = Dataset(_src_fp, "r", format="NETCDF4")
            _initial_abs_time = dataset[geophy_var].initial_time # start of the forecast (absolute time)
            _forecast_time = dataset[geophy_var].forecast_time # forecast time relative to the initial_time
            _forecast_time_units = dataset[geophy_var].forecast_time_units
            if(_forecast_time_units != 'hours'):
                raise Exception('Forecast time units must be of hours')            
            forecast_time_hrs.append(_forecast_time)
        
        # convert the initial_time in UTC (string format) to JDUT1 (example string is: 12/11/2020 (12:00))
        t = astropy.time.Time.strptime(_initial_abs_time, '%m/%d/%Y (%H:%M)') 
        initial_time_UTC = astropy.time.Time(t, scale='utc')
        initial_time_UT1 = initial_time_UTC.ut1 # convert from UTC to UT1
        initial_time_JDUT1 = initial_time_UT1.jd # convert to Julian Date format

        # required forecast time relative to the initial_time
        required_forecast_time_hrs = (time_JDUT1 - initial_time_JDUT1)*24 

        # find dataset closest to the required forecast time
        idx = (np.abs(forecast_time_hrs - required_forecast_time_hrs)).argmin()

        dataset  = Dataset(source_fps[idx], "r", format="NETCDF4")
        lons = dataset.variables['lon_0'][:]
        lats = dataset.variables['lat_0'][:]
        var_data = dataset.variables[geophy_var][:]

        #print(lons)
        #print(lats)
        #print(temperature)

        # interpolate the variable data to the pixel center positions   
        interpl_data = interpl(lons, lats, var_data, pixel_center_pos)

        #print(result)

        '''
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()

        plt.contourf(lons, lats, temperature,60, transform=ccrs.PlateCarree(), cmap=get_cmap("jet"))
        plt.colorbar(ax=ax, shrink=.98)
        plt.show()
        '''
        dataset.close()

        # return the interpolated variable data
        syn_data = namedtuple("syn_data", ["pixel_center_pos", "interpl_data", "geophy_var"])

        return syn_data(
            pixel_center_pos,
            interpl_data,
            geophy_var
        )
    
    def get_id(self):
        """ Get the instrument identifier.

        :returns: instrument identifier.
        :rtype: str

        """
        return self._id
    
    def get_field_of_view(self):
        """ Get the instrument field-of-view.

        :returns: Instrument field-of-view. 
        :rtype: :class:`instrupy.util.ViewGeometry`

        """
        return self.fieldOfView

    def get_scene_field_of_view(self):
        """ Get the scene-field-of-view (sceneFOV).

        :returns: Scene-field-of-view.
        :rtype: :class:`instrupy.util.ViewGeometry`

        """
        return self.sceneFieldOfView

    def get_field_of_regard(self):
        """ Get the instrument field of regard. 

        :returns: Field of regard (list of ``ViewGeometry`` objects). 
        :rtype: list, :class:`instrupy.util.ViewGeometry`

        """
        return self.fieldOfRegard

    def get_orientation(self):
        """ Get the instrument orientation.

        :returns: Instrument orientation.
        :rtype: :class:`instrupy.util.Orientation`

        """
        return self.orientation
    
    def get_pointing_option(self):
        """ Get the list of pointing options.

        :returns: List of pointing options.
        :rtype: list, :class:`instrupy.util.Orientation`

        """
        return self.pointingOption

    def get_pixel_config(self):
        """ Get pixel-configuration.
        
        :returns: Pixel-configuration, i.e. the number of detector rows and columns.
        :rtype: namedtuple, (int, int)
        
        """        
        pixel_config = namedtuple("pixel_config", ["numberDetectorRows", "numberDetectorCols"])

        return pixel_config(self.numberDetectorRows, self.numberDetectorCols)


