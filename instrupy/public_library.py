""" 
.. module:: public_library

:synopsis: *Provides easy interface to the package main functionalities, independent of the instrument type.*

"""
from .util import Entity, MathUtilityFunctions, SensorGeometry, FileUtilityFunctions
from .basic_sensor import BasicSensor
from .passive_optical_scanner import PassiveOpticalScanner
from .synthetic_aperture_radar import SyntheticApertureRadar
import sys
import pandas, csv
import json
import numpy

class InstrumentCoverageParameters():
        """ Data structure to hold instrument coverage related parameters.

        :ivar _id: Payload identifier (Optional)
        :vartype _id: str

        :ivar fov_geom: Geometry of the field-of-view
        :vartype fov_geom: :class:`FOVGeometry`

        :ivar fov_cone: List of cone angles in degrees
        :vartype fov_cone: list, float

        :ivar fov_clock: List of clock angles in degrees
        :vartype fov_clock: list, float

        :ivar fov_at: Along-track FOV in degrees
        :vartype fov_at: float

        :ivar fov_ct: Cross-track FOV in degrees
        :vartype fov_ct: float

        :ivar orien_eu_seq1: Euler sequence 1
        :vartype orien_eu_seq1: int

        :ivar orien_eu_seq2: Euler sequence 2
        :vartype orien_eu_seq2: int

        :ivar orien_eu_seq3: Euler sequence 3
        :vartype orien_eu_seq3: int

        :ivar orien_eu_ang1: Euler angle 1 in degrees
        :vartype orien_eu_ang1: float

        :ivar orien_eu_ang2: Euler angle 2 in degrees
        :vartype orien_eu_ang2: float

        :ivar orien_eu_ang3: Euler angle 3 in degrees
        :vartype orien_eu_ang3: float

        :ivar purely_side_look: Flag to specify if instrument operates in a strictly side-looking viewing geometry.
        :vartype purely_side_look: bool

        :ivar yaw180_flag: Flag applies in case of field-of-regard. If true, it signifies that the field-of-regard includes the field-of-view of payload rotated about nadir by 180 deg. 
        :vartype yaw180_flag: bool
        
        """
        def __init__(self, _id = None, fov_geom = None, fov_cone = None, fov_clock = None, fov_at = None, fov_ct = None, 
                    orien_eu_seq1 = None, orien_eu_seq2 = None, orien_eu_seq3 = None, 
                    orien_eu_ang1 = None, orien_eu_ang2 = None, orien_eu_ang3 = None,
                    purely_side_look = None, yaw180_flag = None):
            
    
            self._id = str(_id) if _id is not None else None
            self.fov_geom = SensorGeometry.get(fov_geom) if fov_geom is not None else None
            self.fov_cone = [float(i) for i in fov_cone] if fov_cone is not None else None
            self.fov_clock = [float(i) for i in fov_clock] if fov_clock is not None else None # clock can be "None" if Conical Sensor
            if(self.fov_clock is None):
                if(self.fov_geom=='CONICAL'):
                    self.fov_clock = [0]
                else:
                    pass
                    #raise RuntimeError("Missing clock angles from instrument coverage specifications for non-conical sensor fov.")
            self.fov_at = float(fov_at) if fov_at is not None else None
            self.fov_ct = float(fov_ct) if fov_ct is not None else None
            self.orien_eu_seq1 = int(orien_eu_seq1) if orien_eu_seq1 is not None else None
            self.orien_eu_seq2 = int(orien_eu_seq2) if orien_eu_seq2 is not None else None
            self.orien_eu_seq3 = int(orien_eu_seq3) if orien_eu_seq3 is not None else None
            self.orien_eu_ang1 = float(orien_eu_ang1) if orien_eu_ang1 is not None else None
            self.orien_eu_ang2 = float(orien_eu_ang2) if orien_eu_ang2 is not None else None
            self.orien_eu_ang3 = float(orien_eu_ang3) if orien_eu_ang3 is not None else None
            self.purely_side_look = bool(purely_side_look) if purely_side_look is not None else None
            self.yaw180_flag = bool(yaw180_flag) if yaw180_flag is not None else None


        def get_as_string(self, param):
            """ Get the necessary instrument coverage specifications in string format so it can be directly passed
                as arguments to the :code:`orbitpropcov` program.
            """
            if(param == 'Geometry'):
                return self.fov_geom.name
            elif(param == 'Clock'):
                clock = [str(i) for i in self.fov_clock] 
                clock = ','.join(clock) 
                return clock
            elif(param == 'Cone'):
                cone = [str(i) for i in self.fov_cone] 
                cone = ','.join(cone) 
                return cone 
            elif(param == 'Orientation'):
                orien = str(self.orien_eu_seq1) + ',' + str(self.orien_eu_seq2) + ',' + str(self.orien_eu_seq3) + ',' + \
                    str(self.orien_eu_ang1) + ',' + str(self.orien_eu_ang2) + ',' + str(self.orien_eu_ang3)
                return orien
            elif(param == 'yaw180_flag'):
                return str(int(self.yaw180_flag))
            else:
                raise RuntimeError("Unknown parameter") 

class Instrument(Entity): 
    """ A class to serve as single interface to different types of instruments.
     
        :ivar _sensor: Object of relevant instrument.
        :vartype _sensor: :class:`instrupy.basic_sensor.BasicSensor` (or) :class:`instrupy.passive_optical_sensor.PassiveOpticalSensor` (or) :class:`instrupy.synthetic_aperture_radar.SyntheticApertureRadar`

        :ivar _instru_id: Instrument id
        :ivartype _instru_id: str

        :ivar _ssen_id: List of subsensor ids
        :ivartype _ssen_id: list, int
    """
    def __init__(self, specs = dict()):
        """ Initialize an instrument.

        :param specs:  Give relevant specifications according to type of instrument.
                       See :ref:`Instrument search-space specifications`. 
        
        :type specs: dict

        """

        switchInstrumentTypes = {"Basic Sensor": BasicSensor, 
                                 "Passive Optical Scanner": PassiveOpticalScanner,
                                 "Synthetic Aperture Radar": SyntheticApertureRadar        
                                 }
        
        try: 
            _type = specs["@type"]            
        except: 
            raise Exception("Type of instrument not defined.")   

        assert(switchInstrumentTypes.get(_type, None)),"Specified instrument type not supported."          
        _instrument = switchInstrumentTypes.get(_type, None)  
        
        try:
            self._sensor = _instrument.from_dict(specs)
        except Exception as e:
            self._sensor = None  
            raise Exception("Something wrong in initializing instrument \n" + str(e))

        self._instru_id = self._sensor.instru_id
        self._ssen_id = self._sensor.ssen_id

        super(Instrument,self).__init__(specs.get("@id", None), "Instrument")

    @staticmethod
    def from_dict(d): 
        """ Parses a instrument from a normalized JSON dictionary. """
        return Instrument(d)

    def to_dict(self):
        return self._sensor.to_dict()
    
    def get_id(self):
        return self._sensor.get_id()

    def get_type(self):
        return self._sensor.get_type()
    
    def get_fov(self, subsensor_id = None):       
        return self._sensor.get_fov(subsensor_id)
    
    def get_pixel_config(self, subsensor_id = None):
        return self._sensor.get_pixel_config(subsensor_id)

    def synthesize_observation(self, time_JDUT1, pixel_center_pos, interpl_method='scipy.interpolate.linear', subsensor_id = None):        
        return self._sensor.synthesize_observation(time_JDUT1=time_JDUT1, pixel_center_pos=pixel_center_pos, ssen_id=subsensor_id)  

    def calc_typ_data_metrics(self, SpacecraftOrbitState, TargetCoords, subsensor_id = None):   
        """ Calculate the typical observation data-metrics.
        
        :param SpacecraftOrbitState: Spacecraft state at the time of observation. This is approximately taken to be the middle (or as close as possible to the middle) of the access interval.

        Dictionary keys are: 
        
        * :code:`Time[JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
        * :code:`x[km]` (:class:`float`), :code:`y[km]` (:class:`float`), :code:`z[km]` (:class:`float`), Cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane at the time of observation.
        * :code:`vx[km/s]` (:class:`float`), :code:`vy[km/s]` (:class:`float`), :code:`vz[km/s]` (:class:`float`), velocity of spacecraft in Earth Centered Inertial frame with equatorial plane at the time of observation.
        
        :paramtype SpacecraftOrbitState: dict
        
        :param TargetCoords: Location of the observation.

                            Dictionary keys are: 
                            
                            * :code:`Lat [deg]` (:class:`float`), :code:`Lon [deg]` (:class:`float`), indicating the corresponding ground-point accessed (latitude, longitude) in degrees.

        :paramtype TargetCoords: dict

        :returns: Typical calculated observation data metrics specific to the instrument type.
        :rtype: dict

        """ 
        return self._sensor.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords, subsensor_id)   

    def get_coverage_specs(self, _ssen_id = None):
        """ Get field-of-view (or scene-field-of-view) and the field-of-regard specifications for the instrument-subsensor. Also returns
            the orientation (common for both FOV, FOR) and if the instrument is constrained to take observations at a strictly 
            side-looking (no squint) target geometry. The :code:`yaw180_flag` indicates if the instrument coverage includes
            orientation of the instrument rotated by 180 deg about the yaw axis from the nominal orientation. This is significant
            for the case of FOR specifications and purely side-looking instruments.
            
            .. todo:: update unit test

            :param _ssen_id: subsensor id
            :paramtype _ssen_id: int

            :returns: JSON string with the coverage specifications            
            :rtype: dict
        
        """
        if _ssen_id is not None:
            _ssen_id_indx = self._ssen_id.index(_ssen_id)
            _ssen = self._sensor.subsensor[_ssen_id_indx]
        else:
            _ssen = self._sensor.subsensor[0]

        # determine if instrument takes observations at a purely side looking geometry
        purely_side_look = False
        instru_type = type(self._sensor).__name__
        if(instru_type == 'SyntheticApertureRadar'):
            purely_side_look = True
        """if(instru_type == 'PassiveOpticalScanner'):
            if(_ssen.scanTechnique == "PUSHBROOM" or _ssen.scanTechnique == "WHISKBROOM"):
                purely_side_look = True """

        orientation_specs = {"eulerAngle1": _ssen.orientation.euler_angle1,
                            "eulerAngle2": _ssen.orientation.euler_angle2, 
                            "eulerAngle3": _ssen.orientation.euler_angle3,
                            "eulerSeq1": _ssen.orientation.euler_seq1,
                            "eulerSeq2": _ssen.orientation.euler_seq2,
                            "eulerSeq3": _ssen.orientation.euler_seq3,
                            }
        
        for_specs = {"geometry": _ssen.fieldOfRegard._geometry,
                    "coneAnglesVector" : _ssen.fieldOfRegard._coneAngleVec_deg,
                    "clockAnglesVector": _ssen.fieldOfRegard._clockAngleVec_deg,
                    "AlongTrackFov": _ssen.fieldOfRegard._AT_fov_deg,
                    "CrossTrackFov": _ssen.fieldOfRegard._CT_fov_deg,
                    "yaw180_flag": _ssen.fieldOfRegard._yaw180_flag
                    }

        if(hasattr(_ssen, 'sceneFieldOfView')):
            if(_ssen.sceneFieldOfView):
                # if there exists a scene FOV, pass the corresponding fov specifications for coverage caclulations
                fov_specs = {"geometry": _ssen.sceneFieldOfView._geometry,
                            "coneAnglesVector" : _ssen.sceneFieldOfView._coneAngleVec_deg,
                            "clockAnglesVector": _ssen.sceneFieldOfView._clockAngleVec_deg,
                            "AlongTrackFov": _ssen.sceneFieldOfView._AT_fov_deg,
                            "CrossTrackFov": _ssen.sceneFieldOfView._CT_fov_deg,
                            "yaw180_flag": _ssen.sceneFieldOfView._yaw180_flag
                            }
                result = {"@id": _ssen._id, "Orientation": orientation_specs, "fieldOfView": fov_specs, "purely_side_look": purely_side_look, "fieldOfRegard": for_specs}
                return json.dumps(result)   
            
        fov_specs = {"geometry": _ssen.fieldOfView._geometry,
                    "coneAnglesVector" : _ssen.fieldOfView._coneAngleVec_deg,
                    "clockAnglesVector": _ssen.fieldOfView._clockAngleVec_deg,
                    "AlongTrackFov": _ssen.fieldOfView._AT_fov_deg,
                    "CrossTrackFov": _ssen.fieldOfView._CT_fov_deg,
                    "yaw180_flag": _ssen.fieldOfView._yaw180_flag
                    }

        result = {"@id": _ssen._id, "Orientation": orientation_specs, "fieldOfView": fov_specs, "purely_side_look": purely_side_look, "fieldOfRegard": for_specs}
        return json.dumps(result)

    def get_FOV_FOR_objs(self, _ssen_id):
        """ Obtain field-of-view (FOV), field-of-regard (FOR) for a given instrument, sub-sensor in the 
            :class:`instrupy.public_library.InstrumentCoverageParameters` object format.

        :param _ssen_id: Subsensor identifier
        :paramtype _ssen_id:  int

        :returns: FOV and FOR related parameters
        :rtype: list, :class:`InstrumentCoverageParameters`
        
        """       
        _ics = FileUtilityFunctions.from_json(self.get_coverage_specs(_ssen_id))
        ics_fldofview = InstrumentCoverageParameters(_ics["@id"], _ics["fieldOfView"]["geometry"], 
                                            _ics["fieldOfView"]["coneAnglesVector"], _ics["fieldOfView"]["clockAnglesVector"],
                                            _ics["fieldOfView"]["AlongTrackFov"], _ics["fieldOfView"]["CrossTrackFov"],
                                            _ics["Orientation"]["eulerSeq1"], _ics["Orientation"]["eulerSeq2"], _ics["Orientation"]["eulerSeq3"],
                                            _ics["Orientation"]["eulerAngle1"], _ics["Orientation"]["eulerAngle2"], _ics["Orientation"]["eulerAngle3"],
                                            _ics["purely_side_look"], _ics["fieldOfView"]["yaw180_flag"]
                                            )           
        ics_fldofreg = InstrumentCoverageParameters(_ics["@id"], _ics["fieldOfRegard"]["geometry"], 
                                            _ics["fieldOfRegard"]["coneAnglesVector"], _ics["fieldOfRegard"]["clockAnglesVector"],
                                            _ics["fieldOfRegard"]["AlongTrackFov"], _ics["fieldOfRegard"]["CrossTrackFov"],
                                            _ics["Orientation"]["eulerSeq1"], _ics["Orientation"]["eulerSeq2"], _ics["Orientation"]["eulerSeq3"],
                                            _ics["Orientation"]["eulerAngle1"], _ics["Orientation"]["eulerAngle2"], _ics["Orientation"]["eulerAngle3"],
                                            _ics["purely_side_look"], _ics["fieldOfRegard"]["yaw180_flag"]
                                            )           
        #print("Instrument + subsensor ID is: ", _ics["@id"])
        #print("fieldOfRegard is: ", _ics["fieldOfRegard"])
        #print("Orientation:", _ics["Orientation"])
        #print("purely_side_look:", _ics["purely_side_look"])
        return [ics_fldofview, ics_fldofreg]   

      

################################################### Legacy functions #################################################################    
    def generate_level0_data_metrics(self, POI_filepath, AccessInfo_filepath, Result_filepath):
 
        #  '''############################### Deprecated ############################### 
        
        #     Generate typical data metrics per access event, per grid-point. 
        #     This function iteratively calls :code:`calc_typ_data_metrics` of the specified instrument type over all access events 
        #     available in the input file.
        #     CSV formatted files are supported to supply the required input data.

        #     :param POI_filepath: Filepath to CSV file containing lat/lon co-ordinates of points of interest along with index.

        #                          First row contains the following header elements, and the following rows contain the data corresponding to these headers.
                                 
        #                          Description of the header elements:
                            
        #                          * :code:`POI`, Index of the point-of-interest.
        #                          * :code:`Lat [deg]`, :code:`Lon [deg]`, latitude, longitude of the point-of-interest.  

        #                          .. note:: Make sure the header titles are as specified above, and the delimiters are commas.
        #     :paramtype POI_filepath: str

             
        #     :param AccessInfo_filepath: Filepath to CSV file containing data of access events and their time-intervals.
        #                        First three rows convey general information. The fourth row conveys the Epoch in Julian Days UT1.
        #                        The fifth row  contains the following header elements, and the following rows contain the data corresponding to these headers.

        #                        Description of the header elements:
                                
        #                        * :code:`AccessFrom[Days]`,  The time at which access starts in [days], referenced to epoch in row-4.
        #                        * :code:`Duration[s]`, Duration of access in [s].
        #                        * :code:`POI` indicating index of ground-point.
        #                        * :code:`EventNum` indicating index of event.
        #                        * :code:`Time[Days]`, Time in [Days] at which the alongside satellite-state is recorded. Referenced to the epoch specified in row-4.
        #                        * :code:`X[km]`, :code:`Y[km]` :code:`Z[km]`, cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
        #                        * :code:`VX[km/s]`, :code:`VY[km/s]`, :code:`VZ[km/s]`, velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.


        #     :paramtype AccessInfo_filepath: str

        #     :param Result_filepath: Filepath to CSV file in which the results are written
        #                         Description of the header elements:
                                
        #                        * :code:`POI` indicating index of ground-point.
        #                        * :code:`Access From [JDUT1]`,  The time at which access starts in [JDUT1]
        #                        * :code:`Duration [s]`, Duration of access in [s].
        #                        * + other header elements specific to the instrument type

        #                        .. note:: this is an **output** file of the function
            
        #  ''':paramtype Result_filepath: str

        
        epoch_JDUT1 = pandas.read_csv(AccessInfo_filepath, skiprows = [0], nrows=1, header=None) # 2nd row contains the epoch
        epoch_JDUT1 = float(epoch_JDUT1[0][0].split()[2])
                
        poi_info_df = pandas.read_csv(POI_filepath)

        # Read the local access.csv file
        access_info_df = pandas.read_csv(AccessInfo_filepath,skiprows = [0,1,2,3]) # read the access times (corresponding to the DSM in which the instrument was used)
        eventIdxArray = access_info_df['EventNum']

        # erase any old file and create new one
        with open(Result_filepath,'w') as f:
            w = csv.writer(f)

        with open(Result_filepath,'a+', newline='') as f:
            w = csv.writer(f)

            # Iterate over all logged access events
            idx = 0
            for eventIdx in eventIdxArray:   
                AccessInfo = dict()
                AccessInfo["Access From [JDUT1]"] = epoch_JDUT1 + access_info_df.loc[idx]["AccessFrom[Days]"] 
                AccessInfo["Access Duration [s]"] = access_info_df.loc[idx]["Duration[s]"] 
                
                poi_indx = access_info_df.loc[idx]["POI"] 
                AccessInfo["Lat [deg]"] = poi_info_df["lat[deg]"][list(poi_info_df["POI"]).index(poi_indx)]
                AccessInfo["Lon [deg]"] = poi_info_df["lon[deg]"][list(poi_info_df["POI"]).index(poi_indx)]

                SpacecraftOrbitState = dict()
                SpacecraftOrbitState["Time[JDUT1]"] = epoch_JDUT1 + access_info_df.loc[idx]["Time[Days]"] 
                SpacecraftOrbitState["x[km]"] = access_info_df.loc[idx]["X[km]"] 
                SpacecraftOrbitState["y[km]"] = access_info_df.loc[idx]["Y[km]"] 
                SpacecraftOrbitState["z[km]"] = access_info_df.loc[idx]["Z[km]"] 
                SpacecraftOrbitState["vx[km/s]"] = access_info_df.loc[idx]["VX[km/s]"] 
                SpacecraftOrbitState["vy[km/s]"] = access_info_df.loc[idx]["VY[km/s]"] 
                SpacecraftOrbitState["vz[km/s]"] = access_info_df.loc[idx]["VZ[km/s]"] 

                obsv_metrics = self._sensor.calc_typ_data_metrics(SpacecraftOrbitState, AccessInfo) # calculate the data metrics specific to the instrument type
                _v = dict({'Access From [JDUT1]':AccessInfo["Access From [JDUT1]"], 'Access Duration [s]':AccessInfo["Access Duration [s]"], 'POI index': access_info_df['POI'][idx]}, **obsv_metrics)
                
                if idx==0: #1st iteration
                    w.writerow(_v.keys())    

                w.writerow(_v.values())

                idx = idx + 1  

    def generate_level1_data_metrics(self, inFilePath, outFilePath):
        ############################### Deprecated ###############################
        
        # Concatenate all the level-0 metrics into a single data-frame
        level0_df = None

        
        if(isinstance(inFilePath, str)):
            inFilePath = [inFilePath]
        
        for _filePath in inFilePath:
            if level0_df is None:
                level0_df = pandas.read_csv(_filePath)       
            else:
                level0_df = level0_df.append(pandas.read_csv(_filePath), ignore_index=True)

        # Sort and Index level0 DF by POI
        level0_df = level0_df.set_index('POI index')
        level0_df = level0_df.sort_index()
     
        # Only consider Valid observation events, not all Access events
        level0_df = level0_df.loc[level0_df['Coverage [T/F]'] == True]

        # Drop some columns
        level0_df = level0_df.drop(columns=['Coverage [T/F]'])
        level0_df = level0_df.drop(columns=['Access From [JDUT1]'])
        level0_df = level0_df.drop(columns=['Access Duration [s]'])

        # Find the mean of observation metrics at each POI
        level1_df_mean = level0_df.groupby('POI index').mean()
        # Rename the column names
        columns_list = list(level1_df_mean.columns.values)
        for col in columns_list:
            level1_df_mean = level1_df_mean.rename(columns={col: 'Mean of ' + col})

        # Find Standard deviation of observation metrics at each POI
        level1_df_sd = level0_df.groupby('POI index').std()
        columns_list = list(level1_df_sd.columns.values)
        for col in columns_list:
            level1_df_sd = level1_df_sd.rename(columns={col: 'SD of ' + col})

        level1_metrics = pandas.concat([level1_df_mean, level1_df_sd], axis=1)  
        level1_metrics.to_csv(outFilePath, index=True, header = True, na_rep='nan')

    def generate_level2_data_metrics(self, inFilePath, outFilePath):
        ############################### Deprecated ###############################

        level1_df = pandas.read_csv(inFilePath)       

        level1_df = level1_df.drop(columns=['POI index'])

        # Find the mean of observation metrics over all POIs
        level2_df_mean = level1_df.mean().to_frame()
        level2_df_mean = level2_df_mean.T
        # Rename the column names
        columns_list = list(level2_df_mean.columns.values)
        
        for col in columns_list:
            level2_df_mean = level2_df_mean.rename(columns={col: 'Mean of ' + col})


        # Find Standard deviation of observation metrics over all POIs
        level2_df_sd = level1_df.std().to_frame().T
        columns_list = list(level2_df_sd.columns.values)
        for col in columns_list:
            level2_df_sd = level2_df_sd.rename(columns={col: 'SD of ' + col})

        level2_metrics = pandas.concat([level2_df_mean, level2_df_sd], axis=1)
        level2_metrics.to_csv(outFilePath,index=False, na_rep='nan')


    def generate_level1_coverage_metrics(self, missionDuration_days, inFilePath, outFilePath):
        ############################### Deprecated ###############################
        # Concatenate all the level-0 metrics into a single data-frame
        level0_df = None

        if(isinstance(inFilePath, str)):
            inFilePath = [inFilePath]

        for _filePath in inFilePath:
            if level0_df is None:
                level0_df = pandas.read_csv(_filePath)       
            else:
                level0_df = level0_df.append(pandas.read_csv(_filePath), ignore_index=True)

        # Keep only columns of interest 
        level0_df = level0_df.filter(items=['POI index', 'Access From [JDUT1]', 'Access Duration [s]', 'Coverage [T/F]'])

        # Sort and Index level0 DF by POI
        level0_df = level0_df.set_index('POI index')
        level0_df = level0_df.sort_index()

        # Only consider Valid observation events, not all Access events
        level0_df = level0_df.loc[level0_df['Coverage [T/F]'] == True]

        # Delete the Coverage [T/F] column
        level0_df = level0_df.drop(columns=['Coverage [T/F]'])    

        # Find the percentage coverage metric for each of the POIs
        percentage_cov_df = level0_df['Access Duration [s]'].groupby('POI index').sum()*100/(missionDuration_days * 86400.0)
        percentage_cov_df = percentage_cov_df.to_frame('Percentage Coverage')

        # Compute gaps at each of the POIs, ignore the gaps at the extremes (starting and ending)
        level0_df = level0_df.groupby('POI index')

        level1_cov_df = pandas.DataFrame(columns=['POI index', 'Max Coverage Gap [Days]', 'Mean Coverage Gap [Days]', 'Time Average Gap [Days]'])

        for name, data_at_poi in level0_df:
            
            # Sort the data for the poi according to Access From [JDUT1]
            data_at_poi = data_at_poi.sort_values(by=['Access From [JDUT1]'])
            
            # Add Access To [JDUT1]
            data_at_poi['Access To [JDUT1]'] = data_at_poi['Access From [JDUT1]'] + data_at_poi['Access Duration [s]']/(60.0*60.0*24.0)
            
            # Delete the Access Duration [s] column
            data_at_poi = data_at_poi.drop(columns=['Access Duration [s]'])  
            
            # get list of access dates From, To alternating at the poi. First entry corresponds to Access From event
            access_events = data_at_poi.values.flatten()
            access_events = pandas.DataFrame({'Time[JDUT1]':access_events}) # convert to Dataframe
            
            # There maybe cases with two or more satellites seeing a poi at same time, in which case their access periods will overlap
            # Current code checks for such a condition, and raises exception if encountered. TODO: Write code to process this situation.
            if(access_events['Time[JDUT1]'].is_monotonic_increasing is False):
                raise Exception("Not prepared for this situation.")
                        
            # remove the 0th entry and odd entries (1st, 3rd, ....) since they correspond to access-duration
            # the remaining data-frame are the gaps at the poi
            gaps = access_events.diff()
            gaps = gaps.iloc[::2]
            gaps = gaps.iloc[1:]
            if gaps.empty: # if no gaps avaliable for the poi, skip to next iteration
                continue            
            
            # compute the gap metrics at the poi
            
            max_coverage_gap_days = gaps.max() # Find the max coverage gap metric at the poi
            
            mean_coverage_gap_days = gaps.mean()

            time_average_gap_days = (gaps**2).sum()/ (missionDuration_days) #if gaps.notnull().any() else numpy.null
            
  
            _cov = {'POI index': [data_at_poi.index.values[0]],
                    'Max Coverage Gap [Days]': [max_coverage_gap_days.values[0]],
                    'Mean Coverage Gap [Days]': [mean_coverage_gap_days.values[0]], 
                    'Time Average Gap [Days]': [time_average_gap_days.values[0]],
                    'Gaps [Days]': list(gaps.transpose().values)
                   }
            level1_cov_df = level1_cov_df.append(pandas.DataFrame(data=_cov), ignore_index=True, sort=True)

        level1_cov_df = level1_cov_df.set_index('POI index')
        level1_cov_df = level1_cov_df.sort_index()
        level1_cov_df = pandas.concat([level1_cov_df, percentage_cov_df], axis=1)  
        level1_cov_df.to_csv(outFilePath, index=True, header = True, na_rep='nan')

    def generate_level2_coverage_metrics(self, inFilePath, outFilePath):
        ############################### Deprecated ###############################

        level1_df = pandas.read_csv(inFilePath)       

        level1_df = level1_df.drop(columns=['POI index'])

        # Find the mean of coverage metrics over all POIs
        level2_df_mean = level1_df.mean().to_frame()
        level2_df_mean = level2_df_mean.T
        # Rename the column names
        columns_list = list(level2_df_mean.columns.values)
        
        for col in columns_list:
            level2_df_mean = level2_df_mean.rename(columns={col: 'Mean of ' + col})


        # Find Standard deviation of coverage metrics over all POIs
        level2_df_sd = level1_df.std().to_frame().T
        columns_list = list(level2_df_sd.columns.values)
        for col in columns_list:
            level2_df_sd = level2_df_sd.rename(columns={col: 'SD of ' + col})

        level2_metrics = pandas.concat([level2_df_mean, level2_df_sd], axis=1)
        level2_metrics.to_csv(outFilePath,index=False, na_rep='nan')



  




