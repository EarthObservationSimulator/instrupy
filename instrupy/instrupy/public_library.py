""" 
.. module:: public_library

:synopsis: *Provides easy interface to the package main functionalities, independent of the instrument type.*

"""
from .util import Entity, MathUtilityFunctions
from .basic_sensor import BasicSensor
from .passive_optical_scanner import PassiveOpticalScanner
from .synthetic_aperture_radar import SyntheticApertureRadar
import sys
import pandas, csv
import json
import numpy

class Instrument(Entity): 
    """ A class to serve as single interface to different types of instruments.
     
        :ivar _sensor: Object of relevant instrument.
        :vartype _sensor: :class:`instrupy.basic_sensor.BasicSensor` (or) :class:`instrupy.passive_optical_sensor.PassiveOpticalSensor` (or) :class:`instrupy.synthetic_aperture_radar.SyntheticApertureRadar`

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

        super(Instrument,self).__init__(specs.get("@id", None), "Instrument")

    @staticmethod
    def from_dict(d): 
        """ Parses a instrument from a normalized JSON dictionary. """
        return Instrument(d)

    def calc_typ_data_metrics(self, epoch_JDUT1, SpacecraftOrbitState, TargetCoords):   
        """ Calculate the typical observation data-metrics.
        
        :param epoch_JDUT1: Epoch of the alongside data in Julian Day UT1.
        :paramtype epoch_JDUT1: float

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
        return self._sensor.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords) 


    def dshield_generate_level0_data_metrics(self, POI_filepath, SatelliteState_filepath, AccessInfo_filepath, Result_filepath): 
        ''' Generate typical data metrics iterating over all access times. The time column of the satellite-state data file and the
            access data file must be referenced to the same epoch and produced at the same time-step size.
            This function iteratively calls :code:`calc_typ_data_metrics` of the specified instrument type over all access events 
            available in the input file.

            :param POI_filepath: Filepath to CSV file containing lat/lon co-ordinates of points of interest along with index.

                                 First row contains the following header elements, and the following rows contain the data corresponding to these headers.
                                 
                                 Description of the header elements:
                            
                                 * :code:`gpi`, Index of the point-of-interest.
                                 * :code:`lat[deg]`, :code:`lon[deg]`, latitude, longitude of the point-of-interest.  

                                 .. note:: Make sure the header titles are as specified above, and the delimiters are commas.

            :paramtype POI_filepath: str

             
            :param SatelliteState_filepath: Filepath to CSV file containing data satellite states at fixed time-step.
                    First four rows convey general information. The second row conveys the Epoch in Julian Days UT1.
                    The fifth row  contains the following header elements, and the following rows contain the data corresponding to these headers.

                    Description of the header elements:
                    
                    * :code:`Time[s]`,  Time referenced to epoch.
                    * :code:`X[km]`, :code:`Y[km]` :code:`Z[km]`, cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
                    * :code:`VX[km/s]`, :code:`VY[km/s]`, :code:`VZ[km/s]`, velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.

            :paramtype SatelliteState_filepath: str

            :param AccessInfo_filepath: Filepath to CSV file containing data of access events and their time-intervals.
                               First three rows convey general information. The second row conveys the Epoch in Julian Days UT1.
                               The fifth row  contains the following header elements, and the following rows contain the data corresponding to these headers.

                               Description of the header elements:
                                
                               * :code:`Time[s]`,  The time at which access starts in seconds, referenced to epoch.
                               * :code:`GP0, GP1, GP2, ....` Columns specific to each ground-point.

            :paramtype AccessInfo_filepath: str

            :param Result_filepath: Filepath to CSV file in which the results are written. First row contains the epoch in Julian Day UT1.
                                The third row contains the column header elements. Description of the header elements:
                                
                               * :code:`observationTime[s]`,  The time at which the observation is made in seconds referenced to the epoch.
                               * :code:`gpi` indicating index of ground-point.
                               * + other header elements containing data-metrics specific to the instrument type

                               .. note:: this is an **output** file of the function
            
            :paramtype Result_filepath: str

        '''
        epoch_JDUT1 = pandas.read_csv(AccessInfo_filepath, skiprows = [0], nrows=1, header=None).astype(str) # 2nd row contains the epoch
        epoch_JDUT1 = float(epoch_JDUT1[0][0].split()[2])
                
        poi_info_df = pandas.read_csv(POI_filepath)
        poi_info_df = poi_info_df.set_index('gpi')

        # Read the access file
        access_info_df = pandas.read_csv(AccessInfo_filepath,skiprows = [0,1,2,3]) # read the access times (corresponding to the DSM in which the instrument was used)
        access_info_df = access_info_df.set_index('Time[s]')

        # Read the satellite state file
        sat_state_df = pandas.read_csv(SatelliteState_filepath,skiprows = [0,1,2,3]) 
        sat_state_df = sat_state_df.set_index('Time[s]')
        sat_state_df = sat_state_df.loc[access_info_df.index] # retain states at only those times in which there are accesses

        # copy second and third row from the original access file
        with open(AccessInfo_filepath, 'r') as f:
            next(f)
            head = [next(f) for x in [1,2]] 

        # erase any old file and create new one
        with open(Result_filepath,'w') as f:
            for r in head:
                f.write(str(r))         

        with open(Result_filepath,'a+', newline='') as f:
            w = csv.writer(f)

            # Iterate over all valid logged access events
            acc_indx = list(access_info_df[access_info_df.notnull()].stack().index) # list of valid access [time, POI]

            idx = 0
            for indx in acc_indx:

                time = float(indx[0])
                poi_indx = int(indx[1][2:])

                TargetCoords = dict()                
                TargetCoords["Lat [deg]"] = poi_info_df.loc[poi_indx]["lat[deg]"]
                TargetCoords["Lon [deg]"] = poi_info_df.loc[poi_indx]["lon[deg]"]

                SpacecraftOrbitState = dict()
                SpacecraftOrbitState["Time[JDUT1]"] = epoch_JDUT1 + time*1.0/86400.0 
                SpacecraftOrbitState["x[km]"] = sat_state_df.loc[time]["X[km]"] 
                SpacecraftOrbitState["y[km]"] = sat_state_df.loc[time]["Y[km]"] 
                SpacecraftOrbitState["z[km]"] = sat_state_df.loc[time]["Z[km]"] 
                SpacecraftOrbitState["vx[km/s]"] = sat_state_df.loc[time]["VX[km/s]"] 
                SpacecraftOrbitState["vy[km/s]"] = sat_state_df.loc[time]["VY[km/s]"] 
                SpacecraftOrbitState["vz[km/s]"] = sat_state_df.loc[time]["VZ[km/s]"] 

                obsv_metrics = self._sensor.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords) # calculate the data metrics specific to the instrument type
                _v = dict({'observationTime[s]':time, 'gpi': poi_indx}, **obsv_metrics)
                if idx==0: #1st iteration
                    w.writerow(_v.keys())    
                w.writerow(_v.values())
                idx = idx + 1
                #print(obsv_metrics)
    

    def get_coverage_specs(self):
        """ Get field-of-view (or scene-field-of-view) and the field-of-regard specifications for the instrument. Also returns
            the orientation (common for both FOV, FOR) and if the instrument is constrained to take observations at a strictly 
            side-looking (no squint) target geometry. The :code:`yaw180_flag` indicates if the instrument coverage includes
            orientation of the instrument rotated by 180 deg about the yaw axis from the nominal orientation. This is significant
            for the case of FOR specifications and purely side-looking instruments.
            
            .. todo:: update unit test

            :returns: JSON string with the coverage specifications                     

            :rtype: dict
        
        """
        # determine if instrument takes observations at a purely side looking geometry
        purely_side_look = False
        instru_type = type(self._sensor).__name__
        if(instru_type == 'SyntheticApertureRadar'):
            purely_side_look = True
        
        orientation_specs = {"eulerAngle1": self._sensor.orientation.x_rot_deg,
                             "eulerAngle2": self._sensor.orientation.y_rot_deg, 
                             "eulerAngle3": self._sensor.orientation.z_rot_deg,
                             "eulerSeq1": 1,
                             "eulerSeq2": 2,
                             "eulerSeq3": 3
                             }
        
        for_specs = {"geometry": self._sensor.fieldOfRegard._geometry,
                     "coneAnglesVector" : self._sensor.fieldOfRegard._coneAngleVec_deg,
                     "clockAnglesVector": self._sensor.fieldOfRegard._clockAngleVec_deg,
                     "AlongTrackFov": self._sensor.fieldOfRegard._AT_fov_deg,
                     "CrossTrackFov": self._sensor.fieldOfRegard._CT_fov_deg,
                     "yaw180_flag": self._sensor.fieldOfRegard._yaw180_flag
                    }

        if(hasattr(self._sensor, 'sceneFieldOfView')):
            if(self._sensor.sceneFieldOfView):
                # if there exists a scene FOV, pass the corresponding fov specifications for coverage caclulations
                fov_specs = {"geometry": self._sensor.sceneFieldOfView._geometry,
                            "coneAnglesVector" : self._sensor.sceneFieldOfView._coneAngleVec_deg,
                            "clockAnglesVector": self._sensor.sceneFieldOfView._clockAngleVec_deg,
                            "AlongTrackFov": self._sensor.sceneFieldOfView._AT_fov_deg,
                            "CrossTrackFov": self._sensor.sceneFieldOfView._CT_fov_deg,
                            "yaw180_flag": self._sensor.sceneFieldOfView._yaw180_flag
                            }
                result = {"Orientation": orientation_specs, "fieldOfView": fov_specs, "purely_side_look": purely_side_look}
                return json.dumps(result)   
               
        fov_specs = {"geometry": self._sensor.fieldOfView._geometry,
                     "coneAnglesVector" : self._sensor.fieldOfView._coneAngleVec_deg,
                     "clockAnglesVector": self._sensor.fieldOfView._clockAngleVec_deg,
                     "AlongTrackFov": self._sensor.fieldOfView._AT_fov_deg,
                     "CrossTrackFov": self._sensor.fieldOfView._CT_fov_deg,
                     "yaw180_flag": self._sensor.fieldOfView._yaw180_flag
                     }


        result = {"Orientation": orientation_specs, "fieldOfView": fov_specs, "purely_side_look": purely_side_look, "fieldOfRegard": for_specs}
        return json.dumps(result)

################################################### Legacy functions #################################################################    
    def generate_level0_data_metrics(self, POI_filepath, AccessInfo_filepath, Result_filepath):
 
        ''' ########################## Legacy function, discontinued ########################## 
        
            Generate typical data metrics per access event, per grid-point. 
            This function iteratively calls :code:`calc_typ_data_metrics` of the specified instrument type over all access events 
            available in the input file.
            CSV formatted files are supported to supply the required input data.

            :param POI_filepath: Filepath to CSV file containing lat/lon co-ordinates of points of interest along with index.

                                 First row contains the following header elements, and the following rows contain the data corresponding to these headers.
                                 
                                 Description of the header elements:
                            
                                 * :code:`POI`, Index of the point-of-interest.
                                 * :code:`Lat [deg]`, :code:`Lon [deg]`, latitude, longitude of the point-of-interest.  

                                 .. note:: Make sure the header titles are as specified above, and the delimiters are commas.
            :paramtype POI_filepath: str

             
            :param AccessInfo_filepath: Filepath to CSV file containing data of access events and their time-intervals.
                               First three rows convey general information. The fourth row conveys the Epoch in Julian Days UT1.
                               The fifth row  contains the following header elements, and the following rows contain the data corresponding to these headers.

                               Description of the header elements:
                                
                               * :code:`AccessFrom[Days]`,  The time at which access starts in [days], referenced to epoch in row-4.
                               * :code:`Duration[s]`, Duration of access in [s].
                               * :code:`POI` indicating index of ground-point.
                               * :code:`EventNum` indicating index of event.
                               * :code:`Time[Days]`, Time in [Days] at which the alongside satellite-state is recorded. Referenced to the epoch specified in row-4.
                               * :code:`X[km]`, :code:`Y[km]` :code:`Z[km]`, cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
                               * :code:`VX[km/s]`, :code:`VY[km/s]`, :code:`VZ[km/s]`, velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.


            :paramtype AccessInfo_filepath: str

            :param Result_filepath: Filepath to CSV file in which the results are written
                                Description of the header elements:
                                
                               * :code:`POI` indicating index of ground-point.
                               * :code:`Access From [JDUT1]`,  The time at which access starts in [JDUT1]
                               * :code:`Duration [s]`, Duration of access in [s].
                               * + other header elements specific to the instrument type

                               .. note:: this is an **output** file of the function
            
            :paramtype Result_filepath: str

        '''
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
        """ ########################## Legacy function, discontinued ########################## """
        
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
        """ ########################## Legacy function, discontinued ########################## """

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
        """ ########################## Legacy function, discontinued ########################## """
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
        """ ########################## Legacy function, discontinued ########################## """

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



  




