""" 
.. module:: public_library

:synopsis: *Contains interface class (of instrupy to external script) for managing space system instruments.*

.. note:: For users who plan to use the InstruPy package without custom modification of the codebase, 
          only this module :mod:`public_library` of the InstruPy package is sufficient to interface 
          external script.
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


    def calc_typ_data_metrics_over_one_access_interval(self, epoch_JDUT1, SpacecraftOrbitState, AccessInfo):    
        return self._sensor.calc_typ_data_metrics_over_one_access_interval(SpacecraftOrbitState, AccessInfo) 

    def generate_level0_data_metrics(self, POI_filepath, AccessInfo_filepath, Result_filepath):
 
        ''' Generate typical data metrics per access event, per grid-point. 
            This function iteratively calls :code:`calc_typ_data_metrics_over_one_access_interval` of the specified instrument type over all access events 
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

            .. seealso::
                * :ref:`poi_file_description`
                * :ref:`access_info_CSV_file_description`
                * :ref:`basic_sensor_csv_output`
                * :ref:`passive_optical_scanner_csv_output`

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

                obsv_metrics = self._sensor.calc_typ_data_metrics_over_one_access_interval(SpacecraftOrbitState, AccessInfo) # calculate the data metrics specific to the instrument type
                _v = dict({'Access From [JDUT1]':AccessInfo["Access From [JDUT1]"], 'Access Duration [s]':AccessInfo["Access Duration [s]"], 'POI index': access_info_df['POI'][idx]}, **obsv_metrics)
                
                if idx==0: #1st iteration
                    w.writerow(_v.keys())    

                w.writerow(_v.values())

                idx = idx + 1
    
    def generate_level1_data_metrics(self, inFilePath, outFilePath):
        
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
    

    def get_coverage_specs(self):
        
        orientation_specs = {"eulerAngle1": self._sensor.orientation.x_rot_deg,
                             "eulerAngle2": self._sensor.orientation.y_rot_deg, 
                             "eulerAngle3": self._sensor.orientation.z_rot_deg,
                             "eulerSeq1": 1,
                             "eulerSeq2": 2,
                             "eulerSeq3": 3
                             }

        fov_specs = { "geometry": self._sensor.fieldOfView._geometry,
                      "coneAnglesVector" : self._sensor.fieldOfView._coneAngleVec_deg,
                      "clockAnglesVector": self._sensor.fieldOfView._clockAngleVec_deg
                    }
                    
        """ Add the along-track fov and cross-track fov. 
            
            TODO: The following code assumes only conical or rectangular as the geometry of the input fov.

            TODO: The following code does NOT take into account possible non-zero yaw rotation 
            (i.e. rotation about the imaging/pointing axis). 
        
        """
        if(self._sensor.fieldOfView._geometry == "CONICAL"):
            fov_specs['AlongTrackFov'] = self._sensor.fieldOfView._coneAngleVec_deg[0]*2
            fov_specs['CrossTrackFov'] = self._sensor.fieldOfView._coneAngleVec_deg[0]*2
        elif(self._sensor.fieldOfView._geometry == "RECTANGULAR"):
            fov_specs['AlongTrackFov'] = self._sensor.fieldOfView.get_rectangular_fov_specs()[0]
            fov_specs['CrossTrackFov'] = self._sensor.fieldOfView.get_rectangular_fov_specs()[1]
        else:
            raise Exception("Do not know along-track, cross-track fovs for the given sensor fov geometry")


        result = {"Orientation": orientation_specs, "fieldOfView": fov_specs}

        return json.dumps(result)




  




