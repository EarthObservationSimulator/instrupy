Output description
******************

Overview
=========================

InstruPy produces three *types* of output files. 

.. csv-table:: Input parameter description 
   :header: Type, Description
   :widths: 10,40

   Level-0 , Contains observation data-metrics  per access-event (recorded by OC) per grid-point.                                            
   Level-1 , Contains statistics of observation data-metrics time  taken over all  access-events over entire mission duration per grid point.
   Level-2 , Contains statistics of observation data-metrics taken over all  mission duration and *then* over all grid-points.                

The level-1 and level-2 data metrics are post-processed from level-0 data metrics.

**Level-0** data metrics
---------------------------

The Level-0 metrics are generated for a hypothetical pixel centered around the ground-point (at which the access-event is calculated). The level-0 metrics depend on the type of instrument used in the Mission. Please refer to the following links for level-0 metric description based on the instrument,
   
    * :ref:`basic_sensor_csv_output`, 
    * :ref:`passive_optical_scanner_csv_output` and 
    * :ref:`synthetic_aperture_radar_csv_output`.


**Level-1** data metrics 
---------------------------

These are produced from the Level-0 data-metrics. There are two kinds of level-1 metrics:

    1. Mean of corresponding level-0 data metric

       :math:`\bar{m}_{p_j}^{l_1} =\dfrac{1}{M}\Sigma_{i=1}^M m_{i,p_j}^{l_0}`

    2. Standard Deviation of corresponding level-0 data metric
    
       :math:`\grave{m}_{p_j}^{l_1}=\sqrt{ \dfrac{1}{M}\Sigma_{i=1}^M (m_{i,p_j}^{l_0} - \bar{m}_{p_j}^{l_1})^2}`

where,

:math:`m_{i,p_j}^{l_0}` is the level-0 metric corresponding to ***valid*** observation :math:`i` at ground-point :math:`p_j`,
The level-0 metric at the corresponding access-event is considered **only if** the corresponding `Coverage [T/F]` flag is true. 

:math:`\bar{m}_{p_j}^{l_1}` is the metric averaged over all ***valid*** observations :math:`1` to :math:`M` at ground-point :math:`p_j`.

:math:`\grave{m}_{p_j}^{l_1}` is the standard-deviation metric over all ***valid*** observations :math:`1` to :math:`M` at ground-point :math:`p_j`.

Note that in general :math:`M \leq Total \hspace{2mm} number \hspace{2mm} of \hspace{2mm} Access \hspace{2mm} Events \hspace{2mm} at \hspace{2mm} p_j`

**Level-2** data metrics 
---------------------------

These are produced from the *Level-1* data-metrics. There are three kinds of level-2 data metrics:

    1. Mean of Mean of corresponding *level-0* data metric

       :math:`\bar{m}^{l_2}=\dfrac{1}{N}\Sigma_{j=1}^N \bar{m}_{p_j}^{l_1}`

    2. Mean of Standard deviation of corresponding *level-0* data metric

       :math:`\widetilde{m}^{l_2}=\dfrac{1}{N}\Sigma_{j=1}^N \grave{m}_{p_j}^{l_1}`

    3. Standard deviation of mean of corresponding *level-0* data metric

       :math:`\grave{m}^{l_2}=\sqrt{\dfrac{1}{N}\Sigma_{i=1}^N (\bar{m}_{p_j}^{l_1} - \bar{m}^{l_2})^2}`

where,

:math:`\bar{m}_{p_j}^{l_1}` is the level-1 mean metric corresponding to ground-point :math:`p_j`,

:math:`\grave{m}_{p_j}^{l_1}` is the level-1 standard-deviation metric corresponding to ground-point :math:`p_j`,

:math:`\bar{m}^{l_2}` is the level-2 mean of mean metric over all observations and ground-points. 

:math:`\widetilde{m}^{l_2}` is the level-2 mean of standard-deviation metric over all observations and ground-points.

:math:`\grave{m}^{l_2}` is the level-2 standard-deviation of mean metric over all observations and ground-points. 

Only Groundpoints :math:`p_1` to :math:`p_N` with atleast one valid observation is considered. 

Note that in general :math:`N \leq Total \hspace{2mm} number \hspace{2mm} of \hspace{2mm} Gridpoints`

Output file-naming convention
==============================

Per architecture there are three CSV files with following names produced:

1. `satID_level0_data_metrics.csv`

A separate file is written corresponding to each satellite in the architecture. The `satID` term corresponds to the 
satellite ID as specified in the respective access info file from which the metrics are generated.

The data is written row-wise, with each row corresponding to a unique access-event. The column names are:

    i. Basic Sensor type instrument

       :code:`POI index`, :code:`Access From [JDUT1]`, :code:`Access Duration [s]`, :code:`Coverage [T/F]`, :code:`Incidence angle [deg]`, :code:`Look angle [deg]`, :code:`Observation Range [km]`, :code:`Solar Zenith [deg]`.
    
    ii. Passive Optical Scanner type instrument

       :code:`POI index`, :code:`Access From [JDUT1]`, :code:`Access Duration [s]`, :code:`Coverage [T/F]`, :code:`Noise-Equivalent delta T [K]`, :code:`DR`, SNR, :code:`Ground Pixel Along-Track Resolution [m]`, :code:`Ground Pixel Cross-Track Resolution [m]`.

    iii. Synthetic Aperture Radar type instrument

        :code:`POI index`, :code:`Access From [JDUT1]`, :code:`Access Duration [s]`, :code:`Coverage [T/F]`, :code:`Incidence Angle [deg]`, :code:`Swath-Width [m]`, :code:`Sigma NEZ Nought [dB]`, :code:`Ground Pixel Along-Track Resolution [m]`, :code:`Ground Pixel Cross-Track Resolution [m]`.

2. `level1_data_metrics.csv`

The data is written row-wise, with each row corresponding to a unique ground-point. The column names are:

    i. Basic Sensor type instrument

       :code:`POI index`, 
       :code:`Mean of Incidence angle [deg]`, :code:`SD of Incidence angle [deg]`, 
       :code:`Mean of Solar Zenith [deg]`, :code:`SD of Solar Zenith [deg]`, 
       :code:`Mean of Look angle [deg]`, :code:`SD of Look angle [deg]`, 
       :code:`Mean of Observation Range [km]`, :code:`SD of Observation Range [km]`.
    
    ii. Passive Optical Scanner type instrument

        :code:`POI index`,
        :code:`Mean of Noise-Equivalent delta T [K]`, :code:`SD of Noise-Equivalent delta T [K]`, 
        :code:`Mean of DR`, :code:`SD of DR`, 
        :code:`Mean of SNR`, :code:`SD of SNR`, 
        :code:`Mean of Ground Pixel Along-Track Resolution [m]`, :code:`SD of Ground Pixel Along-Track Resolution [m]`, 
        :code:`Mean of Ground Pixel Cross-Track Resolution [m]`, :code:`SD of Ground Pixel Cross-Track Resolution [m]`.

    iii. Synthetic Aperture Radar type instrument

         :code:`POI index`, 
         :code:`Mean of Incidence Angle [deg]`, :code:`SD of Incidence Angle [deg]`, 
         :code:`Mean of Swath-Width [m]`, :code:`SD of Swath-Width [m]`, 
         :code:`Mean of Sigma NEZ Nought [dB]`, :code:`SD of Sigma NEZ Nought [dB]`, 
         :code:`Mean of Ground Pixel Along-Track Resolution [m]`, :code:`SD of Ground Pixel Along-Track Resolution [m]`, 
         :code:`Mean of Ground Pixel Cross-Track Resolution [m]`, :code:`SD of Ground Pixel Cross-Track Resolution [m]`.

3. `level2_data_metrics.csv`

The data is one row only. The column names are:

    i. Basic Sensor type instrument

       :code:`Mean of Mean of Incidence angle [deg]`, :code:`SD of Mean of Incidence angle [deg]`, :code:`Mean of SD of Incidence angle [deg]`, :code:`SD of SD of Incidence angle [deg]`, 
       :code:`Mean of Mean of Solar Zenith [deg]`, :code:`SD of Mean of Solar Zenith [deg]`, :code:`Mean of SD of Solar Zenith [deg]`, :code:`SD of SD of Solar Zenith [deg]`, 
       :code:`Mean of Mean of Look angle [deg]`, :code:`SD of Mean of Look angle [deg]`, :code:`Mean of SD of Look angle [deg]`, :code:`SD of SD of Look angle [deg]`, 
       :code:`Mean of Mean of Observation Range [km]`, :code:`SD of Mean of Observation Range [km]`, :code:`Mean of SD of Observation Range [km]`, :code:`SD of SD of Observation Range [km]`
    
    ii. Passive Optical Scanner type instrument

        :code:`Mean of Mean of Noise-Equivalent delta T [K]`, :code:`SD of Mean of Noise-Equivalent delta T [K]`, :code:`Mean of SD of Noise-Equivalent delta T [K]`, :code:`SD of SD of Noise-Equivalent delta T [K]`,    
        :code:`Mean of Mean of DR`, :code:`SD of Mean of DR`, :code:`Mean of SD of DR`, :code:`SD of SD of DR`, 
        :code:`Mean of Mean of SNR`, :code:`SD of Mean of SNR`, :code:`Mean of SD of SNR`, :code:`SD of SD of SNR`,
        :code:`Mean of Mean of Ground Pixel Along-Track Resolution [m]`, :code:`SD of Mean of Ground Pixel Along-Track Resolution [m]`, :code:`Mean of SD of Ground Pixel Along-Track Resolution [m]`, :code:`SD of SD of Ground Pixel Along-Track Resolution [m]`, 
        :code:`Mean of Mean of Ground Pixel Cross-Track Resolution [m]`,  :code:`SD of Mean of Ground Pixel Cross-Track Resolution [m]`,  :code:`Mean of SD of Ground Pixel Cross-Track Resolution [m]`, :code:`SD of SD of Ground Pixel Cross-Track Resolution [m]`

    iii. Synthetic Aperture Radar type instrument

         :code:`Mean of Mean of Incidence Angle [deg]`, :code:`SD of Mean of Incidence Angle [deg]`, :code:`Mean of SD of Incidence Angle [deg]`, :code:`SD of SD of Incidence Angle [deg]`,
         :code:`Mean of Mean of Swath-Width [m]`, :code:`SD of Mean of Swath-Width [m]`, :code:`Mean of SD of Swath-Width [m]`, :code:`SD of SD of Swath-Width [m]`, 
         :code:`Mean of Mean of Sigma NEZ Nought [dB]`, :code:`SD of Mean of Sigma NEZ Nought [dB]`, :code:`Mean of SD of Sigma NEZ Nought [dB]`,  :code:`SD of SD of Sigma NEZ Nought [dB]`,
         :code:`Mean of Mean of Ground Pixel Along-Track Resolution [m]`, :code:`SD of Mean of Ground Pixel Along-Track Resolution [m]`, :code:`Mean of SD Ground Pixel Along-Track Resolution [m]`, :code:`SD of SD Ground Pixel Along-Track Resolution [m]`, 
         :code:`Mean of Mean of Ground Pixel Cross-Track Resolution [m]`,  :code:`SD of Mean of Ground Pixel Cross-Track Resolution [m]`,  :code:`Mean of SD Ground Pixel Cross-Track Resolution [m]`,  :code:`SD of SD Ground Pixel Cross-Track Resolution [m]`

***Note:*** While reading the CSV files, please refer to the columns by their column names and ***not*** the column numbers.

Coverage metrics
=========================

The instrument module also processes the `Coverage [T/F]` level-0 data-output (of any type of instrument) and calculates coverage metrics at each point-of-interest (level-1) and global-region metrics (level-2). Please refer to *SMAD, 3rd edition Section 7.2.3* for description of the level-1 coverage metrics. The level-2 coverage metrics are the mean, and standard-deviation of the level-1 metrics taken over all the POIs. 

Per architecture there are two CSV files with following names produced:

1. `level1_coverage_metrics.csv`

The data is written row-wise, with each row corresponding to a unique point-of-interest. The column names are:

    :code:`POI index`, :code:`Gaps [days]`, :code:`Max Coverage Gap [Days]`, :code:`Mean Coverage Gap [Days]`, :code:`Time Average Gap [Days]`, :code:`Percentage Coverage`, 

2. `level2_coverage_metrics.csv`

The data is one row only. The column names are:

    :code:`Mean of Max Coverage Gap [Days]`, :code:`Mean of Mean Coverage Gap [Days]`, :code:`Mean of Time Average Gap [Days]`, :code:`Mean of Percentage Coverage`,
    :code:`SD of Max Coverage Gap [Days]`, :code:`SD of Mean Coverage Gap [Days]`, :code:`SD of Time Average Gap [Days]`, :code:`SD of Percentage Coverage`