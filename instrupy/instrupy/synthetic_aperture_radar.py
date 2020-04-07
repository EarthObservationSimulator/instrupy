""" 
.. module:: synthetic_aperture_radar

:synopsis: *Module to handle synthetic_aperture_radar type of instrument.*

        Inline code comments contain references to the following articles:

        1. Performance Limits for Synthetic Aperture Radar - second edition SANDIA Report 2006. ----> Main reference.
        2. Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993. ----> Reference for PRF validity calculations, corrections for spaceborne radar.

        
"""

import json
import numpy
import copy
import pandas, csv
from .util import Entity, Orientation, FieldOfView, MathUtilityFunctions, Constants, FileUtilityFunctions

class SyntheticApertureRadar(Entity):
    """A synthetic aperture radar class estimating observation data-metrics from Stripmap mode of imaging.       
      
        A full-swath imaging mode is considered, meaning the swath-width is the area illuminated by the 3-dB beamwidth of the antenna.

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

        :ivar fieldOfView: Field of view specification of instrument. In case of SAR, this is calculated from the user supplied antenna dimensions.
        :vartype fieldOfView: :class:`instrupy.util.FieldOfView` 
     
        :ivar sceneFieldOfView: Field of view corresponding to a "scene" captured by the SAR. A scene is made of multiple concatenated strips.
        :vartype sceneFieldOfView: :class:`instrupy.util.FieldOfView` 

        :ivar dataRate: Rate of data recorded (Mbps) during nominal operations.
        :vartype dataRate: float  

        :ivar pulseWidth: Actual pulse width in (s).
        :vartype pulseWidth: float

        :ivar antennaAlongTrackDim: Antenna size in the along-track direction in (m).
        :vartype antennaAlongTrackDim: float

        :ivar antennaCrossTrackDim: Antenna size in the cross-track direction in (m).
        :vartype antennaCrossTrackDim: float

        :ivar antennaApertureEfficiency: Aperture efficiency of antenna (:math:`0 < \\eta_{ap} < 1`).
        :vartype antennaApertureEfficiency: float

        :ivar operatingFrequency: Operating radar center frequency in (Hz).
        :vartype operatingFrequency: float

        :ivar peakTransmitPower: Peak transmit power in (W).
        :vartype peakTransmitPower: float

        :ivar chirpBandwidth: Bandwidth of radar operation in (Hz).
        :vartype chirpBandwidth: float

        :ivar minimumPRF: The minimum pulse-repetition-frequency of operation in (Hz).
        :vartype minimumPRF: float

        :ivar maximumPRF: The maximum pulse-repetition-frequency of operation in (Hz).
        :vartype maximumPRF: float

        :ivar sceneNoiseTemp: Nominal scene noise temperature in (K).
        :vartype sceneNoiseTemp: float

        :ivar systemNoiseFigure:  (dB) System noise figure for the receiver. The system noise figure includes primarily the noise figure of the front-end Low-Noise
                                  Amplifier (LNA) and the losses between the antenna and the LNA. Typical system noise figures for sub-kilowatt radar systems are 3.0 dB to 3.5 dB
                                  at X-band, 3.5 dB to 4.5 dB at Ku-band, and perhaps 6 dB at Ka-band.
        :vartype systemNoiseFigure: float

        :ivar radarLosses: (dB) These include a variety of losses primarily over the microwave signal path, but doesn't include the atmosphere. Included are a power loss from transmitter power amplifier
                           output to the antenna port, and a two-way loss through the radome. Typical numbers might be 0.5 dB to 2 dB from TX amplifier to the
                           antenna port, and perhaps an additional 0.5 dB to 1.5 dB two-way through the radome.
        :vartype radarLosses: float

        :ivar sigmaNEZ0threshold: The :math:`\\sigma_{NEZ0}` threshold for classification as a valid observation.
        :vartype sigmaNEZ0threshold: float

        :cvar L_r: Reduction in SNR gain due to non-ideal range filtering.
        :vartype L_r: float

        :cvar L_a: Reduction in SNR gain due to non-ideal azimuth filtering.
        :vartype L_a: float

        :cvar a_wa:  Azimuth impulse response broadening factor.
        :vartype a_wa: float

        :cvar a_wr: Range impulse response broadening factor.
        :vartype a_wr: float

        :cvar L_atmos_dB: 2-way atmospheric loss of electromagnetic energy.
        :vartype L_atmos_dB: float

        .. note:: The actual pulse-repetition frequency is taken as the highest PRF within allowed range of PRFs. The highest PRF is chosen since it allows greater :math:`\\sigma_{NEZ0}`
       
   
    """
    L_r = float(1.2)
    L_a = float(1.2)  
    a_wa = float(1.2)  
    a_wr = float(1.2)
    L_atmos_dB = float(2)
    

    def __init__(self, name=None, acronym=None, mass=None,volume=None, power=None,  orientation=None, fieldOfView = None,
            sceneFieldOfView = None, dataRate=None, bitsPerPixel = None, pulseWidth = None, antennaAlongTrackDim= None, 
            antennaCrossTrackDim = None, antennaApertureEfficiency = None, operatingFrequency = None, 
            peakTransmitPower = None, chirpBandwidth = None, minimumPRF = None, maximumPRF = None, 
            radarLosses = None, sceneNoiseTemp = None, systemNoiseFigure = None, sigmaNEZ0threshold = None, _id=None):
        """Initialize a Synthetic Aperture Radar object.

        """          
        self.name = str(name) if name is not None else None
        self.acronym = str(acronym) if acronym is not None else self.name
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None else None
        self.fieldOfView = copy.deepcopy(fieldOfView) if fieldOfView is not None else None
        self.sceneFieldOfView = copy.deepcopy(sceneFieldOfView) if sceneFieldOfView is not None else None
        self.dataRate = float(dataRate) if dataRate is not None else None          
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None 
        self.pulseWidth = float(pulseWidth) if pulseWidth is not None else None
        self.antennaAlongTrackDim = float(antennaAlongTrackDim) if antennaAlongTrackDim is not None else None
        self.antennaCrossTrackDim = float(antennaCrossTrackDim) if antennaCrossTrackDim is not None else None
        self.antennaApertureEfficiency = float(antennaApertureEfficiency) if antennaApertureEfficiency is not None else None
        self.operatingFrequency = float(operatingFrequency) if operatingFrequency is not None else None
        self.peakTransmitPower = float(peakTransmitPower) if peakTransmitPower is not None else None
        self.chirpBandwidth = float(chirpBandwidth) if chirpBandwidth is not None else None        
        self.minimumPRF = float(minimumPRF) if minimumPRF is not None else None   
        self.maximumPRF = float(maximumPRF) if maximumPRF is not None else None 
        self.radarLosses = float(radarLosses) if radarLosses is not None else None
        self.sceneNoiseTemp = float(sceneNoiseTemp) if sceneNoiseTemp is not None else None  
        self.systemNoiseFigure = float(systemNoiseFigure) if systemNoiseFigure is not None else None 
        self.sigmaNEZ0threshold = float(sigmaNEZ0threshold) if sigmaNEZ0threshold is not None else None 
        super(SyntheticApertureRadar,self).__init__(_id, "Synthetic Aperture Radar")
        
    @staticmethod
    def from_dict(d):

        # Only side-looking orientation fo instrument suported for synthetic aperture radar stripmap imaging
        orien_json_str = d.get("orientation", None)
        if(FileUtilityFunctions.from_json(orien_json_str).get("convention",None) != "SIDE_LOOK"):
            raise Exception("Only side-looking orientation of instrument supported for the synthetic aperture radar imaging.")
            
        # check supplied PRF range
        _PRFmin = d.get("minimumPRF", None)
        _PRFmax = d.get("maximumPRF", None)
        if(_PRFmin > _PRFmax):
            raise Exception("PRF minimum must be less than or equal to PRF maximum.")

        # calculate instrument FOV based on antenna dimensions
        D_az_m = d.get("antennaAlongTrackDim", None)
        D_elv_m = d.get("antennaCrossTrackDim", None)
        opWavelength =  Constants.speedOfLight/ d.get("operatingFrequency", None)
        # calculate antenna beamwidth and hence fovs [1] Eqn 41.
        along_track_fov_deg = numpy.rad2deg(opWavelength/ D_az_m)
        cross_track_fov_deg = numpy.rad2deg(opWavelength/ D_elv_m)
        fov_json_str = '{ "sensorGeometry": "RECTANGULAR", "alongTrackFieldOfView":' + str(along_track_fov_deg)+ ',"crossTrackFieldOfView":' + str(cross_track_fov_deg) + '}' 
        
        # initialize "Scene FOV" if required        
        sceneLength2ALtRatio = d.get("sceneLength2AltRatio", None)
        if(sceneLength2ALtRatio):
            sc_AT_fov_deg = numpy.rad2deg(numpy.arctan(sceneLength2ALtRatio)) # approximate along_track_fov_deg
            sc_CT_fov_deg = cross_track_fov_deg
            sc_fov_json_str = '{ "sensorGeometry": "RECTANGULAR", "alongTrackFieldOfView":' + str(sc_AT_fov_deg)+ ',"crossTrackFieldOfView":' + str(sc_CT_fov_deg) + '}' 
            scene_fov = FieldOfView.from_json(sc_fov_json_str)
        else:
            scene_fov = None
 
        return SyntheticApertureRadar(
                        name = d.get("name", None),
                        acronym = d.get("acronym", None),
                        mass = d.get("mass", None),
                        volume = d.get("volume", None),
                        power = d.get("power", None),
                        orientation = Orientation.from_json(d.get("orientation", None)),
                        fieldOfView = FieldOfView.from_json(fov_json_str),
                        sceneFieldOfView = scene_fov,
                        dataRate = d.get("dataRate", None),
                        bitsPerPixel = d.get("bitsPerPixel", None),
                        pulseWidth = d.get("pulseWidth", None),
                        antennaAlongTrackDim = d.get("antennaAlongTrackDim", None),
                        antennaCrossTrackDim = d.get("antennaCrossTrackDim", None),
                        antennaApertureEfficiency = d.get("antennaApertureEfficiency", None),
                        operatingFrequency = d.get("operatingFrequency", None),
                        peakTransmitPower = d.get("peakTransmitPower", None),
                        chirpBandwidth = d.get("chirpBandwidth", None),
                        minimumPRF = d.get("minimumPRF", None),
                        maximumPRF = d.get("maximumPRF", None),
                        radarLosses = d.get("radarLosses", None),
                        sceneNoiseTemp = d.get("sceneNoiseTemp", None),
                        systemNoiseFigure = d.get("systemNoiseFigure", None),
                        sigmaNEZ0threshold = d.get("sigmaNEZ0threshold", None),
                        )

    @staticmethod
    def get_azimuthal_resolution(v_sc_kmps, v_g_kmps, D_az):
        """ Calculate azimuthal resolution taking into consideration difference in platform, footprint velocities.
        See eqn (5.3.6.3) in [2]
        """
        return (D_az/2.0)*(v_g_kmps/ v_sc_kmps)


    def calc_typ_data_metrics_over_one_access_interval(self, SpacecraftOrbitState, AccessInfo):
        ''' Calculate typical observation data metrics during the given access period.

            :param SpacecraftOrbitState: Spacecraft position at the middle (or as close as possible to the middle) of the access interval.

                               Dictionary keys are: 
                               
                               * :code:`Time[JDUT1]` (:class:`float`), Time in Julian Day UT1.
                               * :code:`x[km]` (:class:`float`), :code:`y[km]` (:class:`float`), :code:`z[km]` (:class:`float`), cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
                               * :code:`vx[km/s]` (:class:`float`), :code:`vy[km/s]` (:class:`float`), :code:`vz[km/s]` (:class:`float`), velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.
            :paramtype SpacecraftOrbitState: dict

            
            :param AccessInfo: Access information.

                               Dictionary keys are: 
                                
                               * :code:`Access From [JDUT1]` (:class:`float`) Start absolute time of Access in Julian Day UT1.
                               * :code:`Duration [s]` (:class:`float`): Access duration in [s]. Ignored.
                               * :code:`Lat [deg]` (:class:`float`), :code:`Lon [deg]` (:class:`float`), indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
            :paramtype AccessInfo: dict

            :returns: Typical calculated observation data metrics. It is assumed that the satellite takes *one* observation data sample per access event. Below metrics are 
                      calculated using satellite position data at or near the middle of the access event.
                      Dictionary keys are: 
                    
                      * :code:`Coverage [T/F]` (:class:`bool`) indicating if observation was possible during the access event.
                      * :code:`Sigma NEZ Nought [dB]` (:class:`float`)  The backscatter coefficient :math:`\\sigma_0` of a target for which the signal power level in final
                        image is equal to the noise power level. 
                      * :code:`Ground Pixel Along-Track Resolution [m]` (:class:`float`) Resolution on an hypothetical ground-pixel centered about observation point
                      * :code:`Ground Pixel Cross-Track Resolution [m]` (:class:`float`) Resolution on an hypothetical ground-pixel centered about observation point
                      * :code:`Swath-Width [m]` (:class:`float`) Swath-width of the strip of which the imaged pixel is part off.
                      * :code:`Incidence Angle [deg]` (:class:`float`) Observation incidence angle at the ground-pixel.

            :rtype: dict

            .. note:: It is assumed that the instrument captures **one** *observation/ image* over the entire access interval. 
                      While this is true for stripmap radar imagers, push-broom imagers and whisk-broom imagers, it is not true for all
                      different types of instruments.

            .. note:: We differentiate between **access** and **coverage**. **Access** is when the target location
                      falls under the sensor FOV. **Coverage** is when the target location falls under sensor FOV *and* 
                      can be observed.

            .. todo:: Include frequency dependent atmospheric losses in :math:`\\sigma_{NEZ0}` calculations.
                        
        '''
        
        # Observation time in Julian Day UT1
        tObs_JDUT1 = SpacecraftOrbitState["Time[JDUT1]"]

        # Calculate Target position in ECI frame
        TargetPosition_km = MathUtilityFunctions.geo2eci([AccessInfo["Lat [deg]"], AccessInfo["Lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position in Cartesian coordinates
        SpacecraftPosition_km = numpy.array([SpacecraftOrbitState["x[km]"], SpacecraftOrbitState["y[km]"], SpacecraftOrbitState["z[km]"]])  
        SpacecraftVelocity_kmps = numpy.array([SpacecraftOrbitState["vx[km/s]"], SpacecraftOrbitState["vy[km/s]"], SpacecraftOrbitState["vz[km/s]"]]) 
        v_sc_kmps = numpy.linalg.norm(SpacecraftVelocity_kmps)
        
         #  Calculate range vector between spacecraft and POI (Target)
        range_vector_km = TargetPosition_km - SpacecraftPosition_km

        alt_km = numpy.linalg.norm(SpacecraftPosition_km) - Constants.radiusOfEarthInKM
        look_angle = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(SpacecraftPosition_km)))
        incidence_angle_rad = numpy.arcsin(numpy.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)       

        # Copying values into variables of more code-friendly variables
        Re = Constants.radiusOfEarthInKM * 1e3        
        c = Constants.speedOfLight
        k = Constants.Boltzmann
        tau_p = self.pulseWidth        
        B_T = self.chirpBandwidth
        P_T = self.peakTransmitPower
        D_az = self.antennaAlongTrackDim
        D_elv = self.antennaCrossTrackDim
        fc = self.operatingFrequency
        eta_ap = self.antennaApertureEfficiency               
        PRFmin_Hz = self.minimumPRF    
        PRFmax_Hz = self.maximumPRF   
        L_radar = 10.0**(self.radarLosses/10.0) # convert to linear units
        F_N = 10.0**(self.systemNoiseFigure/10.0) # convert to linear units
        L_atmos = 10.0**(SyntheticApertureRadar.L_atmos_dB/10.0) # convert to linear units
        L_r = SyntheticApertureRadar.L_r
        L_a = SyntheticApertureRadar.L_a
        a_wr = SyntheticApertureRadar.a_wr
        a_wa = SyntheticApertureRadar.a_wa
        T = self.sceneNoiseTemp
              
        v_g_kmps = 1e-3 * MathUtilityFunctions.compute_satellite_footprint_speed(SpacecraftPosition_km*1e3, SpacecraftVelocity_kmps*1e3) # This is approximation, since the image footprint velocity is not necessarily equal to the
                                                # satellite footprint speed. However it is reasonable approximation in case of low-altitudes and small look angles. TBD: Improve the model.

                                  
        instru_look_angle_rad = numpy.abs(numpy.deg2rad(self.orientation.y_rot_deg))

        f_P = SyntheticApertureRadar.find_valid_highest_possible_PRF(PRFmin_Hz, PRFmax_Hz, v_sc_kmps, v_g_kmps, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)
      
        isCovered = False
        rho_a = None
        rho_y = None
        sigma_N_dB = None
        W_gr = None
        theta_i = None

        if (f_P is not None): # Observation is possible at PRF = f_P        
            
            range_km = numpy.linalg.norm(range_vector_km)
            R = range_km*1e3


            h = alt_km * 1e3
            Rs = Re + h

            lamb = c/fc


            # full swath imaging _mode meaning all the swath. My own formulation, not in reference.
            # dimension of antenna along range direction
            # NOTE: While calculating swath width the instrument nominal look angle must be used!!!! Not the Target look angle.            
            theta_elv = lamb/ D_elv # half-power beamwidth illuminating the swath
            gamma_n = instru_look_angle_rad - 0.5*theta_elv
            gamma_f = instru_look_angle_rad + 0.5*theta_elv
            theta_in = numpy.arcsin(numpy.sin(gamma_n)*Rs/Re)
            theta_if = numpy.arcsin(numpy.sin(gamma_f)*Rs/Re)
            alpha_n = theta_in - gamma_n
            alpha_f = theta_if - gamma_f
            alpha_s = alpha_f - alpha_n
            W_gr = Re*alpha_s   

            # look angle to target ground pixel [2] equation 5.1.3.4
            # Note that this is not the same as incidence angle to middle of swath

            theta_i = incidence_angle_rad

                            
            psi_g = numpy.pi/2.0 - theta_i # grazing angle                    
            
            # [1] equation 17, find P_avg, average transmit power
            T_eff = tau_p # approximate effective pulse duration to be actual pulse duration, as in case of matched filter processing
            d = T_eff * f_P # [1] equation 17
            P_avg = d*P_T
            
            # [1] equation 8, find G_A
            A_A = D_elv*D_az
            G_A = 4.0*numpy.pi*eta_ap*A_A/lamb**2  
            
            # [1] equation 37 we can get the sigma_N. Note that the spacecraft speed is used in the equation and not the ground-speed, see [2] for explanation.              
            sigma_N = (265.0*numpy.pi**3*k*T / c)*(R**3 * (v_sc_kmps*1e3) * numpy.cos(psi_g))*(B_T*F_N*L_radar*L_atmos/ (P_avg*G_A**2*lamb**3))*(L_r*L_a/(a_wr*a_wa))
            sigma_N_dB = 10.0*numpy.log10(sigma_N)

            # [1] equations 36, 23 we can get rho_y
            rho_y = a_wr*c/(2*B_T*numpy.cos(psi_g))

            # [1] equation 69 we get minimum possible azimuth resolution (for strip mapping)
            rho_a = SyntheticApertureRadar.get_azimuthal_resolution(v_sc_kmps, v_g_kmps, D_az)

            # check if sigma NEZ nought satisfies uer supplied threshold. Note that lesser is better.
            if(sigma_N_dB < self.sigmaNEZ0threshold):
                isCovered = True

               
        obsv_metrics = {}
        obsv_metrics["Ground Pixel Along-Track Resolution [m]"] = rho_a
        obsv_metrics["Ground Pixel Cross-Track Resolution [m]"] = rho_y
        obsv_metrics["Sigma NEZ Nought [dB]"] = sigma_N_dB
        obsv_metrics["Incidence angle [deg]"] = numpy.rad2deg(theta_i) if theta_i is not None else numpy.nan
        obsv_metrics["Swath-width [km]"] = W_gr/1e3 if W_gr is not None else numpy.nan        
        obsv_metrics["Coverage [T/F]"] = isCovered

        return obsv_metrics
                

    @staticmethod
    def find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc_kmps, v_g_kmps, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc):
        """ Function to find the highest possible pulse repetition frequency within the user supplied range of PRFs, which 
            shall allow observation of target. Not all PRFs are valid and a valid PRF has to be chosen so that it meets all
            the below conditions:

            1. The length of the echo from 3-dB antenna beam illuminated swath is less than inter-pulse period.
            2. The PRFs are high enough to allow for unambiguous detection of doppler shifts.
            3. The echos from target doesn't overlap with a transmit pulse (in the future).
            4. The echo from Nadir (or a previous tranmist pulse) doesn't overlap with the desired echo.

        [2] is the primary reference for this formulation, although some errors have been found (and corrected for the current
        implementation) in the text.

        Of all the available valid PRFs, the highest PRF is chosen since it improves the :math:`\sigma_{NEZ0}` observation data-metric. 

        .. warning:: Be careful to use the instrument look angle and not the target pixel look angle while calculating the swath-width
                     and incidence angle to middle of the swath.

        """
        h = alt_km * 1e3
        Re = Constants.radiusOfEarthInKM * 1e3         
        c = Constants.speedOfLight
        Rs = Re + h # [2]  part of equation 5.1.2.3

        lamb = c/fc        

        # full swath imaging _mode meaning all the swath. My own formulation, not in reference.
        # dimension of antenna along range direction
        # NOTE: While calculating full swath width the instrument look angle must be used!!!! Not the Target look angle.
        theta_elv = lamb/ D_elv # 3-dB beamwidth (full-beamwidth) illuminating the swath
        gamma_n = instru_look_angle_rad - 0.5*theta_elv
        gamma_f = instru_look_angle_rad + 0.5*theta_elv
        theta_in = numpy.arcsin(numpy.sin(gamma_n)*Rs/Re)
        theta_if = numpy.arcsin(numpy.sin(gamma_f)*Rs/Re)
        alpha_n = theta_in - gamma_n
        alpha_f = theta_if - gamma_f

    
        # Preprocessing to check if PRF allows for unambiguous echo detection from swath W_gr        
        Rn = numpy.sqrt(Re**2 + Rs**2 - 2*Re*Rs*numpy.cos(alpha_n)) # [2] equation 5.1.3.9 slant-range to near edge of swath
        Rf = numpy.sqrt(Re**2 + Rs**2 - 2*Re*Rs*numpy.cos(alpha_f)) # [2] equation 5.1.3.10 slant-range to far edge of swath
        
        tau_near = 2*Rn/c # [2] equation 5.1.3.11
        tau_far = 2*Rf/c # [2] equation 5.1.3.12


        PRFmax = 1.0/(2.0*tau_p + tau_far - tau_near) # max allowble PRF [2] equation 5.1.3.13
        PRFmin = v_sc_kmps*1e3/SyntheticApertureRadar.get_azimuthal_resolution(v_sc_kmps, v_g_kmps, D_az) # minimum allowable PRF to satisfy Nyquist sampling criteria [2] equation 5.1.2.1 modified to [2] equation (5.4.4.2)


        f_P = None    
        # Find the highest possible prf within the input prf range which allows for unambiguous echo detection. 
        for _f_P in range(int(f_Pmax), int(f_Pmin), -1): # step down in terms of 1 Hz
            PRFOK = True

            # perform first check of PRF validity
            if((_f_P< PRFmin) or (_f_P>PRFmax)):
                PRFOK = False 
                continue # goto next iteration of _f_P
            
            # perform second check of PRF validity, check that target echo is not eclipsed by transmit pulse
            # inequality [2] 5.1.4.1 
            N = int(_f_P*2.0*Rn/c) + 1
            if(not(((N-1)/(tau_near-tau_p) < _f_P ) and (_f_P < N/(tau_far + tau_p)))):
                PRFOK = False 
                continue # goto next iteration of _f_P
                        
            # perform second check of PRF validity, check that target echo is not eclipsed by naidr echo 
            # from any of the succeeding pulses (from the transmit pulse under consideration to the pulse 
            # just before the echo
            # refer my notes for the nadir interference condition
            # inequality [2] 5.1.5.2 seems wrong. 
            tau_nadir = 2.0*h/c
            M = int(_f_P*2.0*Rf/c) + 1
            for m in range(1,M):
                if(((m/(tau_near - tau_p - tau_nadir)) > _f_P) and (_f_P > (m/(tau_far + tau_p - tau_nadir)))):
                    PRFOK = False # there is nadir echo overlap with desired echo for the mth pulse
                    continue # goto next iteration of _f_P 

            # If control has reached here and if PRFOK is still True, then this is a good PRF
            if(PRFOK == True):
                f_P = _f_P
                break

        return f_P        