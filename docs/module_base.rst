
.. _base_module:

``instrupy.base`` --- Base
********************************

.. _mode_json_obj:

:code:`mode` JSON object format
================================
Several modes (in a list) maybe specified within a single instrument. Each mode corresponds to a specific operating point. For example, 
consider a *Synthetic Aperture Radar* type instrument which operates in both L-band and C-band. Such an instrument is considered
to be made up of two modes with one mode operating at L-band and the other mode at C-band. 
A mode-identifier can be specified by the user with which the corresponding mode can be referenced.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   @id, string,, Unique identifier of mode.

The parameters outside the mode block are used as the common parameters for all the modes, while the parameters specified
within a mode list entry are specific to the particular mode.

Example: The example below is that of a *Basic Sensor* type instrument with two modes. The common parameters for both the modes
are outside the :code:`mode` block. The `NadirObservationMode` has a nadir orientation while the `SideObservationMode`
has an off-nadir orientation.
 
.. code-block:: python

               specs = '{        
                           "@type": "Basic Sensor",
                           "name": "Atom",
                           "@id": "senX",  
                           "mass": 28, 
                           "volume": 0.12, 
                           "power": 32, 
                           "bitsPerPixel": 8, 
                           "fieldOfViewGeometry": {
                                       "sensorGeometry": "CIRCULAR",
                                       "diameter": 35
                                 },
                           "mode":[{
                                    "@id": "NadirObservationMode",                            
                                    "orientation": {
                                          "referenceFrame": "SC_BODY_FIXED",
                                          "convention": "REF_FRAME_ALIGNED"
                                    }      
                                 },
                                 {
                                    "@id": "SideObservationMode",
                                    "orientation": {
                                       "referenceFrame": "SC_BODY_FIXED",
                                       "convention": "SIDE_LOOK",
                                       "sideLookAngle": 30
                                 }       
                                 }
                           ]
                        }'

               x = Instrument.from_json(specs) 

               
.. automodule:: instrupy.base
   :members:
   :undoc-members:
   :show-inheritance: