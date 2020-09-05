Synthetic Aperture Radar Description
*************************************


Types of Operating Modes
==========================
There are two distinct types of modes supported: (1) *Stripmap,* (2) *ScanSAR*. This is to be specified with the :code:`@type` key name 
in the :code:`mode` JSON block. The default :code:`@type` is *Stripmap*. 

The below example declares a *Synthetic Aperture Radar* type instrument with three modes: *stripmap*
(L band), *scansar* (C Band) and *stripmap* (L Band). 


Example:

.. code-block:: python

               {
                  "@id": "lcSAR",
                  "@type": "Synthetic Aperture Radar",
                  "mass": 200,
                  .
                  .
                  .
                  "mode": [
                           {
                              "@type": "stripmap",
                              "@id": "stripmL",
                              "operatingFrequency": 1.2,
                              .
                              .
                              .
                           },
                           {  
                              "@type": "scansar",
                              "@id": "scsC",
                              "operatingFrequency": 5.3,
                              .
                              .
                              .
                           },
                           { 
                              "@type": "stripmap",
                              "@id": "stripmC",
                              "operatingFrequency": 5.3,
                              .
                              .
                              .
                           }
               }

The below example does not have a :code:`mode` JSON specification block, and by default it is assumed to be a single type of mode with
the type of mode as *Stripmap*.

.. code-block:: python

               {
                  "@id": "microLsar",
                  "@type": "Synthetic Aperture Radar",
                  "mass": 200,
                  "operatingFrequency": 1.2,
                  .
                  .
                  .
               }

The below example declares a single mode *Synthetic Aperture Radar* type instrument with the mode type as *ScanSAR*.

.. code-block:: python

               {
                  "@id": "eons",
                  "@type": "Synthetic Aperture Radar",
                  "mass": 200,
                  "mode": [{                           
                              "@type": "scansar",
                              "@id": "scsL",
                              .
                              .
                              .
                           }
               }
.. toctree::
   :maxdepth: 2
   :caption: Detailed description of SAR acquisition modes

   synthetic_aperture_radar_stripmap_desc
   synthetic_aperture_radar_scansar_desc