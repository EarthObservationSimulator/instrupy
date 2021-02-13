""" 
.. module:: base

:synopsis: *Provides the factory class :class:`InstrumentFactory` for dishing out instrument objects
            based on the input dict/json-string specifications.*

"""
from .basic_sensor import BasicSensor
#from .passive_optical_scanner import PassiveOpticalScanner
#from .synthetic_aperture_radar import SyntheticApertureRadar

class InstrumentFactory:
    """ Factory class which allows to register and invoke the appropriate instrument class. 
    
    :class:`BasicSensor`, :class:`PassiveOpticalSensor` and :class:`SyntheticApertureRadar`
    instrument classes are registered in the factory. Additional user-defined instrument classes
    can be registered as shown below: 

    Usage: 
    
    .. code-block:: python
        
        factory = instrupy.InstrumentFactory()
        factory.register_instrument('Doppler Radar 2021', DopplerRadar)
        factory.get_instrument('DopplerRadar')

    :ivar _creators: Dictionary mapping instrument type label to the appropriate instrument class. 
    :vartype _creators: dict

    """

    def __init__(self):
        self._creators = {}
        self.register_instrument('BasicSensor', BasicSensor)
        #self.register_instrument('PassiveOpticalScanner', PassiveOpticalScanner)
        #self.register_instrument('SyntheticApertureRadar', SyntheticApertureRadar)

    def register_instrument(self, _type, creator):
        """ Function to register instruments.

        :var _type: Instrument type (label).
        :vartype _type: str

        :var creator: Instrument class.
        :vartype creator: Instrument class.

        """
        self._creators[_type] = creator

    def get_instrument(self, _type):
        """ Function to get the appropriate instrument instance.

        :var _type: Instrument type (label).
        :vartype _type: str
        
        """
        creator = self._creators.get(_type)
        if not creator:
            raise ValueError(_type)
        return creator()

