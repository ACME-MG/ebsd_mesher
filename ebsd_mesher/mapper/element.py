"""
 Title:         Element
 Description:   Represents an element
 Author:        Janzen Choi

"""

# Libraries
from ebsd_mesher.maths.orientation import deg_to_rad, rad_to_deg

# Element class
class Element:
    
    def __init__(self, phi_1:float, Phi:float, phi_2:float, grain_id:int, degrees:bool=True):
        """
        Contains information about an of a grain
        
        Parameters:
        * `phi_1`:    The phi_1 orientation of the grain
        * `Phi`:      The Phi orientation of the grain
        * `phi_2`:    The phi_2 orientation of the grain
        * `grain_id`: The ID of the grain to which the element belongs
        * `degrees`:  Whether to store the orientations as degrees
        """
        self.set_orientation(phi_1, Phi, phi_2, degrees)
        self.set_grain_id(grain_id)
    
    def set_orientation(self, phi_1:float, Phi:float, phi_2:float, degrees:bool=True) -> None:
        """
        Sets the orientations
        
        Parameters:
        * `phi_1`:   The average phi_1 orientation of the grain
        * `Phi`:     The average Phi orientation of the grain
        * `phi_2`:   The average phi_2 orientation of the grain
        * `degrees`: Whether to store the orientations as degrees
        """
        self.phi_1   = phi_1
        self.Phi     = Phi
        self.phi_2   = phi_2
        self.degrees = degrees

    def get_orientation(self, degrees:bool=True) -> tuple:
        """
        Returns the orientation of the grain
        
        Parameters:
        * `Whether to return the orientations in degrees or radians
        """
        
        # Store as radians/degrees and return as radians/degrees
        if (self.degrees and degrees) or (not self.degrees and not degrees):
            return self.phi_1, self.Phi, self.phi_2
        
        # Store as radians but return as degrees
        elif not self.degrees and degrees:
            return rad_to_deg(self.phi_1), rad_to_deg(self.Phi), rad_to_deg(self.phi_2)
        
        # Store as degrees but return as radians
        elif self.degrees and not degrees:
            return deg_to_rad(self.phi_1), deg_to_rad(self.Phi), deg_to_rad(self.phi_2)

    def set_grain_id(self, grain_id:int) -> int:
        """
        Sets the grain ID
        
        Parameters:
        * `grain_id`: The ID of the grain to which the element belongs
        """
        self.grain_id = grain_id
        
    def get_grain_id(self) -> int:
        """
        Returns the grain ID
        """
        return self.grain_id
