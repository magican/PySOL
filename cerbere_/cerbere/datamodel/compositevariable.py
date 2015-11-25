# encoding: utf-8
from .variable import Variable

__all__ = ['CompositeVariable']

class CompositeVariable(Variable):
    """
    Complex variable that must be described by several quantities or components (such as a vector : (u,v) or (module,direction)
    """
    
    def __init__(self, components=None, shortname=None, description=None, authority=None, standardname=None):
        """
        """
        Variable.__init__(self,shortname, description, authority, standardname)
        self.components = components



