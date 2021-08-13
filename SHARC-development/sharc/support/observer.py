# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 16:30:55 2016

@author: edgar
"""

from abc import ABCMeta, abstractmethod
 
class Observer(object):
    """
    This is the abstract base class of the Observer pattern. A concrete 
    class can extend this class when it wants to be informed of changes in 
    observable objects.
    """
    
    __metaclass__ = ABCMeta
     
    @abstractmethod
    def notify_observer(self, *args, **kwargs):
        """
        This method is called whenever the observed object is changed. An 
        application calls an Observable object's notify_observers method to 
        have all the object's observers notified of the change.
        """
        pass
