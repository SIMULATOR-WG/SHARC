# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 16:32:32 2016

@author: edgar
"""

class Observable(object):
    """
    This class represents an observable object, or "data" in the model-view 
    paradigm. It can be subclassed to represent an object that the application 
    wants to have observed. An observable object can have one or more 
    observers. An observer may be any object that implements interface 
    Observer. After an observable instance changes, an application calling the 
    Observable's notify_observers method causes all of its observers to be 
    notified of the change by a call to their update method.
    
    Attributes
    ----------
        observers : Observer
            The list of observers
    """
 
    def __init__(self):
        self.observers = list()
 
    def add_observer(self, observer):
        """
        Adds a new observer to the current list of observers.
        
        Parameters
        ----------
        observer : Observer
            The observer to be added
        """
        if not observer in self.observers:
            self.observers.append(observer)
             
    def delete_observer(self, observer):
        """
        Deletes an observer from the current list of observers.
        
        Parameters
        ----------
        observer : Observer
            The observer to be deleted
        """
        if observer in self.observers:
            self.observers.remove(observer)
             
    def delete_observers(self):
        """
        Clears the observer list so that this object no longer has any 
        observers.
        """
        if self.observers:
            del self.observers[:]
 
    def notify_observers(self, *args, **kwargs):
        """
        Notifies all observers that this object has changed
        """
        for observer in self.observers:
            observer.notify_observer(*args, **kwargs)
