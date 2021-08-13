# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 16:30:36 2017

@author: edgar
"""

import tkinter.scrolledtext
import queue

class ThreadSafeScrolledText(tkinter.scrolledtext.ScrolledText):
    """
    This is a subclass of ScrolledText that uses the Queue object and make the
    widget ready to be used in a multithreading environment. All the user
    interface code should be run in the main thread; other threads will write
    to the Queue object.
    """
    
    def __init__(self, master, **options):
        """
        Creates the Queue object and starts the update method.
        """
        tkinter.Text.__init__(self, master, **options)
        self.__queue = queue.Queue()
        self.__update()
        
    def write(self, line: str):
        """
        Puts in the queue the line to be displayed.
        """
        self.__queue.put(line)
        
    def __update(self):
        """
        Periodically checks if queue contains something to be printed. If queue
        is empty, exception is thrown and nothing is done.
        """
        try:
            while 1:
                line = self.__queue.get_nowait()
                self.config(state = tkinter.NORMAL)
                self.insert(tkinter.END, str(line))
                self.see(tkinter.END)
                self.config(state = tkinter.DISABLED)
                self.update_idletasks()
        except queue.Empty:
            pass
        # interval is 100 ms
        self.after(100, self.__update)
        