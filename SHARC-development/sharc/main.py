# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:05:49 2016

@author: edgar
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
    
from sharc.model import Model
from sharc.gui.view import View
from sharc.controller import Controller
from sharc.support.logging import Logging

def main():
    Logging.setup_logging()
    
    model = Model()
    view = View()
    controller = Controller()
    
    view.set_controller(controller)
    controller.set_model(model)
    model.add_observer(view)
    
    view.mainloop()


if __name__ == "__main__":
    main()