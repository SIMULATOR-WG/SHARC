# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:05:49 2016

@author: edgar
"""


    
from sharc.model import Model
from sharc.gui.view import View
from sharc.controller import Controller
from sharc.support.logging import Logging

Logging.setup_logging()

model = Model()
view = View()
controller = Controller()

view.set_controller(controller)
controller.set_model(model)
model.add_observer(view)

view.mainloop()
