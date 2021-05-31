# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:05:49 2016

@author: edgar
@modified: Luciano Camilo Fri May 28 15:49:00 2021

SHARC - Version 2.1.05
"""

import sys
import os

from sharc.model import Model
from sharc.gui.view import View
from sharc.controller import Controller
from sharc.support.logging1 import Logging
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))


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
