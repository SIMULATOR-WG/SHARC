# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 18:15:47 2017

@author: edgar
"""

import os
import logging.config

import yaml

class Logging():

    
    @staticmethod
    def setup_logging(default_path='support/logging.yaml', 
                      default_level=logging.INFO, env_key='LOG_CFG'):
        """
        Setup logging configuration
        """
        path = default_path
        value = os.getenv(env_key, None)
        if value:
            path = value
        if os.path.exists(path):
            with open(path, 'rt') as f:
                config = yaml.safe_load(f.read())
            logging.config.dictConfig(config)
        else:
            logging.basicConfig(level=default_level)
