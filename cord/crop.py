from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from .util import *


class District():

  def __init__(self, df, key):
    self.key = key
    for k,v in json.load(open('cord/districts/%s_properties.json' % key)).items():
        setattr(self,k,v)
