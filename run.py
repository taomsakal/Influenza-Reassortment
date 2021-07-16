"""Python script that will run the model when invoked via mesa runserver."""
import random
from random import *
import numpy as np
import itertools
from mesa import Agent, Model
from mesa.time import SimultaneousActivation, StagedActivation
from mesa.datacollection import DataCollector
from parameters import infection_table
from mesa.batchrunner import BatchRunner
import pandas

import model
from model import VirusModel
import matplotlib.pyplot as plt

import model
from model import VirusModel

fixed_params = {"run": "NA",
               "init_pop_size": [1350, 900, 1200, 1500]}
variable_params = {"x":[0]}

batch_run = BatchRunner(VirusModel,
                        fixed_parameters=fixed_params,
                        variable_parameters=variable_params,
                        iterations=1,
                        max_steps=100,
                        display_progress=True)

batch_run.run_all()

agent_data = batch_run.get_collector_agents()
df = agent_data[(0,0)]

df.to_csv("data/data4.csv")