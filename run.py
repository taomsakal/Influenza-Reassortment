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
import pandas as pd

import model
from model import VirusModel
from model import Host
import matplotlib.pyplot as plt

import model
from model import VirusModel


class TestClass:
      def test_recovery(self):
            contactrates = np.array([[[0], [0], [0], [0]],
                        [[0], [0], [0], [0]],
                        [[0], [0], [0.1], [0.03]],
                        [[0], [0], [0.03], [0.1]]])
            initviruses = [(7,3)]
            initpopsize = [100, 100, 100, 100]

            fixed_params = {"run": "NA",
                  "init_pop_size": initpopsize, #default pop: [900, 700, 1000, 850]
                  "recovery_rate":0,
                  "birth_rate":0,
                  "infection_rate":1,
                  "death_rate":0,
                  "init_viruses": initviruses,
                  "contact_rates":contactrates,
                  "fitness_on":False,
                  "seasonal_on":False
                  }
            variable_params = {"it":[0]}
            batch_run = BatchRunner(VirusModel,
                        fixed_parameters=fixed_params,
                        variable_parameters=variable_params,
                        iterations=1,
                        max_steps=100,
                        display_progress=True
                        )
            batch_run.run_all()
            agent_data = list(batch_run.get_collector_agents().values())
            full_data = pd.DataFrame(agent_data[0])
            full_data = full_data.reset_index()
            new_data = full_data.drop(["Species", "Iteration"], axis=1)
            new_data = new_data.groupby(['Step']).sum().reset_index()
            assert (new_data[new_data.Step == 9].iloc[0]["H8N4"] == 200)
      
      def test_recovery2(self):
            contactrates = np.array([[[0], [0], [0], [0]],
                        [[0], [0], [0], [0]],
                        [[0], [0], [0.1], [0.03]],
                        [[0], [0], [0.03], [0.1]]])
            initviruses = [(7,3)]
            initpopsize = [50, 50, 50, 50]

            fixed_params = {"run": "NA",
                  "init_pop_size": initpopsize, #default pop: [900, 700, 1000, 850]
                  "init_viruses": initviruses,
                  "contact_rates":contactrates,
                  "fitness_on":False,
                  "seasonal_on":False
                  }
            variable_params = {"it":[0]}
            batch_run = BatchRunner(VirusModel,
                        fixed_parameters=fixed_params,
                        variable_parameters=variable_params,
                        iterations=1,
                        max_steps=100,
                        display_progress=True
                        )
            batch_run.run_all()
            agent_data = list(batch_run.get_collector_agents().values())
            full_data = pd.DataFrame(agent_data[0])
            full_data = full_data.reset_index()
            return(full_data)
            #new_data = full_data.drop(["Species", "Iteration"], axis=1)
            #new_data = new_data.groupby(['Step']).sum().reset_index()
            #assert (new_data[new_data.Step == 9].iloc[0]["H8N4"] == 0)

      def test_reassortment(self):
            class Ex:
                  def __init__(self):
                        self.viruses ={(0,0), (1,0), (0,1)}
                        self.temp_viruses= {}
                  def is_infectable_by(self, x):
                        return True
            a = Ex()
            Host.recombine(a)
            assert a.viruses == {(0,0), (1,0), (0,1),(1,1)}






fixed_params = {"run": "NA",
               "init_pop_size": [900, 700, 1000, 850]}
variable_params = {"it":[0]}

batch_run = BatchRunner(VirusModel,
                        fixed_parameters=fixed_params,
                        variable_parameters=variable_params,
                        iterations=1,
                        max_steps=10,
                        display_progress=False
                        )

batch_run.run_all()
agent_data = list(batch_run.get_collector_agents().values())
print(agent_data)
full_data = pd.DataFrame(agent_data[0])
full_data = full_data.groupby(['Step','Iteration',"Species"]).sum()

for i in range(1, len(agent_data)-1):
      data = pd.DataFrame(agent_data[i])
      data = data.groupby(['Step','Iteration',"Species"]).sum()
      full_data = pd.concat([full_data, data], axis=0)

print(full_data.head())