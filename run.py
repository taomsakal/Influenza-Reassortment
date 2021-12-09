"""Python script that will run the model when invoked via mesa runserver."""
import numpy as np
import pandas as pd
from mesa.batchrunner import BatchRunner

from model import Host
from model import VirusModel
import time




fixed_params = {"run": "NA",
               "init_pop_size": [250, 250, 250, 250]}
variable_params = {"it":[0]}

batch_run = BatchRunner(VirusModel,
                        fixed_parameters=fixed_params,
                        variable_parameters=variable_params,
                        iterations=1,
                        max_steps=100,
                        display_progress=False
                        )
batch_run.run_all()
agent_data = list(batch_run.get_collector_agents().values())
print(agent_data)
full_data = pd.DataFrame(agent_data[0])
full_data = full_data.groupby(['Step', 'Iteration', "Species"]).sum()

for i in range(1, len(agent_data) - 1):
    data = pd.DataFrame(agent_data[i])
    data = data.groupby(['Step', 'Iteration', "Species"]).sum()
    full_data = pd.concat([full_data, data], axis=0)

print(full_data.head())
print(full_data.tail())
