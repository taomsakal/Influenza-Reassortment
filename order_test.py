import numpy as np
import pandas as pd
from mesa.batchrunner import BatchRunner

from model import VirusModel

contactrates = np.array([[[0], [0], [0], [0]],
                         [[0], [0], [0], [0]],
                         [[0], [0], [0.1], [0.03]],
                         [[0], [0], [0.03], [0.1]]])
initviruses = [(7, 3), (3, 7)]
initpopsize = [0, 0, 200, 200]
fixed_params = {"run": "NA",
                "init_pop_size": initpopsize,  # default pop: [900, 700, 1000, 850]
                "init_viruses": initviruses,
                "immigration_rate": 0,
                "birth_rate": 0,
                "death_rate": 0,
                "contact_rates": contactrates,
                "fitness_on": False,
                "seasonal_on": False
                }
variable_params = {"it": [0]}
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
new_data = full_data.drop(["Iteration", "AgentID"], axis=1)
new_data = new_data.groupby(['Step', "Species"]).sum().reset_index()
print(new_data[['H8N4', 'H4N8', 'H4N4', 'H8N8']])

initviruses = [(8, 4), (4, 8)]
fixed_params = {"run": "NA",
                "init_pop_size": initpopsize,  # default pop: [900, 700, 1000, 850]
                "init_viruses": initviruses,
                "immigration_rate": 0,
                "birth_rate": 0,
                "death_rate": 0,
                "contact_rates": contactrates,
                "fitness_on": False,
                "seasonal_on": False
                }

batch_run2 = BatchRunner(VirusModel,
                         fixed_parameters=fixed_params,
                         variable_parameters=variable_params,
                         iterations=1,
                         max_steps=100,
                         display_progress=True
                         )
batch_run2.run_all()
agent_data = list(batch_run2.get_collector_agents().values())
full_data = pd.DataFrame(agent_data[0])
full_data = full_data.reset_index()
new_data = full_data.drop(["Iteration", "AgentID"], axis=1)
new_data = new_data.groupby(['Step', "Species"]).sum().reset_index()
print(new_data[['H9N5', 'H5N9', 'H5N5', 'H9N9']])
