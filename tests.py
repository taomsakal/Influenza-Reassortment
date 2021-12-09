import numpy as np
import pandas as pd
from mesa.batchrunner import BatchRunner

from model import Host, VirusModel, ONES, ZEROS
from model import VirusModel
import time

testviruses = np.array(
              [ [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.] ])

# Same as testviruses but with a few entries deleted.
testviruses2 = np.array(
              [ [1., 1., 1., 0., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 0., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 0., 0., 0., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 0., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 1., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 0., 1.],
                [1., 1., 1., 1., 0., 0., 1., 1., 0., 0.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 1., 1., 0., 0., 1., 1., 1., 1.] ])

testviruses3  = np.array(
              [ [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.] ])
class TestClass:

    def test_recovery_single(self):
        """Tests if one organism recovers from all viruses."""

        # Setup
        model = VirusModel(init_hosts=False)
        host = Host(model, "Human", viruses=testviruses)

        host.recovery_prob = 1
        host.recover()

        assert np.array_equal(host.viruses, ZEROS.astype(float))  # Make sure recovered from all viruses
        print(host.susceptibility.astype(bool))
        print(np.invert(testviruses.astype(bool)))  # Make sure is immune to those viruses

    def test_recombine(self):
        """Test the recombine function."""

        # Setup
        model = VirusModel(init_hosts=False)
        host = Host(model, "Human")
        host.mutation_prob = 0  # Make sure no mutations happen

        # testviruses already has all possible combinations, so the host's viruses
        # should stay the same after recombination.
        host.viruses = testviruses
        host.recombine()
        assert np.array_equal(host.viruses, testviruses)

        # testviruses2 should recombine into testviruses
        host.viruses = testviruses2
        host.recombine()
        assert np.array_equal(host.viruses, testviruses)

        # A matrix with all possible H1Nx and H*N1
        # This should recombine into all possible viruses
        host.viruses = testviruses3
        host.recombine()
        assert np.array_equal(host.viruses, ONES)



    def test_susceptibility(self):
        """If a host recovers from all viruses they have they should no longer be susceptible to them."""

        model = VirusModel(init_hosts=False)
        host = Host(model, "Human", viruses=testviruses)

        host.recovery_prob = 1
        host.recover()

        assert np.array_equal(host.viruses, ZEROS.astype(float))  # Make sure recovered from all viruses
        print(host.susceptibility.astype(bool))
        print(np.invert(testviruses.astype(bool)))  # Make sure is immune to those viruses


    def test_recovery(self):
        contactrates = np.array([[[0], [0], [0], [0]],
                                 [[0], [0], [0], [0]],
                                 [[0], [0], [0.1], [0.03]],
                                 [[0], [0], [0.03], [0.1]]])
        initviruses = [(7, 3)]
        initpopsize = [100, 100, 100, 100]

        fixed_params = {"run": "NA",
                        "init_pop_size": initpopsize,  # default pop: [900, 700, 1000, 850]
                        "recovery_rate": 0,
                        "birth_rate": 0,
                        "infection_rate": 1,
                        "death_rate": 0,
                        "init_viruses": initviruses,
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
        new_data = full_data.drop(["Species", "Iteration"], axis=1)
        new_data = new_data.groupby(['Step']).sum().reset_index()
        assert (new_data[new_data.Step == 9].iloc[0]["H8N4"] == 200)

    def test_recovery2(self):
        contactrates = np.array([[[0], [0], [0], [0]],
                                 [[0], [0], [0], [0]],
                                 [[0], [0], [0.1], [0.03]],
                                 [[0], [0], [0.03], [0.1]]])
        initviruses = [(7, 3)]
        initpopsize = [50, 50, 50, 50]

        fixed_params = {"run": "NA",
                        "init_pop_size": initpopsize,  # default pop: [900, 700, 1000, 850]
                        "recovery_rate": 1,
                        "birth_rate": 0,
                        "infection_rate": 1,
                        "death_rate": 0,
                        "init_viruses": initviruses,
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
        new_data = full_data.drop(["Species", "Iteration"], axis=1)
        new_data = new_data.groupby(['Step']).sum().reset_index()
        assert (new_data[new_data.Step == 9].iloc[0]["H8N4"] == 0)

    def test_reassortment(self):
        class Ex:
            def __init__(self):
                self.viruses = {(0, 0), (1, 0), (0, 1)}
                self.temp_viruses = {}

            def is_infectable_by(self, x):
                return True

        a = Ex()
        Host.recombine(a)
        assert a.viruses == {(0, 0), (1, 0), (0, 1), (1, 1)}
