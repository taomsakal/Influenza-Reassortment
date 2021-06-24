"""
The agent and simulation class live here.
"""

import random
import numpy as np
from mesa import Agent, Model
from mesa.time import SimultaneousActivation

import parameters


class Host(Agent):

    def __init__(self, model, species, viruses=None):
        """
        Creates a new host organism
        Args:
            model: The model the host is a part of
            species: The speces the host is. Either "Human", "Pig", "Bird", or "Poultry".
            viruses: The set of viruses currently infecting the host.
        """

        super().__init__(model.next_id(), model)  # Initialize basic agent code, assign a unique id
        self.model = model

        self.species = species
        assert species in ["Human", "Pig", "Bird", "Poultry"]

        # Right now the viruses that a host has are listed in a set.
        # todo: We'll explore using a matrix to store this info later to avoid loops.
        if viruses is None:
            self.viruses = set()
        else:
            self.viruses = viruses

    def step(self):
        """
        Steps the organism forward in time. This does the following:

            1. It contacts other hosts
            2. It recieves viruses from the other hosts which can infect it
            3. All reassortments are created in the host
            4. Host birth/death/recovery

        """

        # Contact other hosts
        # todo: make different species have different contact rates
        contacts = self.contacts()

        # Host catches all possible viruses it is susceptible to
        for contact in contacts:
            for virus in contact:
                if self.is_infectable_by(virus):
                    self.viruses.add(virus)  # add virus to the set

        # Existing viruses recombine and best takes over
        self.viruses = self.recombine()

        # Host birth/death
        # Todo: We'll write this part after we get the rest working.

    def is_infectable_by(self, virus):
        """Returns true if the virus can infect the host."""

        pass  # todo: use the H and N tables to check infection
        # Put the table data in parameters.py.
        # You'll have to decide on a good, well performing data format to use.
        # dict and dataframes are two options.

    def recombine(self):
        """Returns the set of all possible combinations of the current viruses in the host."""
        pass  # todo

    def contacts(self):
        """Returns a list of other organism the host has contacted and got viruses from."""

        pass  # todo. self.model.contact_rates should give you the contact rate matrix


class VirusModel(Model):

    def __init__(self, run="NA", init_pop_size=100):
        """

        Args:
            run:
            init_pop_size: The initial population size of each species
        """

        super().__init__()  # Initialize basic agent code, assign a unique id

        self.run = run
        self.running = True  # For batch runs
        self.iteration = 0  # The number of timesteps the simulation has run
        self.schedule = SimultaneousActivation(self)  # Todo: double check how this works.

        # Population sizes
        self.human_pop_size = 0
        self.pig_pop_size = 0
        self.bird_pop_size = 0
        self.poultry_pop_size = 0
        self.total_pop_size = 0

        # Adjacency matrix of contact rates where entry ij is the contact rate species j to species i.
        # 1 -> humans
        # 2 -> pigs
        # 3 -> poultry
        # 4 -> wild birds
        # todo: put more reasonable values in
        self.contact_rates = np.array([[10, 2, 3, 4],
                                       [1, 10, 2, 3],
                                       [1, 2, 30, 4],
                                       [1, 2, 3, 40]])

        # todo: initialize population

    def step(self):
        """Steps the entire model one time step."""

        self.iteration += 1
        self.schedule.step()  # step all agents

        # todo: add data gathering


if __name__ == '__main__':
    pass
