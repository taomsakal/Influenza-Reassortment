"""
The agent and simulation class live here.
"""

import itertools
import json
import math
import random

import numpy as np
from mesa import Agent, Model
from mesa.datacollection import DataCollector
from mesa.time import StagedActivation

from parameters import infection_table
from parameters import viruses as all_viruses

species_infected = np.sum(infection_table, axis=0)
fitness = (np.ones((16, 9)) / species_infected) * np.random.normal(1, 0.2, (16, 9))


class Host(Agent):

    def __init__(self, model, species, viruses=None, immune_virus=None):
        """
        Creates a new host organism.

        Args:
            model: The model the host is a part of.
            species: The species the host is. Either "Human", "Pig", "Bird", or "Poultry".
            viruses: The set of viruses currently infecting the host.
        """

        self.id = model.next_id()
        super().__init__(self.id, model)  # Initialize basic agent code, assign a unique id
        self.model = model
        self.it = self.model.it  # sets iteration for model

        self.species = species  # What species the host is
        assert species in ["Human", "Pig", "Bird", "Poultry"]

        # Assigns each species a number. This is for later adjacency matrix lookups.
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)

        # Temporary holder for viruses after contact but  before all individuals have contacted each other
        self.temp_viruses = set()

        # Right now the viruses that a host has are listed in a set.
        # todo: We'll explore using a matrix to store this info later to avoid loops.
        self.viruses = set()
        self.time_since_infection = 0  # used to measure the time since infection and affects recovery chance
        if (viruses is not None) and self.is_infectable_by(viruses):
            self.viruses.add(viruses)

        # Holds the antigens that the viruses are immune to
        self.h_immune = set()
        self.n_immune = set()
        if (immune_virus is not None) and self.is_infectable_by(immune_virus):
            self.h_immune.add(immune_virus[0])
            self.h_immune.add(immune_virus[1])

        # Dynamically creates variables for agent reporters to measure each virus.
        # For example, self.H1N1 measures the number of H1N1 viruses in a host.
        for i in all_viruses:
            setattr(self, f"H{i[0] + 1}N{i[1] + 1}", int(i in self.viruses))

    def __eq__(self, other):
        return self.id == other.id

    def contract_virus(self):
        """
        Spreads viruses to other agents. Stage 1 of each step.
        """

        # Gets list of contacts
        contacts = self.contacts()

        # Sets infection base rate
        self.base_rate = self.model.infection_rate

        # If there are seasons then change this base rate sinusoidally
        if self.model.seasonal_on:
            self.base_rate = self.base_rate * self.model.amplitude * (math.sin((self.it * math.pi) / self.model.period))

        # Contact other agents and receive viruses
        for contact in contacts:
            for virus in self.viruses:
                if contact.is_infectable_by(virus):

                    self.infection_rate = self.base_rate

                    if self.model.fitness_on:
                        self.infection_rate = self.infection_rate * fitness[virus[0]][virus[1]]

                    # If contact is immune to both antigens, the virus doesnt spread
                    if virus[0] in contact.h_immune and virus[1] in contact.n_immune:
                        self.infection_rate = 0

                    # If contact is immune to one antigen, the virus is less likely to spread
                    elif virus[0] in contact.h_immune:
                        self.infection_rate = self.infection_rate * self.model.cross_immunity_effect
                    elif virus[1] in contact.n_immune:
                        self.infection_rate = self.infection_rate * self.model.cross_immunity_effect

                    if np.random.rand(1, 1)[0, 0] < self.infection_rate:
                        contact.temp_viruses.add(virus)

    def recombine(self):
        """
        Viruses recombine. Stage 2 of each step.
        """
        self.viruses = self.viruses.union(self.temp_viruses)
        self.temp_viruses = set()
        self.h = [item[0] for item in self.viruses]
        self.n = [item[1] for item in self.viruses]

        # takes a long time
        self.viruses = {tuple(x) for x in itertools.product(self.h, self.n) if self.is_infectable_by(x)}

    def recovery(self):
        """
        Agents have a chance of recovery and gain immunity to recovered antigens. Stage 3 of each step.
        """

        # Recovery
        if len(self.viruses) > 0:

            # Increases times since infection
            self.time_since_infection = self.time_since_infection + 1

            # recovery chance increase over time
            self.recovery_chance = self.model.recovery_rate * self.time_since_infection

            # Removes viruses and adds immunity if the host recovers
            if np.random.rand(1, 1)[0, 0] <= self.recovery_chance:
                self.time_since_infection = 0
                self.h_immune = self.h_immune.union({item[0] for item in self.viruses})
                self.n_immune = self.n_immune.union({item[1] for item in self.viruses})
                self.viruses = set()

        # Chance of losing immunity (mimics antigenic drift)
        if len(self.h_immune) > 0:
            self.h_immune = set(
                [i for i in self.h_immune if (np.random.rand(1, 1)[0, 0] < self.model.mutation_rate)])  # 0.95

        if len(self.n_immune) > 0:
            self.n_immune = set(
                [i for i in self.n_immune if (np.random.rand(1, 1)[0, 0] < self.model.mutation_rate)])  # 0.95

        # Sets agent reporters for datacollector
        for i in all_viruses:
            setattr(self, f"H{i[0] + 1}N{i[1] + 1}", int(i in self.viruses))

    def is_infectable_by(self, virus):
        """Returns true if the host can be infected by virus."""
        return infection_table[self.species_id][virus[0]][virus[1]]  # Lookup entry in the infection table

    def birth_death(self):
        """
        Agents have a chance of birth and death. Stage 4 of each step.
        """

        # Chance of birth
        if np.random.rand(1, 1)[0, 0] < self.model.birth_rate:
            init_virus = None
            immune_virus = None
            if random.random() < self.model.immigration_rate:
                init_virus = (random.randint(0, 15), random.randint(0, 8))
            if random.random() < self.model.immigration_rate:
                immune_virus = (random.randint(0, 15), random.randint(0, 8))
            if self.species_id == 0:
                host = Host(self.model, self.species, init_virus, immune_virus)
                self.model.hosts_0.append(host)
            elif self.species_id == 1:
                host = Host(self.model, self.species, init_virus, immune_virus)
                self.model.hosts_1.append(host)
            elif self.species_id == 2:
                host = Host(self.model, self.species, init_virus, immune_virus)
                self.model.hosts_2.append(host)
            else:
                host = Host(self.model, self.species, init_virus, immune_virus)
                self.model.hosts_3.append(host)
            self.model.schedule.add(host)

        # Sets species specific mortality rates (birds and poultry less likely to die from Influenza)
        if self.species_id == 0:
            mortality = 1.25
        elif self.species_id == 1:
            mortality = 1.25
        elif self.species_id == 2:
            mortality = 1.005
        else:
            mortality = 1.005

        # Chance of death increases with number of viruses the host is infected by
        if np.random.rand(1, 1)[0, 0] < self.model.death_rate * (mortality ** len(self.viruses)):
            if self.species_id == 0:
                self.model.hosts_0.remove(self)
            elif self.species_id == 1:
                self.model.hosts_1.remove(self)
            elif self.species_id == 2:
                self.model.hosts_2.remove(self)
            else:
                self.model.hosts_3.remove(self)
            self.model.schedule.remove(self)
            self.viruses = set()
            del self

    def contacts(self):
        """Returns a list of other organism the host has contacted and got viruses from."""
        contacts = []
        num_contacts = int(len(self.model.hosts_0) * self.model.contact_rates[self.species_id][0][0])
        samp = random.sample(self.model.hosts_0, num_contacts)
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_1) * self.model.contact_rates[self.species_id][1][0])
        samp = random.sample(self.model.hosts_1, num_contacts)
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_2) * self.model.contact_rates[self.species_id][2][0])
        samp = random.sample(self.model.hosts_2, num_contacts)
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_3) * self.model.contact_rates[self.species_id][3][0])
        samp = random.sample(self.model.hosts_3, num_contacts)
        contacts = contacts + samp
        return contacts


class VirusModel(Model):

    def __init__(self, run="NA", init_pop_size=[900, 650, 1000, 750], it=0, infection_rate=0.22, recovery_rate=0.2,
                 mutation_rate=0.96, birth_rate=0.04, death_rate=0.039, cross_immunity_effect=0.05, init_viruses=None,
                 immigration_rate=0.9, contact_rates=None, fitness_on=True, seasonal_on=True, period=100,
                 amplitude=0.8):
        """
        The main virus model itself.


        Args:
            init_pop_size: The initial population size of each species [Humans, Pigs, Birds, Poultry]
            x: Batch runner throws an error without a dummy variable to use as a variable parameter, therefore this variable acts as a dummy
            it: Iteration number
        """

        super().__init__()  # Initialize basic agent code, assign a unique id
        self.it = it
        self.run = run
        self.running = True  # For batch runs
        self.model_step = 0  # The number of timesteps the simulation has run
        self.schedule = StagedActivation(self, ["contract_virus", "recombine", "recovery", "birth_death"], True,
                                         True)  # set schedule
        self.infection_rate = infection_rate  # infection rate
        self.recovery_rate = recovery_rate
        self.mutation_rate = mutation_rate
        self.birth_rate = birth_rate
        self.death_rate = death_rate
        self.cross_immunity_effect = cross_immunity_effect
        self.init_viruses = init_viruses
        self.immigration_rate = immigration_rate
        self.fitness_on = fitness_on
        self.seasonal_on = seasonal_on
        self.period = period
        self.amplitude = amplitude

        # Population sizes
        self.human_pop_size = init_pop_size[0]
        self.pig_pop_size = init_pop_size[1]
        self.bird_pop_size = init_pop_size[2]
        self.poultry_pop_size = init_pop_size[3]
        self.total_pop_size = sum(init_pop_size)

        # lists with agents of each species. Used to get contacts
        self.hosts_0 = []
        self.hosts_1 = []
        self.hosts_2 = []
        self.hosts_3 = []

        # It's a weird naming convention but self.model_step actually measures
        # the current step the model is on, not the current model_step
        # todo: could you clarify the difference here? I don't understand it. For now I've renamed this to model_step -Taom
        self.model_step = 0

        # Adjacency matrix of gaussian contact rate distributions where entry ij is 
        # the contact rate species j to species i.
        # 1 -> humans
        # 2 -> pigs
        # 3 -> birds
        # 4 -> poultry
        # todo: put more reasonable values
        if contact_rates is None:
            self.contact_rates = np.array([[[0.01], [0.0045], [0.002], [0.001]],
                                           [[0.0045], [0.0095], [0.003], [0.0045]],
                                           [[0.002], [0.003], [0.09], [0.003]],
                                           [[0.001], [0.0045], [0.003], [0.01]]])
        else:
            self.contact_rates = contact_rates

        # Initialize population
        init_virus = None
        immune_virus = None
        id = 0  # This isn't used? Also has the same name as a built in python function.
        self.all_viruses = list(all_viruses)

        for i in range(self.human_pop_size):
            # if (id % 30== 0):
            #  init_virus = self.all_viruses[int(id/30) % len(self.all_viruses)]
            # id = id +1

            if random.randint(0, 25) == 0:
                if self.init_viruses == None:
                    init_virus = random.choice([(0, 0), (1, 1), (2, 1), (4, 0), (6, 1), (6, 2), (6, 6), (8, 1), (9, 6)])
                else:
                    init_virus = random.choice(init_viruses)

            if random.randint(0, 5) == 0:
                immune_virus = random.choice([(0, 0), (1, 1), (2, 1), (4, 0), (6, 1), (6, 2), (6, 6), (8, 1), (9, 6)])

            host = Host(self, "Human", init_virus, immune_virus)
            self.schedule.add(host)
            init_virus = None
            immune_virus = None
            self.hosts_0.append(host)

        for i in range(self.pig_pop_size):
            # if (id % 30== 0):
            #  init_virus = self.all_viruses[int(id/30) % len(self.all_viruses)]
            # id = id +1

            if random.randint(0, 25) == 0:
                if self.init_viruses == None:
                    init_virus = random.choice([(0, 0), (0, 1), (1, 2), (2, 1), (2, 2), (3, 5), (4, 1), (8, 1)])
                else:
                    init_virus = random.choice(init_viruses)

            if random.randint(0, 5) == 0:
                immune_virus = random.choice([(0, 0), (0, 1), (1, 2), (2, 1), (2, 2), (3, 5), (4, 1), (8, 1)])
            host = Host(self, "Pig", init_virus, immune_virus)
            self.schedule.add(host)
            init_virus = None
            immune_virus = None
            self.hosts_1.append(host)

        for i in range(self.bird_pop_size):
            # if (id % 20== 0):
            #  init_virus = self.all_viruses[int(id/20) % len(self.all_viruses)]
            # id = id +1
            if random.randint(0, 9) == 0:
                if self.init_viruses == None:
                    init_virus = (random.randint(0, 15), random.randint(0, 8))
                else:
                    init_virus = random.choice(init_viruses)

            immune_virus = (random.randint(0, 15), random.randint(0, 8))
            host = Host(self, "Bird", init_virus, immune_virus)
            self.schedule.add(host)
            init_virus = None
            immune_virus = None
            self.hosts_2.append(host)

        for i in range(self.poultry_pop_size):
            # if (id % 20 == 0):
            #  init_virus = self.all_viruses[int(id/20) % len(self.all_viruses)]
            # id = id +1
            if random.randint(0, 9) == 0:
                if self.init_viruses == None:
                    init_virus = (random.randint(0, 12), random.randint(0, 8))
                else:
                    init_virus = random.choice(init_viruses)
            immune_virus = (random.randint(0, 12), random.randint(0, 8))
            host = Host(self, "Poultry", init_virus, immune_virus)
            self.schedule.add(host)
            self.hosts_3.append(host)
            immune_virus = None
            init_virus = None

        # Sets data reporters for each virus
        reporters = "{\"Iteration\":\"it\",\"Species\":\"species\","
        for i in all_viruses:
            name = f"H{i[0] + 1}N{i[1] + 1}"
            reporters = reporters + "\"" + name + "\"" + ":" + "\"" + name + "\"" + ","
        reporters = reporters[:-1] + "}"
        self.datacollector = DataCollector(
            # model_reporters={"Strain_data": Count_Strains}
            agent_reporters=json.loads(reporters)
        )

        # Collect data for initial condition
        self.datacollector.collect(self)

    def step(self):
        """Steps the entire model one time step."""

        # Print current status
        if self.model_step % 10 == 0:
            print(f"Iteration {self.model_step}")

        # Step all agents
        self.schedule.step()
        self.model_step += 1

        # Collect data
        self.datacollector.collect(self)


if __name__ == '__main__':
    pass
