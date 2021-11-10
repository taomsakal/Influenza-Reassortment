"""
The agent and simulation class live here.
"""

import json
import random

import numpy as np
from mesa import Agent, Model
from mesa.datacollection import DataCollector
from mesa.time import StagedActivation


from parameters import infection_table
from parameters import viruses as all_viruses

NUM_H = 17  # Number of rows
NUM_N = 10  # Number of columns

ONES = np.ones([NUM_H, NUM_N])
ZEROS = np.zeros([NUM_H, NUM_N])

rng = np.random.default_rng(seed=2021)

species_infected = np.sum(infection_table, axis=0)
fitness = (np.ones((17, 10)) / species_infected) * rng.normal(1, 0.2, (17, 10))


class Host(Agent):

    def __init__(self, model, species, viruses=ZEROS):
        """
        Creates a new host organism
        Args:
            model: The model the host is a part of
            species: The species the host is. Either "Human", "Pig", "Bird", or "Poultry".
            viruses: The set of viruses currently infecting the host.
        """

        self.rng = np.random.default_rng(seed=2021)
        self.id = model.next_id()
        super().__init__(self.id, model)  # Initialize basic agent code, assign a unique id
        self.model = model
        self.it = 0


        # Viruses will be held in a matrix of size num_H x num_N.
        # The entry i,j represents the probability of having virus HiNj.
        self.viruses = ZEROS + viruses

        self.immunity = ONES  # A matrix. Entry i,j is 1 if susceptible to HiNj, 0 if immune. A value in between represents partial immunity.

        self.species = species  # Species of the host organism
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)  # Species id number for looking up in infection table with

        # holder for viruses after contact before all individuals have contacted each other
        self.temp_viruses = np.array([])

        self.virus_list = [tuple(x) for x in np.asarray(np.where(self.viruses == 1)).T]

        # sets variables for agent reporters to measure each virus
        for i in all_viruses:
            setattr(self, f"H{i[0] + 1}N{i[1] + 1}", int(tuple(i) in self.virus_list))

    def __eq__(self, other):
        return self.id == other.id

    def contract_virus(self):
        """
        STAGE 1
        Contacts other agents and is exposed to viruses.
        """


        contacts = self.contacts()  # Get list of contacts
        exposures = [contact.viruses for contact in contacts]  # Get virus matrices of those exposed to
        transmission_probabilities = [exposure *
                                      self.model.transmission_prob *
                                      self.immunity
                                      for exposure in
                                      exposures]  # Find the transmission probability of each virus, taking into account immunity
        transmitted_viruses = [self.collapse_probabilities(p) for p in
                               transmission_probabilities]  # Decide which viruses actually did infect
        transmitted_viruses = np.sum(transmitted_viruses, axis=0)  # Sum up all the infections from each contact
        transmitted_viruses = infection_table[
                                  self.species_id] * transmitted_viruses  # Filter out the ones species is ressistant to.

        self.temp_viruses = transmitted_viruses  # Store these viruses in the temp viruses variable until after all contacts are finished.

        #
        # for contact in contacts:
        #     viruses_transferred = np.array(self.viruses) * infection_table[contact.species_id] * self.base_rate
        #     viruses_transferred[:,0] = 0
        #     viruses_transferred[0,:] = 0
        #     if(self.model.fitness_on):
        #         viruses_transferred = viruses_transferred * fitness
        #
        #     for i in range (16):
        #         if (contact.viruses[i + 1, 0] == 1):
        #             viruses_transferred[i+1, 1:10] = viruses_transferred[i+1, 1:10] * 0.5
        #     for i in range (9):
        #         if (contact.viruses[0, i+1] == 1):
        #             viruses_transferred[1:16, i+1] = viruses_transferred[1:16, i+1] * 0.5
        #
        #     contact.temp_viruses = helpers.vectransfer(viruses_transferred)

    def reduce_to_one(self, x):
        """
        Returns an array of True an Falses, where True is when an entry is bigger than one.
        This means all values less than one become zero and all above one become one.
        """
        return x > 1

    def recombine(self):
        """
        Transmitted viruses appear in the host and viruses recombine. Stage 2 of each step.
        """

        self.viruses += self.temp_viruses  # Add transmitted viruses to matrix of viruses host has
        self.viruses = self.viruses.astype(float)  # Convert from bool to floats

        H = np.sum(self.viruses, axis=0)  # Sum up rows to get an array that has a positive number in entry i when Hi is present.
        N = np.sum(self.viruses, axis=1)  # Sum up each column to get an array with a positive number in entry i when Ni is present.
        combos = np.outer(N, H)  # Outer product to get every possible combination.
        combos = combos >= ONES  # Ensure any number above 1 becomes 1. All others become zero.
        self.viruses = combos.astype(float)

        # virus_list = np.where(self.viruses[1:16, 1:9] == 1)
        # self.h = set(virus_list[0] + 1)
        # self.n = set(virus_list[1] + 1)
        #
        # for x in self.h:
        #     for y in self.n:
        #         self.viruses[x, y] = 1


    def birth_death(self):
        """
        Agents have a chance of birth and death. Stage 3 of each step.
        """

        #self.virus_list = [tuple(x) for x in np.asarray(np.where(self.viruses == 1)).T]

        # for i in all_viruses:
        #     setattr(self, f"H{i[0] + 1}N{i[1] + 1}", int(tuple(i) in self.virus_list))

        # Todo: move this to model level
        # # Sets species specific mortality rates (birds and poultry less liekly to die from Influenza)
        # if (self.species_id == 0):
        #     mortality = 1.25
        # elif (self.species_id == 1):
        #     mortality = 1.25
        # elif (self.species_id == 2):
        #     mortality = 1.005
        # else:
        #     mortality = 1.005

        # Todo: Chance of death increases with number of viruses the host is infected by. Calculate the number via np.sum(self.viruses)
        # If die just replace host with an empty one. Ie a new organism took the old one's place.
        if self.rng.random(1)[0] < self.model.death_rate:
            self.viruses = ZEROS
            self.immunity = ONES


        #     if (self.species == "Human"):
        #         self.model.hosts_0.remove(self)
        #         self.model.len_hosts_0 = self.model.len_hosts_0 - 1
        #     elif (self.species_id == 1):
        #         self.model.hosts_1.remove(self)
        #         self.model.len_hosts_1 = self.model.len_hosts_1 - 1
        #     elif (self.species_id == 2):
        #         self.model.hosts_2.remove(self)
        #         self.model.len_hosts_2 = self.model.len_hosts_2 - 1
        #     else:
        #         self.model.hosts_3.remove(self)
        #         self.model.len_hosts_3 = self.model.len_hosts_3 - 1
        #
        # self.model.schedule.remove(self)  # Remove current agent
        # self.model.schedule.add(Host(self.model, self.species))  # Replace with new birth

    def contacts(self):
        """Returns a list of other organism the host has contacted and got viruses from."""
        contacts = []

        num_contacts = int(self.model.len_hosts_0 * self.model.contact_rates[self.species_id][0])
        samp = list(self.rng.choice(self.model.hosts_0, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_1 * self.model.contact_rates[self.species_id][1])
        samp = list(self.rng.choice(self.model.hosts_1, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_2 * self.model.contact_rates[self.species_id][2])
        samp = list(self.rng.choice(self.model.hosts_2, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_3 * self.model.contact_rates[self.species_id][3])
        samp = list(self.rng.choice(self.model.hosts_3, num_contacts))
        contacts = contacts + samp

        return contacts

    def collapse_probabilities(self, p):
        """
        Takes a numpy array of probabilities and randomly determine which become zero and which become one.
        We do this by making a uniform random matrix and then comparing the values of the probability matrix to that.
        Those that are less than the random matrix become zero and the rest become one.

        Args:
            p: Numpy probability array

        Returns:
            A probability matrix of Trues and Falses (ie zeros and ones.)
        """

        r = np.random.rand(NUM_H, NUM_N)  # Random matrix with uniform probability between zero and one
        a = r < p
        return a.astype(float)


class VirusModel(Model):

    def __init__(self, run="NA", init_pop_size=[900, 650, 1000, 750], it=0, infection_rate=0.25, recovery_rate=0.2,
                 mutation_rate=0.23, birth_rate=0.04, death_rate=0.03, cross_immunity_effect=0.05, init_viruses=None,
                 immigration_rate=0.02, contact_rates=None, fitness_on=True):
        """
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
        self.schedule = StagedActivation(self, ["contract_virus", "recombine", "birth_death"])  # set schedule
        self.infection_rate = infection_rate
        self.recovery_rate = recovery_rate
        self.mutation_rate = mutation_rate
        self.birth_rate = birth_rate
        self.death_rate = death_rate
        self.cross_immunity_effect = cross_immunity_effect
        self.init_viruses = init_viruses
        self.immigration_rate = immigration_rate
        self.fitness_on = fitness_on
        self.transmission_prob = .5

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

        self.model_step = 0


        # Adjacency matrix of gaussian contact rate distributions where entry ij is 
        # the contact rate species j to species i.
        # 1 -> humans
        # 2 -> pigs
        # 3 -> birds
        # 4 -> poultry
        # todo: put more reasonable values
        if (contact_rates is None):
            self.contact_rates = np.array([[0.01, 0.0045, 0.002, 0.001],
                                           [0.0045, 0.0095, 0.003, 0.0045],
                                           [0.002, 0.003, 0.09, 0.003],
                                           [0.001, 0.0045, 0.003, 0.01]])
        else:
            self.contact_rates = contact_rates

        # initialize population


        # Make a bunch of random organisms for now
        for i in range(np.sum(init_pop_size)):

            species = random.choice(["Human", "Pig", "Bird", "Poultry"])  # Decide species with equal probability.
            init_viruses = np.random.choice([0, 1], size=(NUM_H, NUM_N), p=[.01, .99])  # Randomly decide some viruses it has
            host = Host(self, species, init_viruses)  # Make the host
            self.schedule.add(host)  # Add it to the list of hosts that the model simulates

            # Add host to correct species pool.
            if species == "Human":
                self.hosts_0.append(host)
            if species == "Pig":
                self.hosts_1.append(host)
            if species == "Birds":
                self.hosts_2.append(host)
            if species == "Poultry":
                self.hosts_3.append(host)

            # sets reporters for each virus
        reporters = "{\"Iteration\":\"it\",\"Species\":\"species\","
        for i in all_viruses:
            name = f"H{i[0] + 1}N{i[1] + 1}"
            reporters = reporters + "\"" + name + "\"" + ":" + "\"" + name + "\"" + ","
        reporters = reporters[:-1] + "}"
        self.datacollector = DataCollector(
            # model_reporters={"Strain_data": Count_Strains}
            agent_reporters=json.loads(reporters)
        )

    def step(self):

        print(f"Step {self.model_step}")


        """Steps the entire model one time step."""
        self.model_step += 1

        # collects data after 350 steps
        if (self.model_step >= 0):
            self.datacollector.collect(self)

        self.len_hosts_0 = len(self.hosts_0)
        self.len_hosts_1 = len(self.hosts_1)
        self.len_hosts_2 = len(self.hosts_2)
        self.len_hosts_3 = len(self.hosts_3)

        self.schedule.step()  # step all agents

        # Todo: put this in agent class
        # recovery
        recovering = rng.choice(np.array(self.schedule.agents), int(len(self.schedule.agents) * self.recovery_rate))
        for agent in recovering:
            virus_list = np.where(agent.viruses[1:16, 1:9] == 1)
            h = set(virus_list[0] + 1)
            n = set(virus_list[1] + 1)
            for i in h:
                agent.viruses[i, 0] = 1
            for i in n:
                agent.viruses[0, i] = 1
            agent.viruses[1:16, 1:9] = 0


        # Todo: put this in agent class
        # antigen Drift
        losing_immunity = rng.choice(np.array(self.schedule.agents),
                                     int(len(self.schedule.agents) * self.mutation_rate))
        h = np.ones(17)
        h[:int(17 * self.mutation_rate)] = 1
        rng.shuffle(h)
        n = np.ones(10)
        n[:int(10 * self.mutation_rate)] = 1
        rng.shuffle(n)
        loss_array = np.ones((17, 10))
        loss_array[:, 0] = h
        loss_array[0, :] = n

        for agent in losing_immunity:
            agent.viruses = agent.viruses * loss_array

        # Todo: make it so a new immigant replaces and old one instead of being added to the population?
        #     # Immigration
        # for i in range(int(self.len_hosts_0 * self.immigration_rate)):
        #     init_virus = np.zeros(170)
        #     init_virus[:34] = 1
        #     rng.shuffle(init_virus)
        #     init_virus = np.reshape(init_virus, (17, 10))
        #
        #     host = Host(self, "Human", init_virus)
        #     self.schedule.add(host)
        #     init_virus = None
        #     self.hosts_0.append(host)
        # for i in range(int(self.len_hosts_1 * self.immigration_rate)):
        #     init_virus = np.zeros(170)
        #     init_virus[:34] = 1
        #     rng.shuffle(init_virus)
        #     init_virus = np.reshape(init_virus, (17, 10))
        #
        #     host = Host(self, "Pig", init_virus)
        #     self.schedule.add(host)
        #     init_virus = None
        #     self.hosts_1.append(host)
        # for i in range(int(self.len_hosts_2 * self.immigration_rate)):
        #     init_virus = np.zeros(170)
        #     init_virus[:34] = 1
        #     rng.shuffle(init_virus)
        #     init_virus = np.reshape(init_virus, (17, 10))
        #
        #     host = Host(self, "Bird", init_virus)
        #     self.schedule.add(host)
        #     init_virus = None
        #     self.hosts_2.append(host)
        # for i in range(int(self.len_hosts_3 * self.immigration_rate)):
        #     init_virus = np.zeros(170)
        #     init_virus[:34] = 1
        #     rng.shuffle(init_virus)
        #     init_virus = np.reshape(init_virus, (17, 10))
        #
        #     host = Host(self, "Poultry", init_virus)
        #     self.schedule.add(host)
        #     self.hosts_3.append(host)
        # #     init_virus = None


if __name__ == '__main__':
    pass
