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

# NUM_H by NUM_N matrix of all ones or all zeros. Used in various functions.
ONES = np.ones([NUM_H, NUM_N])
ZEROS = np.zeros([NUM_H, NUM_N])

# Random number generator
rng = np.random.default_rng(seed=2021)


#species_infected = np.sum(infection_table, axis=0)
#fitness = (np.ones((17, 10)) / species_infected) * rng.normal(1, 0.2, (17, 10))


class Host(Agent):

    def __init__(self, model, species, viruses=ZEROS):
        """
        Creates a new host organism
        Args:
            model: The model the host is a part of
            species: The species the host is. Either "Human", "Pig", "Bird", or "Poultry".
            viruses: The set of viruses currently infecting the host.
        """

        # Initialize basic agent code.
        self.id = model.next_id()
        super().__init__(self.id, model)
        self.model = model
        self.it = 0

        # Parameters
        self.mutation_prob = 0.0001  # Probability that one virus mutates
        self.recovery_prob = 0.1  # Probability that host recovers from ALL viruses
        self.death_rate = 0.005  # Probability that host dies


        # Viruses will be held in a matrix of size num_H x num_N.
        # The entry i,j represents the probability of having virus HiNj.
        self.viruses = ZEROS + viruses
        self.H = np.sum(self.viruses,
                        axis=1)  # Sum up rows to get an array that has a positive number in entry i when Hi is present.
        self.N = np.sum(self.viruses,
                        axis=0)  # Sum up each column to get an array with a positive number in entry i when Ni is present.

        # Susceptiblity matrix. Entry i,j is 1 if susceptible to HiNj,
        # 0 if immune. A value in between represents partial susceptibility.
        self.susceptibility = ONES

        self.species = species  # Species of the host organism
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)  # Species id number for looking up in infection table with

        # Sets species specific mortality rates, relative to base rate.
        if self.species == "Human":
            self.death_rate = self.death_rate * .75
        elif self.species == "Pig":
            self.death_rate = self.death_rate * 1.2
        elif self.species == "Bird":
            self.death_rate = self.death_rate * 1
        elif self.species == "Poultry":
            self.death_rate = self.death_rate * 1.5

        # Holder for viruses after contact before all individuals have contacted each other
        self.temp_viruses = ZEROS

        # sets variables for agent reporters to measure each virus
        self.virus_list = [tuple(x) for x in np.asarray(np.where(self.viruses == 1)).T]
        for i in all_viruses:
            setattr(self, f"H{i[0] + 1}N{i[1] + 1}", int(tuple(i) in self.virus_list))

    def __eq__(self, other):
        """Says two agents are equal if they share the same unique id"""
        return self.id == other.id

    def contract_virus(self):
        """
        STAGE 1
        Contacts other agents and get exposed to viruses they have, then decide if those viruses cause infection.

        Todo: Bugfix. It seems the first row and column are always transmitted?
        """

        if rng.random() < .01:
            print(self.species, self.viruses)

        contacts = self.contacts()  # Get list of contacts
        exposures = [contact.viruses for contact in contacts]  # Get virus matrices of those agen was exposed to

        # Find the transmission probability of each virus, taking into account susceptibility
        transmission_probabilities = [exposure *
                                      self.model.transmission_prob *
                                      self.susceptibility
                                      for exposure in
                                      exposures]

        # Decide which viruses successfully infected the agent
        transmitted_viruses = [self.collapse_probabilities(p) for p in
                               transmission_probabilities]

        # Sum up all the infections from each contact
        transmitted_viruses = np.sum(transmitted_viruses, axis=0)

        # Filter out the ones species is resistant to.
        transmitted_viruses = infection_table[
                                  self.species_id] * transmitted_viruses

        # Store these viruses in the temp viruses variable until after all contacts are finished.
        self.temp_viruses = transmitted_viruses




    def recombine(self):
        """
        STAGE 2

        Transmitted viruses are placed inside the host. Then a mutation might happen.
        Finally all viruses within the host recombine.
        """

        self.viruses += self.temp_viruses  # Add transmitted viruses to matrix of viruses host has AFTER all hosts went though the transmission step
        self.viruses = self.viruses.astype(float)  # Convert from bool to floats


        self.H = np.sum(self.viruses, axis=1)  # Sum up rows to get an array that has a positive number in entry i when Hi is present.
        self.N = np.sum(self.viruses, axis=0)  # Sum up each column to get an array with a positive number in entry i when Ni is present.
        self.H, self.N = self.mutate(self.H, self.N)
        combos = np.outer(self.H, self.N)  # Outer product to get every possible combination.
        combos = combos >= ONES  # Ensure any number above 1 becomes 1. All others become zero.
        self.viruses = combos.astype(float)  # Convert from bool to floats


    def mutate(self, H, N):
        """
        Each virus in the host has a chance of mutating into a new virus with either H or N changed.
        This function adds a new H or N protein

        Note: This assumes that a H or N mutation is rare and only happens once in a host. Because of this
        a uniform mutation probability is a close approximation. (Assuming the host always carries
        constant number of viruses, equally distributed between all types it carries.)

        This data is then passed to recombine. Because of how the calculation works
        hosts without any viruses automatically throw away the mutant, which means
        we can skip the check to see if the host has any viruses.

        """

        if rng.random(1) < self.mutation_prob:
            if rng.random(1) < .5:  # Equal chance to mutate into H or N
                # Pick a random index and add to the H list
                new_H_index = rng.integers(NUM_H)
                H[new_H_index] += 1
            else:
                # Pick a random N and add to N list.
                new_N_index = rng.integers(NUM_N)
                N[new_N_index] += 1

        return H, N

    def recover(self):
        """
        STAGE 3

        Recover from all current viruses and become immune to those of that type.
        (This might be a hugely unrealistic approximation with the host recovering from everything at once.)
        """

        if rng.random(1) < self.recovery_prob:

            # ADD IMMUNITY
            # TODO: Outer products may be in wrong order. Need to run tests.

            # Make a matrix where every H row is filled with a positive number if and only if the host has at
            # least one virus with that H protein.
            immune_H = np.outer(self.H, np.ones(NUM_N))

            # Same thing but fill each N column with positive numbers if the host has a virus with that N protein
            immune_N = np.outer(np.ones(NUM_H), self.N)

            # Add these together. All positive entries are the virus the host is immune to.
            immune = immune_N + immune_H
            immune = immune > 0  # Change to boolean matrix with True in entry i,j if the host is immune to HiNj

            # Invert this boolean matrix to figure out what viruses the host is susceptible to.
            self.susceptibility = np.invert(immune).astype(float)  # invert so that is a susceptibility matrix instead of an immunity one.

            # Host recovers from all viruses
            self.viruses = ZEROS


    def birth_death(self):
        """
        STAGE 4

        Agents have a chance to die. If an agent dies a new agent of the same species is born to take its place.

        Todo: we need to decide how having influenza affects this rates.
        (We can calculate the number of viruses via np.sum(self.viruses))
        """

        # If die just replace host with an empty one. Ie a new organism took the old one's place.
        if rng.random(1) < self.death_rate:
            self.viruses = ZEROS
            self.susceptibility = ONES


    def contacts(self):
        """Returns a list of other organism the host has contacted and got viruses from."""

        contacts = []

        num_contacts = int(self.model.len_hosts_0 * self.model.contact_rates[self.species_id][0])
        samp = list(rng.choice(self.model.hosts_0, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_1 * self.model.contact_rates[self.species_id][1])
        samp = list(rng.choice(self.model.hosts_1, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_2 * self.model.contact_rates[self.species_id][2])
        samp = list(rng.choice(self.model.hosts_2, num_contacts))
        contacts = contacts + samp

        num_contacts = int(self.model.len_hosts_3 * self.model.contact_rates[self.species_id][3])
        samp = list(rng.choice(self.model.hosts_3, num_contacts))
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
                 immigration_rate=0.02, contact_rates=None, fitness_on=True, init_hosts=True):
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
        self.schedule = StagedActivation(self, ["contract_virus", "recombine", "recover", "birth_death"])  # set schedule
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
        if contact_rates is None:
            self.contact_rates = np.array([[0.01, 0.0045, 0.002, 0.001],
                                           [0.0045, 0.0095, 0.003, 0.0045],
                                           [0.002, 0.003, 0.09, 0.003],
                                           [0.001, 0.0045, 0.003, 0.01]]) * 1
        else:
            self.contact_rates = contact_rates


        if init_hosts:
            # Initialize population
            # Make a bunch of random organisms for now
            for i in range(np.sum(init_pop_size)):

                species = random.choice(["Human", "Pig", "Bird", "Poultry"])  # Decide species with equal probability.
                init_viruses = np.random.choice([0, 1], size=(NUM_H, NUM_N), p=[.999, .001])  # Randomly decide some viruses it has
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
        if self.model_step >= 0:
            self.datacollector.collect(self)

        self.len_hosts_0 = len(self.hosts_0)
        self.len_hosts_1 = len(self.hosts_1)
        self.len_hosts_2 = len(self.hosts_2)
        self.len_hosts_3 = len(self.hosts_3)

        self.schedule.step()  # step all agents
        self.immigrate(self.immigration_rate)


    def immigrate(self, p):
        """
        A random virus immigrates into a random population with probability p.
        """

        pass

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
