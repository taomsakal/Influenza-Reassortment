"""
The agent and simulation class live here.
"""
 
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

 
def Count_Strains(model):
    
    strains = []
    count = {}
    for agent in model.schedule.agents:
        for virus in agent.viruses:
            name = "H"+str(virus[0]+1)+"N"+str(virus[1]+1)
            count[name] = count.get(name, 0) + 1
    return count
    
    #return(1)

class Host(Agent):
 
    def __init__(self, model, species, viruses=None):
        """
        Creates a new host organism
        Args:
            model: The model the host is a part of
            species: The species the host is. Either "Human", "Pig", "Bird", or "Poultry".
            viruses: The set of viruses currently infecting the host.
        """
 
        super().__init__(model.next_id(), model)  # Initialize basic agent code, assign a unique id
        self.model = model
 
        self.species = species
        assert species in ["Human", "Pig", "Bird", "Poultry"]

        # assigns a species id to use to access values of contact rates adjacency matrix and 
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)

        #holder for viruses after contact before all individuals have contacted each other
        self.temp_viruses = []

        # Right now the viruses that a host has are listed in a set.
        # todo: We'll explore using a matrix to store this info later to avoid loops.
        self.viruses = set()
        if viruses is not None:
            self.viruses.add(viruses)
 
    
    def contract_virus(self):
        contacts = self.contacts
        for contact in contacts():
            for virus in contact.viruses:
                if (self.is_infectable_by(virus)):
                    self.temp_viruses.append(virus)
    
    def recombine(self):
        self.viruses = self.viruses.union(set(self.temp_viruses))
        self.temp_viruses = []
        h = [item[0] for item in self.viruses]
        n = [item[1] for item in self.viruses]
        self.viruses= set(itertools.product(h,n))
    
    # todo: def birth_death()

    def is_infectable_by(self, virus):
        """Returns true if the virus can infect the host."""
        return infection_table[self.species_id][virus[0]][virus[1]]
 
    def contacts(self):
        contacts = []
        """Returns a list of other organism the host has contacted and got viruses from."""
        for i in range(4): 
            contact_rate = self.model.contact_rates[self.species_id][i]
            num_contacts = round(gauss(contact_rate[0], contact_rate[1]))
            if (num_contacts < 0):
                num_contacts=0
            contacts = contacts + sample([x for x in self.model.schedule.agents if x.species_id==i], num_contacts)

        return contacts
 
 
class VirusModel(Model):
 
    def __init__(self, run="NA", init_pop_size=[100, 500, 200, 400], x=0):
        """
 
        Args:
            run:
            init_pop_size: The initial population size of each species [Humans, Pigs, Birds, Poultry]
            x: Batch runner throws an error without a dummy variable to use as a variable parameter
        """
 
        super().__init__()  # Initialize basic agent code, assign a unique id
        self.x = x
        self.run = run
        self.running = True  # For batch runs
        self.iteration = 0  # The number of timesteps the simulation has run
        self.schedule = StagedActivation(self, ["contract_virus", "recombine"], True, True)  # set schedule
 
        # Population sizes
        self.human_pop_size = init_pop_size[0]
        self.pig_pop_size = init_pop_size[1]
        self.bird_pop_size = init_pop_size[2]
        self.poultry_pop_size = init_pop_size[3]
        self.total_pop_size = sum(init_pop_size)
 
        # Adjacency matrix of gaussian contact rate distributions where entry ij is the contact rate species j to species i.
        # 1 -> humans
        # 2 -> pigs
        # 3 -> poultry
        # 4 -> wild birds
        # todo: put more reasonable values and sd
        self.contact_rates = np.array([[[4,1], [2,1], [2,1], [1,1]],
                                       [[2,1], [9,1], [1,1], [2,1]],
                                       [[2,1], [1,1], [8,1], [2,1]],
                                       [[1,1], [2,1], [2,1], [5,1]]])
 
        # initialize population
        init_virus = None
        for i in range(self.human_pop_size):
            if (randint(0,50) == 50):
                init_virus = choice([(0,0),(1,1),(2,1),(4,0),(6,1),(6,2),(6,6),(8,1),(9,6)])
            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.pig_pop_size):
            if (randint(0,50) == 50):
                init_virus = choice([(0,0),(0,1),(1,2),(2,1),(2,2),(3,5),(4,1),(8,1)])
            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.bird_pop_size):
            if (randint(0,50) == 50):
                init_virus = (randint(0,15), randint(0,8))
            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.poultry_pop_size):
            if (randint(0,50) == 50):
                init_virus = (randint(0,12), randint(0,8))
            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            init_virus = None
        
        self.datacollector = DataCollector(
            model_reporters={"Strain_data": Count_Strains})

    def step(self):
        """Steps the entire model one time step."""
 
        self.iteration += 1
        self.datacollector.collect(self)
        self.schedule.step()  # step all agents
 
        # todo: add data gathering
 
 
if __name__ == '__main__':
    pass
