"""
The agent and simulation class live here.
"""
 
from random import *
import numpy as np
import itertools
from mesa import Agent, Model
from mesa.time import StagedActivation
from mesa.datacollection import DataCollector
from parameters import infection_table
from parameters import viruses as all_viruses
import json

 
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
        self.time_since_infection = 0
        if viruses is not None:
            self.viruses.add(viruses)
        self.h_immune = set()
        self.n_immune = set()
        
        for i in all_viruses:
            exec("%s = %d" % ("self." + f"H{i[0]+1}N{i[1]+1}", int(i in self.viruses)))
    
    def contract_virus(self):
        contacts = self.contacts()
        self.temp_viruses = []
        for contact in contacts:
            for virus in contact.viruses:
                if (self.is_infectable_by(virus)):
                    self.infection_rate = self.model.infection_rate
                    if (virus[0] in self.h_immune):
                        self.infection_rate = self.infection_rate * 0.5
                    if (virus[1] in self.n_immune):
                        self.infection_rate = self.infection_rate * 0.5
                    if (rand.random() <= self.infection_rate): #
                        self.temp_viruses.append(virus)
    
    def recombine(self):
        self.viruses = self.viruses.union(set(self.temp_viruses))
        self.h = [item[0] for item in self.viruses]
        self.n = [item[1] for item in self.viruses]
        self.viruses= set([tuple(x) for x in itertools.product(self.h,self.n) if self.is_infectable_by(x)])

    def recovery(self):
        if (len(self.viruses) > 0):
            self.time_since_infection = self.time_since_infection+1
            self.recovery_chance = 0.55 * 0.95**(len(self.viruses)-1) * 1.3**(self.time_since_infection-1)
            if (rand.random() <= self.recovery_chance):
              self.time_since_infection = 0
              self.h_immune = self.h_immune.union(set([item[0] for item in self.viruses]))
              self.n_immune = self.n_immune.union(set([item[1] for item in self.viruses]))
              self.viruses= set()

        if (len(self.h_immune)>0):
            self.h_immune = set([i for i in self.h_immune if (randint(0,19) < 18)])

        if (len(self.n_immune)>0):
            self.n_immune = set([i for i in self.n_immune if (randint(0,19) < 18)])

        for i in all_viruses:
            exec("%s = %d" % ("self."+ f"H{i[0]+1}N{i[1]+1}", int(i in self.viruses)))
    
    def birth_death(self):
        if (randint(0,19) > 16):
            host = Host(self.model, self.species)
            self.model.schedule.add(host)
        if (rand.random() < (0.1 * 1.2**len(self.viruses))):
            self.model.schedule.remove(self)
            self.viruses = set()

    def is_infectable_by(self, virus):
        """Returns true if the virus can infect the host."""
        return infection_table[self.species_id][virus[0]][virus[1]]
 
    def contacts(self):
        contacts = []
        """Returns a list of other organism the host has contacted and got viruses from."""
        for i in range(4): 
            #contact_rate = self.model.contact_rates[self.species_id][i]
            num_contacts = int(self.model.contact_rates[self.species_id][i][0])
            #if (num_contacts < 0):
            #    num_contacts=0
            samp = sample([x for x in self.model.schedule.agents if x.species_id==i], num_contacts)
            contacts = contacts + samp
        return contacts
 
 
class VirusModel(Model):
 
    def __init__(self, run="NA", init_pop_size=[700, 200, 500, 1000], x=0):
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
        self.schedule = StagedActivation(self, ["contract_virus", "recombine", "recovery", "birth_death"], True, True)  # set schedule
        self.infection_rate = 0.11 # infection rate

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
        self.contact_rates = np.array([[[5], [3], [3], [1]],
                                       [[1], [10], [1], [2]],
                                       [[3], [1], [20], [3]],
                                       [[1], [3], [3], [15]]])
 
        # initialize population
        init_virus = None
        for i in range(self.human_pop_size):
            if (randint(0,44) == 0):
                init_virus = choice([(0,0),(1,1),(2,1),(4,0),(6,1),(6,2),(6,6),(8,1),(9,6)])
            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.pig_pop_size):
            if (randint(0,44) == 0):
                init_virus = choice([(0,0),(0,1),(1,2),(2,1),(2,2),(3,5),(4,1),(8,1)])
            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.bird_pop_size):
            if (randint(0,44) == 0):
                init_virus = (randint(0,15), randint(0,8))
            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
        for i in range(self.poultry_pop_size):
            if (randint(0,44) == 0):
                init_virus = (randint(0,12), randint(0,8))
            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            init_virus = None
        

        reporters = "{\"Species\":\"species\","
        for i in all_viruses:
            name = f"H{i[0]+1}N{i[1]+1}"
            reporters = reporters + "\"" + name + "\"" + ":" + "\"" + name + "\"" + ","
        reporters = reporters[:-1] + "}"
        self.datacollector = DataCollector(
            #model_reporters={"Strain_data": Count_Strains}
            agent_reporters = json.loads(reporters) 
        )  

    def step(self):
        """Steps the entire model one time step."""
 
        self.iteration += 1
        self.datacollector.collect(self)
        self.schedule.step()  # step all agents
 
 
if __name__ == '__main__':
    pass