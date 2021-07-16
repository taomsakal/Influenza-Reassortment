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
from parameters import transmission_fitness
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
        self.id = model.next_id()
        super().__init__(self.id, model)  # Initialize basic agent code, assign a unique id
        self.model = model
 
        self.species = species
        assert species in ["Human", "Pig", "Bird", "Poultry"]
 
        # assigns a species id to use to access values of contact rates adjacency matrix and 
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)

        #holder for viruses after contact before all individuals have contacted each other
        self.temp_viruses = set()

        # Right now the viruses that a host has are listed in a set.
        # todo: We'll explore using a matrix to store this info later to avoid loops.
        self.viruses = set()
        self.time_since_infection = 0
        if viruses is not None:
            self.viruses.add(viruses)
        self.h_immune = set()
        self.n_immune = set()
        
        for i in all_viruses:
             setattr(self, f"H{i[0]+1}N{i[1]+1}", int(i in self.viruses))

    def __eq__(self, other):
        return self.id == other.id
    
    def contract_virus(self):
        self.random_numbers = np.random.rand(1,1)[0,0]
        contacts = self.contacts()
        for contact in contacts:
            for virus in self.viruses:
                if (contact.is_infectable_by(virus)):
                    self.infection_rate = self.model.infection_rate
                    if (virus[0] in contact.h_immune and virus[1] in contact.n_immune):
                        self.infection_rate = 0.005
                    elif (virus[0] in contact.h_immune):
                        self.infection_rate = self.infection_rate * 0.5
                    elif (virus[1] in contact.n_immune):
                        self.infection_rate = self.infection_rate * 0.5
                    #self.infection_rate = self.infection_rate * transmission_fitness[self.species_id][virus[0]][virus[1]][contact.species_id]
                    if (np.random.rand(1,1)[0,0] < self.infection_rate): #
                        contact.temp_viruses.add(virus)
    
    def recombine(self):
        self.viruses = self.viruses.union(self.temp_viruses)
        self.temp_viruses = set()
        self.h = [item[0] for item in self.viruses]
        self.n = [item[1] for item in self.viruses]

        # takes a long time
        self.viruses= {tuple(x) for x in itertools.product(self.h,self.n) if self.is_infectable_by(x)}

    def recovery(self):
        if (len(self.viruses) > 0):
            self.time_since_infection = self.time_since_infection+1
            self.recovery_chance = 0.25 * 3**(self.time_since_infection-1) 
            if (np.random.rand(1,1)[0,0] <= self.recovery_chance):
              self.time_since_infection = 0
              self.h_immune = self.h_immune.union({item[0] for item in self.viruses})
              self.n_immune = self.n_immune.union({item[1] for item in self.viruses})
              self.viruses= set()

        if (len(self.h_immune)>0):
            self.h_immune = set([i for i in self.h_immune if (np.random.rand(1,1)[0,0] < 0.96)])

        if (len(self.n_immune)>0):
            self.n_immune = set([i for i in self.n_immune if (np.random.rand(1,1)[0,0] < 0.96)])

        for i in all_viruses:
            setattr(self, f"H{i[0]+1}N{i[1]+1}", int(i in self.viruses))
    
    def birth_death(self):
        if (np.random.rand(1,1)[0,0] < 0.069):
            host = Host(self.model, self.species)
            self.model.schedule.add(host)
            if (self.species_id == 0):
              self.model.hosts_0.append( host )
            elif (self.species_id == 1):
              self.model.hosts_1.append( host )
            elif (self.species_id == 2):
              self.model.hosts_2.append( host )
            else:
              self.model.hosts_3.append( host )
        if (np.random.rand(1,1)[0,0] < 0.06*(1.2**len(self.viruses))):
            if (self.species_id == 0):
              self.model.hosts_0.remove( self )
            elif (self.species_id == 1):
              self.model.hosts_1.remove( self )
            elif (self.species_id == 2):
              self.model.hosts_2.remove( self )
            else:
              self.model.hosts_3.remove( self )
            self.model.schedule.remove(self)
            self.viruses = set()
            del self

    def is_infectable_by(self, virus):
        """Returns true if the virus can infect the host."""
        return infection_table[self.species_id][virus[0]][virus[1]]
 
    def contacts(self):
        contacts = []
        """Returns a list of other organism the host has contacted and got viruses from."""
        for i in range(4): 
            #contact_rate = self.model.contact_rates[self.species_id][i]
            num_contacts = int(len(eval(f"self.model.hosts_{i}")) * self.model.contact_rates[self.species_id][i][0])
            samp = rand.sample(eval(f"self.model.hosts_{i}"), num_contacts)
            contacts = contacts + samp
        return contacts
 
 
class VirusModel(Model):
 
    def __init__(self, run="NA", init_pop_size=[1350, 900, 1200, 1500], x=0):
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
        self.infection_rate = 0.12 # infection rate

        # Population sizes
        self.human_pop_size = init_pop_size[0]
        self.pig_pop_size = init_pop_size[1]
        self.bird_pop_size = init_pop_size[2]
        self.poultry_pop_size = init_pop_size[3]
        self.total_pop_size = sum(init_pop_size)
        self.hosts_0 = []
        self.hosts_1= []
        self.hosts_2= []
        self.hosts_3= []
 
        # Adjacency matrix of gaussian contact rate distributions where entry ij is the contact rate species j to species i.
        # 1 -> humans
        # 2 -> pigs
        # 3 -> birds
        # 4 -> poultry
        # todo: put more reasonable values
        self.contact_rates = np.array([[[0.0065], [0.0044], [0.0044], [0.0022]],
                                       [[0.0033], [0.0083], [0.0022], [0.0033]],
                                       [[0.0033], [0.0033], [0.0117], [0.0044]],
                                       [[0.0044], [0.0033], [0.0033], [0.0107]]])
 
        # initialize population
        init_virus = None
        id = 0
        self.all_viruses = list(all_viruses)
        for i in range(self.human_pop_size):
            if (id % 30== 0):
              init_virus = self.all_viruses[int(id/30) % len(self.all_viruses)]
            id = id +1
            #if (randint(0,44) == 0):
            #    init_virus = choice([(0,0),(1,1),(2,1),(4,0),(6,1),(6,2),(6,6),(8,1),(9,6)])
            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_0.append(host)
        for i in range(self.pig_pop_size):
            if (id % 30== 0):
              init_virus = self.all_viruses[int(id/30) % len(self.all_viruses)]
            id = id +1
            #if (randint(0,44) == 0):
            #    init_virus = choice([(0,0),(0,1),(1,2),(2,1),(2,2),(3,5),(4,1),(8,1)])
            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_1.append(host)
        for i in range(self.bird_pop_size):
            if (id % 20== 0):
              init_virus = self.all_viruses[int(id/20) % len(self.all_viruses)]
            id = id +1
            #if (randint(0,44) == 0):
            #    init_virus = (randint(0,15), randint(0,8))
            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_2.append(host)
        for i in range(self.poultry_pop_size):
            if (id % 20 == 0):
              init_virus = self.all_viruses[int(id/20) % len(self.all_viruses)]
            id = id +1
            #if (randint(0,44) == 0):
            #    init_virus = (randint(0,12), randint(0,8))
            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            self.hosts_3.append(host)
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