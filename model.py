"""
The agent and simulation class live here.
"""
 
import random as rand
import numpy as np
import itertools
from mesa import Agent, Model
from mesa.time import StagedActivation
from mesa.datacollection import DataCollector
from parameters import infection_table
from parameters import viruses as all_viruses
import json
import math

rng = np.random.default_rng(seed=2021)


species_infected = np.sum(infection_table,axis = 0) 
fitness = (np.ones((17,10))/species_infected) * rng.normal(1, 0.2, (17, 10))
class Host(Agent):
  
    def __init__(self, model, species, viruses=None):
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
        self.it = self.model.it #sets iteration from model
        self.species = species
        self.virus_list = []
        assert species in ["Human", "Pig", "Bird", "Poultry"]
        self.death = False
 
        # assigns a species id to use to access values of contact rates adjacency matrix and 
        self.species_id = ["Human", "Pig", "Bird", "Poultry"].index(self.species)

        #holder for viruses after contact before all individuals have contacted each other
        self.temp_viruses = np.array([])
        self.viruses = viruses * infection_table[self.species_id]
        # Right now the viruses that a host has are listed in a set.
        # todo: We'll explore using a matrix to store this info later to avoid loops.
        
        self.virus_list = [tuple(x) for x in np.asarray(np.where(self.viruses==1)).T]
         
        #sets variables for agent reporters to measure each virus
        for i in all_viruses:
            setattr(self, f"H{i[0]+1}N{i[1]+1}", int(tuple(i) in self.virus_list))

    def __eq__(self, other):
        return self.id == other.id
    
    def contract_virus(self):
        """
        Spreads viruses to other agents. Stage 1 of each step.
        """

        #gets list of contacts
        if (not (self.death)):
            contacts = self.contacts()
            self.base_rate = self.model.infection_rate

            for contact in contacts:
                viruses_transferred = np.array(self.viruses) * infection_table[contact.species_id] * self.base_rate
                viruses_transferred[:,0] = 0
                viruses_transferred[0,:] = 0
                if(self.model.fitness_on):
                    viruses_transferred = viruses_transferred * fitness

                for i in range (16):
                    if (contact.viruses[i + 1, 0] == 1):
                        viruses_transferred[i+1, 1:10] = viruses_transferred[i+1, 1:10] * 0.5
                for i in range (9):
                    if (contact.viruses[0, i+1] == 1):
                        viruses_transferred[1:16, i+1] = viruses_transferred[1:16, i+1] * 0.5

                transfer = np.vectorize(self.transfer)
                contact.temp_viruses = transfer(viruses_transferred)


                
    def transfer(self, x):
        return rng.choice(a=[1,0],p=[x,1-x])
    
    def reduce_to_one(self, x):
        if (int(x)>1):
            return(1)
        else:
            return(0)
    
    def recombine(self):
        """
        Viruses recombine. Stage 2 of each step.
        """
        if (not (self.death)):
            if (self.temp_viruses != np.array([])):
                reduce_to_one = np.vectorize(self.reduce_to_one)
                self.viruses = reduce_to_one(self.viruses + self.temp_viruses)
                self.temp_viruses = np.array([])
            virus_list = np.where(self.viruses[1:16,1:9]==1)
            self.h = set(virus_list[0]+1)
            self.n = set(virus_list[1]+1)

            for x in self.h:
                for y in self.n:
                    self.viruses[x,y] = 1
            if (self.death):
                self.model.schedule.remove(self)
             
    def birth_death(self):
        """
        Agents have a chance of birth and death. Stage 4 of each step.
        """

        self.virus_list = [tuple(x) for x in np.asarray(np.where(self.viruses==1)).T]
        
        for i in all_viruses:
            setattr(self, f"H{i[0]+1}N{i[1]+1}", int(tuple(i) in self.virus_list))
             
        if (self.model.iteration>0):
            
            #Sets species specific mortality rates (birds and poultry less liekly to die from Influenza)
            if (self.species_id == 0):
                mortality = 1.25
            elif (self.species_id == 1):
                mortality = 1.25
            elif (self.species_id == 2):
                mortality = 1.005
            else:
                mortality = 1.005
            
            #Chance of death increases with number of viruses the host is infected by
            if (self.rng.random(1)[0] < self.model.death_rate*(mortality**len(self.virus_list))):
                if (self.species == "Human"):
                    self.model.hosts_0.remove( self )
                elif (self.species_id == 1):
                    self.model.hosts_1.remove( self )
                elif (self.species_id == 2):
                    self.model.hosts_2.remove( self )
                else:
                    self.model.hosts_3.remove( self )
                self.viruses = self.viruses * 0
                self.death = True
                

    def contacts(self):
        """Returns a list of other organism the host has contacted and got viruses from."""
        contacts = []
        num_contacts = int(len(self.model.hosts_0) * self.model.contact_rates[self.species_id][0])
        samp = list(self.rng.choice(self.model.hosts_0, num_contacts))
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_1) * self.model.contact_rates[self.species_id][1])
        samp = list(self.rng.choice(self.model.hosts_1, num_contacts))
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_2) * self.model.contact_rates[self.species_id][2])
        samp = list(self.rng.choice(self.model.hosts_2, num_contacts))
        contacts = contacts + samp
        num_contacts = int(len(self.model.hosts_3) * self.model.contact_rates[self.species_id][3])
        samp = list(self.rng.choice(self.model.hosts_3, num_contacts))
        contacts = contacts + samp
        return contacts
 
 
class VirusModel(Model):
 
    def __init__(self, run="NA", init_pop_size=[900, 650, 1000, 750], it=0, infection_rate= 0.25, recovery_rate = 0.2, mutation_rate = 0.23, birth_rate = 0.04, death_rate = 0.03, cross_immunity_effect = 0.05, init_viruses=None, immigration_rate=0.02, contact_rates = None, fitness_on=True):
        """
        Args:
            init_pop_size: The initial population size of each species [Humans, Pigs, Birds, Poultry]
            x: Batch runner throws an error without a dummy variable to use as a variable parameter, therefore this variable acts as a dummy
            it: Iteration number
        """

        rng = np.random.default_rng(seed=2021)
        super().__init__()  # Initialize basic agent code, assign a unique id
        self.it = it
        self.run = run
        self.running = True  # For batch runs
        self.iteration = 0  # The number of timesteps the simulation has run
        self.schedule = StagedActivation(self, ["birth_death", "contract_virus", "recombine"], False, False)  # set schedule 
        self.infection_rate = infection_rate # infection rate
        self.recovery_rate = recovery_rate
        self.mutation_rate = mutation_rate
        self.birth_rate = birth_rate
        self.death_rate = death_rate
        self.cross_immunity_effect = cross_immunity_effect
        self.init_viruses = init_viruses
        self.immigration_rate = immigration_rate
        self.fitness_on = fitness_on

        
        # Population sizes
        self.human_pop_size = init_pop_size[0]
        self.pig_pop_size = init_pop_size[1]
        self.bird_pop_size = init_pop_size[2]
        self.poultry_pop_size = init_pop_size[3]
        self.total_pop_size = sum(init_pop_size)

        #lists with agents of each species. Used to get contacts
        self.hosts_0 = []
        self.hosts_1= []
        self.hosts_2= []
        self.hosts_3= []

        # It's a weird naming convention but self.iteration actually measures 
        # the current step the model is on, not the current iteration
        self.iteration = 0
 
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

        init_virus = None
        self.all_viruses = list(all_viruses)
        
        for i in range(self.human_pop_size):

            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_0.append(host)
        for i in range(self.pig_pop_size):
            
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_1.append(host)
        for i in range(self.bird_pop_size):
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_2.append(host)
        for i in range(self.poultry_pop_size):
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            self.hosts_3.append(host)
            init_virus = None 
        


        #sets reporters for each virus
        reporters = "{\"Iteration\":\"it\",\"Species\":\"species\","
        for i in all_viruses:
            name = f"H{i[0]+1}N{i[1]+1}"
            reporters = reporters + "\"" + name + "\"" + ":" + "\"" + name + "\"" + ","
        reporters = reporters[:-1] + "}"
        self.datacollector = DataCollector(
            #model_reporters={"Strain_data": Count_Strains}
            agent_reporters = json.loads(reporters)
        )

    def step(self):

        rng = np.random.default_rng(seed=2021)
        """Steps the entire model one time step."""
        self.iteration += 1

        #collects data after 350 steps
        if(self.iteration > 0):
            self.datacollector.collect(self)
        
        self.schedule.step()  # step all agents
        print("0")
        #recovery
        recovering = rng.choice(np.array(self.schedule.agents), int(len(self.schedule.agents)*self.recovery_rate))
        for agent in recovering:
            virus_list = np.where(agent.viruses[1:16,1:9]==1)
            h = set(virus_list[0]+1)
            n = set(virus_list[1]+1)
            for i in h:
                agent.viruses[i,0] = 1
            for i in n:
                agent.viruses[0,i] = 1
            agent.viruses[1:16,1:9] = 0
        
        losing_immunity = rng.choice(np.array(self.schedule.agents), int(len(self.schedule.agents)*self.mutation_rate))
        h = np.ones(17)
        h[:int(17*self.mutation_rate)] = 1
        rng.shuffle(h)
        n = np.ones(10)
        n[:int(10*self.mutation_rate)] = 1
        rng.shuffle(n)
        loss_array = np.ones((17,10))
        loss_array[:, 0] = h
        loss_array[0, :] = n

        for agent in losing_immunity:
            agent.viruses = agent.viruses * loss_array
        


        for i in range(int(len(self.hosts_0) * self.birth_rate)):
            init_virus = np.zeros(170)
            init_virus = np.reshape(init_virus, (17,10))
            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_0.append(host)
        for i in range(int(len(self.hosts_1) * self.birth_rate)):
            init_virus = np.zeros(170)
            init_virus = np.reshape(init_virus, (17,10))
            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_1.append(host)
        for i in range(int(len(self.hosts_2) * self.birth_rate)):
            init_virus = np.zeros(170)
            init_virus = np.reshape(init_virus, (17,10))
            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_2.append(host)
        for i in range(int(len(self.hosts_3) * self.birth_rate)):
            init_virus = np.zeros(170)
            init_virus = np.reshape(init_virus, (17,10))
            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            self.hosts_3.append(host)
            init_virus = None 



        for i in range(int(len(self.hosts_0) * self.immigration_rate)):

            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Human", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_0.append(host)
        for i in range(int(len(self.hosts_1) * self.immigration_rate)):
            
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Pig",init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_1.append(host)
        for i in range(int(len(self.hosts_2) * self.immigration_rate)):
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Bird", init_virus)
            self.schedule.add(host)
            init_virus = None
            self.hosts_2.append(host)
        for i in range(int(len(self.hosts_3) * self.immigration_rate)):
            init_virus = np.zeros(170)
            init_virus[:34] = 1
            rng.shuffle(init_virus)
            init_virus = np.reshape(init_virus, (17,10))

            host = Host(self, "Poultry", init_virus)
            self.schedule.add(host)
            self.hosts_3.append(host)
            init_virus = None 


if __name__ == '__main__':
    pass