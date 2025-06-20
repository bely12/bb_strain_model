import numpy as np
import random
import collections
from collections import Counter
import itertools
import math
import copy


class Host:

  def __init__(self, pop_size, species = False):

    if pop_size % 2 != 0:
      raise ValueError("The value of pop_size must be divisible by 2.")
    hosts = []
    existing_ids = []
    
    for i in range(pop_size): # assign a unique id
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)

      if species == True: # assign specific host species
        if i < pop_size//2:
          host_type = 'rodent'
        else:
          host_type = 'bird'
      else:
        host_type = '-'

      # create host
      hosts.append({'id': id, 
                    'infections': [], 
                    'host_type': host_type, 
                    'born': random.choice(range(1,76))})

    self.pop = hosts
    self.pop_size = len(self.pop)
    self.rodent_pop = [d for d in self.pop if d.get('host_type') == 'rodent']
    self.rodent_pop_size = len(self.rodent_pop)
    self.bird_pop = [d for d in self.pop if d.get('host_type') == 'bird']
    self.bird_pop_size = len(self.bird_pop)

  def refresh(self, day): # simulates host birth and death
    for host in self.pop:
      if host['born'] == day:
        host['infections'] = []

################################################################################

class Vector:

  def __init__(self, pop_size, loci, lam):

    if pop_size % 2 != 0:
      raise ValueError("The value of pop_size must be divisible by 2.")

    # define the pathogens and starting infection status
    self.strain_set = [''.join(bits) for bits in itertools.product('01', repeat= loci)]
    nymph_infections = np.random.poisson(lam, size=(pop_size // 2)) # poisson distribution for tick infection status
    
    # initialize the population of ticks
    tick_pop = []
    existing_ids = []
    lin_id = 1

    # create nymphs
    for i in range(len(nymph_infections)):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      if nymph_infections[i] == 0:
        tick_pop.append({'id': id, 'stage': 'nymph', 'strains': []})
      else:
        tick_pop.append({'id': id, 'stage': 'nymph', 'strains': [random.choice(self.strain_set)]}) 
        lin_id += 1

    # create larva
    for i in range(pop_size // 2):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      tick_pop.append({'id': id, 'stage': 'larva', 'strains': []})

    # for use in defining attributes
    nymph_pop = [d for d in tick_pop if d.get('stage') == 'nymph']

    strain_pop = []
    for d in nymph_pop:
      if d['strains'] != []:
        for i in range(len(d['strains'])):
          strain_pop.append(d['strains'][i])

    # define important stuff
    self.pop = tick_pop
    self.pop_size = len(self.pop)
    self.nymph_pop = nymph_pop # probably won't use this
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.year = 0
    self.poisson_infection_status = nymph_infections
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)

    # get strain frequncies
    counts = dict(Counter(strain_pop))
    # total_counts = sum(counts.values())
    # frequencies = {key: round((value / total_counts), 2) for key, value in counts.items()}
    # self.current_strain_count = len(frequencies) # to monitor number of strains in pop
    self.current_strain_count = 0
    for genotype in self.strain_set:
      if genotype not in counts:
        counts[genotype] = 0.0
      else:
        self.current_strain_count += 1
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts), 2) for key, value in counts.items()}
    self.strain_frequencies = frequencies.copy()

  
  ## function molts larva to nymphs and hatches new nymphs; current nymphs age out to adults and are not included in sim
  def update_pop(self):
    # remove old nymphs
    filtered_list = [d for d in self.pop if d.get('stage') != 'nymph']
    updated_pop = filtered_list

    # molt larva to nymphs
    for tick in updated_pop:
      if tick['stage'] == 'larva':
        tick['stage'] = 'nymph' 
    existing_ids = []
    for i in range(len(updated_pop)):
      existing_ids.append(updated_pop[i]['id'])

    # create new larva
    for i in range(self.pop_size // 2):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      updated_pop.append({'id': id,'stage': 'larva', 'strains': []})

    # redefine pop
    self.pop = updated_pop

    # for use in defining attributes
    nymph_pop = [d for d in self.pop if d.get('stage') == 'nymph']
    strain_pop = []
    for d in nymph_pop:
      if d['strains'] != []:
        for i in range(len(d['strains'])):
          strain_pop.append(d['strains'][i])

    # redefine attributes
    self.nymph_pop = nymph_pop 
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.year += 1

    # get strain frequncy data
    counts = dict(Counter(strain_pop))
    # total_counts = sum(counts.values())
    # frequencies = {key: round((value / total_counts), 2) for key, value in counts.items()}
    # self.current_strain_count = len(frequencies) # to monitor number of strains in pop
    self.current_strain_count = 0
    for genotype in self.strain_set:
      if genotype not in counts:
        counts[genotype] = 0.0
      else:
        self.current_strain_count += 1
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts), 2) for key, value in counts.items()}
    self.strain_frequencies = frequencies.copy() # master list of freqs even if 0 

  # function assigns a host for each tick to bite and the day/s it bites on
  def interaction(self, hosts):
    interaction_list = []

    # assign tick-host interactions
    for i in range(len(self.pop)):
      # larva bites
      if self.pop[i]['stage'] == 'larva':
        z = random.randint(1, len(hosts)-1)
        chosen_host = hosts[z]['id']
        bite_day = random.randint(75, 150)
      # nymph bites
      if self.pop[i]['stage'] == 'nymph':
        z = random.randint(1, len(hosts)-1)
        chosen_host = hosts[z]['id']
        bite_day = random.randint(1, 75)

      # add info to interaction list
      interaction_list.append({'tick_id': self.pop[i]['id'],
                              'host_id': chosen_host,
                              'bite_day': bite_day})
    return interaction_list


  def recombination(self, rate = 0.01):
    break_point = len(self.strain_set[0]) // 2
    for tick in self.pop:
      # only move forward if there are multiple strains being carried by the tick
      if tick['strains'] != [] and len(tick['strains']) > 1:
        strains = tick['strains']
      else:
        continue
      # cycle through each of the strains carried by the tick and decide if recombination will happen
      for j in range(len(strains)):
        if random.random() < rate:
          # if yes - define the receiving strain; make a temp copy of all strains minus the receiving strain and assign a donor
          strain = strains[j]
          strains_copy = copy.deepcopy(strains)
          strains_copy.remove(strain)
          donor = random.choice(strains_copy)
          
          # using a 50:50 chance of whether recomb is first half of seq or last; induce recombination
          if random.random() < 0.5:
            new_strain = strain[:break_point] + donor[break_point:]
          else:
            new_strain = donor[:break_point] + strain[break_point:]
          tick['strains'][j] = new_strain # recombined strain replaces parent strain
          
################################################################################

def tick2host_transmission(tick, host, cross_protection = 'high'):

  # create the transmission community from tick
  if tick['strains'] == []:
    tick_transmission_community = []
  if random.random() < 0.7 and tick['strains'] != []:
    number_of_strains = random.randint(1, len(tick['strains']))
    tick_transmission_community = random.sample(tick['strains'], number_of_strains)
  else:
    tick_transmission_community = []

  # transmission
  if tick_transmission_community != []:
    transmitted_strains = []

    # if the host has no infections
    if host['infections'] == []:
      for strain in tick_transmission_community:
        if random.random() < 0.7:
          transmitted_strains.append(strain)
      return transmitted_strains

    # if the host has infections
    if host['infections'] != []:
      seen_allele_1 = []
      seen_allele_2 = []
      for infection_strain in host['infections']:
        if infection_strain[0] in seen_allele_1:
          continue
        else:
          seen_allele_1.append(infection_strain[0])
        
        if infection_strain[1] in seen_allele_2:
          continue
        else:
          seen_allele_2.append(infection_strain[1])
      
      for strain in tick_transmission_community:
        matches = 0
        allele_1 = strain[0]
        allele_2 = strain[1]
        if allele_1 in seen_allele_1:
          matches += 1
        if allele_2 in seen_allele_2:
          matches += 1
      
      if matches == 2:
        fitness = 0
      if matches == 0: 
        fitness = 0.7 # transmission probability same as when there are no current infections

      if cross_protection == 'high' and matches == 1:
        fitness = 0.2 
      if cross_protection == 'medium' and matches == 1:
        fitness = 0.4
      if cross_protection == 'low' and matches == 1:
        fitness = 0.6
      if cross_protection == 'none' and matches == 1:
        fitness = 0.7
      
      if random.random() < fitness:
        transmitted_strains.append(strain)
    return transmitted_strains
  
  else:
    return []

def host2tick_transmission(host):
  # create transmission community from host
  if host['infections'] == []:
    return []

  if random.random() < 0.7 and host['infections'] != []:
    number_of_strains = random.randint(1, len(host['infections']))
    host_transmission_community = random.sample(host['infections'], number_of_strains)
  else:
    return []
  
  # transmission
  transmitted_strains = []
  for strain in host_transmission_community:
    if random.random() < 0.7:
      transmitted_strains.append(strain)
  return transmitted_strains
