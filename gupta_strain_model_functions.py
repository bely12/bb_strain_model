import numpy as np
import random
import collections
from collections import Counter
import itertools
import math
import copy


class Host:

  def __init__(self, pop_size):

    hosts = []
    existing_ids = []
    
    for i in range(pop_size): # assign a unique id
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)

      # create host
      hosts.append({'id': id, 
                    'infections': [],
                    'born': random.choice(range(1,76))})

    # define some attributes
    self.pop = hosts
    self.pop_size = len(self.pop)

  def refresh(self, day): # simulates host birth and death
    for host in self.pop:
      if host['born'] == day:
        host['infections'] = []

################################################################################

class Pathogen:

  def __init__(self):

    # define pathogen alleles for 10 loci 
    allele_A = ['A', 'a']
    allele_B = ['B', 'b']
    allele_C = ['C', 'c']
    allele_D = ['D', 'd']
    allele_E = ['E', 'e']
    allele_F = ['F', 'f']
    allele_G = ['G', 'g']
    allele_H = ['H', 'h']
    allele_I = ['I', 'i']
    allele_J = ['J', 'j']

    allele_sets = {'A': allele_A,
                  'B': allele_B,
                  'C': allele_C,
                  'D': allele_D,
                  'E': allele_E,
                  'F': allele_F,
                  'G': allele_G,
                  'H': allele_H,
                  'I': allele_I,
                  'J': allele_J}

    self.genotypes = [''.join(allele) for allele in itertools.product(*allele_sets)]
  
class Vector:

  def __init__(self, pop_size, pathogens, lam):

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
        starter_infections = random.sample(pathogens, k=nymph_infections[i])
        tick_pop.append({'id': id, 'stage': 'nymph', 'strains': starter_infections}) 
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
    self.nymph_pop = nymph_pop
    self.year = 0
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)

    # get strain frequncies
    counts = dict(Counter(strain_pop))
    self.current_strain_count = 0
    for genotype in pathogens:
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


  def recombination(self, lam = 0.01):
    for tick in self.pop:
      # only move forward if there are multiple strains being carried by the tick
      if tick['strains'] != [] and len(tick['strains']) > 1:
        strains = tick['strains']
      else:
        continue

      for strain in strains: # cycle through variants
        recombination_status = (np.random.poisson(lam, 10) > 0).astype(int) # determine recombination for each loci
          
        if (recombination_status == 1).any(): # check to make sure at least one recombination event will happen
          strains_copy = copy.deepcopy(strains)
          strains_copy.remove(strain)
          donor = random.choice(strains_copy)
          
          loci_list = list(strain) # convert string to a list so I can mutate it
          loci_index = 0
          for decision in recombination_status:
            if decision == 1:
              loci_list[loci_index] = donor[loci_index]
            loci_index +=1
          
          new_strain = "".join(loci_list)
          if new_strain not in strains:
            strains[strains.index(strain)] = new_strain
          
################################################################################

def tick2host_transmission(tick, host):

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
      for transmission_strain in tick_transmission_community:
        most_similar = []
        for host_strain in host['infections']:
          matches = sum(loci_host == loci_vector for loci_host, loci_vector in zip(host_strain, transmission_strain))
          most_similar.append(matches)
        immune_val = max(most_similar)
        fitness = (10-immune_val) * 0.1
        if random.random() < fitness:
          transmitted_strains.append(transmission_strain)
    
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
