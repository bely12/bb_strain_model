# vector-host-pathogen strain structure model functions and classes
import numpy as np
import random
import collections
from collections import Counter
import itertools
import math
import copy

# calculate hamming distance for 2 strains
def hamming_distance(string1, string2):
  distance = 0
  for i in range(len(string1)):
      if string1[i] != string2[i]:
          distance += 1
  return distance

# assigns a probability of transmission based on antigenic distances; output is a list of dictionaries [{distance: N, probability: n},..]
def distance_probabilities(values, curve, prob_min=0.0, prob_max=1.0):
  if curve == 'sigmoid':
    midpoint = (max(values) + min(values)) / 2
    def sigmoid(x):
      return round(prob_min + (prob_max - prob_min) / (1 + math.exp(-1.0 * (x - midpoint))),4)
    probabilities = [sigmoid(val) for val in values]
  if curve == 'linear':
    max_val = max(values)
    min_val = min(values)
    probabilities = [(prob_min + (prob_max - prob_min) * (val - min_val) / max_val) for val in values]
  return dict(zip(values, probabilities))

# exponential function for determining transmission probability based on MNP scheme; 0 = bird specialist, 1 = rodent specialist, 0.5 = generalist
def mnp_probabilities(strains):
  mnp_values = dict(zip(strains, np.round(np.random.uniform(0.0, 1.0, len(strains)),4)))
  res = []
  for strain in mnp_values:
    res.append({'strain': strain, 'mnp_value': mnp_values[strain], 'rodent': round((mnp_values[strain]**2),4), 'bird': round((1-mnp_values[strain])**2,4)})
  return res


# class for managing host population
class Host:

  def __init__(self, pop_size, MNP = False):

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

      if MNP == True: # assign specific host species
        if i < pop_size//2:
          host_type = 'rodent'
        else:
          host_type = 'bird'
      else:
        host_type = '-'

      hosts.append({'id': id, 'infection': [], 'host_type': host_type}) # create host

    self.pop = hosts
    self.pop_size = len(self.pop)
    self.rodent_pop = [d for d in self.pop if d.get('host_type') == 'rodent']
    self.rodent_pop_size = len(self.rodent_pop)
    self.bird_pop = [d for d in self.pop if d.get('host_type') == 'bird']
    self.bird_pop_size = len(self.bird_pop)

# class for managing tick population while tracking borrelia population data
class Vector:

  def __init__(self, pop_size, n_strains, strain_length, lam):

    if pop_size % 2 != 0:
      raise ValueError("The value of pop_size must be divisible by 2.")

    # define the pathogens and starting infection status
    self.strain_set = [''.join(bits) for bits in itertools.product('01', repeat= strain_length)]
    self.current_strains = random.sample(self.strain_set, n_strains)
    nymph_infections = np.random.poisson(lam, size=(pop_size // 2)) # poisson distribution for tick infection status

    # initialize the population of ticks
    tick_pop = []
    existing_ids = []

    # create nymphs
    for i in range(len(nymph_infections)):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      tick_pop.append({'id': id, 'stage': 'nymph', 'strains': random.sample(self.current_strains, min(nymph_infections[i], len(self.current_strains)))})

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
      strain_pop.extend(d['strains'])

    # define important stuff
    self.pop = tick_pop
    self.current_strain_count = len(self.current_strains)
    self.nymph_pop = nymph_pop # probably won't use this
    self.avg_carried = round(np.mean([len(tick['strains']) for tick in self.nymph_pop if len(tick['strains']) > 0]),4)
    self.pop_size = len(self.pop)
    self.year = 0
    self.poisson_infection_status = nymph_infections
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.diversity = round(1 - (sum((count / len(strain_pop)) ** 2 for count in Counter(strain_pop).values())), 3)

    # get strain frequncy data
    totals = [item for d in self.pop for item in d['strains']]
    counts = dict(Counter(totals))

    # calculate frequencies for current strains, define as attribute
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts) * 100, 2) for key, value in counts.items()}
    self.current_strain_frequencies = frequencies.copy()

    # add strains not present and assign freq as 0, define as attribute
    for item in self.strain_set:
      if item not in frequencies:
        frequencies[item] = 0.0
    self.strain_set_frequencies = frequencies

    self.avg_gen_dist = [{'year': self.year, 'avg_dist': 0.0}]
    self.antigen_distances = [] # create empty list upon intitializing pop, will append with ant dists during sim

  ## function molts larva to nymphs and hatches new nymphs; current nymphs age out to adults and are not included in sim
  def update_pop(self):
    # remove old nymphs
    filtered_list = [d for d in self.pop if d.get('stage') != 'nymph']
    updated_pop = filtered_list

    # molt larva to nymphs
    for tick in updated_pop:
      if tick['stage'] == 'larva':
        tick['stage'] = 'nymph' ##### I can probably just use this without the conditional statement above it #####
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
      strain_pop.extend(d['strains'])

    # redefine attributes
    self.nymph_pop = nymph_pop # I probably don't need this but I'll leave it here for now
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.avg_carried = round(np.mean([len(tick['strains']) for tick in self.nymph_pop if len(tick['strains']) > 0]),4)
    self.diversity = round(1 - (sum((count / len(strain_pop)) ** 2 for count in Counter(strain_pop).values())),3)
    self.year += 1

    # get strain frequncy data
    totals = [item for d in self.pop for item in d['strains']]
    counts = dict(Counter(totals))

    # calculate frequencies for current strains, define as attribute
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts) * 100, 2) for key, value in counts.items()}

    self.current_strains = list(frequencies.keys())
    self.current_strain_count = len(self.current_strains)
    self.current_strain_frequencies = frequencies.copy()

    # add strains not present and assign freq as 0, define as attribute
    for item in self.strain_set:
      if item not in frequencies:
        frequencies[item] = 0.0
    self.strain_set_frequencies = frequencies

    totals = [item for d in self.pop for item in d['strains']]
    generation_distances = []
    for i in range(len(totals)):
      for j in range(len(totals)):
        if i != j:
          dist = hamming_distance(totals[i], totals[j])
          self.antigen_distances.append(dist)
          generation_distances.append(dist)
    self.avg_gen_dist = round(np.mean(generation_distances),3)

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

  def mnp_recs(self, values):
    mnp_vals = []
    for strain in self.current_strains:
      mnp_vals.append(values.get(strain)['mnp_value'])
    return mnp_vals
  
  def recombination(self, rate = 0.5):
    break_point = len(self.strain_set[0]) // 2
    for i in range(len(self.pop)):
      strains = self.pop[i]['strains']
      if len(strains) <= 1:
        continue
      else:
        for j in range(len(strains)):
          if random.random() < rate:
            strain = strains[j]
            strains_copy = copy.deepcopy(strains)
            strains_copy.remove(strain)
            donor = random.choice(strains_copy)
            if random.random() < 0.5:
              new_strain = strain[:break_point] + donor[break_point:]
            else:
              new_strain = strain[break_point:] + donor[:break_point]

            self.pop[i]['strains'][j] = new_strain

  def mutate(self, rate = 0.01):
    for i in range(len(self.pop)):
      temp_strain_set = copy.deepcopy(self.pop[i]['strains'])
      for j in range(len(temp_strain_set)):
        mut_pois_dist = np.random.poisson(rate, len(temp_strain_set[j]))
        new_string = []
        for k in range(len(temp_strain_set[j])):
          if mut_pois_dist[k] > 0:
            new_string.append(str(random.randint(0,1)))
          else:
            new_string.append(temp_strain_set[j][k]) #if no, select same bit
        mutated_string = ''.join(new_string) #join bits together
        temp_strain_set[j] = mutated_string #change in the temp list
      self.pop[i]['strains'] = temp_strain_set #update the pop list

# class for managing Borrelia populations during infection
class Pathogen:
  def __init__(self, tick_strains ,pop_size = 100):
    # create transmission community
    if tick_strains == []:
      transmission_community = []

    # decide if transmission community forms, create transmission community
    if random.random() < 0.7 and tick_strains != []:
      number_of_strains = random.randint(1, len(tick_strains))
      transmission_community = random.sample(tick_strains, number_of_strains)
    else:
      transmission_community = []

    # create the populations for each pathogen strain in transmission community
    bb_pop = []
    for i in transmission_community:
      temp_pop = []
      for n in range(pop_size):
        temp_pop.append(i)
      bb_pop.append({'starting_strain': i,'strain_pop': temp_pop})

    self.pop = bb_pop

  # this function determines which strains will be transferred from the vector to the host and returns a list of those strains
  def tick2host_transmission(self, host_strains, values, host_type=None, cross_reactivity = False, MNP = False, uniform_fitness = False):
    # grab most dominant variant from each strain community, all all variants if no selection
    winning_strains = []
    for i in range(len(self.pop)):
      count = Counter(self.pop[i].get("strain_pop", []))
      most_common = count.most_common(1)
      most_frequent_strain = most_common[0][0] if most_common else None
      winning_strains.append(most_frequent_strain)

    # determine transmission
    transmitted_strains = []

    # if fitness is randomly assigned
    if cross_reactivity == False and MNP == False:
      if uniform_fitness == True:
        transmitted_strains = random.sample(winning_strains, random.randint(0, len(winning_strains)))
      else: 
        for strain in winning_strains:
          if random.random() < values[strain]:
            transmitted_strains.append(strain)
      return transmitted_strains

    # if cross reactivity is on
    if host_strains == [] and cross_reactivity == True:
      for strain in winning_strains:
        if random.random() < 0.7:
          transmitted_strains.append(strain)
      return transmitted_strains
    if cross_reactivity == True and host_strains != []:
      for strain in winning_strains:
        distances = []
        for strain2 in host_strains:
          distances.append(hamming_distance(strain, strain2)) # to avoid using antigen distance matrix
        x = min(distances)
        if random.random() < values[x]:
          transmitted_strains.append(strain)
      return transmitted_strains

    # if host specialization is on
    if MNP == True:
      for strain in winning_strains:
        if random.random() < values.get(strain)[host_type]:
        #z = [d for d in values if d['strain'] == strain]
        #if random.random() < z[0][host_type]:
          transmitted_strains.append(strain)
      return transmitted_strains

  # host to tick
  def host2tick_transmission(self, strains):
    if strains == []:
      return []

    if len(strains) > 5:
      strains = random.sample(strains, random.randint(1, 5))

    if random.random() < 0.7:
      number_of_strains = random.randint(1, len(strains))
      transmission_community = random.sample(strains, number_of_strains)
      transmitted_strains = []
      #all strains in transmission community have 0.7 probability of getting transmitted
      for strain in transmission_community:
        if random.random() < 0.7:
          transmitted_strains.append(strain)
      return transmitted_strains
    else:
      return []
    
