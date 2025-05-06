import numpy as np
import random
import collections
from collections import Counter
import itertools
import math
import copy

# functions
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
def mnp_probability(strains):
  mnp_values = dict(zip(strains, np.round(np.random.uniform(0.0, 1.0, len(strains)),4)))
  res = []
  for strain in mnp_values:
    res.append({'strain': strain, 'mnp_value': mnp_values[strain], 'rodent': round((mnp_values[strain]**2),4), 'bird': round((1-mnp_values[strain])**2,4)})
  return res


################################################################################
# classes for managing hosts and vectors

class Host:

  def __init__(self, pop_size, host_specialization = False):

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

      if host_specialization == True: # assign specific host species
        if i < pop_size//2:
          host_type = 'rodent'
        else:
          host_type = 'bird'
      else:
        host_type = '-'

      # create host
      hosts.append({'id': id, 'infections': [], 'host_type': host_type, 'born': random.choice(range(1,76))})

    self.pop = hosts
    self.pop_size = len(self.pop)
    self.rodent_pop = [d for d in self.pop if d.get('host_type') == 'rodent']
    self.rodent_pop_size = len(self.rodent_pop)
    self.bird_pop = [d for d in self.pop if d.get('host_type') == 'bird']
    self.bird_pop_size = len(self.bird_pop)

  def refresh(self, day):
    for host in self.pop:
      if host['born'] == day:
        host['infections'] = []

#####

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
        variant = random.choice(self.current_strains)
        tick_pop.append({'id': id, 'stage': 'nymph', 'strains': [{'lineage_id': lin_id,'variant':variant,'history': [variant], 'time': [0]}]})
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
          strain_pop.append(d['strains'][i]['variant'])

    # define important stuff
    self.pop = tick_pop
    self.current_strain_count = len(self.current_strains)
    self.nymph_pop = nymph_pop # probably won't use this
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.pop_size = len(self.pop)
    self.year = 0
    self.poisson_infection_status = nymph_infections
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.diversity = round(1 - (sum((count / len(strain_pop)) ** 2 for count in Counter(strain_pop).values())), 3)

    # get strain frequncy data
    counts = dict(Counter(strain_pop))

    # calculate frequencies for current strains, define as attribute
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts) * 100, 2) for key, value in counts.items()}
    self.current_strain_frequencies = frequencies.copy()

    # add strains not present and assign freq as 0, define as attribute
    for item in self.strain_set:
      if item not in frequencies:
        frequencies[item] = 0.0
    self.strain_set_frequencies = frequencies

    self.avg_gen_dist = 0.0

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
          strain_pop.append(d['strains'][i]['variant'])

    # redefine attributes
    self.nymph_pop = nymph_pop 
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.diversity = round(1 - (sum((count / len(strain_pop)) ** 2 for count in Counter(strain_pop).values())),3)
    self.year += 1

    # get strain frequncy data
    counts = dict(Counter(strain_pop))

    # calculate frequencies for current strains, define as attribute
    total_counts = sum(counts.values())
    frequencies = {key: round((value / total_counts) * 100, 2) for key, value in counts.items()}

    self.current_strains = list(frequencies.keys())
    self.current_strain_count = len(self.current_strains)
    self.current_strain_frequencies = frequencies.copy()

#******************************** turned off *********************************
    # add strains not present and assign freq as 0, define as attribute
    # for item in self.strain_set:
    #   if item not in frequencies:
    #     frequencies[item] = 0.0
    # self.strain_set_frequencies = frequencies

  def sample(self):
    # step 1 - get 100 infected ticks
    infected_ticks = [d for d in self.nymph_pop if len(d['strains']) != 0]
    
    # step 2 - sample the ticks and get all strains they carry
    if len(infected_ticks) < 100:
      sampled_ticks = infected_ticks
    else:
      sampled_ticks = random.sample(infected_ticks, 100)
    self.sampled_ticks = sampled_ticks
    sampled_strains = []
    for d in sampled_ticks:
      for i in range(len(d['strains'])):
        sampled_strains.append(d['strains'][i]['variant'])
    self.sampled_strains = sampled_strains
    
    # step 3 compute antigen distances and save
    generation_distances = []
    for i in range(len(sampled_strains)):
      for j in range(len(sampled_strains)):
        if i != j:
          dist = hamming_distance(sampled_strains[i], sampled_strains[j])
          #self.antigen_distances.append(dist)
          generation_distances.append(dist)
    self.avg_gen_dist = round(np.mean(generation_distances),3)
    self.generation_distances = generation_distances

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
            new_strain = strain['variant'][:break_point] + donor['variant'][break_point:]
          else:
            new_strain = donor['variant'][:break_point] + strain['variant'][break_point:]
          # update ticks strains to reflect recombination events
          if new_strain != strain['variant']:
            tick['strains'][j]['history'].append(new_strain) # update history
            tick['strains'][j]['time'].append(self.year)
          tick['strains'][j]['variant'] = new_strain # recombined strain replaces parent strain


  def mutate(self, rate = 0.01):
    for i in range(len(self.pop)):
      if self.pop[i]['strains'] != []:
        temp_strain_set = copy.deepcopy(self.pop[i]['strains'])
        for j in range(len(temp_strain_set)):
          variant = temp_strain_set[j]['variant']
          mut_pois_dist = np.random.poisson(rate, len(variant))
          new_string = []
          for k in range(len(variant)):
            if mut_pois_dist[k] > 0:
              new_bit = str(random.randint(0,1))
              new_string.append(new_bit)
            else:
              new_string.append(variant[k])
          mutated_string = ''.join(new_string)
          if mutated_string != variant:
            temp_strain_set[j]['history'].append(mutated_string)
            temp_strain_set[j]['time'].append(self.year)
          temp_strain_set[j]['variant'] = mutated_string
        self.pop[i]['strains'] = temp_strain_set
      else:
        continue

################################################################################
# transmission functions

def tick2host_transmission(tick, host, transmission_probabilities, host_type=None, cross_reactivity = False, host_specialization = False, uniform_fitness = False):

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

    # for cross reactivity
    if host['infections'] == [] and cross_reactivity == True:
      for strain in tick_transmission_community:
        if random.random() < 0.7:
          transmitted_strains.append(strain)
      return transmitted_strains

    if cross_reactivity == True and host['infections'] != []:
      for strain in tick_transmission_community:
        distances = []
        for strain2 in host['infections']:
          distances.append(hamming_distance(strain['variant'], strain2['variant']))
        x = min(distances)
        if random.random() < transmission_probabilities[x]:
          transmitted_strains.append(strain)
      return transmitted_strains

    # for host specialization
    if host_specialization == True:
      for strain in tick_transmission_community:
        if random.random() < transmission_probabilities.get(strain['variant'])[host_type]:
          transmitted_strains.append(strain)
      return transmitted_strains

    # for uniform fitness
    if uniform_fitness == True:
      for strain in tick_transmission_community:
        #if random.random() < transmission_probabilities.get(strain['variant']):
        if random.random() < 0.5:
          transmitted_strains.append(strain)
      return transmitted_strains

  # if transmission community was originally empty
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
