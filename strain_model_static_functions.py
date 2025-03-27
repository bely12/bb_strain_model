import numpy as np
import random
import itertools
from itertools import combinations
import math
#import copy

# reboot for static strains

# general functions
def hamming_distance(string1, string2):
  distance = 0
  for i in range(len(string1)):
      if string1[i] != string2[i]:
          distance += 1
  return distance

def assign_hs(spec_type): # assigns a single hs value
  if spec_type == 'generalist':
    hs_val = 0.5
  if spec_type == 'rodent':
    hs_val = random.uniform(0.51, 1.0)
  if spec_type == 'bird':
    hs_val = random.uniform(0.0, 0.49)
  return round(hs_val, 2)

def assign_hs_v2(spec_types): # takes a list of strain hs types; can be used to create their values within Pathogen class
  hs_vals = []
  for spec in spec_types:
    if spec == 'any':
      hs_val = round(random.uniform(0.0,1.0),2)
    if spec == 'generalist':
      hs_val = 0.5
    if spec == 'rodent':
      hs_val = round(random.uniform(0.51, 1.0),2)
    if spec == 'bird':
      hs_val = round(random.uniform(0.0, 0.49),2)
    
    if spec == 'spec_bias': 
      if random.random() <= 0.5:
        hs_val = round(random.uniform(0.0,0.2),2)
      else:
        hs_val = round(random.uniform(0.8,1.0),2)
    if spec == 'gen_bias':
      hs_val = round(random.uniform(0.4, 0.6),2)
    
    hs_vals.append(hs_val)
  return hs_vals

def assign_antigen(n): # assigns antigens randomly to each strain being created; use when creating Pathogen class
  antigen_vals = []
  for i in range(n):
    antigen_vals.append(round(random.uniform(0.0, 2.0),2))
  return antigen_vals

def antigen_distance(antigen1, antigen2):
  lrg = max(antigen1, antigen2)
  sml = min(antigen1, antigen2)
  dist1= sml + (2.0-lrg)
  dist2 = lrg - sml
  return min(dist1, dist2)

def dispersion(antigen_value_list): # this only works correctly for 3 strains!
  all_dists = []
  combos = list(itertools.combinations(antigen_value_list, 2))
  for pair in combos:
    all_dists.append(antigen_distance(pair[0], pair[1]))
  return sum(all_dists)

def spec_weight(pathogens):
  vals = []
  for i in range(len(pathogens)):
    vals.append(((pathogens[i]['hs'] - 0.5)**2)/3)
  return round(math.sqrt(sum(vals)),2)

def transmission_probability(strain, host):
  dists = []
  for i in range(len(host['infections'])):
    dists.append(antigen_distance(strain['antigen'], host['infections'][i]['antigen']))

  if host['host_type'] == 'rodent':
    return (strain['hs']**2) * (1 - math.exp(-2.0 * min(dists)))
  
  if host['host_type'] == 'bird':
    return (1-(strain['hs']**2)) * (1 - math.exp(-2.0 * min(dists)))

def tick2host_transmission(tick, host):
  # create the transmission community from tick
  if tick['strains'] == []:
    return []
  else:
    tick_transmission_community = random.sample(tick['strains'], k = random.randint(1,len(tick['strains'])))

  # transmission
  transmitted_strains = []
  if host['infections'] == []:
    for strain in tick_transmission_community:
      
      if host['host_type'] == 'rodent':
        transmission_prob = (strain['hs']**2) * (1 - math.exp(-2.0))
      
      if host['host_type'] == 'bird':
        transmission_prob = (1-(strain['hs']**2)) * (1 - math.exp(-2.0))
      
      if random.random() < transmission_prob:
        transmitted_strains.append(strain)
    return transmitted_strains

  else:
    for strain in tick_transmission_community:
      if random.random() < transmission_probability(strain, host):
        transmitted_strains.append(strain)
    return transmitted_strains


def host2tick_transmission(host):
  # create transmission community from host
  if host['infections'] == []:
    return []
  else:
    if random.random() < 0.7:
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

################################################################################
# host population
class Host:
  def __init__(self, n_rodents, n_birds):

    hosts = []
    existing_ids = []
    for i in range(n_rodents):
      while True: # assign a unique id
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      hosts.append({'id': id, 'infections': [], 'host_type': 'rodent', 'born': random.choice(range(1,76))})

    for i in range(n_birds):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      hosts.append({'id': id, 'infections': [], 'host_type': 'bird', 'born': random.choice(range(1,76))})

    # define pops and sizes
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

################################################################################
# pathogen population
class Pathogen:
  def __init__(self, n, hs, antigen): # hs and antigen args should be lists equal to length of n

    self.strain_community = []
    for i in range(n):
      self.strain_community.append({'strain': str(i), 'hs': hs[i], 'antigen': antigen[i]})


################################################################################
# vector population
class Vector:
  def __init__(self, pop_size, strains, lam):

    if pop_size % 2 != 0:
      raise ValueError("The value of pop_size must be divisible by 2.")

    # define the pathogens and starting infection status
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
      if nymph_infections[i] == 0:
        tick_pop.append({'id': id, 'stage': 'nymph', 'strains': []})
      else:
        if nymph_infections[i] >len(strains):
          variants = random.sample(strains, k=len(strains))
          tick_pop.append({'id': id, 'stage': 'nymph', 'strains': variants})
        if nymph_infections[i] <len(strains):
          variants = random.sample(strains, k=nymph_infections[i])
          tick_pop.append({'id': id, 'stage': 'nymph', 'strains': variants})

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
    self.nymph_pop = nymph_pop
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.pop_size = len(self.pop)
    self.year = 0
    self.poisson_infection_status = nymph_infections
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)


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
    self.nymph_pop = nymph_pop # I probably don't need this but I'll leave it here for now
    self.num_carried = sum(len(tick['strains']) for tick in self.nymph_pop) / len(self.nymph_pop)
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.year += 1


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
