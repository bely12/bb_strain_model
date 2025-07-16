import numpy as np
import random
import collections
from collections import Counter
import itertools
import math
import copy

### functions ###
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

def adaptive_fitness(val, host, gen_fit = 'low'):
  
  if gen_fit == 'high':
    if host == 'bird':
      probability = 1 - abs(val)**1.5
    elif host == 'rodent':
      probability = 1 - abs(val - 1)**1.5
    return round(probability, 3)
  
  elif gen_fit == 'low':
    if host == 'bird':
      probability = (1-val)**2
    elif host == 'rodent':
      probability = val**2
    return round(probability, 3)
  
  elif gen_fit == 'linear':
    if host == 'bird':
      probability = val
    elif host == 'rodent':
      probability = 1 - val
    return round(probability,3)
  

def adaptive_trait_val(adaptive_gene):
  catch = []
  for i in range(len(adaptive_gene)):
    catch.append(int(adaptive_gene[i]))
  return sum(catch)/len(adaptive_gene)

def spec_weight(adaptive_vals_list):
  weights = []
  for val in adaptive_vals_list:
    weights.append(((val - 0.5)**2)/len(adaptive_vals_list))
  return round(math.sqrt(sum(weights)),2)

def calc_pop_cycle(rodents, birds, intensity):
  return abs(rodents-birds) / ((rodents + birds) * intensity)

################################################################################

class Host:

  def __init__(self, rodents, birds):
    
    pop_size = rodents + birds
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

      if i < rodents:
        host_type = 'rodent'
      else:
        host_type = 'bird'

      # create host
      hosts.append({'id': id, 'infections': [], 'host_type': host_type, 'born': random.choice(range(1,76))})

    self.pop = hosts
    self.pop_size = len(self.pop)
    self.rodent_pop = [d for d in self.pop if d.get('host_type') == 'rodent']
    self.rodent_pop_size = len(self.rodent_pop)
    self.bird_pop = [d for d in self.pop if d.get('host_type') == 'bird']
    self.bird_pop_size = len(self.bird_pop)

  def refresh(self, day, switch_ids, dynamic, growing_species):
    for host in self.pop:
      if host['born'] == day:
        host['infections'] = [] # clear infection status 
        
        if dynamic: # for use when fluctuating species numbers
          if host['id'] in switch_ids:
            host['host_type'] = growing_species
    
    # redefine sub_pops
    self.rodent_pop = [d for d in self.pop if d.get('host_type') == 'rodent']
    self.rodent_pop_size = len(self.rodent_pop)
    self.bird_pop = [d for d in self.pop if d.get('host_type') == 'bird']
    self.bird_pop_size = len(self.bird_pop)
  
  def host_switch(self, intensity, declining_species):
    n = int(intensity * len(self.pop)) # how many hosts to switch ### will need a check built in to ensure even number

    if declining_species == 'rodent':
      hosts_to_switch = random.sample(self.rodent_pop, n) # pick n random members of rodent pop
    elif declining_species == 'bird':
      hosts_to_switch = random.sample(self.bird_pop, n) # pick n random members of bird pop

    self.switch_ids = []
    for host in hosts_to_switch:
      self.switch_ids.append(host['id'])

#####

class Vector:

  def __init__(self, pop_size, n_strains, strain_length, gene, lam, adpt_sel = False):

    if pop_size % 2 != 0:
      raise ValueError("The value of pop_size must be divisible by 2.")
    if strain_length % 2 != 0:
      raise ValueError("The length of antigen sequence needs to be divisible by 2.")

    # create pathogen strains
    seq_bits = ['0','1']
    self.ancestral_strains = []
    for i in range(n_strains):
      if gene != 'modular':
        while True:
          variant = ''.join(random.choices(seq_bits, k=strain_length))
          if variant not in self.ancestral_strains:
            self.ancestral_strains.append(variant)
            break
      else: 
        while True:
          first_half = ''.join(random.choices(seq_bits, k=strain_length//2)) # random 
          second_half = '01' * ((strain_length //2)//2) # yields generalist start value
          variant = first_half + second_half
          if variant not in self.ancestral_strains:
            self.ancestral_strains.append(variant)
            break
    
    if adpt_sel or gene == 'multi':
      adpt_strain_set = ['01'*(strain_length//2), '10'*(strain_length//2), '0'*(strain_length//2) + '1'*(strain_length//2), '1'*(strain_length//2) + '0'*(strain_length//2)]
      adaptive_gene = random.choice(adpt_strain_set)

    # infection status for nymphs in starting tick pop
    nymph_infections = np.random.poisson(lam, size=(pop_size // 2))
    
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
        if gene != 'multi': 
          if adpt_sel == False:
            variant = random.choice(self.ancestral_strains) # pick a variant (random used in case I ever start with more than one strain)
          else: 
            variant = adaptive_gene
          tick_pop.append({'id': id, 'stage': 'nymph', 'strains': [{'lineage_id': lin_id, 'variant':variant, 'history': [variant]}]}) 
          lin_id += 1

        else: 
          variant = random.choice(self.ancestral_strains)
          tick_pop.append({'id': id, 'stage': 'nymph', 'strains': [{'lineage_id': lin_id, 'variant':variant, 'adaptive_gene': adaptive_gene, 'history': [variant], 'adp_history': [adaptive_gene]}]}) 
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

    # define pop
    self.pop = tick_pop
    
    ### POPULATION MONITORING ###
    nymph_pop = [d for d in tick_pop if d.get('stage') == 'nymph']
    strain_pop = []
    for d in nymph_pop:
      if d['strains'] != []:
        for i in range(len(d['strains'])):
          strain_pop.append(d['strains'][i]['variant'])

    self.year = 0
    self.poisson_infection_status = nymph_infections
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.num_carried = round(sum(len(tick['strains']) for tick in nymph_pop) / sum(1 for d in nymph_pop if d.get("strains")),3)
    self.avg_gen_dist = 0.0 # set to 0 at start of sim; will change if sampling is performed

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
    for i in range(len(self.pop) // 2):
      while True:
        new_id = random.randint(10000, 99999)
        if new_id not in existing_ids:
          id = new_id
          break
      existing_ids.append(id)
      updated_pop.append({'id': id,'stage': 'larva', 'strains': []})

    # redefine pop
    self.pop = updated_pop
    self.year += 1

  # updates infection rate and avg strains carried 
  def pop_stats(self):
    nymph_pop = [d for d in self.pop if d.get('stage') == 'nymph']
    strain_pop = []
    for d in nymph_pop:
      if d['strains'] != []:
        for i in range(len(d['strains'])):
          strain_pop.append(d['strains'][i]['variant'])

    # self.nymph_pop = nymph_pop 
    self.infection_rate = round(sum(1 for d in nymph_pop if d.get("strains")) / len(nymph_pop) * 100, 3)
    self.num_carried = round(sum(len(tick['strains']) for tick in nymph_pop) / sum(1 for d in nymph_pop if d.get("strains")),3)

  ## function samples ticks to get a whole bunch of statistics
  def sample(self, gene):
    
    # step 1 - sample 100 infected ticks to get strains they carry
    nymph_pop = [d for d in self.pop if d.get('stage') == 'nymph'] 
    infected_ticks = [d for d in nymph_pop if len(d['strains']) != 0] 
    if len(infected_ticks) < 100:
      sampled_ticks = infected_ticks
    else:
      sampled_ticks = random.sample(infected_ticks, 100)
    self.sampled_ticks = sampled_ticks
    
    sampled_strains = []
    sampled_adaptive_seqs = []
    sampled_combined_seqs = []
    for d in sampled_ticks:
      for i in range(len(d['strains'])):
        sampled_strains.append(d['strains'][i]['variant'])
        if gene == 'multi':
          sampled_adaptive_seqs.append(d['strains'][i]['adaptive_gene'])
          sampled_combined_seqs.append(d['strains'][i]['variant'] + d['strains'][i]['adaptive_gene'])
    
    self.sampled_strains = sampled_strains
    self.sampled_adaptive = sampled_adaptive_seqs
    self.sampled_combined = sampled_combined_seqs

    # step 2 compute antigen distances and adaptive vals
    generation_distances = []
    generation_adaptive_vals = []

    for i in range(len(sampled_strains)):
      # first get all adaptive trait values 
      if gene == 'single':
        generation_adaptive_vals.append(adaptive_trait_val(sampled_strains[i]))
      if gene == 'modular':
        generation_adaptive_vals.append(adaptive_trait_val(sampled_strains[i][len(sampled_strains[i])//2:]))
      if gene == 'multi': 
        generation_adaptive_vals.append(adaptive_trait_val(sampled_adaptive_seqs[i]))

      # then get pairwise distances 
      for j in range(i + 1, len(sampled_strains)):
        dist = hamming_distance(sampled_strains[i], sampled_strains[j])
        generation_distances.append(dist)

    # step 3 - define values for data collection during sim
    self.avg_gen_dist = round(np.mean(generation_distances),3)
    self.sampled_distances = generation_distances
    self.sampled_adaptive_vals = generation_adaptive_vals
    self.sampled_spec_weight = spec_weight(self.sampled_adaptive_vals)

  ## function assigns a host for each tick to bite and the day/s it bites on
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

  def recombination(self, replace, gene = 'NA', rate = 0.01):
    for tick in self.pop:
      # only move forward if there are multiple strains being carried by the tick
      if tick['strains'] != [] and len(tick['strains']) > 1:
        strains = copy.deepcopy(tick['strains'])
      else:
        continue
      
      # cycle through each of the strains carried by the tick and decide if recombination will happen
      updated_strains = []
      for j in range(len(strains)):
        strain = strains[j]
        if random.random() < rate:
          break_start = random.randint(0, len(strains[j]['variant'])-1)
          break_end = random.randint(break_start, len(strains[j]['variant']))
          strains_copy = copy.deepcopy(strains) # doing this so I can remove the recieving strain when choosing donor
          strains_copy.remove(strain)
          donor = random.choice(strains_copy)
          recombinant = strain['variant'][:break_start] + donor['variant'][break_start:break_end] + strain['variant'][break_end:]
        else: 
          recombinant = strain['variant']
          
        # reapeat for adaptive gene if applicable 
        if random.random() < rate and gene == 'multi':
          break_start = random.randint(0, len(strains[j]['adaptive_gene'])-1)
          break_end = random.randint(break_start, len(strains[j]['adaptive_gene']))
          strains_copy = copy.deepcopy(strains)
          strains_copy.remove(strain)
          donor = random.choice(strains_copy)
          recombinant2 = strain['adaptive_gene'][:break_start] + donor['adaptive_gene'][break_start:break_end] + strain['adaptive_gene'][break_end:]
        elif gene == 'multi':
          recombinant2 = strain['adaptive_gene']

        # update with recombinant
        if replace: 
          if recombinant != strain['variant']:
            strains[j]['history'].append(recombinant) 
            strains[j]['variant'] = recombinant
          if gene == 'multi' and recombinant2 != strain['adaptive_gene']:
            strains[j]['adp_history'].append(recombinant2)
            strains[j]['adaptive_gene'] = recombinant2

        elif not replace:  
          if gene != 'multi':
            if recombinant != strain['variant']:
              updated_strains.append({'lineage_id': strain['lineage_id'], 'variant': recombinant, 'history': strain['history']+[recombinant]})

          elif recombinant != strain['variant'] and recombinant2 != strain['adaptive_gene']:
            tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': recombinant, 'adaptive_gene': recombinant2, 'history': strain['history']+[recombinant], 'adp_history': strain['adp_history']+[recombinant2]})
          elif recombinant != strain['variant']:
            tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': recombinant, 'adaptive_gene': strain['adaptive_gene'], 'history': strain['history']+[recombinant], 'adp_history': strain['adp_history']})
          elif recombinant2 != strain['adaptive_gene']:
            tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': strain['variant'], 'adaptive_gene': recombinant2, 'history': strain['history'], 'adp_history': strain['adp_history']+[recombinant2]})

      # update the actual pop
      if replace:
        tick['strains'] = strains
      else:
        tick['strains'].extend(updated_strains)

  def mutate(self, replace, gene = 'NA', rate = 0.01): 
    for tick in self.pop:
      if tick['strains'] != []:
        strains = copy.deepcopy(tick['strains'])
        if not replace:
          strains_to_add = []
        for strain in strains:
          mutated_variant = None
          mutated_adp_gene = None
          
          variant = strain['variant']
          mut_pois_dist = np.random.poisson(rate, len(variant))
          if mut_pois_dist.sum() > 0:
            new_string = []
            for k in range(len(variant)):
              if mut_pois_dist[k] > 0:
                new_bit = str(random.randint(0,1))
                new_string.append(new_bit)
              else:
                new_string.append(variant[k])
            mutated_variant = ''.join(new_string)
          
          if gene == 'multi':
            adp_gene = strain['adaptive_gene']
            mut_pois_dist2 = np.random.poisson(rate, len(adp_gene))
            if mut_pois_dist2.sum() > 0:
              new_string2 = []
              for k in range(len(adp_gene)):
                if mut_pois_dist2[k] > 0:
                  new_bit2 = str(random.randint(0,1))
                  new_string2.append(new_bit2)
                else:
                  new_string2.append(adp_gene[k])
              mutated_adp_gene = ''.join(new_string2)
          
          if replace:
            if mutated_variant is not None and mutated_variant != variant:
              strain['variant'] = mutated_variant
              strain['history'].append(mutated_variant)
            if gene == 'multi' and mutated_adp_gene is not None and mutated_adp_gene != adp_gene:
              strain['adaptive_gene'] = mutated_adp_gene
              strain['adp_history'].append(mutated_adp_gene)
          
          if not replace: 
            if gene == 'multi': 
              if mutated_variant is not None and mutated_adp_gene is not None: 
                strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'adaptive_gene': mutated_adp_gene, 'history': strain['history']+[mutated_variant], 'adp_history': strain['adp_history']+[mutated_adp_gene]})
              elif mutated_variant is not None:
                strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'adaptive_gene': strain['adaptive_gene'], 'history': strain['history']+[mutated_variant], 'adp_history': strain['adp_history']})
              elif  mutated_adp_gene is not None:
                strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': strain['variant'], 'adaptive_gene': mutated_adp_gene, 'history': strain['history'], 'adp_history': strain['adp_history']+[mutated_adp_gene]})
            elif mutated_variant is not None:
              strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'history': strain['history']+[mutated_variant]})
        
        if replace:
          tick['strains'] = strains
        elif strains_to_add != []:
          tick['strains'].extend(strains_to_add)


  # def recombination(self, replace, gene = 'NA', rate = 0.01):
  #   for tick in self.pop:
  #     # only move forward if there are multiple strains being carried by the tick
  #     if tick['strains'] != [] and len(tick['strains']) > 1:
  #       strains = tick['strains']
  #     else:
  #       continue
      
  #     # cycle through each of the strains carried by the tick and decide if recombination will happen
  #     for j in range(len(strains)):
  #       if random.random() < rate:
  #         break_start = random.randint(0, len(strains[j]['variant'])-1)
  #         break_end = random.randint(break_start, len(strains[j]['variant']))
  #         strain = strains[j]
  #         strains_copy = copy.deepcopy(strains) # doing this so I can remove the recieving strain when choosing donor
  #         strains_copy.remove(strain)
  #         donor = random.choice(strains_copy)
  #         recombinant = strain['variant'][:break_start] + donor['variant'][break_start:break_end] + strain['variant'][break_end:]

  #         # update with recombinant
  #         if recombinant != strain['variant'] and replace == True:
  #           tick['strains'][j]['history'].append(recombinant) 
  #           tick['strains'][j]['variant'] = recombinant 
  #         if recombinant != strain['variant'] and replace == False:
  #           if gene != 'multi': # accounts for different dicitonary structure of multi gene mode
  #             tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': recombinant, 'history': strain['history']+[recombinant]})
  #           else:
  #             tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': recombinant, 'adaptive_gene': strain['adaptive_gene'], 'history': strain['history']+[recombinant], 'adp_history': strain['adp_history']})

  #       # reapeat for adaptive gene if applicable 
  #       if random.random() < rate and gene == 'multi':
  #         break_start = random.randint(0, len(strains[j]['adaptive_gene'])-1)
  #         break_end = random.randint(break_start, len(strains[j]['adaptive_gene']))
  #         strain = strains[j]
  #         strains_copy = copy.deepcopy(strains)
  #         strains_copy.remove(strain)
  #         donor = random.choice(strains_copy)
  #         recombinant2 = strain['adaptive_gene'][:break_start] + donor['adaptive_gene'][break_start:break_end] + strain['adaptive_gene'][break_end:]

  #         # update ticks strains seqs and history; decide to replace or add to community
  #         if recombinant2 != strain['adaptive_gene'] and replace == True:
  #           tick['strains'][j]['adp_history'].append(recombinant2) 
  #           tick['strains'][j]['adaptive_gene'] = recombinant2 
  #         if recombinant2 != strain['adaptive_gene'] and replace == False:
  #           tick['strains'].append({'lineage_id': strain['lineage_id'], 'variant': strain['variant'], 'adaptive_gene': recombinant2, 'history': strain['history'], 'adp_history': strain['adp_history']+[recombinant2]})

  # def mutate(self, replace, gene = 'NA', rate = 0.01): 
  #   for tick in self.pop:
  #     if tick['strains'] != []:
  #       if not replace:
  #         strains_to_add = []
  #       for strain in tick['strains']:
  #         mutated_variant = None
  #         mutated_adp_gene = None
          
  #         variant = strain['variant']
  #         mut_pois_dist = np.random.poisson(rate, len(variant))
  #         if mut_pois_dist.sum() > 0:
  #           new_string = []
  #           for k in range(len(variant)):
  #             if mut_pois_dist[k] > 0:
  #               new_bit = str(random.randint(0,1))
  #               new_string.append(new_bit)
  #             else:
  #               new_string.append(variant[k])
  #           mutated_variant = ''.join(new_string)
          
  #         if gene == 'multi':
  #           adp_gene = strain['adaptive_gene']
  #           mut_pois_dist2 = np.random.poisson(rate, len(adp_gene))
  #           if mut_pois_dist2.sum() > 0:
  #             new_string2 = []
  #             for k in range(len(adp_gene)):
  #               if mut_pois_dist2[k] > 0:
  #                 new_bit2 = str(random.randint(0,1))
  #                 new_string2.append(new_bit2)
  #               else:
  #                 new_string2.append(adp_gene[k])
  #             mutated_adp_gene = ''.join(new_string2)
          
  #         if replace:
  #           if mutated_variant is not None and mutated_variant != variant:
  #             strain['variant'] = mutated_variant
  #             strain['history'].append(mutated_variant)
  #           if gene == 'multi' and mutated_adp_gene is not None and mutated_adp_gene != adp_gene:
  #             strain['adaptive_gene'] = mutated_adp_gene
  #             strain['adp_history'].append(mutated_adp_gene)
          
  #         if not replace: 
  #           if gene == 'multi': 
  #             if mutated_variant is not None and mutated_adp_gene is not None: 
  #               strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'adaptive_gene': mutated_adp_gene, 'history': strain['history']+[mutated_variant], 'adp_history': strain['adp_history']+[mutated_adp_gene]})
  #             elif mutated_variant is not None:
  #               strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'adaptive_gene': strain['adaptive_gene'], 'history': strain['history']+[mutated_variant], 'adp_history': strain['adp_history']})
  #             elif  mutated_adp_gene is not None:
  #               strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': strain['variant'], 'adaptive_gene': mutated_adp_gene, 'history': strain['history'], 'adp_history': strain['adp_history']+[mutated_adp_gene]})
  #           elif mutated_variant is not None:
  #             strains_to_add.append({'lineage_id': strain['lineage_id'], 'variant': mutated_variant, 'history': strain['history']+[mutated_variant]})
        
  #       if not replace and strains_to_add != []:
  #         tick['strains'].extend(strains_to_add)

  # def old_mutate(self, rate = 0.01):
  #   for i in range(len(self.pop)):
  #     if self.pop[i]['strains'] != []:
  #       #temp_strain_set = copy.deepcopy(self.pop[i]['strains'])
  #       #for j in range(len(temp_strain_set)):
  #       for j in range(len(self.pop[i]['strains'])): # new
  #         #variant = temp_strain_set[j]['variant']
  #         variant = self.pop[i]['strains'][j]['variant'] # new
  #         mut_pois_dist = np.random.poisson(rate, len(variant))
  #         new_string = []
  #         for k in range(len(variant)):
  #           if mut_pois_dist[k] > 0:
  #             new_bit = str(random.randint(0,1))
  #             new_string.append(new_bit)
  #           else:
  #             new_string.append(variant[k])
  #         mutated_string = ''.join(new_string)
  #         if mutated_string != variant:
  #           #temp_strain_set[j]['history'].append(mutated_string)
  #           self.pop[i]['strains'][j]['history'].append(mutated_string) # new
  #         #temp_strain_set[j]['variant'] = mutated_string
  #         self.pop[i]['strains'][j]['variant'] = mutated_string # new
  #       #self.pop[i]['strains'] = temp_strain_set
  #     else:
  #       continue

################################################################################
# transmission functions

def tick2host_transmission(tick, host, selection, gene, transmission_probabilities=None, host_type=None, gen_fit = 'high'):
  # selection: immune, adaptive, neutral, hybrid
  # gene: modular, multi 
  # gen_fit: high, low; how fit are generalists with adaptive selection
  
  # CREATE TRANSMISSION COMMUNITY #
  if tick['strains'] == []:
    tick_transmission_community = []
  if random.random() < 0.7 and tick['strains'] != []:
    number_of_strains = random.randint(1, len(tick['strains']))
    tick_transmission_community = random.sample(tick['strains'], number_of_strains)
  else:
    tick_transmission_community = []

  if tick_transmission_community != []:
    transmitted_strains = []

    # IMMUNE SELECTION #
    if selection == 'immune' and host['infections'] == []:
    #if host['infections'] == [] and cross_reactivity == True and not host_specialization:
      for strain in tick_transmission_community:
        if random.random() < 0.7:
          transmitted_strains.append(strain)
      return transmitted_strains

    if selection == 'immune' and host['infections']:
      #if cross_reactivity and host['infections'] and not host_specialization:
        for strain in tick_transmission_community:
          distances = []
          for strain2 in host['infections']:
            distances.append(hamming_distance(strain['variant'], strain2['variant']))
          x = min(distances)
          if random.random() < transmission_probabilities[x]:
            transmitted_strains.append(strain)
        return transmitted_strains

    # ADAPTIVER SELECTION #
    if selection == 'adaptive':
    #if host_specialization == True and not cross_reactivity:
      for strain in tick_transmission_community:
        hs_prob = adaptive_fitness(adaptive_trait_val(strain['variant']), host=host_type, gen_fit=gen_fit)
        if random.random() < hs_prob:
          transmitted_strains.append(strain)
      return transmitted_strains

    # NEUTRAL SELECTION #
    if selection == 'neutral':
      for strain in tick_transmission_community:
        if random.random() < 0.5:
          transmitted_strains.append(strain)
      return transmitted_strains
    
    # HYBRID SELECTION #
    if selection == 'hybrid':
      for strain in tick_transmission_community:
        if gene == 'modular':
          hs_prob = adaptive_fitness(adaptive_trait_val(strain['variant'][len(strain['variant'])//2:]), host=host_type, gen_fit=gen_fit)
        if gene == 'multi':
          hs_prob = adaptive_fitness(adaptive_trait_val(strain['adaptive_gene']), host=host_type, gen_fit=gen_fit)
        
        if host['infections'] == []:
          if random.random() < hs_prob:
            transmitted_strains.append(strain)
        else:
          distances = []
          for strain2 in host['infections']:
            if gene == 'modular':
              distances.append(hamming_distance(strain['variant'][:len(strain['variant'])//2], strain2['variant'][:len(strain2['variant'])//2]))
            if gene == 'multi':
              distances.append(hamming_distance(strain['variant'], strain2['variant']))
          x = min(distances)
          immune_prob = transmission_probabilities[x]
          total_prob = hs_prob * immune_prob
          if random.random() < total_prob:
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
