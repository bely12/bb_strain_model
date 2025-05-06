import numpy as np
import random
import test_strain_model_functions as be
from tqdm import tqdm
import argparse
from argparse import RawTextHelpFormatter

##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='Vector borne pathogen evolution simulation'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-mode', type=int, help='1: cross reactivity, 2: host specialization, 3: uniform fitness\n')
parser.add_argument('-antigen_length', type=int, help='size of bit string representing antigen')
parser.add_argument('-vector_pop_size', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-years', type=int, help='number of years to simulate')
parser.add_argument('-mut', type=float, default=0.01, help='per site mutation rate in antigen')
parser.add_argument('-recomb_rate', type=float, default = 0.01, help = 'recombination rate')
parser.add_argument('-out', default=None, help='prefix for output file')
parser.add_argument('-run_tag', default = 1, help='unique id for each sim run when doing batch')
args = parser.parse_args()

##### set key parameters #####
vector_pop_size = args.vector_pop_size # needs to be an even number, at least 50
antigen_length = args.antigen_length
mutation_rate = args.mut
recombination_rate = args.recomb_rate
transmission_curve = 'sigmoid' # choose sigmoid or linear
sim_years = args.years

if args.mode == 1:
  cross_reactivity = True
  host_specialization = False
  uniform_fitness = False
if args.mode == 2:
  cross_reactivity = False
  host_specialization = True
  uniform_fitness = False
if args.mode == 3:
  cross_reactivity = False
  host_specialization = False
  uniform_fitness = True

# print sim parameters
print('Parameters','\n',
      'Vector pop size: ',vector_pop_size,'\n',
      'host pop size: ', round(vector_pop_size/50),'\n',
      'antigen length: ',antigen_length,'\n',
      'per site mutation rate: ', mutation_rate,'\n'
      'cross reactivity: ', cross_reactivity,'\n',
      'host specialization: ', host_specialization,'\n',
      'uniform fitness value: ', uniform_fitness, '\n',
      'simulated years: ', sim_years,'\n',
      'batch progress (current run): ', args.run_tag)

##### initialize populations #####
ticks = be.Vector(pop_size= vector_pop_size, n_strains= 1, strain_length= antigen_length, lam= 0.5)
host_pop_size = round(len(ticks.pop)/50)
hosts = be.Host(host_pop_size, host_specialization)

##### create transmission probability tables #####
if cross_reactivity == True:
  probabilities = be.distance_probabilities(list(range(antigen_length+1)), curve= transmission_curve)
if host_specialization == True:
  mnp_values = be.mnp_probability(ticks.strain_set)
  mnp_values2 = {item['strain']: {'mnp_value': item['mnp_value'], 'rodent': item['rodent'], 'bird': item['bird']} for item in mnp_values}
if uniform_fitness == True:
  fitness_values = dict(zip(ticks.strain_set, np.round(np.random.uniform(0.0, 1.0, len(ticks.strain_set)),4)))

##### initialize data collection #####
unique_lineage_ids = set()
for tick in ticks.nymph_pop:
  for strain in tick['strains']:
    unique_lineage_ids.add(strain['lineage_id'])

all_data = [{'run_tag': args.run_tag,
             'year': ticks.year,
             'infection_rate': ticks.infection_rate,
             'diversity': ticks.diversity,
             'avg_antigen_distance': ticks.avg_gen_dist,
             'active_strain_count': ticks.current_strain_count,
             'avg_strains_carried': ticks.num_carried,
             'unique_lineages': len(unique_lineage_ids)}]

if host_specialization == True:
  mnp_container = [ticks.mnp_recs(mnp_values2)]

########## sim ##########

#for year in range(sim_years):
for year in tqdm(range(sim_years)):
  
  # set tick/host interactions for the year
  interactions = ticks.interaction(hosts.pop)

  # sim through the days of the year
  for day in range(1, 151):

    # refresh the host pop for births and deaths
    hosts.refresh(day)

    # iterate through interactions to see if there is a bite for current day
    for j in range(len(interactions)):
      if interactions[j]['bite_day'] == day:

        # find the host and tick involved in interaction
        current_host = next((item for item in hosts.pop if item["id"] == interactions[j]['host_id']), None)
        current_tick = next((item for item in ticks.pop if item["id"] == interactions[j]['tick_id']), None)
        
        # determine strains to be transmitted from tick to host
        if cross_reactivity == True:
          strains_from_tick = be.tick2host_transmission(tick=current_tick, host=current_host, transmission_probabilities=probabilities, cross_reactivity= True)
  
        if host_specialization == True:
          strains_from_tick = be.tick2host_transmission(tick=current_tick, host=current_host, host_type= current_host['host_type'], transmission_probabilities=mnp_values2, host_specialization=True)

        if uniform_fitness == True:
          strains_from_tick = be.tick2host_transmission(tick=current_tick, host=current_host, transmission_probabilities=fitness_values ,uniform_fitness=True)
        
        # determine strains to be transmitted from host to tick
        strains_from_host = be.host2tick_transmission(current_host)
        
        # update host infection community
        for item in hosts.pop:
          if item["id"] == current_host['id']:
            item["infections"] = current_host['infections'] + [item for item in strains_from_tick if (item not in current_host['infections'])]
            break
        
        # update tick infection community
        for item in ticks.pop:
          if item['id'] == current_tick['id']:
            item['strains'] = current_tick['strains'] + [item for item in strains_from_host if item not in current_tick['strains']]
            break

  # update populations
  ticks.update_pop()
  
  # sampling for pairwise antigen distances and variant frequencies
  ticks.sample()

  # collect data
  unique_lineage_ids = set()
  for tick in ticks.nymph_pop:
    for strain in tick['strains']:
      unique_lineage_ids.add(strain['lineage_id'])

  all_data.append({'run_tag': args.run_tag,
                   'year': ticks.year,
                   'infection_rate': ticks.infection_rate,
                   'diversity': ticks.diversity,
                   'avg_antigen_distance': ticks.avg_gen_dist,
                   'active_strain_count': ticks.current_strain_count,
                   'avg_strains_carried': ticks.num_carried,
                   'unique_lineages': len(unique_lineage_ids)})

  if host_specialization == True:
    mnp_container.append(ticks.mnp_recs(mnp_values2))

  # mutate and recombination
  if year+1 != sim_years:
    ticks.mutate(rate = args.mut)
    ticks.recombination(rate = args.recomb_rate)


##### data output #####
if args.out != None:
  import csv
  import pandas as pd
  from collections import Counter
  
  run = 'run_' + str(args.run_tag)

  ##### save main simulation data to output file #####
  df = pd.DataFrame(all_data)
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_sim_data.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_sim_data.tsv', mode= 'a', sep='\t', index=False, header=False)

  ##### save end of sim variant frequencies to output file #####
  # get variant counts in populaton 
  strain_pop = []
  for d in ticks.nymph_pop:
    if d['strains'] != []:
      for i in range(len(d['strains'])):
        strain_pop.append(d['strains'][i]['variant'])
  counts = dict(Counter(strain_pop))

  # store in dataframe
  df_strain_pop = pd.DataFrame(counts, index=[0])
  df_strain_pop = df_strain_pop.transpose()
  df_strain_pop.rename(columns={0: 'counts'}, inplace=True)
  df_strain_pop['frequency'] = df_strain_pop.index.map(ticks.current_strain_frequencies).fillna('')
  df_strain_pop.reset_index(inplace=True)
  df_strain_pop.rename(columns={'index': 'variant'}, inplace=True)
  df_strain_pop['run_id'] = run

  # write file
  if str(args.run_tag) == "1":
    df_strain_pop.to_csv(args.out+'_variant_frequencies.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df_strain_pop.to_csv(args.out+'_variant_frequencies.tsv', mode= 'a', sep='\t', index=False, header=False)


  ##### save SAMPLED strain counts and frequencies to output file #####
  strain_counts = dict(Counter(ticks.sampled_strains))
  my_dict = []
  for strain in strain_counts:
    my_dict.append({'strain': strain, 'count': strain_counts[strain], 'freq': round(strain_counts[strain]/len(ticks.sampled_strains)*100,2)})
  df_sampled_strains = pd.DataFrame(my_dict)
  df_sampled_strains['run_id'] = run
  if str(args.run_tag) == "1":
    df_sampled_strains.to_csv(args.out+'_sampled_variant_frequencies.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df_sampled_strains.to_csv(args.out+'_sampled_variant_frequencies.tsv', mode= 'a', sep='\t', index=False, header=False)
  

  ##### save sampled pairwise antigen distances #####
  with open(args.out + '_sampled_pairwise_antigen_dists.tsv', 'a', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for value in ticks.generation_distances:
      writer.writerow([run, value])

  ##### save lineage histories to output file #####
  # grab records
  records = []
  for tick in ticks.pop:
    if tick['strains'] != []:
      for i in range(len(tick['strains'])):
        records.append({'lineage_id': tick['strains'][i]['lineage_id'], 'history': tick['strains'][i]['history'], 'time': tick['strains'][i]['time']})
  new_records = []
  for item in records:
    new_records.append(dict(zip(item['time'], item['history'])))

  # put into a dataframe
  df = pd.DataFrame(new_records)
  df = df.reindex(columns=range(sim_years)) # create complete set of columns for years
  df = df.ffill(axis=1) # fill in empty rows with state of variant for that year
  df.drop_duplicates(inplace=True) # drop all duplicate rows
  df['run_id'] = run # create a column for run id and set to index
  df.set_index('run_id', inplace=True)

  # write file
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_lineage_history.tsv', mode='a', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_lineage_history.tsv', mode='a', sep='\t', index=False, header=False)
  

  ##### mnp values in final population to output file ##### 
  if host_specialization == True:
    with open(args.out + '_sim_mnp_pop_values.tsv', 'a', newline='') as file:
      writer = csv.writer(file, delimiter='\t')
      for value in mnp_container[-1]:
        writer.writerow([run, value])