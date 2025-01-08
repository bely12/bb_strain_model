import numpy as np
import random
import model_functions_lineage_tracking as be
from tqdm import tqdm
import argparse
from argparse import RawTextHelpFormatter

##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='Vector borne pathogen evolution simulation'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-mode', type=int, help='1: cross reactivity, 2: host specialization, 3: uniform fitness\n')
parser.add_argument('-antigen_length', type=int, help='size of bit string representing antigen')
parser.add_argument('-mut', type=float, default=0.01, help='per site mutation rate in antigen')
parser.add_argument('-vector_pop_size', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-years', type=int, help='number of years to simulate')
parser.add_argument('-plots', default=None, help='set to True if you want plots as output; specify file namein args.plot_file')
parser.add_argument('-out', default=None, help='prefix for output file')
parser.add_argument('-run_tag', default = 1, help='unique id for each sim run when doing batch')
parser.add_argument('-recomb_rate', type=float, default = 0.01, help = 'recombination rate')
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

# print sim parameters; comment out if running batches
print('Parameters','\n',
      'Vector pop size: ',vector_pop_size,'\n',
      'host pop size: ', round(vector_pop_size/50),'\n',
      'antigen length: ',antigen_length,'\n',
      'per site mutation rate: ', mutation_rate,'\n'
      'cross reactivity: ', cross_reactivity,'\n',
      'host specialization: ', host_specialization,'\n',
      'uniform fitness value: ', uniform_fitness, '\n',
      'simulated years: ', sim_years)


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
  #hosts = be.Host(host_pop_size, host_specialization)

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
  
  # list of final strains at end of sim and their frequencies
  # with open(args.out + '_sim_final_antigens.tsv', 'a', newline='') as file:
  #   writer = csv.writer(file, delimiter='\t')
  #   for strain, frequency in ticks.current_strain_frequencies.items():
  #     writer.writerow([run, strain, frequency])
  
  # get final variants and their counts and frequency in the population HAVE TO ADD A COLUMN WITH RUN ID AND EITHER APPEND OR CONCAT DFS FOR RUNNING IN BATCHES!!!
  strain_pop = []
  for d in ticks.nymph_pop:
    if d['strains'] != []:
      for i in range(len(d['strains'])):
        strain_pop.append(d['strains'][i]['variant'])
  counts = dict(Counter(strain_pop))
  df = pd.DataFrame(counts, index=[0])
  df = df.transpose()
  df.rename(columns={0:'counts'}, inplace=True)
  df['frequency'] = df.index.map(ticks.current_strain_frequencies).fillna('')
  df.reset_index(inplace=True)
  df.rename(columns={'index': 'variant'}, inplace=True)
  df['run_id'] = run
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_final_variant_pop.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_final_variant_pop.tsv', mode= 'a', sep='\t', index=False, header=False)

  ### get lineage histories ### 
  # collect histories from population
  records = []
  for tick in ticks.pop:
    if tick['strains'] != []:
      for i in range(len(tick['strains'])):
        records.append({'lineage_id': tick['strains'][i]['lineage_id'], 'history': tick['strains'][i]['history']})
  history_list = []
  for item in records:
    history_list.append(item['history'])
  # get unique histories only
  temp_history = list(map(tuple, history_list))
  unique_history = set(temp_history)
  # build into a dataframe and save
  df = pd.DataFrame({'trace': list(unique_history)})
  df2 = df['trace'].apply(pd.Series)
  df2['run_id'] = run
  df2.to_csv(args.out+'_final_variant_histories.tsv', mode= 'a', sep='\t', index=False, header=False)

  ### sampled pairwise antigen distances ###
  sampled_distances = random.sample(ticks.antigen_distances, 1000)
  with open(args.out + '_sampled_pairwise_antigen_dists.tsv', 'a', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for value in sampled_distances:
      writer.writerow([run, value])
  
  ### mnp values in final population ###  not relevant for batch
  if host_specialization == True:
    with open(args.out + '_sim_mnp_pop_values.tsv', 'w', newline='') as file:
      writer = csv.writer(file, delimiter='\t')
      for value in mnp_container[-1]:
        writer.writerow([run, value])

  # main data
  df = pd.DataFrame(all_data)
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_sim_data.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_sim_data.tsv', mode= 'a', sep='\t', index=False, header=False)