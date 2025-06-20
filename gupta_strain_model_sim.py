import numpy as np
import random
import gupta_strain_model_functions as be
from tqdm import tqdm
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
    description='simulation of 2 loci with 2 alleles each with recombination and different levels of immune cross protection'
    ' \n', formatter_class=RawTextHelpFormatter)
# parameters
parser.add_argument('-loci', type=int, default=2, help='number of 2 allele loci')
parser.add_argument('-cp', default='high', help='degree of immune cross protection [high, medium, low]')
parser.add_argument('-vec', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-yrs', type=int, help='number of years to simulate')
parser.add_argument('-rec', type=float, default = 0.01, help = 'recombination rate')
# record keeping
parser.add_argument('-out', default=None, help='prefix for output file')
parser.add_argument('-run_tag', default = 1, help='unique id for each sim run when doing batch')
args = parser.parse_args()

print('Parameters','\n',
      'Vector pop size: ',args.vec,'\n',
      'host pop size: ', round(args.vec/50),'\n',
      'number of loci: ',args.loci,'\n',
      'recombination rate: ', args.rec,'\n',
      'cross protection level: ', args.cp, '\n',
      'simulated years: ', args.yrs,'\n',
      'batch progress (current run): ', args.run_tag)

##### INITIALIZE POPULATIONS #####
ticks = be.Vector(pop_size= args.vec, loci= args.loci, lam= 0.5)
hosts = be.Host(round(len(ticks.pop)/50))

##### INITIALIZE DATA COLLECTION #####
all_data = [{'run_tag': args.run_tag,
             'year': ticks.year,
             'infection_rate': ticks.infection_rate,
             'distinct_strains': ticks.current_strain_count}]

freq_history = [{'year': ticks.year, 'pop': ticks.strain_frequencies}]

##### START SIMULATION #####
for year in tqdm(range(args.yrs)):
  
  interactions = ticks.interaction(hosts.pop) # determine vector-host interactions for the year

  for day in range(1, 151):

    hosts.refresh(day) # simulates births and deaths of hosts

    for j in range(len(interactions)): # check for vector-host interactions on current day
      if interactions[j]['bite_day'] == day:
        current_host = next((item for item in hosts.pop if item["id"] == interactions[j]['host_id']), None)
        current_tick = next((item for item in ticks.pop if item["id"] == interactions[j]['tick_id']), None)
        
        strains_from_tick = be.tick2host_transmission(current_tick, current_host, args.cp) # transmission
        strains_from_host = be.host2tick_transmission(current_host)
        
        for item in hosts.pop: # update infections
          if item["id"] == current_host['id']:
            item["infections"] = current_host['infections'] + [item for item in strains_from_tick if (item not in current_host['infections'])]
            break
        for item in ticks.pop:
          if item['id'] == current_tick['id']:
            item['strains'] = current_tick['strains'] + [item for item in strains_from_host if item not in current_tick['strains']]
            break

  ticks.update_pop() # molt ticks

  # collect data
  all_data.append({'run_tag': args.run_tag,
                   'year': ticks.year,
                   'infection_rate': ticks.infection_rate,
                   'distinct_strains': ticks.current_strain_count})

  freq_history.append({'year': ticks.year, 'pop': ticks.strain_frequencies})

  if year+1 != args.yrs:
    ticks.recombination(rate = args.rec) # recombination


##### DATA OUTPUT #####
if args.out != None:
  import csv
  import pandas as pd
  from collections import Counter
  
  run = 'run_' + str(args.run_tag)

  # GENERAL SIM DATA 
  df = pd.DataFrame(all_data)
  df['cp'] = args.cp
  df['rec_rate'] = args.rec
  df['vector_pop_size'] = args.vec
  
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_sim_data.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_sim_data.tsv', mode= 'a', sep='\t', index=False, header=False)

  # ALLELE FREQUENCIES
  rows = []
  for entry in freq_history:
    year = entry['year']
    freqs = entry['pop']
    for genotype, freq in freqs.items():
      rows.append({'run_id': run,'year': year, 'genotype': genotype, 'frequency': freq})
  df = pd.DataFrame(rows)
  df['cp'] = args.cp
  df['rec_rate'] = args.rec
  df['vector_pop_size'] = args.vec
  
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_variant_frequencies.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_variant_frequencies.tsv', mode='a', sep='\t', index=False, header=False)
