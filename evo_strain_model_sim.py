import numpy as np
import random
import evo_strain_model_functions as be
from tqdm import tqdm
import argparse
from argparse import RawTextHelpFormatter

##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='Vector borne pathogen evolution simulation'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-selection', help='immune, adaptive, neutral, hybrid\n')
parser.add_argument('-gene', help='single, modular, multi')
parser.add_argument('-len', type=int, help='size of bit string representing gene')
parser.add_argument('-vec', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-yrs', type=int, help='number of years to simulate')
parser.add_argument('-gen_fit', type=str, default='high', help='should generalists have high or low fitness vals?; [low, high]')
parser.add_argument('-mut', type=float, default=0.01, help='per site mutation rate in antigen')
parser.add_argument('-rec', type=float, default = 0.01, help = 'recombination rate')
parser.add_argument('-replace', type=lambda x: (str(x).lower() == 'true'), default=False, help='should recombinants or mutated variants replace (defualt) the variant or be added')
parser.add_argument('-out', default=None, help='prefix for output file')
parser.add_argument('-run_tag', default = 1, help='unique id for each sim run when doing batch')
args = parser.parse_args()


# print sim parameters
print('Parameters','\n',
      'selection: ', args.selection,'\n',
      'gene type: ', args.gene,'\n',
      'gene length: ',args.len,'\n',
      'per site mutation rate: ', args.mut,'\n',
      'recombination rate: ', args.rec,'\n',
      'Vector pop size: ',args.vec,'\n',
      'host pop size: ', round(args.vec/50),'\n',
      'sim years: ', args.yrs,'\n',
      'batch progress (current run): ', args.run_tag)

##### initialize populations #####
ticks = be.Vector(pop_size= args.vec, n_strains= 1, strain_length= args.len, gene=args.gene, lam= 0.5)
host_pop_size = round(len(ticks.pop)/50)
hosts = be.Host(host_pop_size)

##### create immune selection probability tables #####
if args.selection == 'hybrid' and args.gene == 'modular':
  probabilities = be.distance_probabilities(list(range((args.len//2)+1)), curve= 'sigmoid')
elif args.selection in ['immune', 'hybrid'] and args.gene in ['single', 'multi']:
  probabilities = be.distance_probabilities(list(range(args.len+1)), curve= 'sigmoid')

##### initialize data collection #####
unique_lineage_ids = set()
temp_nymph_pop = [d for d in ticks.pop if d.get('stage') == 'nymph']
for tick in temp_nymph_pop:
  for strain in tick['strains']:
    unique_lineage_ids.add(strain['lineage_id'])

all_data = [{'year': ticks.year,
             'infection_rate': ticks.infection_rate, # only available with ticks.pop_stats()
             'avg_strains_carried': ticks.num_carried, # ^
             'unique_lineages': len(unique_lineage_ids),
             'avg_antigen_distance': ticks.avg_gen_dist, # only available with ticks.sample()
             'spec_weight': 0.0}] # ^


########## sim ##########

#for year in range(sim_years):
for year in tqdm(range(args.yrs)):
  
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
        if args.selection == 'immune':
          strains_from_tick = be.tick2host_transmission(tick=current_tick, 
                                                        host=current_host,
                                                        selection=args.selection,
                                                        gene=args.gene,
                                                        transmission_probabilities=probabilities)

        elif args.selection == 'adaptive':
          strains_from_tick = be.tick2host_transmission(tick=current_tick, 
                                                        host=current_host,
                                                        selection=args.selection,
                                                        gene=args.gene, 
                                                        host_type= current_host['host_type'], 
                                                        gen_fit=args.gen_fit)

        elif args.selection == 'neutral':
          strains_from_tick = be.tick2host_transmission(tick=current_tick, 
                                                        host=current_host, 
                                                        selection=args.selection,
                                                        gene=args.gene)
        
        elif args.selection == 'hybrid':
          strains_from_tick = be.tick2host_transmission(tick=current_tick, 
                                                        host=current_host,
                                                        selection=args.selection,
                                                        gene=args.gene,  
                                                        transmission_probabilities=probabilities,
                                                        host_type= current_host['host_type'], 
                                                        gen_fit=args.gen_fit)
        
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

  # update population
  ticks.update_pop()
  ticks.pop_stats()
  ticks.sample(args.gene)
  
  # collect data
  unique_lineage_ids = set()
  temp_nymph_pop = [d for d in ticks.pop if d.get('stage') == 'nymph']
  for tick in temp_nymph_pop:
    for strain in tick['strains']:
      unique_lineage_ids.add(strain['lineage_id'])

  all_data.append({'year': ticks.year,
                   'infection_rate': ticks.infection_rate, # only available with ticks.pop_stats()
                   'avg_strains_carried': ticks.num_carried, # ^
                   'unique_lineages': len(unique_lineage_ids),
                   'avg_antigen_distance': ticks.avg_gen_dist, # only available with ticks.sample()
                   'spec_weight': ticks.sampled_spec_weight}) # ^

  # mutate and recombination
  if year+1 != args.yrs:
    ticks.mutate(replace=args.replace, gene= args.gene, rate = args.mut)
    #ticks.recombination(replace=args.replace, gene=args.gene, rate = args.rec)

# sample the final population 
ticks.sample(args.gene)

##### data output #####
if args.out != None:
  import csv
  import pandas as pd
  from collections import Counter
  
  run = 'run_' + str(args.run_tag)

  ##### YEARLY SYSTEM STATS #####
  df = pd.DataFrame(all_data)
  df['run_id'] = run
  df['selection'] = args.selection
  df['gene_type'] = args.gene
  df['rec_rate'] = args.rec
  df['mut_rate'] = args.mut
  df['gen_fitness'] = args.gen_fit
  
  if str(args.run_tag) == "1":
    df.to_csv(args.out+'_sim_run_stats.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df.to_csv(args.out+'_sim_run_stats.tsv', mode= 'a', sep='\t', index=False, header=False)

  ##### FINAL VARIANT FREQUENCIES *SAMPLED* #####
  strain_counts = dict(Counter(ticks.sampled_strains))
  my_dict = []
  for strain in strain_counts:
    my_dict.append({'strain': strain, 'count': strain_counts[strain], 'freq': round(strain_counts[strain]/len(ticks.sampled_strains)*100,2)})
  df_sampled_strains = pd.DataFrame(my_dict)
  df_sampled_strains['run_id'] = run
  df_sampled_strains['selection'] = args.selection
  df_sampled_strains['gene_type'] = args.gene
  df_sampled_strains['rec_rate'] = args.rec
  df_sampled_strains['mut_rate'] = args.mut
  df_sampled_strains['gen_fitness'] = args.gen_fit
  
  if str(args.run_tag) == "1":
    df_sampled_strains.to_csv(args.out+'_variant_freqs_sampled.tsv', mode= 'w', sep='\t', index=False, header=True)
  else:
    df_sampled_strains.to_csv(args.out+'_variant_freqs_sampled.tsv', mode= 'a', sep='\t', index=False, header=False)
  
  ##### FINAL ADAPTIVE GENE FREQUENCIES *SAMPLED* #####
  if args.gene in ['multi']:
    strain_counts2 = dict(Counter(ticks.sampled_adaptive))
    my_dict = []
    for strain in strain_counts2:
      my_dict.append({'strain': strain, 'count': strain_counts2[strain], 'freq': round(strain_counts2[strain]/len(ticks.sampled_adaptive)*100,2)})
    df_sampled_adaptive_genes = pd.DataFrame(my_dict)
    df_sampled_adaptive_genes['run_id'] = run
    df_sampled_adaptive_genes['selection'] = args.selection
    df_sampled_adaptive_genes['gene_type'] = args.gene
    df_sampled_adaptive_genes['rec_rate'] = args.rec
    df_sampled_adaptive_genes['mut_rate'] = args.mut
    df_sampled_adaptive_genes['gen_fitness'] = args.gen_fit
    
    if str(args.run_tag) == "1":
      df_sampled_adaptive_genes.to_csv(args.out+'_adpt_gene_freqs_sampled.tsv', mode= 'w', sep='\t', index=False, header=True)
    else:
      df_sampled_adaptive_genes.to_csv(args.out+'_adpt_gene_freqs_sampled.tsv', mode= 'a', sep='\t', index=False, header=False)

    strain_counts3 = dict(Counter(ticks.sampled_combined))
    my_dict = []
    for strain in strain_counts3:
      my_dict.append({'strain': strain, 'count': strain_counts3[strain], 'freq': round(strain_counts3[strain]/len(ticks.sampled_adaptive)*100,2)})
    df_sampled_combined_genes = pd.DataFrame(my_dict)
    df_sampled_combined_genes['antigen'] = df_sampled_combined_genes['strain'].str[:args.len]
    df_sampled_combined_genes['adpt_gene'] = df_sampled_combined_genes['strain'].str[args.len:]    
    df_sampled_combined_genes['run_id'] = run
    df_sampled_combined_genes['selection'] = args.selection
    df_sampled_combined_genes['gene_type'] = args.gene
    df_sampled_combined_genes['rec_rate'] = args.rec
    df_sampled_combined_genes['mut_rate'] = args.mut
    df_sampled_combined_genes['gen_fitness'] = args.gen_fit
    
    if str(args.run_tag) == "1":
      df_sampled_combined_genes.to_csv(args.out+'_combined_gene_freqs_sampled.tsv', mode= 'w', sep='\t', index=False, header=True)
    else:
      df_sampled_combined_genes.to_csv(args.out+'_combined_gene_freqs_sampled.tsv', mode= 'a', sep='\t', index=False, header=False)


  ##### PAIRWISE ANTIGEN DISTANCES *SAMPLED* #####
  with open(args.out + '_pw_dists_sampled.tsv', 'a', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for value in ticks.sampled_distances:
      writer.writerow([run, value, args.selection, args.gene, args.rec, args.mut, args.gen_fit])
 
  ##### ADAPTIVE TRAIT VALUES *SAMPLED* #####
  if args.selection in ['adaptive', 'hybrid']:
    with open(args.out + '_adpt_vals_sampled.tsv', 'a', newline='') as file:
      writer = csv.writer(file, delimiter='\t')
      for value in ticks.sampled_adaptive_vals:
        writer.writerow([run, value, args.selection, args.gene, args.rec, args.mut, args.gen_fit])

  
  ##### FINAL VARIANT FREQUENCIES *ENTIRE POPULATION* #####
  # get variant counts in populaton 
  # strain_pop = []
  # adaptive_pop = []
  # for d in ticks.nymph_pop:
  #   if d['strains'] != []:
  #     for i in range(len(d['strains'])):
  #       strain_pop.append(d['strains'][i]['variant'])
  #       if args.gene == 'multi':
  #         adaptive_pop.append(d['strains'][i]['adaptive_gene'])
  
  # counts = dict(Counter(strain_pop))
  # total_counts = sum(counts.values())
  # frequencies = {key: round((value / total_counts) * 100, 2) for key, value in counts.items()}
  # df_strain_pop = pd.DataFrame(counts, index=[0])
  # df_strain_pop = df_strain_pop.transpose()
  # df_strain_pop.rename(columns={0: 'counts'}, inplace=True)
  # df_strain_pop['frequency'] = df_strain_pop.index.map(frequencies).fillna('')
  # df_strain_pop.reset_index(inplace=True)
  # df_strain_pop.rename(columns={'index': 'variant'}, inplace=True)
  # df_strain_pop['run_id'] = run

  # write file
  # if str(args.run_tag) == "1":
  #   df_strain_pop.to_csv(args.out+'_variant_freqs.tsv', mode= 'w', sep='\t', index=False, header=True)
  # else:
  #   df_strain_pop.to_csv(args.out+'_variant_freqs.tsv', mode= 'a', sep='\t', index=False, header=False)
  
  # if args.gene == 'multi':
  #   counts2 = dict(Counter(adaptive_pop))
  #   total_counts2 = sum(counts2.values())
  #   frequencies = {key: round((value / total_counts2) * 100, 2) for key, value in counts2.items()}
  #   df_strain_pop = pd.DataFrame(counts2, index=[0])
  #   df_strain_pop = df_strain_pop.transpose()
  #   df_strain_pop.rename(columns={0: 'counts'}, inplace=True)
  #   df_strain_pop['frequency'] = df_strain_pop.index.map(ticks.adaptive_strain_frequencies).fillna('')
  #   df_strain_pop.reset_index(inplace=True)
  #   df_strain_pop.rename(columns={'index': 'adaptive_gene'}, inplace=True)
  #   df_strain_pop['run_id'] = run

  #   write file
  #   if str(args.run_tag) == "1":
  #     df_strain_pop.to_csv(args.out+'_adpt_gene_freqs.tsv', mode= 'w', sep='\t', index=False, header=True)
  #   else:
  #     df_strain_pop.to_csv(args.out+'_adpt_gene_freqs.tsv', mode= 'a', sep='\t', index=False, header=False)

  ##### save lineage histories to output file #####
  # grab records
  # records = []
  # records2 = []
  # for tick in ticks.pop:
  #   if tick['strains'] != []:
  #     for i in range(len(tick['strains'])):
  #       records.append({'lineage_id': tick['strains'][i]['lineage_id'], 'history': tick['strains'][i]['imm_history'], 'history2': tick['strains'][i]['adp_history']})
  # new_records = []
  # for item in records:
  #   new_records.append(dict(zip(item['time'], item['history'])))

  # put into a dataframe
  # df = pd.DataFrame(new_records)
  # df = df.reindex(columns=range(sim_years)) # create complete set of columns for years
  # df = df.ffill(axis=1) # fill in empty rows with state of variant for that year
  # df.drop_duplicates(inplace=True) # drop all duplicate rows
  # df['run_id'] = run # create a column for run id and set to index
  # df.set_index('run_id', inplace=True)

  # # write file
  # if str(args.run_tag) == "1":
  #   df.to_csv(args.out+'_lineage_history.tsv', mode='a', sep='\t', index=False, header=True)
  # else:
  #   df.to_csv(args.out+'_lineage_history.tsv', mode='a', sep='\t', index=False, header=False)
  