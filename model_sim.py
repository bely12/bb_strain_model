# vector-host-pathogen simulation
import numpy as np
import random
import model_functions as be
from tqdm import tqdm
import argparse
from argparse import RawTextHelpFormatter

##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='Vector borne pathogen evolution simulation'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-mode', type=int, help='1: cross reactivity, 2: host specialization, 3:random fitness, 4: uniform fitness value\n')
parser.add_argument('-antigen_length', type=int, help='size of bit string representing antigen')
parser.add_argument('-mut', type=float, default=0.01, help='per site mutation rate in antigen')
parser.add_argument('-vector_pop_size', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-years', type=int, help='number of years to simulate')
parser.add_argument('-plots', default=None, help='set to True if you want plots as output; specify file namein args.plot_file')
parser.add_argument('-out', default=None, help='prefix for output file')
parser.add_argument('-run_tag', default = 1, help='unique id for each sim run when doing batch')
args = parser.parse_args()

##### set key parameters #####
vector_pop_size = args.vector_pop_size # needs to be an even number, at least 50
antigen_length = args.antigen_length
mutation_rate = args.mut
transmission_curve = 'sigmoid' # choose sigmoid or linear
sim_years = args.years
if args.mode == 1:
  cross_reactivity = True
  host_specialization = False
  uniform_fitness_value = False
if args.mode == 2:
  cross_reactivity = False
  host_specialization = True
  uniform_fitness_value = False
if args.mode == 3:
  cross_reactivity = False
  host_specialization = False
  uniform_fitness_value = False
if args.mode == 4:
  cross_reactivity = False
  host_specialization = False
  uniform_fitness_value = True

# print sim parameters; comment out if running batches
print('Parameters','\n',
      'Vector pop size: ',vector_pop_size,'\n',
      'host pop size: ', round(vector_pop_size/50),'\n',
      'antigen length: ',antigen_length,'\n',
      'per site mutation rate: ', mutation_rate,'\n'
      'cross reactivity: ', cross_reactivity,'\n',
      'host specialization: ', host_specialization,'\n',
      'uniform fitness value: ', uniform_fitness_value, '\n',
      'simulated years: ', sim_years)



##### initialize populations #####
ticks = be.Vector(pop_size= vector_pop_size, n_strains= 1, strain_length= antigen_length, lam= 0.5)
host_pop_size = round(len(ticks.pop)/50)
hosts = be.Host(host_pop_size, host_specialization)



##### create transmission probability tables #####
if cross_reactivity == True:
  probabilities = be.distance_probabilities(list(range(antigen_length+1)), curve= transmission_curve)
if host_specialization == True:
  mnp_values = be.mnp_probabilities(ticks.strain_set)
  mnp_values2 = {item['strain']: {'mnp_value': item['mnp_value'], 'rodent': item['rodent'], 'bird': item['bird']} for item in mnp_values}
if cross_reactivity == False and host_specialization == False:
  fitness_values = dict(zip(ticks.strain_set, np.round(np.random.uniform(0.0, 1.0, len(ticks.strain_set)),4)))



##### initialize data containers ##### * for batch running, use pairwise_antigen_distances rather than big_data
#strain_freqs = [{'year': ticks.year,'strain_community': ticks.strain_set_frequencies}] #don't use this, it blows the sim up
big_data = [{'year': ticks.year,
             'infection_rate': ticks.infection_rate,
             'diversity': ticks.diversity,
             'antigenic_distance': 0,
             'active_strains': ticks.current_strain_count,
             'avg_carried': ticks.avg_carried}]
pairwise_antigen_distances = []

if host_specialization == True:
  mnp_pop_values = [ticks.mnp_recs(mnp_values2)]



##### sim ##### don't use tqdm if doing batches? 
#for year in range(sim_years):
for i in tqdm(range(sim_years)):
  interactions = ticks.interaction(hosts.pop) #schedule current years vector-host interactions 

  for day in range(1, 151): # go through the days in the year
    for j in range(len(interactions)): # for any bites that day...
      if interactions[j]['bite_day'] == day:

        current_host = next((item for item in hosts.pop if item["id"] == interactions[j]['host_id']), None) # invovled host
        current_tick = next((item for item in ticks.pop if item["id"] == interactions[j]['tick_id']), None) # involved vector
        bb = be.Pathogen(current_tick['strains']) # initialize pathogen transmission community 

        for i in range(15): # genetic algorithm for evolving pathogens; n rounds of mutation and selection
          bb.mutate(rate= mutation_rate)
          if cross_reactivity == True:
            bb.selection(host_infections=current_host['infection'], values=probabilities, cross_reactivity= True)
          if host_specialization == True:
            bb.selection(host_infections=current_host['infection'], host_type= current_host['host_type'], values=mnp_values2, MNP=True)
          if cross_reactivity == False and host_specialization == False and uniform_fitness_value == False:
            bb.selection(host_infections=current_host['infection'], values=fitness_values)
          if uniform_fitness_value == True:
            bb.selection(host_infections=current_host['infection'], values=fitness_values, uniform_fitness= True)

        # determine strains to be transmitted from tick to host
        if cross_reactivity == True:
          strains_from_tick = bb.tick2host_transmission(current_host['infection'], values=probabilities, cross_reactivity= True)
        if cross_reactivity == False and host_specialization == False and uniform_fitness_value == False:
          strains_from_tick = bb.tick2host_transmission(current_host['infection'], values=fitness_values)
        if host_specialization == True:
          strains_from_tick = bb.tick2host_transmission(current_host['infection'], host_type= current_host['host_type'], values=mnp_values2, MNP=True)
        if uniform_fitness_value == True:
          strains_from_tick = bb.tick2host_transmission(current_host['infection'], values=fitness_values, uniform_fitness=True)

        strains_from_host = bb.host2tick_transmission(current_host['infection']) # determine strains to be transmitted from host to tick

        # update host and vector info
        for item in hosts.pop: 
          if item["id"] == current_host['id']:
            item["infection"] = current_host['infection'] + [item for item in strains_from_tick if (item not in current_host['infection'])]
            break
        for item in ticks.pop:
          if item['id'] == current_tick['id']:
            item['strains'] = current_tick['strains'] + [item for item in strains_from_host if item not in current_tick['strains']]
            break

  # update populations
  ticks.update_pop()
  hosts = be.Host(host_pop_size, host_specialization)

  # collect data
  #strain_freqs.append({'year': ticks.year,'strain_community': ticks.strain_set_frequencies})
  big_data.append({'year': ticks.year,
                   'infection_rate': ticks.infection_rate,
                   'diversity': ticks.diversity,
                   'antigenic_distance': ticks.avg_gen_dist,
                   'active_strains': ticks.current_strain_count,
                   'avg_carried': ticks.avg_carried})
  pairwise_antigen_distances.append(ticks.antigen_distances) 
  if host_specialization == True:
    mnp_pop_values.append(ticks.mnp_recs(mnp_values2))


##### data output #####
if args.out != None:
  import csv
  import pandas as pd
  run_tag = 'run_' + str(args.run_tag)
  
  # list of final strains at end of sim 
  with open(args.out + '_sim_final_antigens.tsv', 'a', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for strain in ticks.current_strains:
      writer.writerow([run_tag, strain])
  
  # pairwise antigen distances for histogram
  all_distances = [value for sublist in pairwise_antigen_distances for value in sublist]
  sampled_distances = random.sample(all_distances, 1000)
  with open(args.out + '_sampled_pairwise_antigen_dists.tsv', 'a', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    #for value in all_distances:
    for value in sampled_distances:
      writer.writerow([run_tag,value])
  
  # all mnp values for histogram *turn on if doing single run and/or testing; not relevant for batch
  if host_specialization == True:
    all_mnp = [value for sublist in mnp_pop_values for value in sublist]
    with open(args.out + '_sim_mnp_pop_values.tsv', 'w', newline='') as file:
      writer = csv.writer(file, delimiter='\t')
      for value in all_mnp:
        writer.writerow([value])

  # main data
  df = pd.DataFrame(big_data)
  df.to_csv(args.out+'_sim_data.tsv', sep='\t', index=False, header=True)
  
  # strain frequencies by year *** blows up sim, keep shut off ***
  # rows = []
  # for entry in strain_freqs:
  #   year = entry['year']
  #   for strain, frequency in entry['strain_community'].items():
  #     rows.append({'year': year, 'strain': strain, 'frequency': frequency})
  # df = pd.DataFrame(rows)
  # df.to_csv(args.out+'_sim_strain_freqs.tsv', sep='\t', index=False, header=True)


##### plots ##### this is mainly for testing purposes 
if args.plots == 'True':
  
  # import packages for plotting
  import pandas as pd
  import matplotlib.pyplot as plt
  from scipy.cluster.hierarchy import dendrogram, linkage
  from scipy.spatial.distance import pdist
  from matplotlib.backends.backend_pdf import PdfPages
  
  with PdfPages(args.out + '_sim_plots.pdf') as pdf:
    # dummy plot for parameters text
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('off')
    sim_parameters = (
        f"Parameters:\n"
        f"Vector pop size: {vector_pop_size}\n"
        f"Host pop size: {host_pop_size}\n"
        f"Antigen length: {antigen_length}\n"
        f"Mutation rate: {mutation_rate}\n"
        f"Cross reactivity: {cross_reactivity}\n"
        f"Host specialization: {host_specialization}\n"
        f"Unifrom fitness value: {uniform_fitness_value}\n"
        f"Simulated years: {sim_years}"
    )
    ax.text(0.1, 0.9, sim_parameters, fontsize=12, verticalalignment='top', horizontalalignment='left')
    pdf.savefig(fig)
    plt.close(fig)

    # strain frequencies by year
    # if args.out == None:
    #   rows = []
    #   for entry in strain_freqs:
    #     year = entry['year']
    #     for strain, frequency in entry['strain_community'].items():
    #       rows.append({'year': year, 'strain': strain, 'frequency': frequency})
    #   df = pd.DataFrame(rows)
    # strain_frequencies = df.groupby(['year', 'strain'])['frequency'].mean().unstack()
    # strain_frequencies.plot(kind='line', legend=False)
    # plt.xlabel('Year')
    # plt.ylabel('Frequency')
    # plt.title('Strain Frequencies')
    # pdf.savefig()
    # plt.close()

    # plot changes in avg antigen distance in population
    years = [d['year'] for d in big_data]
    avg_pop_dists = [d['antigenic_distance'] for d in big_data]
    plt.plot(years, avg_pop_dists)
    plt.xlabel('Year')
    plt.ylabel('Avg distance')
    plt.title('Average pairwise distance of antigens')
    plt.grid(True)
    pdf.savefig()
    plt.close()

    # plot histogram of pairwise antigen distances
    all_distances = [value for sublist in pairwise_antigen_distances for value in sublist]
    plt.hist(all_distances, bins=antigen_length, density=True, alpha=0.7, color='blue')
    plt.title('Antigenic distances')
    plt.xlabel('pairwise distance')
    plt.ylabel('density')
    pdf.savefig()
    plt.close()

    # cluster the final strains into groups
    data = np.array([[int(bit) for bit in string] for string in ticks.current_strains])
    distance_matrix = pdist(data, metric='hamming')
    linked = linkage(distance_matrix, method='ward')
    plt.figure(figsize=(10, 7))
    dendrogram(linked, orientation='top', leaf_rotation=90)
    plt.title('Antigen Clusters')
    plt.xlabel('Antigens')
    plt.ylabel('Distance')
    pdf.savefig() 
    plt.close()

    # plot simpson diversity index
    diversity_values = [d['diversity'] for d in big_data]
    plt.plot(years, diversity_values)
    plt.xlabel('Year')
    plt.ylabel('Simpson Diversity Index')
    plt.title('Simpson Diversity Index Over Time')
    plt.grid(True)
    pdf.savefig()
    plt.close()

    # plot histogram of mnp values
    if host_specialization == True:
      all_mnp = [value for sublist in mnp_pop_values for value in sublist]
      plt.hist(all_mnp, bins=10, density=True)
      plt.title('Host Specialization Distribution')
      plt.xlabel('mnp value')
      plt.ylabel('density')
      pdf.savefig()
      plt.close()

    # plot yearly active strain counts
    active_strains = [entry['active_strains'] for entry in big_data]
    plt.plot(years, active_strains)
    plt.xlabel('Year')
    plt.ylabel('Variants')
    plt.title('Number of variants in population')
    pdf.savefig() 
    plt.close()

    # average number of strains carried by an infected tick
    n = [d['avg_carried'] for d in big_data]
    plt.plot(years, n)
    plt.xlabel('Year')
    plt.ylabel('number of strains')
    plt.title('Average strains carried by infected nymphs')
    plt.grid(True)
    pdf.savefig()
    plt.close()

    # plot infection rate
    infection_rates = [d['infection_rate'] for d in big_data]
    plt.plot(years, infection_rates)
    plt.xlabel('Year')
    plt.ylabel('Infection Rate')
    plt.title('Nymph Infection Rate Over Time')
    plt.grid(True)
    pdf.savefig()
    plt.close()