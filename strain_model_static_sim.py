import strain_model_static_functions as bb
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter

### ARGUMENTS ###
parser = argparse.ArgumentParser(
    description='Vector borne pathogen evolution simulation; takes defined sets of 2-3 strains as a community and tests stability'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-n_strains', type=int, help='how many strains in the community; 2 or 3')

def list_of_strings(arg):
    return arg.split(',')
parser.add_argument('-spec_types', type=list_of_strings)

parser.add_argument('-vector_pop_size', type=int, help='size of vector(tick) population; needs to be at least 50 and even number')
parser.add_argument('-rodent_pop_size', type=int, help='number of rodents in host pop')
parser.add_argument('-bird_pop_size', type=int, help='number of birds in host pop')
parser.add_argument('-lam', default=0.7, help='infection rate lambda val for poisson distribution')

parser.add_argument('-years', type=int, help='number of years to simulate')
#parser.add_argument('-replicates', type=int, default=1, help='number of times to run sim with specified conditions')
parser.add_argument('-run_tag', default=1, help='for running in batches')
parser.add_argument('-out', default='results', help='prefix for output file')

args = parser.parse_args()



### define parameters ###
n_strains = args.n_strains
spec_types = args.spec_types
vector_pop_size = args.vector_pop_size
rodent_pop_size = args.rodent_pop_size
bird_pop_size = args.bird_pop_size
sim_years = args.years
#repeat = args.replicates
#run_tag = args.run_tag
print('working on replicate simulation ',args.run_tag)

### start new replication ###
# for i in range(repeat):
#     run_tag = run_tag + 1
#     print('working on replicate simulation ',run_tag,' out of ',repeat)

### initialize pathogen, vector, host populations ###
patho = bb.Pathogen(n = n_strains,
            hs = bb.assign_hs_v2(spec_types),
            antigen = bb.assign_antigen(n_strains))

ticks = bb.Vector(pop_size=vector_pop_size,
                strains=patho.strain_community,
                lam=args.lam)

hosts = bb.Host(n_rodents= rodent_pop_size,
            n_birds= bird_pop_size)


### data collection ###
sim_hs_vals = []
sim_antigen_vals = []
for i in range(len(patho.strain_community)):
    sim_hs_vals.append(patho.strain_community[i]['hs'])
    sim_antigen_vals.append(patho.strain_community[i]['antigen'])

if n_strains == 2:
    antigen_variation = round(bb.antigen_distance(sim_antigen_vals[0], sim_antigen_vals[1]),2)
if n_strains == 3:
    antigen_variation = round(bb.dispersion(sim_antigen_vals),2)

spec_weight_val = bb.spec_weight(patho.strain_community)

run_data = {'run_tag': args.run_tag,
        'v_pop': vector_pop_size,
        'rodent_pop': rodent_pop_size,
        'bird_pop': bird_pop_size,
        'hs_vals': sim_hs_vals,
        'antigen_vals': sim_antigen_vals,
        'antigen_variation': antigen_variation,
        'spec_weight': spec_weight_val}


### start sim ###
for year in range(sim_years):
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
                strains_from_tick = bb.tick2host_transmission(tick=current_tick, host=current_host)

                # determine strains to be transmitted from host to tick
                strains_from_host = bb.host2tick_transmission(current_host)

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
    
    #kill the sim once pop is deemed unstable
    seen = []
    for i in range(len(ticks.nymph_pop)):
        if len(seen) == len(patho.strain_community):
            break
        if ticks.nymph_pop[i]['strains'] != []:
            for j in range(len(ticks.nymph_pop[i]['strains'])):
                if ticks.nymph_pop[i]['strains'][j] not in seen:
                    seen.append(ticks.nymph_pop[i]['strains'][j])
    if len(seen) != len(patho.strain_community):
        break


# identify what strains are left after sim finishes
seen = []
for i in range(len(ticks.nymph_pop)):
    if ticks.nymph_pop[i]['strains'] != []:
        for j in range(len(ticks.nymph_pop[i]['strains'])):
            if ticks.nymph_pop[i]['strains'][j] not in seen:
                seen.append(ticks.nymph_pop[i]['strains'][j])

# identify the trait values of remaining strains
if len(seen) != 0:
    final_sim_hs_vals = []
    final_sim_antigen_vals = []
    for i in range(len(seen)):
        final_sim_hs_vals.append(seen[i]['hs'])
        final_sim_antigen_vals.append(seen[i]['antigen'])
else:
    final_sim_antigen_vals = 'none'
    final_sim_hs_vals = 'none'

# label sim outcome as stable or unstable
if len(seen) == len(patho.strain_community):
    result = "stable"
else:
    result = "unstable"

# add results data to dictionary
run_data.update({'final_hs': final_sim_hs_vals,
                'final_ant': final_sim_antigen_vals,
                'result': result})

# write results to file
df = pd.DataFrame([run_data])
if args.run_tag == "1":
    df.to_csv(args.out+'_static_sim_results.tsv', mode= 'w', sep='\t', index=False, header=True)
else:
    df.to_csv(args.out+'_static_sim_results.tsv', mode= 'a', sep='\t', index=False, header=False)
