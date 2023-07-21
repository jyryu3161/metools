'''
Created on 2014. 12. 8.

@author: user
'''

import sys

import pp

import Multi_Objective_GA
from gurobipy import *

def do_work(title, crossover_rate, mutation_rate, num_target):
    # TODO: parameter optimization.
    obj = Multi_Objective_GA.Multi_Obj_GA()
    obj.read_model('./20150216_CacMBELv160_SimulationModel.xml')
    obj.set_ratio_of_minimum_target_production(0.0)  # 0.0 ~ 1.0 (%)
    obj.set_gene_manipulation_list([0.0, 1.0])  # For KO simulation : [0.0, 1.0] , For KD simulation : [0.5, 1.0]
    obj.set_survival_growth_thrshold(0.05)

    # Genetic Algorithm Parameters
    obj.set_population_size(120)
    obj.set_crossover_rate(crossover_rate)  # 0 ~ 0.5
    obj.set_mutation_rate(mutation_rate)  #
    obj.set_max_mutation_limit(num_target)
    obj.set_max_generation_limit(10)
    obj.set_multi_obj(('EX_Biomass', 'EX_BUOH_Ext_'), (1.0, 1.0))

    obj.set_time_limit(60 * 60 * 60)
    # obj.read_essential_genes('./RawData/Not_Removed_genes.txt') # essential genes or excluded genes as target candidates
    obj.prepare_init_condition(GrowthRxn='EX_Biomass', TargetRxn='EX_BUOH_Ext_')
    obj.run_GA(method="moma")
    obj.save_result('%s_GA_moma_KDKO_C%s_M%s_N%s' % (title, crossover_rate, mutation_rate, num_target))

    return


if __name__ == '__main__':
    ppservers = ()
    if len(sys.argv) > 1:
        ncpus = int(sys.argv[1])
        # Creates jobserver with ncpus workers
        job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"

    # The following submits 8 jobs and then retrieves the results
    jobs = []
    for num_target in [3, 4]:
        for crossover_rate in [0.4, 0.5, 0.6]:
            for mutation_rate in [0.2, 0.3]:
                title = '20150227'
                description = 'GA_%s_N%s_C%s_M%s' % (title, num_target, crossover_rate, mutation_rate)
                jobs.append((description, job_server.submit(func=do_work, args=(
                    description, crossover_rate, mutation_rate, num_target), depfuncs=(),
                                                            modules=("Multi_Objective_GA",))))

    for input, job in jobs:
        print "Multi-objective", input, "is", job()

    job_server.print_stats()
