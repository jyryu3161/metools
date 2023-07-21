'''
Created on 2014. 12. 8.

@author: user
'''

import sys

import pp

import Single_Objective_GA
from gurobipy import *

def do_work(title, crossover_rate, mutation_rate, num_target):
    filename = 'ref_flux.txt'
    
    obj = Single_Objective_GA.Single_Obj_GA()  # Genetic algorithm
    obj.read_model('./iJO1366_Modify_GPR_Simple_GPR_RM_BLOCKED_15dap_V2.xml')
    obj.set_reference_flux(filename)
    
    obj.set_gene_manipulation_list(
        [0.0, 0.5, 1.0])  # For KO simulation : [0.0, 1.0] , For KD & KO simulation : [0.0, 0.5, 1.0]
    obj.set_survival_growth_thrshold(0.05)
    # set_wild_type_flux_for_moma

    # Set Genetic Algorithm Parameters
    obj.set_population_size(120)
    obj.set_crossover_rate(crossover_rate)  #
    obj.set_mutation_rate(mutation_rate)  #
    obj.set_max_mutation_limit(num_target)
    obj.set_max_generation_limit(1000)
    obj.set_time_limit(60 * 60 * 60)

    # Run GA    
    obj.prepare_init_condition(GrowthRxn='Ec_biomass_iJO1366_core_53p95M', TargetRxn='EX_15dap_LPAREN_e_RPAREN_')
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
    for each_iter in range(1):
        for num_target in [1, 3]:
            for crossover_rate in [0.6]:
                for mutation_rate in [0.2, 0.3]:
                    title = '20150302_moma_'
                    description = 'Iteration_%s_GA_%s_N%s_C%s_M%s' % (
                        each_iter, title, num_target, crossover_rate, mutation_rate)

                    jobs.append((description, job_server.submit(func=do_work, args=(
                        description, crossover_rate, mutation_rate, num_target), depfuncs=(),
                                                                modules=("Single_Objective_GA",))))

    for input, job in jobs:
        print "Sum of primes below", input, "is", job()

    job_server.print_stats()
