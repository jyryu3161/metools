'''
Created on 2014. 12. 8.

@author: user
'''

from multiprocessing import Process

import Single_Objective_GA
from gurobipy import *

def do_work(crossover_rate, mutation_rate, num_target):
    obj = Single_Objective_GA.Single_Obj_GA()  # Genetic algorithm
    obj.read_model('./RawData/Jaeho_iAF1260.xml')
    obj.set_ratio_of_minimum_target_production(0.1)  # 0.0 ~ 1.0
    obj.set_elite_offspring_conservation_ratio(0.0)  # 0.0 ~ 1.0
    obj.set_gene_manipulation_list(
        [0.0, 0.5, 1.0])  # For KO simulation : [0.0, 1.0] , For KD & KO simulation : [0.0, 0.5, 1.0]
    obj.set_survival_growth_thrshold(0.05)
    # set_wild_type_flux_for_moma

    # Set Genetic Algorithm Parameters
    obj.set_population_size(100)
    obj.set_crossover_rate(crossover_rate)  #
    obj.set_mutation_rate(mutation_rate)  #
    obj.set_max_mutation_limit(num_target)
    obj.set_max_generation_limit(5000)
    obj.set_time_limit(60 * 60 * 48)

    # Run GA
    obj.read_no_target_genes(
        './RawData/Not_Removed_genes.txt')  # essential genes or excluded genes as target candidates
    obj.prepare_init_condition(GrowthRxn='Ec_biomass_iAF1260_core_59p81M', TargetRxn='EX_4AMINOBUTANOL')
    obj.run_GA(method="moma")
    obj.save_result('150226_Result_moma_KDKO_C%s_M%s_N%s' % (crossover_rate, mutation_rate, num_target))
    return


def do_work2(crossover_rate, mutation_rate, num_target):
    obj = Single_Objective_GA.Single_Obj_GA()  # Genetic algorithm
    obj.read_model('./RawData/Jaeho_iAF1260.xml')
    obj.set_ratio_of_minimum_target_production(0.1)  # 0.0 ~ 1.0
    obj.set_elite_offspring_conservation_ratio(0.0)  # 0.0 ~ 1.0
    obj.set_gene_manipulation_list(
        [0.0, 0.5, 1.0])  # For KO simulation : [0.0, 1.0] , For KD & KO simulation : [0.0, 0.5, 1.0]
    obj.set_survival_growth_thrshold(0.05)
    # set_wild_type_flux_for_moma

    # Set Genetic Algorithm Parameters
    obj.set_population_size(10)
    obj.set_crossover_rate(crossover_rate)  #
    obj.set_mutation_rate(mutation_rate)  #
    obj.set_max_mutation_limit(num_target)
    obj.set_max_generation_limit(10)
    obj.set_time_limit(60 * 60 * 48)

    # Run GA
    obj.read_no_target_genes(
        './RawData/Not_Removed_genes.txt')  # essential genes or excluded genes as target candidates
    obj.prepare_init_condition(GrowthRxn='Ec_biomass_iAF1260_core_59p81M', TargetRxn='EX_4AMINOBUTANOL')
    obj.run_GA(method="moma")
    obj.save_result('150227_Result_moma_KDKO_C%s_M%s_N%s' % (crossover_rate, mutation_rate, num_target))
    return


if __name__ == '__main__':
    # Multi-processing.
    print 'START'
    import time

    s = time.time()
    process_cnt = 0
    process_list = []
    for num_target in [3, 4]:
        for crossover_rate in [0.5, 0.6, 0.7]:
            for mutation_rate in [0.2, 0.3]:
                each_process = Process(target=do_work, args=(crossover_rate, mutation_rate, num_target))
                process_list.append(each_process)
            #                 process_cnt+=1
            #                 if process_cnt == 4:
            #                     for each_process in process_list:
            #                         each_process.start()
            #                     for each_process in process_list:
            #                         each_process.join()
            #                     process_list=[]
            #                     process_cnt=0

    for each_process in process_list:
        each_process.start()
    for each_process in process_list:
        each_process.join()

    print 'END'
    e = time.time()
    print e - s
