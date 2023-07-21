'''
Created on 2014. 12. 8.

@author: user
'''
import os
import time

from deap import base, creator, tools
import numpy
from pylab import *

import Single_Objective_GA
from gurobipy import *

class Multi_Obj_GA(Single_Objective_GA.Single_Obj_GA):
    '''
    classdocs
    OptGene.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        print 'Multi objective GA init'

        Single_Objective_GA.Single_Obj_GA.__init__(self)

        self.obj_rxns = ''
        self.weights = ''

    # self.negative_data_value
    # self.small_value_tuple

    def set_multi_obj(self, objrxn, weights):
        self.obj_rxns = objrxn
        self.weights = weights
        negative_data_set = []
        for each_weight in weights:
            if each_weight > 0:
                negative_data_set.append(0.0)
            elif each_weight < 0:
                negative_data_set.append(10000)

        self.negative_data_value = tuple(negative_data_set)

        return

    def eval_func_fba(self, individual):
        individual = individual[0]
        mutation_constraint_info = {}

        ko_index_list = list(numpy.where(individual != 1.0)[0])

        for idx in ko_index_list:
            target_gene = self.metabolic_gene_mapping[idx]
            mutation_constraint_info[target_gene] = individual[idx]

        if len(ko_index_list) > self.max_mutations:
            return self.negative_data_value

        pre_ko_reactions = self.get_reactions_from_multiple_gene_deletion(self.pre_knockout_genes)

        # For KO, KD
        target_reaction_info = {}
        for each_gene in mutation_constraint_info:  # get_reactions_from_multiple_genes
            if mutation_constraint_info[each_gene] < 1.0:
                target_reactions = self.get_reactions_from_multiple_gene_deletion([each_gene])
                if len(target_reactions) > 0:
                    for each_target_reaction in target_reactions:
                        target_reaction_info[each_target_reaction] = mutation_constraint_info[each_gene]
            else:
                target_reactions = self.get_reactions_from_multiple_genes_amp([each_gene])  # amplification
                if len(target_reactions) > 0:
                    for each_target_reaction in target_reactions:
                        target_reaction_info[each_target_reaction] = mutation_constraint_info[each_gene]


        ## KO constraint
        pre_ko_reactions = list(set(pre_ko_reactions))
        constraint_info = {}
        for each_reaction in pre_ko_reactions:
            constraint_info[each_reaction] = [0.0, 0.0]  # KO

        ## Flux Constraint
        for each_reaction in target_reaction_info:
            if target_reaction_info[each_reaction] == 0.0:
                constraint_info[each_reaction] = [0.0, 0.0]  # For KO
            elif target_reaction_info[each_reaction] < 1.0:
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction],
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction]]  # For Knockdown
            elif target_reaction_info[each_reaction] > 1.0:
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction],
                    1000.0]  # For amplification

        Model_Stat, ObjVal, FluxDist = self.run_FBA(flux_constraints=constraint_info)  # Growth MAX

        if Model_Stat == 2:  # self.GrowthRxn
            Mutant_Growth = FluxDist[self.GrowthRxn]
            Mutant_Production = FluxDist[self.target_reaction]
            if Mutant_Growth < self.wild_growth * self.growth_threshold:
                return self.negative_data_value
            else:
                target_obj_rxn_list = []
                for each_target_rxn in self.obj_rxns:
                    target_obj_rxn_list.append(abs(FluxDist[each_target_rxn]))
                return tuple(target_obj_rxn_list)
        else:
            return self.negative_data_value

    def eval_func_moma(self, individual):
        individual = individual[0]
        mutation_constraint_info = {}

        ko_index_list = list(numpy.where(individual != 1.0)[0])

        for idx in ko_index_list:
            target_gene = self.metabolic_gene_mapping[idx]
            mutation_constraint_info[target_gene] = individual[idx]

        if len(ko_index_list) > self.max_mutations:
            return self.negative_data_value

        pre_ko_reactions = self.get_reactions_from_multiple_gene_deletion(self.pre_knockout_genes)

        # For KO, KD
        target_reaction_info = {}
        for each_gene in mutation_constraint_info:  # get_reactions_from_multiple_genes
            if mutation_constraint_info[each_gene] < 1.0:
                target_reactions = self.get_reactions_from_multiple_gene_deletion([each_gene])
                if len(target_reactions) > 0:
                    for each_target_reaction in target_reactions:
                        target_reaction_info[each_target_reaction] = mutation_constraint_info[each_gene]
            else:
                target_reactions = self.get_reactions_from_multiple_genes_amp([each_gene])  # amplification
                if len(target_reactions) > 0:
                    for each_target_reaction in target_reactions:
                        target_reaction_info[each_target_reaction] = mutation_constraint_info[each_gene]


        ## KO constraint
        pre_ko_reactions = list(set(pre_ko_reactions))
        constraint_info = {}
        for each_reaction in pre_ko_reactions:
            constraint_info[each_reaction] = [0.0, 0.0]  # KO

        ## Flux Constraint
        for each_reaction in target_reaction_info:
            if target_reaction_info[each_reaction] == 0.0:
                constraint_info[each_reaction] = [0.0, 0.0]  # For KO
            elif target_reaction_info[each_reaction] < 1.0:
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction],
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction]]  # For Knockdown
            elif target_reaction_info[each_reaction] > 1.0:
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction],
                    1000.0]  # For amplification

        Model_Stat, ObjVal, FluxDist = self.run_MOMA(wild_flux=self.wild_flux_dist, flux_constraints=constraint_info)  # Growth MAX
        if Model_Stat == 2:  # self.GrowthRxn
            Mutant_Growth = FluxDist[self.GrowthRxn]
            Mutant_Production = FluxDist[self.target_reaction]
            if Mutant_Growth < self.wild_growth * self.growth_threshold:
                return self.negative_data_value
            else:
                target_obj_rxn_list = []
                for each_target_rxn in self.obj_rxns:
                    target_obj_rxn_list.append(abs(FluxDist[each_target_rxn]))
                return tuple(target_obj_rxn_list)
        else:
            return self.negative_data_value

    def run_GA(self, method='moma'):
        import random

        random.seed()
        GA_Result_Container = {}

        creator.create("Fitness", base.Fitness, weights=self.weights)  # Multi objective
        creator.create("Individual", numpy.ndarray, fitness=creator.Fitness)

        toolbox = base.Toolbox()
        toolbox.register("attr_bool", self.set_init_individual)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        if method == 'moma':
            toolbox.register("evaluate", self.eval_func_moma)  # moma
        elif method == 'fba':
            toolbox.register("evaluate", self.eval_func_fba)  # fba evaluation. single target max.
        else:
            toolbox.register("evaluate", self.eval_func_fba)  # fba evaluation. single target max.

        toolbox.register("mate", tools.cxTwoPoints)
        toolbox.register("kd_mutate", self.mutate_gene)
        toolbox.register("select", tools.selNSGA2)  # # Multi objective

        CXPB = self.cross_rate
        MUTPB = self.mut_rate

        pop = toolbox.population(n=self.pop_size * 5)
        fitnesses = list(map(toolbox.evaluate, pop))

        fitness_sum = 0
        for each_fitness in fitnesses:
            for each_value in each_fitness:
                fitness_sum = fitness_sum + each_value

        for ind, fit in zip(pop, fitnesses):
            if fitness_sum == 0:
                ind.fitness.values = self.negative_data_value
            else:
                ind.fitness.values = fit

        pop = toolbox.select(pop, self.pop_size)

        for i in xrange(self.max_generation):
            print '%d Generation' % (i + 1)
            GA_Result_Container[i + 1] = {}

            #offspring = tools.selTournamentDCD(pop, len(pop))
            #offspring = list(map(toolbox.clone, offspring))

            offspring = [toolbox.clone(ind) for ind in pop]

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB:
                    toolbox.mate(child1[0], child2[0])
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                if random.random() < MUTPB:
                    toolbox.kd_mutate(mutant[0])
                    del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            fitness_sum = 0

            for each_fitness in fitnesses:
                for each_value in each_fitness:
                    fitness_sum = fitness_sum + each_value

            for ind, fit in zip(invalid_ind, fitnesses):
                if fitness_sum == 0:
                    ind.fitness.values = self.negative_data_value
                else:
                    ind.fitness.values = fit

            pop = toolbox.select(pop + offspring, self.pop_size)

            pop_optimal_objective_list = []
            for j in xrange(len(self.weights)):
                temp_fitness_list = []
                for ind in pop:
                    temp_fitness_list.append(ind.fitness.values[j])
                if self.weights[j] < 0.0:
                    pop_optimal_objective_list.append(sum(temp_fitness_list) / float(self.pop_size))
                elif self.weights[j] > 0.0:
                    pop_optimal_objective_list.append(sum(temp_fitness_list) / float(self.pop_size))

            print 'Average fitness values : ', pop_optimal_objective_list

            for j in xrange(len(self.weights)):
                GA_Result_Container[i + 1]['%d_Objective_AVG\t' % (j + 1)] = pop_optimal_objective_list[j]

            self.end_time = time.time()
            if self.check_time_limit() == True:
                break

        self.final_population = pop
        self.GA_Result_Container = GA_Result_Container
        self.save_final_population_result(pop)

    def save_final_population_result(self, pop):
        self.final_population_result = {}
        cnt = 1
        for ind in pop:
            Genes = ind[0]
            ko_index_list = list(numpy.where(Genes != 1.0)[0])
            rm_gene_list = []
            target_genes = []
            for index in ko_index_list:
                ko_gene = self.metabolic_gene_mapping[int(index)]
                ko_gene2 = '%s(%s)' % (ko_gene, Genes[index])
                rm_gene_list.append(ko_gene)
                target_genes.append(ko_gene2)

            rm_gene_list = list(set(rm_gene_list))
            ko_reactions = self.get_reactions_from_multiple_gene_deletion(rm_gene_list)
            ko_reactions = list(set(ko_reactions))

            self.final_population_result[cnt] = {}
            self.final_population_result[cnt]['Fitness'] = ind.fitness.values
            self.final_population_result[cnt]['Genes'] = ','.join(target_genes)
            self.final_population_result[cnt]['Reactions'] = ','.join(ko_reactions)
            cnt += 1
        return

    def save_result(self, outputdir):
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        self.write_simulation_condition(outputdir)
        self.write_final_population(outputdir)
        return

    def write_final_population(self, outputdir):
        fp = open(outputdir + '/Final_Result.txt', 'w')
        for each_idx in self.final_population_result:
            fitness = self.final_population_result[each_idx]['Fitness']
            genes = self.final_population_result[each_idx]['Genes']
            reactions = self.final_population_result[each_idx]['Reactions']
            print >> fp, '%d \t Fitness \t %s \t genes \t %s \t reactions \t %s' % (each_idx, fitness, genes, reactions)
        fp.close()
        return

    def write_simulation_condition(self, outputdir):
        fp = open(outputdir + '/Simulation_Condition.txt', 'w')
        print >> fp, 'Population size :', self.pop_size
        print >> fp, 'No. of total candidate targets : ', len(self.ME_candidate_genes)
        print >> fp, 'Mutation rate :', self.mut_rate
        print >> fp, 'Limit of mutations :', self.max_mutations
        print >> fp, 'Objective reactions :', self.obj_rxns
        print >> fp, 'Objective weights :', self.weights
        print >> fp, 'Elapsed time :', self.end_time - self.start_time

        for each_key in self.GA_Result_Container:
            print >> fp, '%d\tgeneration\t' % (each_key),
            for each_obj_key in self.GA_Result_Container[each_key]:
                print >> fp, '%s\t%s\t' % (each_obj_key, self.GA_Result_Container[each_key][each_obj_key]),
            print >> fp
        fp.close()


if __name__ == '__main__':
    # TODO: parameter optimization.
    obj = Multi_Obj_GA()
    obj.read_model('./Koma_EFICAz_target_model_eco_v3.xml')
    obj.set_ratio_of_minimum_target_production(0.01)  # 0.0 ~ 1.0 (%)
    obj.set_gene_manipulation_list([0.0, 1.0])  # For KO simulation : [0.0, 1.0] , For KD simulation : [0.5, 1.0]
    obj.set_survival_growth_thrshold(0.05)

    # Genetic Algorithm Parameters
    obj.set_population_size(120)
    obj.set_crossover_rate(0.7)  # 0 ~ 0.5
    obj.set_mutation_rate(0.2)  #
    obj.set_max_mutation_limit(3)
    obj.set_max_generation_limit(5000)
    obj.set_multi_obj(('Ec_biomass_iAF1260_core_59p81M', 'EX_cellulose_LPAREN_e_RPAREN_'), (1.0, 1.0))

    obj.set_time_limit(60 * 60 * 60)
    # obj.read_essential_genes('./RawData/Not_Removed_genes.txt') # essential genes or excluded genes as target candidates
    obj.prepare_init_condition(GrowthRxn='Ec_biomass_iAF1260_core_59p81M', TargetRxn='EX_cellulose_LPAREN_e_RPAREN_')
    obj.run_GA(method="moma")
    obj.save_result('20160604_multi_optgene')


