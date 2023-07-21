'''
Created on 2014. 12. 8.

@author: user
'''
import os
import time

from cobra.io import read_sbml_model, write_sbml_model
from deap import base, creator, tools
import numpy
from pylab import *
import pylab
import matplotlib.pyplot as plt
from me_targeting.flux_analysis import Simulator
from gurobipy import *

class Single_Obj_GA(Simulator.Simulator):
    '''
    classdocs
    OptGene.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        print 'Single objective GA init'

        self.pop_size = 10  # Population size
        self.mut_rate = 0.2  # Mutation rate
        self.cross_rate = 0.5  # Crossover rate
        self.max_mutations = 5  # Maximum number of mutations
        self.max_generation = 100  # Maximum number of generations
        self.growth_threshold = 0.01  # Minimum cell survival threshold, wild type growth * 0.01
        self.best_offspring_num = 0
        self.timelimit = 60 * 60 * 48  # Limit of computation time

        self.wild_growth = 0.0
        self.wild_target_production = 0.0

        self.target_reaction = ''
        self.GrowthRxn = ''

        self.pre_knockout_genes = []
        self.essential_genes = []

        self.timelimit = 60
        self.target_lb = 0.0

        self.start_time = time.time()
        self.end_time = time.time()

        self.final_population_result = {}
        self.final_population = ''

        self.kd_ratio_list = []
        self.experimental_flux_data = {}

    def check_time_limit(self):
        print 'Elapsed time : ', self.end_time - self.start_time
        if self.end_time - self.start_time > self.timelimit:
            return True
        else:
            return False

    def read_no_target_genes(self, filename):
        self.essential_genes = []
        fp = open(filename)
        for line in fp:
            self.essential_genes.append(line.strip())
        fp.close()

    def set_elite_offspring_conservation_ratio(self, p=0.0):
        self.best_offspring_num = int(self.pop_size * p)

    def set_gene_manipulation_list(self, kd_ratio_list):
        self.kd_ratio_list = kd_ratio_list
        return

    def set_crossover_rate(self, cross_rate):
        self.cross_rate = cross_rate

    def set_population_size(self, pop_size):
        self.pop_size = pop_size

    def set_target_reaction(self, target_reaction):
        self.target_reaction = target_reaction
        self.target_idx = 1

    def set_mutation_rate(self, mut_rate):
        self.mut_rate = mut_rate

    def set_max_mutation_limit(self, max_mutation=5):
        self.max_mutations = max_mutation

    def set_max_generation_limit(self, max_generation):
        self.max_generation = max_generation

    def set_time_limit(self, timelimit):
        self.timelimit = timelimit
        return

    def set_ratio_of_minimum_target_production(self, lb=0.0):
        self.target_lb = lb
        return

    def set_survival_growth_thrshold(self, growth_threshold=0.01):
        self.growth_threshold = growth_threshold
        return

    def set_init_individual(self):
        import random
        random.seed()
        GeneMaxLength = len(self.ME_candidate_genes)
        rand_threshold = 1.0 / float(GeneMaxLength) * 0.1
        individual_gene = numpy.ones((1, GeneMaxLength), dtype=numpy.float64)
        individual_gene = individual_gene[0]
        idx_list = []
        mutatino_num = random.randint(0, self.max_mutations)
        for i in xrange(mutatino_num):
            idx = random.randint(0, len(individual_gene) - 1)
            idx_list.append(idx)

        for each_idx in idx_list:
            individual_gene[each_idx] = random.choice(self.kd_ratio_list)  # Gene KO
        return individual_gene

    def calc_wild_dist(self):
        pre_ko_reactions = self.get_reactions_from_multiple_gene_deletion(self.pre_knockout_genes)
        target_reactions = list(set(pre_ko_reactions))

        constraint_info = {}
        for each_reaction in target_reactions:
            constraint_info[each_reaction] = [0.0, 0.0]  # KO

        if self.experimental_flux_data == {}:
            Model_Stat, ObjVal, FluxDist = self.run_pFBA(new_objective=self.target_reaction,
                                                         flux_constraints=constraint_info)
            WildTargetProductionMin = FluxDist[self.target_reaction] * self.target_lb
            constraint_info[self.target_reaction] = [WildTargetProductionMin, 1000.0]
            # Calc wild flux distribution
            Model_Stat, ObjVal, FluxDist = self.run_pFBA(flux_constraints=constraint_info)
            self.wild_flux_dist = FluxDist
            WildGrowth = ObjVal
            WildTargetProduction = FluxDist[self.target_reaction]
            constraint_info[self.GrowthRxn] = [WildGrowth, WildGrowth]

            # Calc theorectical max production of target
            Model_Stat, ObjVal, FluxDist = self.run_pFBA(new_objective=self.target_reaction,
                                                         flux_constraints=constraint_info)
            self.wild_growth = WildGrowth
            self.wild_target_production = WildTargetProduction

        else:
            model_stat, objective_value, fba_fluxdist = self.run_pFBA(flux_constraints=constraint_info)
            for each_reaction in self.experimental_flux_data:
                constraint_info[each_reaction] = [self.experimental_flux_data[each_reaction],
                                                  self.experimental_flux_data[each_reaction]]
            model_stat, objective_value, fitted_fluxdist = self.run_LP_MOMA(constraint_info, fba_fluxdist)
            self.wild_flux_dist = fitted_fluxdist
            WildGrowth = fitted_fluxdist[self.GrowthRxn]
            WildTargetProduction = fitted_fluxdist[self.target_reaction]

        self.wild_growth = WildGrowth
        self.wild_target_production = WildTargetProduction

        print "Wild Growth : ", WildGrowth
        print "Wild Target Production : ", WildTargetProduction

    def remove_transporter_genes(self):
        Transporter_gene_list = []
        cobra_model = self.cobra_model
        for each_reaction in cobra_model.reactions:
            compartment_list = []
            for each_metabolite in each_reaction.get_reactants() + each_reaction.get_products():
                compartment_list.append(each_metabolite.compartment)
            if len(set(compartment_list)) > 1:
                for each_gene in each_reaction.get_gene():
                    Transporter_gene_list.append(each_gene.id)

        Transporter_gene_list = list(set(Transporter_gene_list))
        print 'No. of transporter related genes : ', len(Transporter_gene_list)
        return Transporter_gene_list

    def prepare_init_condition(self, GrowthRxn, TargetRxn, NotRemovedGenes=[], pre_knockout_genes=[]):
        self.GrowthRxn = GrowthRxn
        self.target_reaction = TargetRxn

        # Remove transporter related genes.
        cobra_model_genes = []
        self.pre_knockout_genes = pre_knockout_genes
        for each_gene in self.cobra_model.genes:
            cobra_model_genes.append(each_gene.id)

        Transporter_gene_list = self.remove_transporter_genes()
        multitarget_gene_list = self.remove_multitarget_genes()

        Rmoved_gene_list = multitarget_gene_list + Transporter_gene_list
        Rmoved_gene_list = list(set(Rmoved_gene_list))

        ME_candidate_genes = list(set(cobra_model_genes).difference(set(Rmoved_gene_list)))
        ME_candidate_genes = list(set(ME_candidate_genes).difference(set(pre_knockout_genes)))
        ME_candidate_genes = list(set(ME_candidate_genes).difference(set(self.essential_genes)))

        self.metabolic_gene_mapping = {}
        for i in xrange(len(ME_candidate_genes)):
            self.metabolic_gene_mapping[i] = ME_candidate_genes[i]
        self.ME_candidate_genes = ME_candidate_genes
        print 'No. of candidate genes : ', len(ME_candidate_genes)

        #
        self.calc_wild_dist()

    def eval_func_fba(self, individual):
        individual = individual[0]
        mutation_constraint_info = {}

        ko_index_list = list(numpy.where(individual != 1.0)[0])

        for idx in ko_index_list:
            target_gene = self.metabolic_gene_mapping[idx]
            mutation_constraint_info[target_gene] = individual[idx]

        if len(ko_index_list) > self.max_mutations:
            return (0.0,)

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
            if target_reaction_info[each_reaction] == 0.0:  # For KO
                constraint_info[each_reaction] = [0.0, 0.0]
            elif target_reaction_info[each_reaction] < 1.0:  # For KD
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction],
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction]]  # For Knockdown
            elif target_reaction_info[each_reaction] > 1.0:  # For Amplification
                constraint_info[each_reaction] = [
                    self.wild_flux_dist[each_reaction] * target_reaction_info[each_reaction], 1000.0]

        Model_Stat, ObjVal, FluxDist = self.run_pFBA(flux_constraints=constraint_info)  # Growth MAX

        if Model_Stat == 2:  # self.GrowthRxn
            Mutant_Growth = FluxDist[self.GrowthRxn]
            Mutant_Production = FluxDist[self.target_reaction]
            if Mutant_Growth < self.wild_growth * self.growth_threshold:
                return (0.0,)
            else:
                return (Mutant_Production,)
        else:
            return (0.0,)

    def eval_func_moma(self, individual):
        individual = individual[0]
        mutation_constraint_info = {}

        ko_index_list = list(numpy.where(individual != 1.0)[0])

        for idx in ko_index_list:
            target_gene = self.metabolic_gene_mapping[idx]
            mutation_constraint_info[target_gene] = individual[idx]

        if len(ko_index_list) > self.max_mutations:
            return (0.0,)

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
        #
        Model_Stat, ObjVal, FluxDist = self.run_LP_MOMA(flux_constraints=constraint_info,
                                                        wild_opt_flux=self.wild_flux_dist)  # Growth MAX
        if Model_Stat == 2:  # self.GrowthRxn
            Mutant_Growth = FluxDist[self.GrowthRxn]
            Mutant_Production = abs(FluxDist[self.target_reaction])
            if Mutant_Growth < self.wild_growth * self.growth_threshold:
                return (0.0,)
            else:
                return (Mutant_Production,)
        else:
            return (0.0,)

    def mutate_gene(self, individual):
        import random
        idx = random.randint(0, len(individual) - 1)
        original_value = individual[idx]

        individual[idx] = random.choice(list(set(self.kd_ratio_list).difference(set([original_value]))))  # Gene KO
        return individual

    def run_GA(self, method='moma'):
        import random
        random.seed()
        GA_Result_Container = {}

        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

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
        toolbox.register("select", tools.selRoulette)
        toolbox.register("bestselect", tools.selBest)

        CXPB = self.cross_rate
        MUTPB = self.mut_rate

        pop = toolbox.population(n=self.pop_size * 5)

        fitnesses = list(map(toolbox.evaluate, pop))

        fitness_sum = 0

        for each_fitness in fitnesses:
            fitness_sum = fitness_sum + each_fitness[0]

        for ind, fit in zip(pop, fitnesses):
            if fitness_sum == 0:
                ind.fitness.values = (0.0001,)
            else:
                ind.fitness.values = fit

        for i in xrange(self.max_generation):
            print '%d Generation' % (i + 1)
            GA_Result_Container[i + 1] = {}

            offspring = toolbox.select(pop, self.pop_size - self.best_offspring_num)
            bestoffspring = toolbox.bestselect(pop, self.best_offspring_num)

            offspring = offspring + bestoffspring
            offspring = list(map(toolbox.clone, offspring))

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
            fitness_sum = 0
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for each_fitness in fitnesses:
                fitness_sum = fitness_sum + each_fitness[0]

            for ind, fit in zip(invalid_ind, fitnesses):
                if fitness_sum == 0:
                    ind.fitness.values = (0.0001,)
                else:
                    ind.fitness.values = fit

            pop[:] = offspring

            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x * x for x in fits)
            std = abs(sum2 / length - mean ** 2) ** 0.5

            print("Maximum fitness value : %s" % max(fits))

            GA_Result_Container[i + 1]['MAX'] = max(fits)
            GA_Result_Container[i + 1]['MIN'] = min(fits)
            GA_Result_Container[i + 1]['AVG'] = mean
            GA_Result_Container[i + 1]['STD'] = std

            self.end_time = time.time()
            if self.check_time_limit() == True:
                break

        self.final_population = pop
        self.GA_Result_Container = GA_Result_Container
        self.save_final_population_result(pop)
        return

    def get_Fitness_info(self):
        return self.GA_Result_Container

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
            self.final_population_result[cnt]['Fitness'] = ind.fitness.values[0]
            self.final_population_result[cnt]['Genes'] = ','.join(target_genes)
            self.final_population_result[cnt]['Reactions'] = ','.join(ko_reactions)
            cnt += 1
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
        print >> fp, 'Wild type growth :', self.wild_growth
        print >> fp, 'Target reaction :', self.target_reaction
        print >> fp, 'Elapsed time :', self.end_time - self.start_time

        for each_key in self.GA_Result_Container:
            print >> fp, '%d\tgeneration\t MAX\t%s\tMIN\t%s\tMEAN\t%s\tSTD\t%s\t' % (
                each_key, self.GA_Result_Container[each_key]['MAX'], self.GA_Result_Container[each_key]['AVG'],
                self.GA_Result_Container[each_key]['MAX'], self.GA_Result_Container[each_key]['STD'])
        fp.close()

    def draw_convergence_curve(self, outputdir):
        fig = plt.figure()
        font = {'family': 'sans-serif',
                'color': 'black',
                'weight': 'normal',
                'size': 16
                }
        x = self.GA_Result_Container.keys()
        y = []
        for each_key in self.GA_Result_Container:
            y.append(self.GA_Result_Container[each_key]['MAX'])

        plt.plot(x, y, 'k')
        plt.title('Convergence curve', fontdict=font)
        plt.xlabel('Generation number', fontdict=font)
        plt.ylabel('Objective function', fontdict=font)

        plt.tight_layout()
        fig.savefig(outputdir + '/Convergence_Curve.png', dpi=80)

        pylab.close(fig)
        plt.clf()
        return

    def save_result(self, outputdir):
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        self.write_simulation_condition(outputdir)
        self.draw_convergence_curve(outputdir)
        self.write_final_population(outputdir)
        return

    def write_best_cobra_model(self, ModelNum=1, outputdir=''):
        result_list = []
        for each_key in self.final_population_result:
            fitness_score = self.final_population_result[each_key]['Fitness']
            gene_list = self.final_population_result[each_key]['KO Genes'].split(',')
            reaction_list = self.final_population_result[each_key]['KO Reactions'].split(',')
            result_list.append([fitness_score, reaction_list, gene_list])
        result_list.sort()
        result_list.reverse()
        fp = open(outputdir + '/Model Info.txt', 'w')
        cnt = 1
        for each_result in result_list[:ModelNum]:
            new_model = copy.deepcopy(self.cobra_model)
            fitness = str(each_result[0]).strip()

            rm_rxns_candidates = each_result[1]
            rm_genes_candidates = each_result[2]

            rm_rxns = []
            for each_reaction in rm_rxns_candidates:
                idx = new_model.reactions.index(each_reaction)
                rm_rxns.append(new_model.reactions[idx])
            new_model.remove_reactions(rm_rxns)
            write_sbml_model(new_model, outputdir + '/Output_%d_Fitness_%s_model.xml' % (cnt, fitness), use_fbc_package=False)
            print >> fp, 'Output_%d_Fitness_%s_model.xml\t%s\t%s' % (
                cnt, fitness, rm_rxns_candidates, rm_genes_candidates)
            cnt += 1
        fp.close()
        return

    def validation_kd_target(self, kd_gene_list=[]):
        pre_ko_reactions = self.get_reactions_from_multiple_gene_deletion(self.pre_knockout_genes)
        ko_target_reactions = pre_ko_reactions

        target_reactions = self.get_reactions_from_multiple_gene_deletion(kd_gene_list)
        target_reactions = list(set(target_reactions))

        constraint_info = {}
        for each_reaction in pre_ko_reactions:
            constraint_info[each_reaction] = [0, 0]  # KO

        target_info_dic = {}
        target_info_dic[tuple(kd_gene_list)] = {'reaction': target_reactions, 'mode': 'a'}

        df, Result_Dic = self.run_flux_tuning(target_info_dic=target_info_dic, wild_flux=self.wild_flux_dist,
                                              flux_constraints=constraint_info, manipulation_level=0.2)  # Growth MAX

        target_production_list = []

        for each_key in Result_Dic:
            Mutant_Growth = Result_Dic[each_key][self.GrowthRxn]
            Target = Result_Dic[each_key][self.target_reaction]
            target_production_list.append(Target)
            print  'Condition:\t%s\tGrowth :\t%s\tProduct\t%s' % (each_key, Mutant_Growth, Target)

    def set_experimental_fluxdata(self, experimental_flux_dist):
        self.experimental_flux_data = experimental_flux_dist

    def remove_multitarget_genes(self):
        multitarget_gene_list = []
        cobra_model = self.cobra_model
        for each_gene in cobra_model.genes:
            if len(each_gene.get_reaction()) >= 5:  # 5, 9.5%
                multitarget_gene_list.append(each_gene.id)

        multitarget_gene_list = list(set(multitarget_gene_list))
        print 'No. of multitarget related genes : ', len(multitarget_gene_list)
        return multitarget_gene_list


if __name__ == '__main__':
    experimental_flux_dist = {}

    obj = Single_Obj_GA()  # Genetic algorithm
    obj.read_model('./RawData/Jungeun_iAF1260_4hpl.xml')
    obj.set_ratio_of_minimum_target_production(0.01)  # 0.0 ~ 1.0
    obj.set_elite_offspring_conservation_ratio(0.0)  # 0.0 ~ 1.0
    obj.set_gene_manipulation_list(
        [0.5, 1.0])  # For KO simulation : [0.0, 1.0] , For KD & KO simulation : [0.0, 0.5, 1.0]
    obj.set_survival_growth_thrshold(0.01)
    obj.set_experimental_fluxdata(experimental_flux_dist)

    # Set Genetic Algorithm Parameters
    obj.set_population_size(100)
    obj.set_crossover_rate(0.7)  #
    obj.set_mutation_rate(0.2)  #
    obj.set_max_mutation_limit(5)
    obj.set_max_generation_limit(10)
    obj.set_time_limit(60 * 60 * 48)

    # Run GA
    obj.read_no_target_genes(
        './RawData/Not_Removed_genes.txt')  # essential genes or excluded genes as target candidates
    obj.prepare_init_condition(GrowthRxn='Ec_biomass_iAF1260_core_59p81M', TargetRxn='EX_4hpl_LPAREN_e_RPAREN_')
    obj.run_GA(method="moma")
    obj.save_result('150305_Result_Moma_3COMPLEX_1')
