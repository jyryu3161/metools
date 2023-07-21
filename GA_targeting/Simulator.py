'''
Created on 2014. 6. 4.

@author: user
'''

import copy

import numpy as np
from pandas import DataFrame
import cobra.io as io

from gurobipy import *

class Simulator(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def get_simulation_condition(self):

        return

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

    def set_target_reaction(self, target_reaction):
        self.target_reaction = target_reaction
        self.target_idx = 1

    def set_max_KD(self, max_KD=3):
        self.max_KD = max_KD

    def prepare_init_condition(self, TargetRxn, GrowthRxn, NotRemovedGenes=[], pre_knockout_genes=[]):
        self.growth_reaction = GrowthRxn
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
        ME_candidate_genes = list(set(ME_candidate_genes).difference(set(self.task_genes)))
        if len(self.sRNA_list) > 0:
            ME_candidate_genes = list(set(ME_candidate_genes) & set(self.sRNA_list))

        self.metabolic_gene_mapping = {}
        for i in xrange(len(ME_candidate_genes)):
            self.metabolic_gene_mapping[i] = ME_candidate_genes[i]

        self.ME_candidate_genes = ME_candidate_genes
        print 'No. of sRNA genes : ', len(self.sRNA_list)
        print 'No. of candidate genes : ', len(ME_candidate_genes)

        self.calc_wild_dist()

    def set_ratio_of_minimum_target_production(self, lb=0.0):
        self.target_lb = lb
        return

    def calc_wild_dist(self):
        pre_ko_reactions = self.get_reactions_from_multiple_gene_deletion(self.pre_knockout_genes)
        target_reactions = list(set(pre_ko_reactions))

        constraint_info = {}
        for each_reaction in target_reactions:
            constraint_info[each_reaction] = [0.0, 0.0]  # KO

        if self.experimental_flux_data == {}:
            Model_Stat, ObjVal, FluxDist = self.run_FBA(new_objective=self.target_reaction,
                                                        flux_constraints=constraint_info)
            WildTargetProductionMin = FluxDist[self.target_reaction] * self.target_lb
            constraint_info[self.target_reaction] = [WildTargetProductionMin, 1000.0]
            # Calc wild flux distribution
            Model_Stat, ObjVal, FluxDist = self.run_FBA(flux_constraints=constraint_info)
            self.wild_flux_dist = FluxDist
            WildGrowth = ObjVal
            WildTargetProduction = FluxDist[self.target_reaction]
            constraint_info[self.growth_reaction] = [WildGrowth, WildGrowth]

            # Calc theorectical max production of target
            Model_Stat, ObjVal, FluxDist = self.run_FBA(new_objective=self.target_reaction,
                                                        flux_constraints=constraint_info)
            self.wild_growth = WildGrowth
            self.wild_target_production = WildTargetProduction

        else:
            model_stat, objective_value, fba_fluxdist = self.run_FBA(flux_constraints=constraint_info)
            for each_reaction in self.experimental_flux_data:
                constraint_info[each_reaction] = [self.experimental_flux_data[each_reaction],
                                                  self.experimental_flux_data[each_reaction]]
            model_stat, objective_value, fitted_fluxdist = self.run_MOMA(constraint_info, fba_fluxdist)
            self.wild_flux_dist = fitted_fluxdist
            WildGrowth = fitted_fluxdist[self.growth_reaction]
            WildTargetProduction = fitted_fluxdist[self.target_reaction]

        Model_Stat, ObjVal, FluxDist = self.run_FBA(new_objective=self.target_reaction,
                                                    flux_constraints=constraint_info)
        WildTheorecticalMaxProduction = FluxDist[self.target_reaction]

        self.wild_growth = WildGrowth
        self.wild_target_production = WildTargetProduction

        print "Wild Growth : ", WildGrowth
        print "Wild Target Production : ", WildTargetProduction
        print "Wild Theoretical Target Max Production : ", WildTheorecticalMaxProduction

    def read_essential_genes(self, filename):
        self.task_genes = []
        fp = open(filename)
        for line in fp:
            self.task_genes.append(line.strip())
        fp.close()

    def read_sRNA_list(self, filename):
        self.sRNA_list = []
        fp = open(filename)
        for line in fp:
            self.sRNA_list.append(line.strip())
        fp.close()

    def run_flux_tuning(self, target_info_dic, wild_flux, flux_constraints={}, manipulation_level=0.2,
                        method='moma', inf_flag=False):
        import itertools
        print 'Start Manipulate_flux ... '
        targets = target_info_dic.keys()

        each_target_info = {}
        manipulation_combination = []
        all_target_reactions = []
        index_cnt = 0
        for each_target in targets:
            target_reactions = target_info_dic[each_target]['reaction']
            all_target_reactions += target_reactions
            mode = target_info_dic[each_target]['mode']
            each_target_info[each_target] = [target_reactions]
            if mode == 'a':
                manipulation_ratio_list = list(np.arange(manipulation_level, 1, manipulation_level))
                manipulation_ratio_list = np.around(manipulation_ratio_list, 3)
            elif mode == 'i':
                manipulation_ratio_list = list(np.arange(1, 2 + manipulation_level, manipulation_level))
                manipulation_ratio_list = np.around(manipulation_ratio_list, 3)
            print manipulation_ratio_list
            manipulation_combination.append(manipulation_ratio_list)
            each_target_info[each_target] = [target_reactions, index_cnt]
            index_cnt += 1

        all_target_reactions = list(set(all_target_reactions))
        all_manipulation_candidates = list(itertools.product(*manipulation_combination))

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        objective = self.objective

        Smatrix = self.Smatrix

        LowerFlux = self.Lower_Boundary_Constraint
        UpperFlux = self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux.keys():
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux.keys():
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('Manipulate_flux_sum')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        v = {}
        for each_reaction in model_reactions:
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
                                            ub=flux_constraints[each_reaction][1], name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
                                            name=each_reaction)
        m.update()

        # Mass balance
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
                                 pairs.select(each_metabolite, '*')) == 0)
        m.update()

        Result_Dic = {}

        for manipulation_candidate in all_manipulation_candidates:
            reaction_flux_ratio_info = {}
            for i in range(len(each_target_info)):
                target = targets[i]
                index = each_target_info[target][1]
                k_info = manipulation_candidate[index]
                target_reactions = each_target_info[target][0]
                for each_reaction in target_reactions:
                    reaction_flux_ratio_info[each_reaction] = k_info

            const_reaction_info_dic = {}
            addr_list = []
            for each_target_reaction in all_target_reactions:
                const_reaction_info_dic[each_target_reaction] = reaction_flux_ratio_info[each_target_reaction] * \
                                                                wild_flux[each_target_reaction]

            m.update()
            # m.setObjective(quicksum( (v[each_target_reaction] - const_reaction_info_dic[each_target_reaction])*(v[each_target_reaction] - const_reaction_info_dic[each_target_reaction]) for each_target_reaction in all_target_reactions ), GRB.MINIMIZE)
            # m.optimize()
            # moma_status = m.status
            # print m.status

            moma_const_dic = {}
            moma_status = 2
            if moma_status == 2:
                for each_reaction in all_target_reactions:
                    moma_const_dic[each_reaction] = const_reaction_info_dic[each_reaction]

                for each_addr in addr_list:
                    removeConstraintIndex = m.getConstrs().index(each_addr)
                    m.remove(m.getConstrs()[removeConstraintIndex])
                m.update()

                for each_target_reaction in all_target_reactions:
                    addr = m.addConstr(v[each_target_reaction] == (moma_const_dic[each_target_reaction]))
                    addr_list.append(addr)
                m.update()

                if method == 'fba':
                    m.setObjective(v[objective], GRB.MAXIMIZE)
                elif method == 'moma':
                    m.setObjective(quicksum((v[each_reaction] - wild_flux[each_reaction]) * (
                        v[each_reaction] - wild_flux[each_reaction]) for each_reaction in model_reactions),
                                   GRB.MINIMIZE)
                else:
                    m.setObjective(v[objective], GRB.MAXIMIZE)

                m.update()
                m.optimize()
                if m.status == 2:
                    simulation_obj = m.ObjVal
                    print 'Objective : ', simulation_obj
                    each_flux_dist = {}
                    for each_reaction in model_reactions:
                        each_flux_dist[each_reaction] = v[each_reaction].x

                    # make key
                    Item_Key = ''
                    for each_key in reaction_flux_ratio_info:
                        each_value = reaction_flux_ratio_info[each_key]
                        Item_Key = Item_Key + '[%s_%.2f]#' % (each_key, each_value)

                    Result_Dic[Item_Key] = each_flux_dist

            for each_addr in addr_list:
                removeConstraintIndex = m.getConstrs().index(each_addr)
                m.remove(m.getConstrs()[removeConstraintIndex])

            m.update()
        df = DataFrame.from_dict(Result_Dic)
        return df.T, Result_Dic

    def run_pFBA(self, new_objective='', flux_constraints={}, inf_flag=False, mode='max'):
        internal_min_flag = True
        # print 'Start simple FBA simulation ... '

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        if new_objective == '':
            objective = self.objective
        else:
            objective = new_objective

        Smatrix = self.Smatrix

        LowerFlux = self.Lower_Boundary_Constraint
        UpperFlux = self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux.keys():
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux.keys():
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('FBA')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        v = {}
        fplus = {}
        fminus = {}

        m.update()

        for each_reaction in model_reactions:
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
                                            ub=flux_constraints[each_reaction][1], name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
                                            name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
        m.update()

        # Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
                                 pairs.select(each_metabolite, '*')) == 0)

        m.update()
        if mode == 'max':
            m.setObjective(v[objective], GRB.MAXIMIZE)
        elif mode == 'min':
            m.setObjective(v[objective], GRB.MINIMIZE)

        m.optimize()
        print m.status
        if m.status == 2:
            obj_value = m.ObjVal
            print 'Objective value : ', m.ObjVal
            if internal_min_flag == True:
                m.addConstr(fplus[objective] - fminus[objective] == obj_value)

                m.addConstr(quicksum(
                    (fplus[reaction] - fminus[reaction]) * coffvalue[metabolite, reaction] for metabolite, reaction in
                    pairs.select(each_metabolite, '*')) == 0)

                for each_reaction in model_reactions:
                    m.addConstr(fplus[each_reaction] - fminus[each_reaction] == v[each_reaction])

                m.update()
                m.setObjective(
                    quicksum((fplus[each_reaction] + fminus[each_reaction]) for each_reaction in model_reactions),
                    GRB.MINIMIZE)

                m.optimize()
                if m.status == 2:
                    obj_value = m.ObjVal
                    print 'Flux minimization objective value : ', m.ObjVal
                    ReactionFlux = {}
                    for reaction in model_reactions:
                        ReactionFlux[reaction] = float(v[reaction].x)
                        if abs(float(v[reaction].x)) <= 1e-6:
                            ReactionFlux[reaction] = 0.0
                    return m.status, obj_value, ReactionFlux
            else:
                ReactionFlux = {}
                for reaction in model_reactions:
                    ReactionFlux[reaction] = float(v[reaction].x)
                    if abs(float(v[reaction].x)) <= 1e-6:
                        ReactionFlux[reaction] = 0.0

                return m.status, obj_value, ReactionFlux
        return m.status, False, False

    def run_MOMA(self, wild_flux={}, flux_constraints={}, objective_reaction_list=[], inf_flag=False,
                 norm='Euclidean'):
        print 'run moma'
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        LowerFlux = self.Lower_Boundary_Constraint
        UpperFlux = self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux:
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux:
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('MOMA')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        v = {}

        for each_reaction in model_reactions:
            if each_reaction in flux_constraints:
                v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
                                            ub=flux_constraints[each_reaction][1], name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
                                            name=each_reaction)

        m.update()

        # Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
                                 pairs.select(each_metabolite, '*')) == 0)

        m.update()

        if norm == 'Euclidean':
            m.setObjective(quicksum(
                (v[each_reaction] - wild_flux[each_reaction]) * (v[each_reaction] - wild_flux[each_reaction])
                for each_reaction in model_reactions), GRB.MINIMIZE)

            m.optimize()


        if m.status == 2:

            ReactionFlux = {}
            for reaction in model_reactions:
                ReactionFlux[reaction] = float(v[reaction].x)
                if abs(float(v[reaction].x)) <= 1e-6:
                    ReactionFlux[reaction] = 0.0
                # print reaction, float(v[reaction].x)

            return m.status, m.ObjVal, ReactionFlux
        else:
            return m.status, False, False

    def run_LP_MOMA(self, wild_flux={}, flux_constraints={}, inf_flag=False):
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        LowerFlux = self.Lower_Boundary_Constraint
        UpperFlux = self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux.keys():
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux.keys():
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('Flux Fitting')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        target_reactions = wild_flux.keys()
        v = {}
        fplus = {}
        fminus = {}

        for each_reaction in model_reactions:
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
                                            ub=flux_constraints[each_reaction][1], name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
                                            name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)

        m.update()

        for each_reaction in model_reactions:
            m.addConstr(v[each_reaction] == (fplus[each_reaction] - fminus[each_reaction]))
            m.addConstr(fplus[each_reaction], GRB.GREATER_EQUAL, v[each_reaction] - wild_flux[each_reaction],
                        name=each_reaction)
            m.addConstr(fminus[each_reaction], GRB.GREATER_EQUAL, wild_flux[each_reaction] - v[each_reaction],
                        name=each_reaction)

        m.update()

        # Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(
                (fplus[reaction] - fminus[reaction]) * coffvalue[metabolite, reaction] for metabolite, reaction in
                pairs.select(each_metabolite, '*')) == 0)

        m.update()

        m.setObjective(quicksum(
            ((fplus[each_reaction] + fminus[each_reaction]) - wild_flux[each_reaction]) for each_reaction in
            target_reactions), GRB.MINIMIZE)

        m.optimize()

        # print 'Model Status : ', m.status
        if m.status == 2:
            # print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in model_reactions:
                ReactionFlux[reaction] = float(v[reaction].x)
                if abs(float(v[reaction].x)) <= 1e-6:
                    ReactionFlux[reaction] = 0.0

            return m.status, m.ObjVal, ReactionFlux
        else:
            return m.status, False, False

    def run_FBA(self, new_objective='', flux_constraints={}, inf_flag=False, internal_min_flag=False, mode='max'):
        # print 'Start simple FBA simulation ... '

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        if new_objective == '':
            objective = self.objective
        else:
            objective = new_objective

        Smatrix = self.Smatrix

        LowerFlux = self.Lower_Boundary_Constraint
        UpperFlux = self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux.keys():
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux.keys():
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = Model('FBA')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        v = {}
        fplus = {}
        fminus = {}

        m.update()

        for each_reaction in model_reactions:
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
                                            ub=flux_constraints[each_reaction][1], name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
                                            name=each_reaction)
            fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
            fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
        m.update()

        # Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
                                 pairs.select(each_metabolite, '*')) == 0)

        m.update()
        if mode == 'max':
            m.setObjective(v[objective], GRB.MAXIMIZE)
        elif mode == 'min':
            m.setObjective(v[objective], GRB.MINIMIZE)

        m.optimize()
        print m.status
        if m.status == 2:
            obj_value = m.ObjVal
            print 'Objective value : ', m.ObjVal
            if internal_min_flag == True:
                m.addConstr(fplus[objective] - fminus[objective] == obj_value)

                m.addConstr(quicksum(
                    (fplus[reaction] - fminus[reaction]) * coffvalue[metabolite, reaction] for metabolite, reaction in
                    pairs.select(each_metabolite, '*')) == 0)

                for each_reaction in model_reactions:
                    m.addConstr(fplus[each_reaction] - fminus[each_reaction] == v[each_reaction])

                m.update()
                m.setObjective(
                    quicksum((fplus[each_reaction] + fminus[each_reaction]) for each_reaction in model_reactions),
                    GRB.MINIMIZE)
                m.optimize()
                if m.status == 2:
                    obj_value = m.ObjVal
                    print 'Flux minimization objective value : ', m.ObjVal
                    ReactionFlux = {}
                    for reaction in model_reactions:
                        ReactionFlux[reaction] = float(v[reaction].x)
                        if abs(float(v[reaction].x)) <= 1e-6:
                            ReactionFlux[reaction] = 0.0
                    return m.status, obj_value, ReactionFlux
            else:
                ReactionFlux = {}
                for reaction in model_reactions:
                    ReactionFlux[reaction] = float(v[reaction].x)
                    if abs(float(v[reaction].x)) <= 1e-6:
                        ReactionFlux[reaction] = 0.0

                return m.status, obj_value, ReactionFlux
        return m.status, False, False

    def read_model(self, filename, model_calc_bool=True):
        ## save model information
        model = io.sbml.create_cobra_model_from_sbml_file(filename)
        return self.load_cobra_model(model, model_calc_bool)

    def load_cobra_model(self, cobra_model, model_calc_bool=True):
        self.cobra_model = cobra_model
        model = cobra_model
        model.optimize()
        # print('\nSimulated growth rate is %1.3f' % model.solution.f)
        model_metabolites = []
        model_reactions = []
        Lower_Boundary_Constraint = {}
        Upper_Boundary_Constraint = {}
        objective_reaction = ''
        for each_metabolite in model.metabolites:
            model_metabolites.append(str(each_metabolite.id))

        Smatrix = {}

        for each_reaction in model.reactions:
            if each_reaction.objective_coefficient == 1.0:
                objective_reaction = str(each_reaction.id)

            reactant_list = each_reaction.get_reactants()
            reactant_coff_list = each_reaction.get_coefficients(reactant_list)
            product_list = each_reaction.get_products()
            product_coff_list = each_reaction.get_coefficients(product_list)

            for i in range(len(reactant_list)):
                Smatrix[(str(reactant_list[i].id), str(each_reaction.id))] = reactant_coff_list[i]

            for i in range(len(product_list)):
                Smatrix[(str(product_list[i].id), str(each_reaction.id))] = product_coff_list[i]

            model_reactions.append(str(each_reaction.id))
            lb = each_reaction.lower_bound
            ub = each_reaction.upper_bound
            if lb < -1000.0:
                lb = float('-inf')
            if ub > 1000.0:
                ub = float('inf')
            Lower_Boundary_Constraint[str(each_reaction.id)] = lb
            Upper_Boundary_Constraint[str(each_reaction.id)] = ub

        self.model_metabolites = model_metabolites
        self.model_reactions = model_reactions
        self.Smatrix = Smatrix
        self.Lower_Boundary_Constraint = Lower_Boundary_Constraint
        self.Upper_Boundary_Constraint = Upper_Boundary_Constraint
        self.objective = objective_reaction

        if model_calc_bool == True:
            self.update_model_info(model)

        return (model_metabolites, model_reactions, Smatrix, Lower_Boundary_Constraint, Upper_Boundary_Constraint,
                objective_reaction)

    def update_model_info(self, model):
        self.model_reaction_info = {}
        self.model_genes = []
        self.model_gene_info = {}
        self.model_metabolite_target_Reaction_info = {}

        self.subsystem = []
        self.subsystem_reaction_info = {}

        self.model_genes = model.genes
        self.metabolite_list = []

        for reaction in model.reactions:
            # print reaction
            self.model_reaction_info[str(reaction.id)] = {}
            self.model_reaction_info[str(reaction.id)]['GPR_str'] = reaction.gene_reaction_rule
            self.model_reaction_info[str(reaction.id)]['GPR_list'] = self.convertGPRstringToListFormat(
                reaction.gene_reaction_rule)
            self.model_reaction_info[str(reaction.id)]['genes'] = reaction.get_gene()
            self.model_reaction_info[str(reaction.id)]['products'] = reaction.get_products()
            self.model_reaction_info[str(reaction.id)]['reactants'] = reaction.get_reactants()

            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.metabolite_list:
                    self.metabolite_list.append(str(each_reactant))

            for each_product in self.model_reaction_info[str(reaction.id)]['products']:
                if each_product not in self.metabolite_list:
                    self.metabolite_list.append(str(each_product))

            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.model_metabolite_target_Reaction_info.keys():
                    self.model_metabolite_target_Reaction_info[str(each_reactant)] = [str(reaction.id)]
                else:
                    self.model_metabolite_target_Reaction_info[str(each_reactant)].append(str(reaction.id))

            self.model_reaction_info[str(reaction.id)]['subsystem'] = str(reaction.subsystem)

            if str(reaction.subsystem) not in self.subsystem_reaction_info.keys():
                self.subsystem_reaction_info[str(reaction.subsystem)] = [str(reaction.id)]
            else:
                self.subsystem_reaction_info[str(reaction.subsystem)].append(str(reaction.id))

            self.subsystem.append(str(reaction.subsystem))
            for gene in reaction.get_gene():
                str_gene = str(gene)
                if str_gene not in self.model_gene_info.keys():
                    self.model_gene_info[str_gene] = [str(reaction.id)]
                else:
                    self.model_gene_info[str_gene].append(str(reaction.id))

        self.metabolite_list = list(set(self.metabolite_list))

        self.model = model.to_array_based_model()
        self.subsystem = list(set(self.subsystem))
        return model

    def convertGPRstringToListFormat(self, line):
        calcNewList = []
        geneList = []
        line = line.strip()
        calcList = line.split('or')
        for c in calcList:
            c = c.replace('(', '')
            c = c.replace(')', '')
            c = c.replace(' ', '')
            c = c.strip()
            if 'and' in c:
                newlist = c.split('and')
                newlist = list(set(newlist))
                newlist.sort()
                calcNewList.append(newlist)
            else:
                geneid = c.strip()
                if geneid not in calcNewList:
                    calcNewList.append(geneid)
        return calcNewList

    def get_reactions_from_multiple_genes_amp(self, gene_list):
        reaction_list = []
        reaction_candidates = []
        for each_gene in gene_list:
            gene_name = str(each_gene)
            tmp_reaction_candidates = self.model_gene_info[gene_name]
            reaction_candidates+=tmp_reaction_candidates

        reaction_candidates = list(set(reaction_candidates))

        for reaction in reaction_candidates:
            boolean_list = []
            GPR_list = self.model_reaction_info[reaction]['GPR_list']
            for each_item in GPR_list:
                if type(each_item) == list:
                    tmp_boolean_list = []
                    for each_gene in each_item:
                        if each_gene in gene_list:
                            tmp_boolean_list.append(2.0)
                        else:
                            tmp_boolean_list.append(1.0)
                    value=1
                    for j in xrange(len(tmp_boolean_list)):
                        value = value * tmp_boolean_list[j]
                    boolean_list.append(value)
                else:
                    if each_item in gene_list:
                        boolean_list.append(2.0)
                    else:
                        boolean_list.append(1.0)
            value=0
            for i in xrange(len(boolean_list)):
                if boolean_list[i] == 2.0:
                    reaction_list.append(reaction)
        reaction_list = list(set(reaction_list))
        return reaction_list

    def get_reactions_from_multiple_gene_deletion(self, gene_list, GPR_calc=True):
        reaction_list = []
        if GPR_calc == True:
            reaction_candidates = []
            for each_gene in gene_list:
                gene_name = str(each_gene)
                tmp_reaction_candidates = self.model_gene_info[gene_name]
                reaction_candidates += tmp_reaction_candidates

            reaction_candidates = list(set(reaction_candidates))

            for reaction in reaction_candidates:
                boolean_list = []
                GPR_list = self.model_reaction_info[reaction]['GPR_list']
                for each_item in GPR_list:
                    if type(each_item) == list:
                        tmp_boolean_list = []
                        for each_gene in each_item:
                            if each_gene in gene_list:
                                tmp_boolean_list.append(0.0)
                            else:
                                tmp_boolean_list.append(1.0)
                        value = 1
                        for j in range(len(tmp_boolean_list)):
                            value = value * tmp_boolean_list[j]
                        boolean_list.append(value)
                    else:
                        if each_item in gene_list:
                            boolean_list.append(0.0)
                        else:
                            boolean_list.append(1.0)
                value = 0
                for i in range(len(boolean_list)):
                    value = value + boolean_list[i]
                if value == 0:
                    reaction_list.append(reaction)
        else:
            temp_reaction_list = []
            for each_gene in gene_list:
                idx = self.cobra_model.genes.index(each_gene)
                for each_reaction in self.cobra_model.genes[idx].get_reaction():
                    temp_reaction_list.append(each_reaction.id)

            reaction_list = list(set(temp_reaction_list))

        return reaction_list


    def get_metabolic_models_from_active_reactions(self, FluxInfo):
        cobra_model = copy.deepcopy(self.cobra_model)
        remove_reaction_candidate_list = []

        for each_reaction in cobra_model.reactions:
            if each_reaction.id in FluxInfo:
                if abs(FluxInfo[each_reaction.id]) <= 10.0e-6:
                    remove_reaction_candidate_list.append(each_reaction)
        cobra_model.remove_reactions(remove_reaction_candidate_list)
        cobra_model.optimize()
        print 'Simulated Growth : ', cobra_model.solution.f
        return cobra_model

    def get_cobra_model(self):
        return copy.deepcopy(self.cobra_model)


if __name__ == '__main__':
    obj = Simulator()
    obj.read_model('./Ecoli_MFA.xml')
