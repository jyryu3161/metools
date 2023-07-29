#Copyright 2014-2016 BioInformatics Research Center, KAIST

import glob
import logging
import os
import sys
import time
import warnings
import itertools
import numpy as np
import pandas as pd
import tqdm

from cobra.io import read_sbml_model, write_sbml_model
from me_targeting.flux_analysis import FSEOF
from me_targeting.flux_analysis import FVSEOF
from me_targeting.flux_analysis import Simulator

def run_MOMA_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, minimum_production_rate, constraints, targeting_mode='reaction', target_num=1):

    candidate_reactions = []    
    for each_reaction in cobra_model.reactions:
        metabolites = each_reaction.reactants + each_reaction.products
        compartments = [each_metabolite.compartment for each_metabolite in metabolites]
        compartments = list(set(compartments))
        if len(compartments) == 1:
            if len(each_reaction.genes) > 0:
                candidate_reactions.append(each_reaction.id)
    
    combinatorial_reaction_sets = []
    for each_reaction_set in itertools.combinations(candidate_reactions, target_num):        
        combinatorial_reaction_sets.append(each_reaction_set)        
    
    fp = open(output_dir+'/MOMA_result_target_num_%s.txt'%(target_num), 'w')
    fp.write('reaction(s)\tbiomass flux\ttarget flux\n')
    
    flux_constraints = {}
    obj = Simulator.Simulator()
    obj.load_cobra_model(cobra_model)
    model_status, objective, flux = obj.run_FBA(new_objective = target_reaction, mode='max')

    max_target_production_rate = flux[target_reaction]
    lb_target_production_rate = flux[target_reaction] * minimum_production_rate
    flux_constraints[target_reaction] = [lb_target_production_rate, 1000.0]
    model_status, objective, wild_flux = obj.run_FBA(new_objective = biomass_reaction, flux_constraints = flux_constraints, mode='max', internal_flux_minimization=True)
    
    for each_reaction_set in tqdm.tqdm(combinatorial_reaction_sets): 
        key_string = None
        for each_set in each_reaction_set:  
            key_string = ';'.join(each_reaction_set)
            
        simulation_flux_constraints = {}
        
        for each_reaction in each_reaction_set:
            simulation_flux_constraints[each_reaction] = [0.0, 0.0]
        
        model_status, objective, perturbed_flux = obj.run_MOMA(wild_flux=wild_flux, flux_constraints=simulation_flux_constraints)
        if objective == False:
            continue
            
        biomass_flux = perturbed_flux[biomass_reaction]
        target_flux = perturbed_flux[target_reaction]   
        # print('%s\t%s\t%s'%(key_string, biomass_flux, target_flux))  
        
        fp.write('%s\t%s\t%s\n'%(key_string, biomass_flux, target_flux))                

    fp.close()
    return

def run_MOMA_targeting_gene(output_dir, cobra_model, biomass_reaction, target_reaction, minimum_production_rate, constraints, targeting_mode='gene', target_num=1):

    candidate_genes = []    
    for each_gene in cobra_model.genes:
        candidate_genes.append(each_gene.id)
    
    combinatorial_gene_sets = []
    for each_set in itertools.combinations(candidate_genes, target_num):        
        combinatorial_gene_sets.append(each_set)        
    
    fp = open(output_dir+'/MOMA_result_target_num_%s.txt'%(target_num), 'w')
    fp.write('reaction(s)\tbiomass flux\ttarget flux\n')
    
    flux_constraints = {}
    obj = Simulator.Simulator()
    obj.load_cobra_model(cobra_model)
    model_status, objective, flux = obj.run_FBA(new_objective = target_reaction, mode='max')

    max_target_production_rate = flux[target_reaction]
    lb_target_production_rate = flux[target_reaction] * minimum_production_rate
    flux_constraints[target_reaction] = [lb_target_production_rate, 1000.0]
    model_status, objective, wild_flux = obj.run_FBA(new_objective = biomass_reaction, flux_constraints = flux_constraints, mode='max', internal_flux_minimization=True)
    for each_gene_set in tqdm.tqdm(combinatorial_gene_sets): 
        key_string = None
        key_list = []
        each_reaction_set = []
        for each_gene in each_gene_set:
            key_list.append(each_gene)
            cobra_gene = cobra_model.genes.get_by_id(each_gene)
            for each_reaction in cobra_gene.reactions:
                each_reaction_set.append(each_reaction.id)
                
        if len(each_reaction_set) == 0:
            continue
            
        each_reaction_set = list(set(each_reaction_set))
        key_string = ';'.join(key_list)
        key_reaction_string = ';'.join(each_reaction_set)

        simulation_flux_constraints = {}
        
        for each_reaction in each_reaction_set:
            simulation_flux_constraints[each_reaction] = [0.0, 0.0]
        
        model_status, objective, perturbed_flux = obj.run_MOMA(wild_flux=wild_flux, flux_constraints=simulation_flux_constraints)
        
        if objective == False:
            continue
            
        biomass_flux = perturbed_flux[biomass_reaction]
        target_flux = perturbed_flux[target_reaction]   
        fp.write('%s(%s)\t%s\t%s\n'%(key_string, key_reaction_string, biomass_flux, target_flux))                

    fp.close()
    return

def run_FSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction):    
    obj = FSEOF.FSEOF()
    obj.load_cobra_model(cobra_model)
    df = obj.run_FSEOF(biomass_reaction, target_reaction)
    df.to_csv(output_dir+'/fseof.csv')
    result_df = obj.result_summary()
    result_df.to_csv(output_dir+'/fseof_summary_result.csv')
    return
    
def run_FVSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction): 
    obj = FVSEOF.FVSEOF()
    obj.load_cobra_model(cobra_model)
    ResultDF = obj.run_FVSEOF(biomass_reaction, target_reaction, {})
    ResultDF.to_csv(output_dir+'/fvseof_result_df.csv')   
    df = obj.result_summary()
    df.to_csv(output_dir+'/fvseof_result_summary.csv')
    return

def run_heterogenous_reaction_addition_simulation(cobra_model):
    
    return

def main():
    start = time.time()    
    warnings.filterwarnings("ignore")
    
    model_file = './raw_data/iML1515_octa_exporter.xml' # iML1515.xml' # model file
    output_dir = './output_octa' # output directory
    mode = 'gene' # or reaction
    biomass_reaction = 'BIOMASS_Ec_iML1515_core_75p37M' #
    target_reaction = 'EX_octa_e' # target reaction
    
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    cobra_model = read_sbml_model(model_file)
    
    constraints = {}
    if mode == 'gene':
        run_MOMA_targeting_gene(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'gene', 1)            
        run_MOMA_targeting_gene(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'gene', 2)
    else:
        run_MOMA_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'reaction', 1)            
        run_MOMA_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'reaction', 2)
    run_FSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction)
    run_FVSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction)
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

if __name__ == '__main__':
    main()
