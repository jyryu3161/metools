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
import matplotlib.pyplot as plt
import copy
import umap
from sklearn.decomposition import PCA
import seaborn as sns

from cobra.io import read_sbml_model, write_sbml_model

from me_targeting.flux_analysis import FSEOF
from me_targeting.flux_analysis import FVSEOF
from me_targeting.flux_analysis import Simulator

from adjustText import adjust_text

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

def target_summary(output_dir, cobra_model):
    moma_result_file = output_dir + '/MOMA_result_target_num_1.txt'
    df = pd.read_csv(moma_result_file, sep='\t')
    
    df = df[df['biomass flux'] > 0.05]
    target_genes = []
    for each_row, each_df in df.iterrows():
        target_str = each_df['reaction(s)']
        gene = target_str.split('(')[0].strip()
        target_genes.append(gene)

    candidate_reactions = []
    for each_gene in cobra_model.genes:
        if each_gene.id not in target_genes:
            continue
        
        if len(each_gene.reactions) == 1:
            for each_rxn in each_gene.reactions:
                candidate_reactions.append(each_rxn.id)
                
    candidate_reactions = list(set(candidate_reactions))
    fseof_result_file = output_dir + '/fseof_summary_result.csv'
    df2 = pd.read_csv(fseof_result_file, sep=',', index_col=0)
    df2 = df2[abs(df2['pearson'])>0.4]
    
    candidate_reactions = list(set(candidate_reactions) & set(df2.index))
    df2.loc[candidate_reactions].to_csv(output_dir + '/fseof_summary_result_filtered.csv')
    
    ###############
#     fvseof_result_file = output_dir + '/fvseof_result_summary.csv'
#     df2 = pd.read_csv(fvseof_result_file, sep=',', index_col=0)
#     df2 = df2[abs(df2['AVG_R'])>0.8]
    
#     candidate_reactions = list(set(candidate_reactions) & set(df2.index))
#     df2.loc[candidate_reactions].to_csv(output_dir + '/fvseof_summary_result_filtered.csv')
    return

def flux_response_analysis(output_dir, model, biomass_reaction, prod_target_reaction):
    
    plt.style.use('ggplot')

    viz_output_dir = output_dir +'/flux_response_analysis_results/'
    try:
        os.mkdir(viz_output_dir)
    except:
        pass
    
    df1 = pd.read_csv(output_dir + '/fseof_summary_result_filtered.csv', index_col=0)
    
    obj = Simulator.Simulator()
    obj.load_cobra_model(model)
    
    for status in ['DOWN', 'UP']:
        df = df1[df1['status']==status]
        target_reactions = df.index
        
        simulation_results = {}
        for reaction_id in tqdm.tqdm(target_reactions):
            simulation_results[reaction_id] = {}
            target_reaction = model.reactions.get_by_id(reaction_id)
            production_reaction = model.reactions.get_by_id(prod_target_reaction)

            # Target 반응의 최소 및 최대 플럭스를 계산합니다.
            a,b,flux1 = obj.run_FBA(new_objective=target_reaction.id, mode='min')
            a,b,flux2 = obj.run_FBA(new_objective=target_reaction.id, mode='max')

            min_flux = flux1[reaction_id]
            max_flux = flux2[reaction_id]

            target_fluxes = np.linspace(min_flux, max_flux, 10)
            production_fluxes = []

            # 각 Target 플럭스 값에 대해 생산 반응의 플럭스를 계산합니다.
            i = 1
            tmp_results = {}
            for flux in target_fluxes:
                constraints={}
                constraints[reaction_id] = [flux, flux]
                a,b,flux = obj.run_FBA(new_objective=production_reaction.id, flux_constraints=constraints, mode='max')
                production_fluxes.append(flux[production_reaction.id])
                tmp_results[i] = flux[production_reaction.id]
                i+=1
            simulation_results[reaction_id] = tmp_results

            # 결과를 시각화합니다.
            plt.figure(figsize=(10, 6))
            plt.plot(target_fluxes, production_fluxes, label=production_reaction.id)
            plt.xlabel('%s (mmol/gDCW/h)'%(target_reaction.id))
            plt.ylabel('%s (mmol/gDCW/h)'%('Target product'))

            plt.legend()
            plt.savefig(viz_output_dir + '/%s_%s.png'%(status, reaction_id), format='png')
            plt.show()

        result_df = pd.DataFrame.from_dict(simulation_results)
        result_df.to_csv(output_dir+'/viz_summary_%s.csv'%(status))
        
    for status in ['DOWN', 'UP']:
        df = pd.read_csv(output_dir + '/viz_summary_%s.csv' % (status), index_col=0)
        df = df.T

        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df)

        plt.figure(figsize=(10, 7))
        plt.scatter(pca_result[:, 0], pca_result[:, 1], color='skyblue', label='Gene', alpha=0.8)

        texts = []
        for i, label in enumerate(df.index):
            texts.append(plt.text(pca_result[i, 0], pca_result[i, 1], label, ha='center', va='center'))

        adjust_text(texts)

        plt.title('PCA Visualization')
        plt.xlabel('PC_1')
        plt.ylabel('PC_2')
        plt.legend()
        plt.savefig(output_dir + '/PCA_Visualization_%s.png' % (status), format='png')
        plt.show()

        # 클러스터링 히트맵 그리기
        plt.figure(figsize=(10, 7))
        sns.clustermap(df, method='average', cmap='vlag', col_cluster=False, figsize=(10, 7))
        plt.savefig(output_dir + '/Clustering_%s.png' % (status), format='png')
        plt.show()
    return

def sampling_analysis(output_dir, model, biomass_reaction, prod_target_reaction):
    plt.style.use('ggplot')
    viz_output_dir = output_dir +'/flux_sampling_results/'
    try:
        os.mkdir(viz_output_dir)
    except:
        pass

    return

def main():
    # 0.5.11
    
    start = time.time()    
    warnings.filterwarnings("ignore")
    
    model_file = './raw_data/iML1515_octa_exporter_fade_del.xml' # iML1515.xml' # model file
    output_dir = './output_octa_fade' # output directory
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
        # run_MOMA_targeting_gene(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'gene', 2)
    else:
        run_MOMA_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'reaction', 1)            
        # run_MOMA_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, 0.1, constraints, 'reaction', 2)
    run_FSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction)
    run_FVSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction)
    
#     # summary
    target_summary(output_dir, cobra_model)
    flux_response_analysis(output_dir, cobra_model, biomass_reaction, target_reaction)
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

if __name__ == '__main__':
    main()
