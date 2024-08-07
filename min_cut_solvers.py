import numpy as np

import networkx as nx

from qiskit_optimization.applications import Maxcut
# from qiskit.algorithms import NumPyMinimumEigensolver
# from qiskit_optimization.algorithms import MinimumEigenOptimizer


#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings('ignore')

# Qiskit 
# import qiskit
# from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization import QuadraticProgram
# qiskit.utils.algorithm_globals.massive=True

# Dwaves
import dimod
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave.system import LeapHybridSampler
from dwave.samplers import SimulatedAnnealingSampler
import dwave.inspector

import numpy as np
import pandas as pd
from sympy import *
import re
import time
import os
import gurobipy as gp
import utils
import json
from collections import Counter
from gurobipy import GRB

# from qiskit.algorithms.optimizers import COBYLA

################### ########
#Different distributions data generator functions


def create_dir(path, log=False):
    if not os.path.exists(path):
        if log:
            print('The directory', path, 'does not exist and will be created')
        os.makedirs(path)
    else:
        if log:
            print('The directory', path, ' already exists')




#################################### SOLVER

def numpy_for_qubo(qubo, p=None):                      # Classical solver for QUBO
    """
    numpy_for_qubo solves the given QUBO using Numpy library functions
    :param
    qubo: CSG problem instance reduced to the form of qubo.

    return:
    result: An array of binary digits which denotes the solutionn of the input qubo problem
    """
    exact_mes = NumPyMinimumEigensolver()
    exact = MinimumEigenOptimizer(exact_mes)  # using the exact classical numpy minimum eigen solver
    result = exact.solve(qubo)
    return result




def solve_QUBO(linear, quadratic, algo, p=1):
    """
    solve_QUBO is a higher order function to solve QUBO using the given algo parameter function
    :param
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.
    algo: a callback function for qaoa_for_qubo or numpy_for_qubo

    return:
    result: An array of binary digits which denotes the solutionn of the input qubo problem


    """
    keys = list(linear.keys())
    keys.sort(key=natural_keys)
    qubo = QuadraticProgram()
    for key in keys:
        qubo.binary_var(key)                                                             # initialize the binary variables
    qubo.minimize(linear=linear, quadratic=quadratic)                                  # initialize the QUBO maximization problem instance
    op, offset = qubo.to_ising()
    qp=QuadraticProgram()
    qp.from_ising(op, offset, linear=True)
    result = algo(qubo, p)
    return result

###############################################

def reduce_quadratic(in_dict):
    
    quadratic = in_dict.copy()
    for k,v in in_dict.items():
        # the dict needs to be non-empty, so the last entry can stay
        if v == 0 and len(quadratic) > 1:
            del quadratic[k]
    return quadratic


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.ht
    For example: Built-in function ['x_8','x_10','x_1'].sort() will sort as ['x_1', 'x_10', 'x_8']
    But using natural_keys as callback function for sort() will sort as ['x_1','x_8','x_10']
    param:
    text: a list of strings ending with numerical characters

    return:
    sorted list in a human way
    """
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def atoi(text):
    """
    Function returns the corresponding value of a numerical string as integer datatype
    param:
    text: string conntaining only numerical charcaters

    return:
    integer value corresponding to the input text
    """
    return int(text) if text.isdigit() else tex                
                 
                 

def exact_solver(linear, quadratic, offset = 0.0):
    """
    Solve Ising hamiltonian or qubo problem instance using dimod API.
    dimod is a shared API for samplers.It provides:
    - classes for quadratic models—such as the binary quadratic model (BQM) class that contains Ising and QUBO models used by samplers such as the D-Wave system—and higher-order (non-quadratic) models.
    - reference examples of samplers and composed samplers.
    - abstract base classes for constructing new samplers and composed samplers.

    :params
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.
    offset: Constant energy offset associated with the Binary Quadratic Model.

    :return
    sample_set: Samples and any other data returned by dimod samplers.
    """
    vartype = dimod.BINARY

    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype)
    sampler = dimod.ExactSolver()
    sample_set = sampler.sample(bqm)
    return sample_set



def dwave_solver(linear, quadratic, offset = 0.0, runs=1000, chain_str = 1, anneal_time=1, **kwargs):
    """
    Solve Ising hamiltonian or qubo problem instance using dimod API for using dwave system.

    :params
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.
    runs: Number of repeated executions

    :return
    sample_set: Samples and any other data returned by dimod samplers.
    """
      
    vartype = dimod.BINARY

    strength = abs(chain_str * max(max(linear.values(),key=abs), max(quadratic.values(),key=abs)))
    
    quadratic = reduce_quadratic(quadratic)
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype)
    dw_2Q =DWaveSampler(token = "jul-be2e8cd3e7ad9063f86b6f11a0f6e025d01e485c", solver={'topology__type': 'chimera'})
     
    def_time = dw_2Q.properties["default_annealing_time"]
    
    
    
    sampler = EmbeddingComposite(DWaveSampler(token = "jul-be2e8cd3e7ad9063f86b6f11a0f6e025d01e485c", solver={'topology__type': 'chimera'}))
    sample_set = sampler.sample(bqm, num_reads=runs, annealing_time=anneal_time*def_time, chain_strength=strength, label=kwargs['label']+"_Q2k")
    
    return sample_set

def dwave_advantage_solver(linear, quadratic, offset = 0.0, runs=100, chain_str = 0.3, anneal_time=1, **kwargs):
    """
    Solve Ising hamiltonian or qubo problem instance using dimod API for using dwave system.

    :params
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.
    runs: Number of repeated executions

    :return
    sample_set: Samples and any other data returned by dimod samplers.
    """
    null_time = time.time()
    vartype = dimod.BINARY
       
    quadratic = reduce_quadratic(quadratic)
    strength = abs(chain_str* max(max(linear.values(),key=abs), max(quadratic.values(),key=abs))) 
    first_time = time.time()
    
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype)
    second_time = time.time()
    dw_ad = DWaveSampler(token = "jul-5215c8008132f206eba40344f90b0cfe8a586600", solver={'topology__type': 'pegasus'})
    third_time = time.time()
    def_time = dw_ad.properties["default_annealing_time"]   
    sampler = EmbeddingComposite(dw_ad)
    
    ct_start = time.time()
    sample_set = sampler.sample(bqm, num_reads=runs, chain_strength=strength, annealing_time=anneal_time*def_time)
    call_time = time.time()
    
    if 'ext_label' in kwargs:
        label = kwargs['ext_label']
        t_id = label.split('_')[-1]
        sparsity = label.split('_')[-2]
        dist = label.split('_')[-3]
        agents = label.split('_')[-4]
        seed = label.split('_')[-5]
        
        path = f"runtime/{seed}/{agents}/{dist}/{sparsity}/{t_id}/"
        utils.create_dir(path)
        try:
            with open(path+'log.json', 'r') as f:
                data = json.load(f)
        except:
            data = {'qpu_sampling_time': 0,
                    'qpu_anneal_time_per_sample': 0,
                    'qpu_readout_time_per_sample': 0,
                    'qpu_access_time': 0,
                    'qpu_access_overhead_time': 0,
                    'qpu_programming_time': 0,
                    'qpu_delay_time_per_sample': 0,
                    'post_processing_overhead_time': 0,
                    'total_post_processing_time': 0,
                    'call_time' :0}
            
        time_dict = sample_set.info['timing']
        
        time_dict['call_time'] = call_time - ct_start
        
        
        
        dic = dict(Counter(data)+Counter(time_dict))     
        with open(path+'log.json', 'w') as f:
            json.dump(dic, f) 
      
        end_time = time.time()
#         print("\nTiming:")
#         print(f"reduce,strengh,vartype: {first_time-null_time}")
#         print(f"bqm: {second_time-first_time}")
#         print(f"sampler:{third_time-second_time}")
#         print(f"embedding def: {ct_start-third_time}")
#         print(f"Call time: {call_time-ct_start}")
#         print(f"dict, os: {end_time-call_time}")
    return sample_set

########################################
def dwave_hybrid_solver(linear, quadratic, offset = 0.0, runs=1000, **kwargs):
    """
    Solve Ising hamiltonian or qubo problem instance using dimod API for using dwave system.

    :params
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.
    runs: Number of repeated executions

    :return
    sample_set: Samples and any other data returned by dimod samplers.
    """
    
    vartype = dimod.BINARY
    quadratic = reduce_quadratic(quadratic)
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, vartype)
    dw_hyb = LeapHybridSampler(token = "jul-be2e8cd3e7ad9063f86b6f11a0f6e025d01e485c", solver={'category': 'hybrid'})
    
    try:
        sample_set = dw_hyb.sample(bqm, time_limit=kwargs['tl'], label=kwargs['label']+"_hyb")
    except:
        sample_set = dw_hyb.sample(bqm, label="_hyb")
    return sample_set

######################################



def extract_best_result(df):
    """
    A function to fetch the binary string with least energy of the input hamiltonian operator

    :params
    df: 

    :return
    x: an array of binary digits specifies the solution of the input problem instance
    fval: value of the operator corresponding to the binary digits in x
    """
    
    row_min = df[df.energy == df.energy.min()]

    cols = []
    for col in df.columns:
        if 'x_' in col:
            cols.append(col)

    x = []

    for col in cols:
        x.append(row_min.iloc[0][col])

    fval = row_min.energy.iloc[0]
    return x, fval

def from_bin_to_var(x, dictionary):
    """
    function to convert binary string to coalition structure

    :params
    x: an array of binary digits (specifies the solution of the input problem instance)
    dictionary: dictionary with coalitions as keys and coalition values as values (CSG problem instance)
    
    :return
    solution: list of lists. coalition structure.
    """
    solution = []
    for i in range(len(x)):
        if x[i] == 1:
            print(list(dictionary.keys())[i])
            solution.append(list(dictionary.keys())[i])
    return solution


def create_QUBO(linear_dict, quadratic_dict):    
    """
    create a QUBO problem instance using the linear and quadratic coefficients.

    :params
    linear: dictionary of linear coefficient terms in the QUBO formulation of the CSG problem.
    quadratic: dictionary of quadratic coefficient terms in the QUBO formulation of the CSG problem.

    :return
    Object of QuadraticProgram class corresponding to the input linear and quadratic coefficients.
    """
    qubo = QuadraticProgram()
    
    keys = list(linear_dict.keys())
    keys.sort(key=natural_keys)
    
    for key in keys:
        qubo.binary_var(key) 


    qubo.minimize(linear=linear_dict, quadratic=quadratic_dict)
    return qubo


def from_columns_to_string(df):

    cols = []
    for col in df.columns:
        if 'x_' in col:
            cols.append(col)

    df['x'] = 'x'
    for index, row in df.iterrows():
        x = ''
        for col in cols:
            x = x + str(row[col])
        df.loc[index, 'x'] = x
    return df[['x', 'num_occurrences', 'energy']]


def get_ordered_solution(dictionary):
    """
    Reordering the (key,value) pairs in the dictionary to fetch only the values in order.
    param:
    dictionary: input dictionary.

    return:
    solution: list of values after reordering the dictionary elements.
    """
    sortedDict = dict(sorted(dictionary.items(), key=lambda x: x[0].lower()))
    solution = []
    for k, v in sortedDict.items():
        solution.append(v)

    return solution



def results_from_dwave(sample_set, exact=False):
    """
    Fetch the details of the output_ from D-Wave system (Quantum Annealing).

    :params
    sample_set: Samples and any other data returned by dimod samplers.

    :return
    solution: a list of binary values corresponding to the solution provided by the output_ of D-Wave device.
    fval: The function value (operator value) of the input hamiltonian corresponding to the solution.
    prob: Probability of the solution.
    rank: rank of the solution out of all possible binary arrays.
    time: time taken by the D-Wave device to compute the solution.
    """
    df = sample_set.to_pandas_dataframe()
    row_min = df[df.energy == df.energy.min()]

    cols = []
    for col in df.columns:
        if 'x_' in col:
            cols.append(col)

    cols.sort(key=natural_keys)
    
    solution = []

    for col in cols:
        solution.append(row_min.iloc[0][col])
    
    
    fval = row_min.energy.iloc[0]

    if not exact:
        occ_min_fval = row_min.num_occurrences.to_list()[0]
        occurences = df.num_occurrences.to_list()

        occurences = sorted(occurences, reverse=True)
        
        time = sample_set.info['timing']['qpu_sampling_time']/1000

        rank = occurences.index(occ_min_fval)+1
        prob = occ_min_fval / sum(df.num_occurrences)
    else:
        rank = 1
        prob = 1
        time = 1

    return solution, fval, prob, rank, time



def min_cut_brute_force(n_agents, induced_subgraph_game, **kwargs):
    #print("Received n_agents, induced_subgraph_game",n_agents, induced_subgraph_game)
    G=nx.Graph()
    G.add_nodes_from(np.arange(0,n_agents,1))
    elist = [tuple((int(x)-1 for x in key.split(',')))+tuple([induced_subgraph_game[key]*-1]) for key in induced_subgraph_game]
    G.add_weighted_edges_from(elist)
    w = [[G.get_edge_data(i,j,default = {'weight': 0})['weight']  for j in range(n_agents)] for i in range(n_agents)]

    x = [0] * n_agents
    cost = [[w[i][j]*x[i]*(1-x[j]) for j in range(n_agents)] for i in range(n_agents)]
    best_cost_brute = sum(induced_subgraph_game.values())
    xbest_cut_brute = x
  
    for b in range(1, 2**(n_agents-1)):

        x = [int(t) for t in reversed(list(bin(b)[2:].zfill(n_agents)))]
        cost = [[w[i][j]*x[i]*(1-x[j]) for j in range(n_agents)] for i in range(n_agents)]
        cost = sum([sum(i) for i in cost])+sum(induced_subgraph_game.values())
        if cost > best_cost_brute:
            best_cost_brute = cost
            xbest_cut_brute = x
    return np.array(xbest_cut_brute), (abs(best_cost_brute) + best_cost_brute) / 2




def min_cut_qiskit_classical_eigensolver(n_agents, induced_subgraph_game, **kwargs):
    G=nx.Graph()
    G.add_nodes_from(np.arange(0,n_agents,1))
    elist = [tuple((int(x)-1 for x in key.split(',')))+tuple([induced_subgraph_game[key]*-1]) for key in induced_subgraph_game]
    G.add_weighted_edges_from(elist)
    w = [[G.get_edge_data(i,j,default = {'weight': 0})['weight']  for j in range(n_agents)] for i in range(n_agents)]
    w = np.array([np.array(row) for row in w])
    max_cut = Maxcut(w)
    qp = max_cut.to_quadratic_program()
    #qubitOp, offset = qp.to_ising()
    exact = MinimumEigenOptimizer(NumPyMinimumEigensolver())
    result = exact.solve(qp)
    return result.x, (abs(result.fval)+result.fval)/2



def gurobi_solver(n, induced_subgraph_game, **kwargs):
    linear, quadratic = get_linear_quadratic_coeffs(n, induced_subgraph_game)
    
#     if 'ext_label' in kwargs:
#         label = kwargs['ext_label']
#         t_id = label.split('_')[-1]
#         sparsity = label.split('_')[-2]
#         dist = label.split('_')[-3]
#         agents = label.split('_')[-4]
#         seed = label.split('_')[-5]
        
#         path = f"runtime/{seed}/{agents}/{dist}/{sparsity}/{t_id}/"
#         utils.create_dir(path)
#         try:
#             with open(path+'log_g.json', 'r') as f:
#                 gr_data = json.load(f)
#         except:
#             gr_data = {'runs' : 0}
            
    qubo_matrix = np.zeros([n,n])
    linear_enum = {}
    for i, (key,value) in enumerate(linear.items()):
        linear_enum[key] = i
        qubo_matrix[i,i] = value
     
   
    for key,value in quadratic.items():
        qubo_matrix[linear_enum[key[0]], linear_enum[key[1]]] = value
        qubo_matrix[linear_enum[key[1]], linear_enum[key[0]]] = value
    model = gp.Model()
   
    #initialize binary variables
    x = model.addVars(n, vtype=GRB.BINARY)
    #objective function to minimize the QUBO matrix
    obj_expr = gp.quicksum(qubo_matrix[i, j] * x[i] * x[j] for i in range(n) for j in range(n))
    model.setObjective(obj_expr)
    model.setParam('OutputFlag', 0)
    model.optimize()
    
    
#     gr_data['runs'] += 1     
#     with open(path+'log_g.json', 'w') as f:
#         json.dump(gr_data, f)
        
    if model.status == GRB.OPTIMAL:
        solution = [int(x[i].X) for i in range(n)]
        binary_string = ''.join(str(bit) for bit in solution)
        int_string = [int(k) for k in list(binary_string)]
        return int_string, None
    else:
        return None, None



def min_cut_dwave_annealer(n_agents, induced_subgraph_game, save_log=False, name_folder='distribution', n_samples= 1000, n_run=1, **kwargs):
    
    linear, quadratic = get_linear_quadratic_coeffs(n_agents, induced_subgraph_game)
    
    try:
        sample = dwave_solver(linear, quadratic, offset = kwargs['offset'], runs=n_samples, chain_str=kwargs['chain_str'], label=kwargs['label'])
    except:
        sample = dwave_solver(linear, quadratic, offset = kwargs['offset'], runs=n_samples, label=kwargs['label'])
    
    
    
    
    # if save_log:
    #     path = os.path.join('QA_results', name_folder, str(n), 'run_'+str(run))
    #     create_dir(path)
    #     try:
    #         sample.to_pandas_dataframe().to_csv(os.path.join(path, 'solutions.csv'))
    #         save_json(os.path.join(path, 'log'), sample.info)
    #     except:
    #         print("\n *** Warning: results for",  name_folder, "with", n_agents, "agents not saved** \n")
    dwave_annealer_solution=[]
    
    d = sample.first[0]
    p = {}
    
    for k,v in d.items():
        p[int(k.split('_')[-1])]=v
    
    d = dict(sorted(p.items()))
    for key, value in d.items():
        dwave_annealer_solution.append(value)
    dwave_annealer_solution = np.array(dwave_annealer_solution)
    dwave_annealer_value = from_columns_to_string(sample.to_pandas_dataframe()).loc[0,'energy']
    # print("s: ", n_agents, " - time: ", sample.info['timing']['qpu_sampling_time']/10**6)
    #dwave_annealer_tte = sample.info['timing']['qpu_sampling_time']/10**6
    return dwave_annealer_solution, dwave_annealer_value


def min_cut_dwave_advantage(n_agents, induced_subgraph_game, save_log=False, name_folder='distribution', n_samples= 100, n_run=1, **kwargs):
    
    linear, quadratic = get_linear_quadratic_coeffs(n_agents, induced_subgraph_game)
    
    

    sample = dwave_advantage_solver(linear, quadratic, runs=n_samples, **kwargs)
    # if save_log:
    #     path = os.path.join('QA_results', name_folder, str(n), 'run_'+str(run))
    #     create_dir(path)
    #     try:
    #         sample.to_pandas_dataframe().to_csv(os.path.join(path, 'solutions.csv'))
    #         save_json(os.path.join(path, 'log'), sample.info)
    #     except:
    #         print("\n *** Warning: results for",  name_folder, "with", n_agents, "agents not saved** \n")
    dwave_annealer_solution=[]
    d = sample.first[0]
    p = {}
    
    for k,v in d.items():
        p[int(k.split('_')[-1])]=v
    
    d = dict(sorted(p.items()))
    
    for key, value in d.items():
        dwave_annealer_solution.append(value)
    dwave_annealer_solution = np.array(dwave_annealer_solution)
    dwave_annealer_value = from_columns_to_string(sample.to_pandas_dataframe()).loc[0,'energy']
    # print("s: ", n_agents, " - time: ", sample.info['timing']['qpu_sampling_time']/10**6)
#     print("ST:", sample.info['timing']['qpu_sampling_time']/10**6)
    return dwave_annealer_solution, dwave_annealer_value

###########################


def min_cut_dwave_hybrid(n_agents, induced_subgraph_game, save_log=False, name_folder='distribution', n_samples= 2000, n_run=1, **kwargs):
    
    linear, quadratic = get_linear_quadratic_coeffs(n_agents, induced_subgraph_game)
    sample = dwave_hybrid_solver(linear, quadratic, runs=n_samples)
    # if save_log:
    #     path = os.path.join('QA_results', name_folder, str(n), 'run_'+str(run))
    #     create_dir(path)
    #     try:
    #         sample.to_pandas_dataframe().to_csv(os.path.join(path, 'solutions.csv'))
    #         save_json(os.path.join(path, 'log'), sample.info)
    #     except:
    #         print("\n *** Warning: results for",  name_folder, "with", n_agents, "agents not saved** \n")
    dwave_annealer_solution=[]
    d = sample.first[0]
    p = {}
    
    for k,v in d.items():
        p[int(k.split('_')[-1])]=v
    
    d = dict(sorted(p.items()))
    
    
    for key, value in d.items():
        dwave_annealer_solution.append(value)
    dwave_annealer_solution = np.array(dwave_annealer_solution)
    dwave_annealer_value = from_columns_to_string(sample.to_pandas_dataframe()).loc[0,'energy']
    # print("s: ", n_agents, " - time: ", sample.info['timing']['qpu_sampling_time']/10**6)
    #dwave_annealer_tte = sample.info['timing']['qpu_sampling_time']/10**6
    return dwave_annealer_solution, dwave_annealer_value


def min_cut_dwave_exact(n_agents, induced_subgraph_game, save_log=False, name_folder='distribution', n_samples= 1000, n_run=1, **kwargs):
    
    linear, quadratic = get_linear_quadratic_coeffs(n_agents, induced_subgraph_game)
    
    sample = exact_solver(linear, quadratic, offset = 0.0)
    
    
    
    
    # if save_log:
    #     path = os.path.join('QA_results', name_folder, str(n), 'run_'+str(run))
    #     create_dir(path)
    #     try:
    #         sample.to_pandas_dataframe().to_csv(os.path.join(path, 'solutions.csv'))
    #         save_json(os.path.join(path, 'log'), sample.info)
    #     except:
    #         print("\n *** Warning: results for",  name_folder, "with", n_agents, "agents not saved** \n")
    dwave_annealer_solution=[]
    
    d = sample.first[0]
    p = {}
    
    for k,v in d.items():
        p[int(k.split('_')[-1])]=v
    
    d = dict(sorted(p.items()))
    
    for key, value in d.items():
        dwave_annealer_solution.append(value)
    dwave_annealer_solution = np.array(dwave_annealer_solution)
    dwave_annealer_value = from_columns_to_string(sample.to_pandas_dataframe()).loc[0,'energy']
    # print("s: ", n_agents, " - time: ", sample.info['timing']['qpu_sampling_time']/10**6)
    #dwave_annealer_tte = sample.info['timing']['qpu_sampling_time']/10**6
    return dwave_annealer_solution, dwave_annealer_value


############################




def get_linear_quadratic_coeffs(n_agents, induced_subgraph_game):
    
    G=nx.Graph()
    G.add_nodes_from(np.arange(0,n_agents,1))
    elist = [tuple((int(x)-1 for x in key.split(',')))+tuple([induced_subgraph_game[key]*-1]) for key in induced_subgraph_game]
    G.add_weighted_edges_from(elist)
    w = [[G.get_edge_data(i,j,default = {'weight': 0})['weight']  for j in range(n_agents)] for i in range(n_agents)]
    w = np.array([np.array(row) for row in w])
    max_cut = Maxcut(w)
    
    qp = max_cut.to_quadratic_program()
    linear = qp.objective.linear.coefficients.toarray(order=None, out=None)
    quadratic = qp.objective.quadratic.coefficients.toarray(order=None, out=None)
    linear = {'x_'+str(idx):-round(value,2) for idx,value in enumerate(linear[0])}
    quadratic = {('x_'+str(iy),'x_'+str(ix)):-quadratic[iy, ix] for iy, ix in np.ndindex(quadratic.shape) if iy<ix and abs(quadratic[iy, ix])!=0}
    return linear, quadratic