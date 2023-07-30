# TODO: could use revision

import numpy as np
from copy import deepcopy
import random
import time
import math

from ft_learn.ft.fault_tree import FaultTree, str2ft, mcsInit
from ft_learn.ft.fault_tree import save_results
from ft_learn.ft.ft_elements import BE, AND, OR
import ft_learn.helper as helper
import ft_learn.moea.genetic_operators as genetic_operators
import ft_learn.moea.fitness as fitness
import os, os.path

def generate_initial_population(basic_events):
    ft1 = FaultTree(OR(deepcopy(basic_events)))
    ft2 = FaultTree(AND(deepcopy(basic_events)))
    return [ft1, ft2]

def upfold(saving_folder):
    # Save FT results:
    if os.path.isdir(saving_folder) is False:
        run_ = 0
        os.makedirs(saving_folder + '/run_' + str(run_))
    elif os.path.isdir(saving_folder) is True:
        files_in_folder = []
        for name in os.listdir(saving_folder):
            files_in_folder.append(name)
        run_ = int(len(files_in_folder))
        os.makedirs(saving_folder + '/run_' + str(run_))
    # Update:
    saving_folder = saving_folder + '/run_' + str(run_) 
    return saving_folder




def perform_genetic_ftmoea(dataset=[], MCSs=[], bes=[], population_size=100, ft_as_input='', generations=100, convergence_criterion=10, multi_objective_function=[1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           config_gen_op=None, selection_strategy='elitist', saving_results_=0, path_save_results="", debugging=False, seg_size=4):
    """
    Learns a FT consistent with a given dataset, using genetic operations.
    :param dataset: Matrix containing all cut sets.
    :param MCSs: Cut sets.
    :param bes: List of all BEs which occur in the data.
    :param population_size: The population size which is maintained after each iteration of genetic operations.
    :param ft_as_input:
    :param generations: The number of iterations performed before the algorithm stops.
    :param convergence_criterion: The number of iterations without score improvement after which the algorithm stops.
    :param multi_objective_function:
    :param config_gen_op: Configuration for genetic operations.
    :param selection_strategy:
    :param saving_results_:
    :param path_save_results:
    :param debugging: Whether debugging should be enabled. This option ensures reproducible random generations.
    :return: Consistent fault tree.
    """

    print('... FT-MOEA initialized ...')
    
    # Initialize random seed
    if debugging:
        random.seed(42)
    else:
        random.seed()

    if not config_gen_op:
        config_gen_op = genetic_operators.GenOpConfig(1.0)

    if len(MCSs)==0:
        multi_objective_function[0] = 0
    
    if path_save_results != '':
        # Check if the saving folder exist, otherwise create it.
        path_save_results = upfold(path_save_results)

    multi_objective_function = list(-1*np.array(multi_objective_function)) # Minimize the size of the FT
    if MCSs != []:
        ft_from_MCSs = MCSs
        
    t = [time.time()]
    dict_iterations = []

    def index_of(a,list):
        for i in range(0,len(list)):
            if list[i] == a:
                return i
        return -1

    def sort_by_values(list1, values):
        sorted_list = []
        while (len(sorted_list) != len(list1)):
            if index_of(min(values), values) in list1:
                sorted_list.append(index_of(min(values), values))
            values[index_of(min(values), values)] = math.inf
        return sorted_list

    # ------------------------------------
    # Create BEs
    basic_events = [BE(be) for be in bes['all']]

    # Create the initial population:
    print("create initial population")
    ft_as_input1 = 'OR(BE12,OR(BE12,OR(AND(BE7,BE8,BE10),AND(BE5,OR(BE8,BE1),BE2,BE12)),AND(BE5,BE1,BE3), OR(BE12, OR(AND(BE2, BE5, BE4)))), OR(BE12, AND(BE2, BE5, BE4), OR(AND(OR(OR(BE12, AND(BE5, OR(BE2, OR(BE1))), AND(BE6, OR(BE10, BE11))), BE12, AND(BE11, BE9, AND(BE7)), OR(BE12, AND(BE7, OR(BE11, BE10)), AND(BE6, BE10))), OR(AND(BE5, OR(BE2, BE1), BE4), AND(OR(BE12, OR(BE11, BE10), BE5), OR(BE12, AND(OR(BE6, BE12, AND(BE7, OR(BE11, BE10))), AND(OR(BE8, BE9, BE12), AND(OR(BE7, BE6), OR(BE11, BE10)))), OR(BE12, AND(BE2, BE5, BE3)))))), BE12))'
    ft_as_input2 = 'OR(BE8, BE4, AND(BE9, OR(BE7, BE10, BE6, BE5)), OR(BE8, BE4, AND(OR(BE6, BE3, OR(BE10, BE6, BE8, BE1)), AND(OR(BE8, BE9, BE5, OR(BE5, AND(BE7))), OR(BE6, BE2, OR(BE10))))))'
    # print("scoop", mcsInit(ft_as_input2, bes)[0])
    # print("scoop", mcsInit(ft_as_input2, bes)[1])
    # print("scoop", mcsInit(ft_as_input2, bes)[2])
    if ft_as_input == '':
        # Linard's (2019)
        initial_population = generate_initial_population(basic_events)
    else:
        initial_population = [str2ft(ft_as_input)]
        
    # initial_population = mcsInit(ft_as_input2, bes, dataset)
    print("pop", initial_population)
    # ------------------------------------
    while len(initial_population) < population_size:
        initial_population = genetic_operators.apply_genetic_operators(initial_population, basic_events, config_gen_op, debugging)

    t.append(time.time())
    start_t = time.time()
    print("start fitness function")
    raw_fts = fitness.cost_function(initial_population, dataset, bes, population_size, ft_from_MCSs, multi_objective_function, seg_size)

    # sortedPeople,fitnesses,fitness_dict

    if path_save_results != '':
        #Saving dataset
        save_results(raw_fts,t[-1]-t[0],path_save_results,dataset,ft_from_MCSs,multi_objective_function)

    print('Gen. \t Fitness Pop.             \t Fitness best \t                   Best individual')
    print('0','\t    ϕ_c=',"{:.4f}".format(np.mean(raw_fts[1][:,0])),', ϕ_d=',"{:.4f}".format(np.mean(raw_fts[1][:,2])),', ϕ_r=',"{:.4f}".format(np.mean(raw_fts[1][:,3])),', ϕ_im=',"{:.4f}".format(np.mean(raw_fts[1][:,4])),'\t /  ϕ_c=',"{:.4f}".format(raw_fts[1][-1,0]),', ϕ_acc=',"{:.4f}".format(np.mean(raw_fts[1][-1,11])),', ϕ_d=',"{:.4f}".format(raw_fts[1][-1,2]),', ϕ_r=',"{:.4f}".format(raw_fts[1][-1,3]),', ϕ_im=',"{:.4f}".format(raw_fts[1][-1,4]),', ϕ_prec=',"{:.2f}".format(raw_fts[1][-1,5]), ', ϕ_spec=',"{:.4f}".format(raw_fts[1][-1,6]), ', ϕ_sens=',"{:.4f}".format(raw_fts[1][-1,7]),', ϕ_npv=',"{:.4f}".format(raw_fts[1][-1,8]),', ϕ_fnr=',"{:.4f}".format(raw_fts[1][-1,9]),', ϕ_fpr=',"{:.4f}".format(raw_fts[1][-1,10]), ', ϕ_s=',"{:.2f}".format(raw_fts[1][-1,1]),'elapsed_time=',"{:.2f}".format(time.time()-start_t),'\t',raw_fts[0][-1])

    dict_iterations.append([str(raw_fts[0][-1])] + np.mean(raw_fts[1], axis=0).tolist() + raw_fts[1][-1].tolist())
    population = raw_fts[0]
    conv = 0

    for i in range(1, generations):
        start_t = time.time()
        t.append(time.time())
        st = 0

        while st == 0:
            new_population = []
            st = 1
            while len(new_population) < population_size:
                new_population = genetic_operators.apply_genetic_operators(population, basic_events, config_gen_op, debugging)
            raw_fts = fitness.cost_function(new_population, dataset, bes, population_size, ft_from_MCSs, multi_objective_function)

        if path_save_results != '':
            #Saving dataset
            save_results(raw_fts,t[-1]-t[0],path_save_results,dataset,ft_from_MCSs,multi_objective_function)
        
        dict_iterations.append([str(raw_fts[0][-1])] + np.mean(raw_fts[1],axis=0).tolist() + raw_fts[1][-1].tolist()  )
        print(str(i),'\t    ϕ_c=',"{:.4f}".format(np.mean(raw_fts[1][:,0])),', ϕ_d=',"{:.4f}".format(np.mean(raw_fts[1][:,2])),', ϕ_r=',"{:.4f}".format(np.mean(raw_fts[1][:,3])),', ϕ_im=',"{:.4f}".format(np.mean(raw_fts[1][:,4])),'\t /  ϕ_c=',"{:.4f}".format(raw_fts[1][-1,0]),', ϕ_acc=',"{:.4f}".format(np.mean(raw_fts[1][-1,11])),', ϕ_d=',"{:.4f}".format(raw_fts[1][-1,2]),', ϕ_r=',"{:.4f}".format(raw_fts[1][-1,3]),', ϕ_im=',"{:.4f}".format(raw_fts[1][-1,4]), ', ϕ_prec=',"{:.4f}".format(raw_fts[1][-1,5]), ', ϕ_spec=',"{:.4f}".format(raw_fts[1][-1,6]),', ϕ_sens=',"{:.4f}".format(raw_fts[1][-1,7]),', ϕ_npv=',"{:.4f}".format(raw_fts[1][-1,8]),', ϕ_fnr=',"{:.4f}".format(raw_fts[1][-1,9]),', ϕ_fpr=',"{:.4f}".format(raw_fts[1][-1,10]), ', ϕ_s=',"{:.2f}".format(raw_fts[1][-1,1]),'elapsed_time=',"{:.2f}".format(time.time()-start_t),'\t',raw_fts[0][-1])

        # ------------------------------
        # Convergence criteria:

        objectives = [i for i in range(len(multi_objective_function)) if multi_objective_function[i] != 0]
        conv += 1
        for objective_type in objectives:
            if (dict_iterations[-2][objective_type + 13] != dict_iterations[-1][objective_type + 13]):
                conv = 0       
        if conv >= convergence_criterion-1: #or ( dict_iterations[-1][4] == 1.0 and dict_iterations[-1][6] == 1.0 ):
            print('... FT-MOEA finalized ...')
            return raw_fts[0], t, raw_fts[1]
        # if multi_objective_function == [-1,-1,-1]:
            
        #     if dict_iterations[-2][4] == dict_iterations[-1][4] and dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][6] == dict_iterations[-2][6]:
        #         conv += 1
        #     else:
        #         conv = 0

        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1: #or ( dict_iterations[-1][4] == 1.0 and dict_iterations[-1][6] == 1.0 ):
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]
        
        # elif multi_objective_function == [-1,0,-1]:
            
        #     if dict_iterations[-2][4] == dict_iterations[-1][4] and dict_iterations[-1][6] == dict_iterations[-2][6]:
        #         conv += 1
        #     else:
        #         conv = 0

        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1 or ( dict_iterations[-1][4] == 0.0 and dict_iterations[-1][6] == 0.0 ):
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]
            
        # elif multi_objective_function == [0,0,-1]: 
        
        #     if dict_iterations[-1][6] == dict_iterations[-2][6]:
        #         conv += 1
        #     else:
        #         conv = 0
                
        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1 or np.max(raw_fts[1][:,2]) == 0.0:
        #         if np.max(raw_fts[1][:,2]) == 0.0:
        #             print('Max. Accuracy achieved.')
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]
        
        # elif multi_objective_function == [0,-1,0]: 
        
        #     if dict_iterations[-1][5] == dict_iterations[-2][5]:
        #         conv += 1
        #     else:
        #         conv = 0
                
        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1:
        #         if np.max(raw_fts[1][:,2]) == 0.0:
        #             print('Max. Accuracy achieved.')
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]
        
        
        # elif multi_objective_function == [0,-1,-1]: 
        
        #     if dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][6] == dict_iterations[-2][6]:
        #         conv += 1
        #     else:
        #         conv = 0
                
        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1:
        #         if np.max(raw_fts[1][:,2]) == 0.0:
        #             print('Max. Accuracy achieved.')
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]
        
        # elif multi_objective_function == [-1,-1,0]: 

        #     if dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][4] == dict_iterations[-2][4]:
        #         conv += 1
        #     else:
        #         conv = 0
                
        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1:
        #         if np.max(raw_fts[1][:,2]) == 0.0:
        #         #if np.max(trimmed_fts[1][:,2]) == 1.0:
        #             print('Max. Accuracy achieved.')
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]


        # elif multi_objective_function == [-1,0,0]: 
            
        #     if dict_iterations[-1][4] == dict_iterations[-2][4]:
        #         conv += 1
        #     else:
        #         conv = 0
                
        #     #return best FT is convergence or best fitness reached
        #     if conv >= convergence_criterion-1 or np.max(raw_fts[1][:,0]) == 0.0:
        #         if np.max(raw_fts[1][:,0]) == 0.0:
        #             print('Max. Accuracy achieved.')
        #         print('... FT-MOEA finalized ...')
        #         return raw_fts[0], t, raw_fts[1]



        """
            SELECTION STRATEGY
        """

        def roulette_wheel(popul, population_size, fitness_dict):
            chosen = []
            for n in range(population_size):
                r = random.random()
                for ft in popul:
                    if r <= fitness_dict[str(ft)]:
                        if ft not in chosen:
                            chosen.append(ft)
                            break
            return chosen

        def sus(popul, population_size, fitness_dict):

            total_fitness = 0.0
            for fit in fitness_dict:
                total_fitness += fitness_dict[fit]
            point_distance = total_fitness / population_size
            start_point = random.uniform(0, point_distance)
            points = [start_point + i * point_distance for i in range(population_size)]

            parents = []
            while len(parents) < population_size:
                random.shuffle(popul)
                i = 0
                while i < len(points) and len(parents) < population_size:
                    j = 0
                    while j < len(popul):
                        subset_sum = 0.0
                        for k in range(0, j):
                            subset_sum += fitness_dict[str(popul[k])]
                        if subset_sum > points[i]:
                            parents.append(popul[j])
                            break
                        j += 1
                    i += 1

            return list(parents)

        def tournament(popul, population_size, fitness_dict):
            tournament_size = 4
            chosen = []
            while len(chosen) < population_size:
                best = None
                subpop = random.sample(list(popul), tournament_size)
                for f in subpop:
                    if best == None:
                        best = f
                        continue
                    if fitness_dict[str(f)] > fitness_dict[str(best)]:
                        best = f

                chosen.append(best)

            return chosen

        def random_select(popul, population_size, fitness_dict):
            return random.sample(list(popul), population_size)

        def elitism(popul, population_size, fitness_dict):
            return popul[-population_size:]

        # SELECTION
        sortedPeople, fitness_dict = raw_fts[0], raw_fts[2]

        if selection_strategy == 'elitist':
            population = elitism(sortedPeople, population_size, fitness_dict)
        elif selection_strategy == 'roulette_wheel':
            population = roulette_wheel(sortedPeople, population_size, fitness_dict)
        elif selection_strategy == 'sus':
            population = sus(sortedPeople, population_size, fitness_dict)
        elif selection_strategy == 'tournament':
            population = tournament(sortedPeople, population_size, fitness_dict)
        elif selection_strategy == 'random_select':
            population = random_select(sortedPeople, population_size, fitness_dict)
        else:
            population = elitism(sortedPeople, population_size, fitness_dict)

    print('... FT-MOEA finalized ...')
    return raw_fts[0], t, raw_fts[1]
