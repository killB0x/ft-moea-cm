# TODO: could use revision

import numpy as np
import itertools

import ft_learn.helper as helper
from ft_learn.ft.ft_elements import BE, AND, OR


def get_gate_statistics(ft):  # Marijn
    result = []
    gates = ft.get_all_gates()
    i_AND = 0
    i_OR = 0

    for gate in gates:

        if isinstance(gate, AND):
            gate_name = 'AND_' + str(i_AND)
            i_AND += 1
        elif isinstance(gate, OR):
            gate_name = 'OR_' + str(i_OR)
            i_OR += 1
        else:
            gate_name = ""
            assert False

        no_be = 0
        for child in gate.children:
            if isinstance(child, BE):
                no_be += 1

        result.append([gate_name, no_be])

    return result


def fitness(ft, dataset):
    count_true = 0
    for data in dataset:
        if data['T'] == ft.evaluate(data):
            count_true += 1
    return count_true / len(dataset)


def fitness_count(ft, dataset):
    count_true = 0
    total_counts = 0
    for data in dataset:
        if data['T'] == ft.evaluate(data):
            count_true += data['N']
        total_counts += data['N']
    return count_true / total_counts


# Fitness function:
def fitness_cutsets(ft1, cs_base, bes):
    # We obtain the cut sets of both trees:
    c1 = helper.cutsets_from_ft(ft1, bes)
    # c2 = helper.cutsets_from_ft(ft2)
    c2 = cs_base

    if len(c1.shape) == 1:
        c1 = c1.reshape((1, c1.shape[0]))

    if len(c2.shape) == 1:
        c2 = c2.reshape((1, c2.shape[0]))

    # Make cut set matrices the same size:
    if c1.shape[1] > c2.shape[1]:
        t1 = np.zeros((c2.shape[0], c1.shape[1]), dtype=int)
        for i in range(0, c2.shape[1]):
            if sum(c2[:, i]) > 0:
                t1[:, i] = c2[:, i]
        c2 = t1
    elif c1.shape[1] < c2.shape[1]:
        t1 = np.zeros((c1.shape[0], c2.shape[1]), dtype=int)
        for i in range(0, c1.shape[1]):
            if sum(c1[:, i]) > 0:
                t1[:, i] = c1[:, i]
        c1 = t1

    # RV-Correlation:
    c_cos_coef = np.trace(np.matmul(np.matmul(c1.T, c1).T, np.matmul(c2.T, c2))) / np.sqrt(
        np.trace(np.matmul(np.matmul(c1.T, c1).T, np.matmul(c1.T, c1))) * np.trace(np.matmul(np.matmul(c2.T, c2).T, np.matmul(c2.T, c2))))

    return c_cos_coef


#%% Compute metrics
def compute_metrics_fts(initial_population,dataset,ft_from_MCSs,multi_objective_function,bes, seg_size=4):
    
    #show us initial population
    fitnesses = []
    fitness_dict = {}

    #Compute accuracy based on the MCS:
    if multi_objective_function[0] != 0:
        acc_mcs = []
        for ft in initial_population:
            acc_mcs.append([ft.phi_c(ft_from_MCSs,bes['all'])])
        acc_mcs = np.array(acc_mcs)
    else:
        acc_mcs = np.ones((len(initial_population),2))*-1
        
    #Compute the size:
    if multi_objective_function[1] != 0 :
        tree_size = []
        for ft in initial_population:
            #tree_size.append([len(ft.get_all_bes()) + len(get_gate_statistics(ft))])
            tree_size.append([ft.phi_s()])
        tree_size = np.array(tree_size)
    else:
        tree_size = np.ones((len(initial_population),2))*-1
    
    #Compute accuracy based on the dataset:
    if multi_objective_function[2] != 0:
        acc_data = []
        for ft in initial_population:
            #start_point_time = time.time()
            acc_data.append([ft.phi_d(dataset)])
            #print("Calculated confusion matrix for ft in", time.time() - start_point_time, "seconds")
        acc_data = np.array(acc_data)
    else:
        acc_data = np.ones((len(initial_population),2))*-1
    if multi_objective_function[3] != 0:
        rand_seg_acc = []
        for ft in initial_population:
            rand_seg_acc.append([ft.phi_r(dataset, seg_size)])
        rand_seg_acc = np.array(rand_seg_acc)
    else:
        rand_seg_acc = np.ones((len(initial_population),2))*-1

    if multi_objective_function[4] != 0:
        im_data = []
        for ft in initial_population:
            im_data.append([ft.phi_im(dataset)])
        im_data = np.array(im_data)
    else:
        im_data = np.ones((len(initial_population),2))*-1

    if multi_objective_function[5] != 0:
        prec_data = []
        for ft in initial_population:
            prec_data.append([ft.phi_prec(dataset)])
        prec_data = np.array(prec_data)
    else:
        prec_data = np.ones((len(initial_population),2))*-1

    if multi_objective_function[6] != 0:
        spec_data = []
        for ft in initial_population:
            spec_data.append([ft.phi_spec(dataset)])
        spec_data = np.array(spec_data)
    else:
        spec_data = np.ones((len(initial_population),2))*-1

    if multi_objective_function[7] != 0:
        sens_data = []
        for ft in initial_population:
            sens_data.append([ft.phi_sens(dataset)])
        sens_data = np.array(sens_data)
    else:
        sens_data = np.ones((len(initial_population),2))*-1

    if multi_objective_function[8] != 0:
        npv_data = []
        for ft in initial_population:
            npv_data.append([ft.phi_npv(dataset)])
        npv_data = np.array(npv_data)
    else:
        npv_data = np.ones((len(initial_population),2))*-1
    
    if multi_objective_function[9] != 0:
        fnr_data = []
        for ft in initial_population:
            fnr_data.append([ft.phi_fnr(dataset)])
        fnr_data = np.array(fnr_data)
    else:
        fnr_data = np.ones((len(initial_population),2))*-1

    if multi_objective_function[10] != 0:
        fpr_data = []
        for ft in initial_population:
            fpr_data.append([ft.phi_fpr(dataset)])
        fpr_data = np.array(fpr_data)
    else:
        fpr_data = np.ones((len(initial_population),2))*-1
    
    if multi_objective_function[11] != 0:
        acc_single_data = []
        for ft in initial_population:
            acc_single_data.append([ft.phi_acc(dataset)])
        acc_single_data = np.array(acc_single_data)
    else:
        acc_single_data = np.ones((len(initial_population),2))*-1

    fitnesses = np.column_stack((acc_mcs[:,0],tree_size[:,0],acc_data[:,0], rand_seg_acc[:,0], im_data[:, 0], prec_data[:, 0], spec_data[:, 0], sens_data[:, 0], npv_data[:, 0], fnr_data[:, 0], fpr_data[:, 0], acc_single_data[:, 0]))
    for ft,fi1 in zip(initial_population,fitnesses):
        fitness_dict[str(ft)]  = fi1
    
    return fitnesses, fitness_dict


# %% Fitness function
def cost_function(initial_population, dataset, bes, population_size, ft_from_MCSs, multi_objective_function, seg_size=4):
    
    
    # Compute the metrics from the FTs
    fitnesses,fitness_dict = compute_metrics_fts(initial_population,dataset,ft_from_MCSs,multi_objective_function,bes, seg_size)

    # ------------------------------------------------
    # Pareto efficiency:
    a = np.diag(multi_objective_function)
    if multi_objective_function[1] == 0:#Does not optimize size
        b = np.array(fitnesses)
        b[:,1] = -1
    else:   
        b = np.array(fitnesses)
    c = np.matmul(a, b.T).T
    # Perform Non-dominated sorting genetic algorithm II (NSGA-II)
    inputPoints = helper.adjust_array(list(c))
    pt,pt_idx = nsgaII(inputPoints,population_size)
    f=[]
    for i in pt:
        if len(i)>1:
            nds = [fitnesses[i,0] + fitnesses[i,2] + fitnesses[i,11]][0].argsort()
            f.append(np.array(i)[nds].tolist())
        else:
            f.append(i)
    # Update the order leaving the best accuracy always on top
    pt = f
    pt_idx = list(itertools.chain.from_iterable(f))
    raw_fts = [np.flip(np.array(initial_population)[pt_idx],0).tolist(),np.flip(fitnesses[pt_idx],0),fitness_dict]
    # ------------------------------------------------

    return raw_fts


def nsgaII(data, pt_size):
    # Identify all the fronts:
    all_fronts = []
    data_copy = data[:]
    idx_data = list(range(0, len(data)))
    while data:
        _, _, idx = simple_cull(inputPoints=data.copy())  # Identify the pareto front
        c = [e for i, e in enumerate(idx_data) if i in idx]
        c.reverse()  # In this way, we give more priority to those FTs with higher accuracy.
        all_fronts.append(c)
        # all_fronts.append(idx[:])
        idx_data = [j for i, j in enumerate(idx_data) if i not in idx]
        data = [j for i, j in enumerate(data) if i not in idx]

    size_per_front = []
    for i in all_fronts:
        size_per_front.append(len(i))

    cumsum_size_per_front = np.cumsum(np.array(size_per_front))

    final_pt_idx = []  # These are the fronts that will pass to the next generation.
    pt_out = []
    cont = 0
    for i in cumsum_size_per_front <= pt_size:
        if i:
            final_pt_idx.append(cont)
            pt_out.append(all_fronts[cont])
        cont += 1

    # Check whether we need to include part of the next front:
    fronts_left_out = [i for i, x in enumerate(cumsum_size_per_front >= pt_size) if x]
    if cumsum_size_per_front[fronts_left_out[0]] > pt_size:
        # Crowding Distance Sorting.
        cws_data = []
        cws_dx = []
        for i in all_fronts[fronts_left_out[0]]:
            cws_data.append(data_copy[i])
            cws_dx.append(i)
        a = calc_crowding_distance(np.array(cws_data))
        a = np.flip(a.argsort())
        a = a[: - cumsum_size_per_front[fronts_left_out[0]] + pt_size]
        b = []
        for i in a:
            b.append(cws_dx[i])
        pt_out.append(b)

    pt_out_index = np.array([])
    for i in pt_out:
        pt_out_index = np.concatenate((pt_out_index, i))

    return pt_out, pt_out_index.astype(int)


# %% Pareto front
'''
Adapted from this: # https://code.activestate.com/recipes/578287-multidimensional-pareto-front/
'''


def simple_cull(inputPoints):
    inputPoints_ori = inputPoints[:]
    paretoPoints = set()
    candidateRowNr = 0
    dominatedPoints = set()
    idx = []
    while True:
        candidateRow = inputPoints[candidateRowNr]
        inputPoints.remove(candidateRow)
        rowNr = 0
        nonDominated = True
        while len(inputPoints) != 0 and rowNr < len(inputPoints):
            row = inputPoints[rowNr]
            if dominates(candidateRow, row):
                # If it is worse on all features remove the row from the array
                inputPoints.remove(row)
                dominatedPoints.add(tuple(row))
            elif dominates(row, candidateRow):
                nonDominated = False
                dominatedPoints.add(tuple(candidateRow))
                rowNr += 1
            else:
                rowNr += 1

        if nonDominated:
            # add the non-dominated point to the Pareto frontier
            paretoPoints.add(tuple(candidateRow))
            idx.append(inputPoints_ori.index(candidateRow))

        if len(inputPoints) == 0:
            break
    return paretoPoints, dominatedPoints, idx


def dominates(row, candidateRow):
    return sum([row[x] >= candidateRow[x] for x in range(len(row))]) == len(row)


# %% Crowding Distance Sorting
'''
The first Crowding Distance Sorting was taken from:
https://github.com/msu-coinlab/pymoo/blob/20abef1ade71915352217400c11ece4c2f35163e/pymoo/algorithms/nsga2.py

'''


def calc_crowding_distance(F):
    infinity = 1e+14

    n_points = F.shape[0]
    n_obj = F.shape[1]

    if n_points <= 2:
        return np.full(n_points, infinity)
    else:

        # sort each column and get index
        I = np.argsort(F, axis=0, kind='mergesort')

        # now really sort the whole array
        F = F[I, np.arange(n_obj)]

        # get the distance to the last element in sorted list and replace zeros with actual values
        dist = np.concatenate([F, np.full((1, n_obj), np.inf)]) \
               - np.concatenate([np.full((1, n_obj), -np.inf), F])

        index_dist_is_zero = np.where(dist == 0)

        dist_to_last = np.copy(dist)
        for i, j in zip(*index_dist_is_zero):
            dist_to_last[i, j] = dist_to_last[i - 1, j]

        dist_to_next = np.copy(dist)
        for i, j in reversed(list(zip(*index_dist_is_zero))):
            dist_to_next[i, j] = dist_to_next[i + 1, j]

        # normalize all the distances
        norm = np.max(F, axis=0) - np.min(F, axis=0)
        norm[norm == 0] = np.nan
        dist_to_last, dist_to_next = dist_to_last[:-1] / norm, dist_to_next[1:] / norm

        # if we divided by zero because all values in one columns are equal replace by none
        dist_to_last[np.isnan(dist_to_last)] = 0.0
        dist_to_next[np.isnan(dist_to_next)] = 0.0

        # sum up the distance to next and last and norm by objectives - also reorder from sorted list
        J = np.argsort(I, axis=0)
        crowding = np.sum(dist_to_last[J, np.arange(n_obj)] + dist_to_next[J, np.arange(n_obj)], axis=1) / n_obj

    # replace infinity with a large number
    crowding[np.isinf(crowding)] = infinity

    return crowding
