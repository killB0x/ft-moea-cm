import numpy as np
import re
from copy import deepcopy
from sympy.logic.boolalg import to_dnf


from ft_learn.ft.mcs import MinCutSets
from ft_learn.ft.ft_elements import BE, AND, OR
import ft_learn.logic.boolean_logic as boolean_logic
import os


def cutsets_from_ft(ft, bes):
    """
    Get cut set matrix from fault tree.
    :param ft: Fault tree.
    :param bes: List of BEs.
    :return: Cut-set matrix corresponding to given fault tree.
    """
    formula = boolean_logic.fault_tree_to_sympy_formula(ft)
    min_cut_sets = boolean_logic.sympy_formula_to_mcs(formula, bes)
    be_indices = [i for i in range(len(bes))]
    return min_cut_sets.get_matrix(be_indices)


# TODO: could use revision

# %%
def del_list_indexes(l, id_to_del):
    somelist = [i for j, i in enumerate(l) if j not in id_to_del]
    return somelist


# %%
def compare(a, b):
    """
    Compare two base strings, disregarding whitespace
    """
    return re.sub("\s*", "", a) == re.sub("\s*", "", b)


# %%
def flatten(lis):
    flat_list = [item for sublist in lis for item in sublist]
    return flat_list


# %%
def check_empty_objects(ft):
    if str(ft).find("()") != -1:
        idx = True
    else:
        idx = False
    return idx


# %%
def cs2str(A):
    # A = np.array([[0, 1, 0, 1, 0, 0, 0, 1, 0],
    #        [1, 0, 1, 1, 0, 0, 0, 1, 0],
    #        [1, 0, 0, 1, 1, 0, 0, 1, 0],
    #        [1, 0, 0, 1, 0, 1, 0, 1, 0],
    #        [1, 0, 0, 1, 0, 0, 1, 1, 0],
    #        [1, 0, 0, 1, 0, 0, 0, 1, 1]]).tolist()

    all_IE = []
    for i in A:
        j = (np.where(np.array(i) == 1)[0] + 1).astype(int).tolist()
        IE = '('
        for k in j:
            IE += 'BE' + str(k) + ' & '
        IE = IE[0:-2] + ')'
        all_IE.append(IE)
    str_ = ''
    for i in all_IE:
        str_ += i + ' | '
    str_ = str_[0:-2]
    # str_cutsets  = to_dnf(str_)

    return str_


# %%
def get_ft_structure_paths(ft1):
    op = [m.start() for m in re.finditer('\(', ft1)]
    cp = [m.start() for m in re.finditer('\)', ft1)]
    cm = np.array([m.start() for m in re.finditer(',', ft1)])

    m = []
    for i in op:
        m.append(np.array(cp) - i)
    m = np.array(m)

    k = []
    cp_, op_ = deepcopy(cp), deepcopy(op)
    while len(m) > 0:
        idx = np.where(m == np.min(m[m > 0]))
        k.append([op_[idx[0][0]], cp_[idx[1][0]]])
        m = np.delete(np.delete(m, idx[0][0], 0), idx[1][0], 1)
        cp_.pop(idx[1][0]), op_.pop(idx[0][0])
    del (cp_, op_)

    k = np.array(k)
    m = []
    for i in range(0, len(k)):
        m.append(np.logical_and(list(k[i, 0] < k[:, 0]), list(k[i, 1] > k[:, 1])))
    m.reverse()

    ft_str = {}
    cont_and, cont_or = 0, 0
    for i in k:
        if any([any(i[0] > cm), any(i[0] > op)]):
            # Looks for the closest comma:
            if len(cm[i[0] > cm]) > 0:
                a = [max(cm[i[0] > cm])]
            else:
                a = []
            # Looks for the closest parenthesis:
            b = [max(np.array(op)[i[0] > op])]
            c = max(a + b)
        else:
            c = -1
        if ft1[c + 1:i[0]] == 'AND':
            cont_and += 1
            ft_str[ft1[c + 1:i[1] + 1]] = 'AND_' + str(cont_and)
        elif ft1[c + 1:i[0]] == 'OR':
            cont_or += 1
            ft_str[ft1[c + 1:i[1] + 1]] = 'OR_' + str(cont_or)

    ft_str_final = {}
    cont = 1
    while cont == 1:
        a = list(ft_str.keys())
        cont = 0

        for i in a:
            keys_all = []
            for j in a:
                if j.find(i) != -1 and i != j:
                    keys_all.append(j)
                    cont = 1

            for j in keys_all:
                old_j = j
                ft_str[j.replace(i, ft_str[i])] = ft_str[j]
                del ft_str[old_j]

            if len(keys_all) > 0:
                # To update the keys
                break

    ft_str_final = dict(map(reversed, ft_str.items()))
    ft_str_final_list = [(k, v) for k, v in ft_str_final.items()]
    ft_str_final_list.reverse()

    # %% Check the paths per BE:
    table_structure = []
    for i in ft_str_final_list:
        temp1 = re.split(',|\(|\)', i[1])
        temp1 = del_list_indexes(temp1, [i for i, x in enumerate(temp1) if x == ''])
        del temp1[0]

        # Look for the children BEs:
        be_s = []
        for j in temp1:
            if j.find("BE") != -1:
                # There is a BE under the gate:
                be_s.append(int(j.strip()[2:]))

        # Look for the children gates:
        gates = []
        if len(temp1) > 0:
            for j in temp1:
                if j.find("BE") == -1:
                    # There is at least one gate:
                    gates.append(j)

        table_structure.append([i[0][0:i[0].index('_')], i[0], gates, be_s])

    for i in table_structure:

        if i[1][0:i[1].index('_')] == 'AND':
            i[1] = "AND_" + str(- int(i[1][i[1].index('_') + 1:]) + cont_and + 1)
        elif i[1][0:i[1].index('_')] == 'OR':
            i[1] = "OR_" + str(- int(i[1][i[1].index('_') + 1:]) + cont_or + 1)

        if len(i[2]) > 0:
            new_var = []
            for j in i[2]:
                if j[0:j.index('_')] == 'AND':
                    new_var.append("AND_" + str(- int(j[j.index('_') + 1:]) + cont_and + 1))
                elif j[0:j.index('_')] == 'OR':
                    new_var.append("OR_" + str(- int(j[j.index('_') + 1:]) + cont_or + 1))
            i[2] = new_var

    list_table_str = []
    for i in table_structure:
        for j in i[2]:
            list_table_str.append([i[1], j])

        for j in i[3]:
            list_table_str.append([i[1], "BE" + str(j)])

    paths = []
    for j in list_table_str:
        # Check that the first node is a root node:
        if any(np.array(list_table_str, dtype='object')[:, 1] == j[0]) == False:
            master_var = [j]
            cont = 1
            while cont == 1:
                cont = 0
                new_master_var = []
                for k in master_var:
                    indx = [i for i, x in enumerate(np.array(list_table_str, dtype='object')[:, 0].tolist()) if x == k[-1]]
                    if len(indx) > 0:
                        cont = 1
                        for p in indx:
                            new_master_var.append(k + [list_table_str[p][1]])
                    else:
                        new_master_var.append(k)

                if cont == 1:
                    # Here you can update the master variable.
                    master_var = new_master_var
                else:
                    # You add the results to the path:
                    paths = paths + master_var

    # Make paths into strings:
    path_str = []
    path_str_num = []
    for i in paths:
        s = '-'
        s2 = '-'
        i.reverse()
        for j in i:
            s = s + j + '-'
            if str(j).find('_') != -1:
                if j[0:j.index('_')] == 'AND':
                    s2 = s2 + 'AND' + '-'
                elif j[0:j.index('_')] == 'OR':
                    s2 = s2 + 'OR' + '-'
            else:
                s2 = s2 + j + '-'

        path_str.append(s2[1:-1])
        path_str_num.append(s[1:-1])

    return path_str, path_str_num, table_structure


#%% Get indexes from string
def strfind(string,s):
    lst= []
    for i in range(len(string)):
        if (string[i] == s):
            lst.append(i)
    return lst
#%%
def delete_str_index(string, indices):
    z = bytearray(string.encode())
    for i in sorted(indices, key=abs, reverse=True):
        del z[i]
    return z.decode()
#%%

def cutsets_from_ft_string(ft_string):
    
    # Delete the empty spaces:
    ft_string = delete_str_index(ft_string,strfind(ft_string,' '))
    
    _, _, tab = get_ft_structure_paths(ft_string)

    cond = {}
    for i in tab:
        str_ = ''

        # What is the type of gate?
        if i[0] == 'AND':
            g = '&'
        elif i[0] == 'OR':
            g = '|'
        else:
            g = ""
            assert False

        # Add the basic events in that gate:
        for j in i[-1]:
            str_ = str_ + 'BE' + str(j) + ' ' + g + ' '

        # Add the other gates within the parent gate.
        for j in i[-2]:
            str_ = str_ + ' ' + j + ' ' + g

        cond[i[1]] = str_[0:-2]

    tab.reverse()

    for i in tab:
        if len(i[2]) > 0:
            # Parent string:
            ps = cond[i[1]]
            ps2 = re.split('[ ]', ps)
            # Children string:
            for j in i[2]:
                cs = cond[j]
                # idxs = [i for i, x in enumerate(ps2) if x == j]
                idxs = [k for k, x in enumerate(ps2) if compare(x, j)]
                for k in idxs:
                    ps2[k] = ps2[k].replace(ps2[k], '(' + cs + ')')
            ps = ' '.join(ps2)
            cond[i[1]] = ps

    output_string = cond[tab[-1][1]]
    cutsets = to_dnf(output_string)

    # Compute the length of each cutset
    a = re.split('\|', str(cutsets))
    b = []
    for i in a:
        b.append(len(re.split('&', i)))

    # Construct the cut set matrix:
    cutset_length = b
    be_all = []
    for i in a:
        j = re.split('[&()BE]', i)
        be = []
        for k in j:
            try:
                be.append(int(k))
            except:
                continue
        be_all.append(be)
    max_be = max(list(flatten(be_all)))
    del a
    cutset_matrix = []
    for i in be_all:
        a = np.zeros(max_be)
        a[np.array(i) - 1] = 1
        cutset_matrix.append(a.astype(int))

    if len(cutset_matrix) == 1:
        cutset_matrix = cutset_matrix[0]

    mincutsets = MinCutSets()
    mincutsets.compute_from_cut_sets(np.array(cutset_matrix))
    be_indices = [i for i in range(len(cutset_matrix[0]))]
    cutset_matrix = mincutsets.get_matrix(be_indices)

    return cs2str(cutset_matrix), np.sum(cutset_matrix, axis=1), cutset_matrix


# %%
def saving_results(raw_fts, trimmed_fts, time, saving_folder):
    str_fts_sorted = []
    str_fts_sorted_trimmed = []
    for i, j in zip(raw_fts[0], trimmed_fts[0]):
        str_fts_sorted.append(str(i))
        str_fts_sorted_trimmed.append(str(j))

    files_in_folder = []
    for name in os.listdir(saving_folder + '/'):
        if os.path.isfile(os.path.join(saving_folder + '/', name)):
            files_in_folder.append(name)

    num_files_dir = int(len(files_in_folder))
    saving_dir = saving_folder + '/' + 'gen_' + str(num_files_dir)
    mdic = {'parameters': raw_fts[1], 'parameters_trimmed': trimmed_fts[1], 'fts_sorted': str_fts_sorted, 'fts_sorted_trimmed': str_fts_sorted_trimmed, 'time': time}
    savemat(saving_dir + '.mat', mdic)




#%%
def cutsets2ft(cut_sets,bes):
    bes = np.array(bes)
    fts = 'OR('
    for i in cut_sets.tolist():
        j = 'AND('
        for k in bes[np.array(i)==1]:
            j = j + k + ','
        j = j[0:-1] + '),'
        fts = fts + j
    fts = fts[0:-1] + ')'
    return fts
#%%
def getMCSs(all_comb_failure):
    if len(all_comb_failure.shape)>1:
        cut_sets = [];
        while all_comb_failure.shape[0]>0:
            # Find minimal combinations:
            min_comb = sum(all_comb_failure.T)
            min_cut = np.where(min_comb == min_comb.min())[0][0]
            ci = all_comb_failure[min_cut,:]
            to_delete = np.where(all_comb_failure[:,ci==1].all(1))[0]
            all_comb_failure = np.delete(all_comb_failure,to_delete,0)
            if len(to_delete)>0:
                cut_sets.append(ci) 
        cut_sets = np.array(cut_sets).astype(int)
    else:
        cut_sets = all_comb_failure.reshape(1,len(all_comb_failure))    
    return cut_sets

#%%
def get_gate_statistics(ft):#Marijn
    result = []
    gates = ft.get_all_gates()
    i_AND = 0
    i_OR = 0
    
    for gate in gates:
        
        if isinstance(gate, AND):
            gate_name = 'AND_' + str(i_AND)
            i_AND+=1
        if isinstance(gate, OR):
            gate_name = 'OR_' + str(i_OR)
            i_OR+=1
        
        N_BE = 0
        for child in gate.children:
            if isinstance(child, BE):
                N_BE +=1
        
        result.append([gate_name, N_BE])
    
    return result
#%%
def adjust_array(a): 
    b = []
    for i in a:
        b.append([float(j) for j in i])
    return b
#%%

    


#%%
def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list