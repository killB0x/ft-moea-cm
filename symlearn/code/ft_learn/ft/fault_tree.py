from ft_learn.ft.ft_elements import BE, AND, OR, VOT
import ft_learn.helper as helper
import numpy as np
from copy import deepcopy
import ft_learn.moea.fitness as fitness
import os
from scipy.io import savemat


class FaultTree:
    """
    Simple representation of fault trees.
    """

    def __init__(self, top_event, name=""):
        """
        Constructor.
        :param top_event: Top-level event (TLE).
        :param name: Name
        """
        self.top_event = top_event
        self.name = name

    def copy(self, allow_duplicates=False):
        """
        Create copy of this fault tree.
        This acts as deepcopy.
        :param allow_duplicates: Whether duplicate children are allowed for gates.
        :return: Copy of fault tree.
        """

        def construct_iter(element, elements, allow_duplicates):
            if isinstance(element, BE):
                if element.name in elements:
                    # Get existing BE
                    return elements[element.name]
                else:
                    # Create new BE
                    be = BE(element.name)
                    elements[element.name] = be
                    return be
            else:
                duplicate = False
                # if not allow_duplicates and str(element) in elements:
                if str(element) in elements:
                    # Use existing gate (only when no duplicates are allowed)
                    # Get existing gate
                    gate = elements[str(element)]
                    # Check whether children coincide
                    # If not, create a new duplicate gate
                    for child in element.children:
                        if child not in gate.children:
                            if allow_duplicates:
                                duplicate = True
                                break
                            else:
                                assert False
                    if not duplicate:
                        # Use existing gate
                        return gate

                # Create new gate (element does not exist yet or duplicate elements can occur)
                # Start by creating name
                name = element.name
                if duplicate:
                    # Create new name
                    name += "a"
                # Create children
                children = [construct_iter(child, elements, allow_duplicates) for child in element.children]  # elements is updated due to call-by-reference
                # Create gate
                if isinstance(element, OR):
                    gate = OR(children, name=name)
                elif isinstance(element, AND):
                    gate = AND(children, name=name)
                elif isinstance(element, VOT):
                    gate = VOT(children, element.threshold, name=name)
                else:
                    assert False
                elements[str(element)] = gate
                return gate

        top_event = construct_iter(self.top_event, dict(), allow_duplicates=allow_duplicates)
        return FaultTree(top_event, name=self.name)

    @classmethod
    def create_from_flat_representation(cls, flat_representation, basic_events):
        # TODO: find better representation
        def construct_children(element):
            if "BE" in element:
                return basic_events[element]
            else:
                children = []
                for child in element[1]:
                    children.append(construct_children(child))
                if element[0] == "OR":
                    return OR(children)
                elif element[0] == "AND":
                    return AND(children)
                else:
                    assert False

        return FaultTree(construct_children(flat_representation))

    def remove_gate(self, gate):
        """
        Remove gate from fault tree and update parent/children relations.
        :param gate: Gate to remove.
        """
        # Remove the gate as parent of its children.
        # We do not use gate.remove_child() as we still might want to access the children of the gate
        for child in gate.children:
            child.parents.remove(gate)
        # Remove gate from parents
        if gate != self.top_event:
            assert len(gate.parents) >= 1
            for parent in gate.parents:
                parent.remove_child(gate)

    def evaluate(self, be_values):
        """
        Evaluate the fault tree with respect to the given BE values.
        :param be_values: Boolean values for each BE.
        :return: Evaluation result of the top-level event.
        """
        return self.top_event.evaluate(be_values)

    def get_all_bes(self, sort=False):
        # TODO: avoid re-computation
        def iter_all_bes(elem):
            if isinstance(elem, BE):
                return {elem}
            else:
                return set().union(*[iter_all_bes(child) for child in elem.children])

        all_bes = list(iter_all_bes(self.top_event))
        if sort:
            all_bes.sort(key=lambda x: x.name)
        return all_bes

    def get_all_gates(self, sort=False):
        # TODO: avoid re-computation
        # TODO: avoid duplicates
        def iter_all_gates(elem):
            if isinstance(elem, BE):
                return []
            else:
                s = [elem]
                for child in elem.children:
                    s.extend(iter_all_gates(child))
                return s

        all_gates = iter_all_gates(self.top_event)
        if sort:
            all_gates.sort(key=lambda x: str(x))
        return all_gates

    def count_connections(self):
        def iter_connections(elem):
            if isinstance(elem, BE):
                return 1
            else:
                no = 1
                for child in elem.children:
                    no += iter_connections(child)
                return no

        return iter_connections(self.top_event)

    def to_string(self, include_names=True):
        """
        Return string representation of the fault tree.
        :param include_names: Whether to include gate names into string output.
        :return: String representation.
        """
        return self.top_event.to_string(include_names=include_names)

    def print_parents(self):

        def iter_parents(elem):
            s = "{}: {}".format(elem, elem.parents)

            if isinstance(elem, BE):
                return s
            else:
                children = [iter_parents(child) for child in elem.children]
                return s + "\n" + "\n".join(children)

        return iter_parents(self.top_event)

    def __eq__(self, other):
        if isinstance(other, FaultTree):
            # Only compare the structure but not the naming of the gates
            return self.to_string(include_names=False) == other.to_string(include_names=False)
        else:
            return False

    def __str__(self):
        return self.to_string(include_names=False)

    def simplify(self):
        """
        Simplify the fault tree by applying simple rewrite rules:
        1. Remove empty gates without children.
        2. Deleting gates with a single successors.
        3. Removing consecutive OR- or AND-gates
        :return: True iff the fault tree was simplified.
        """

        def iter_simplify(element):
            """
            Iterator to apply simplifications on the fault-tree structure.
            Needs iterative application until a fixpoint is reached.
            :param element: Current element.
            :return: True iff some simplification took place. This does not mean that all possible simplifications took place.
            """
            if isinstance(element, BE):
                # Cannot simplify BEs
                return False

            # Delete emtpy gates without children
            # ----------
            if len(element.children) == 0:
                for parent in list(element.parents):
                    parent.remove_child(element)
                return True

            # Delete gates with a single successor
            # ----------
            if len(element.children) == 1:
                child = element.children[0]
                for parent in list(element.parents):
                    parent.remove_child(element)
                    if child not in parent.children:
                        parent.add_child(child)
                return True

            # Remove consecutive OR- or AND-gates
            # ----------
            simplified = False
            for child in list(element.children):
                if type(child) == type(element):
                    # Remove child and add children to element
                    element.remove_child(child)
                    for grandchild in list(child.children):
                        child.remove_child(grandchild)
                        if grandchild not in element.children:
                            element.add_child(grandchild)
                    simplified = True

            if simplified:
                return True

            # Repeated elements under the same gate cannot occur because add_child() does not allow this
            # TODO: check for differently named gates?

            # Recursively apply simplifications
            for child in element.children:
                if iter_simplify(child):
                    return True

            return False

        changed = True
        while changed:
            changed = iter_simplify(self.top_event)
    

    def phi_d(self,dataset):
        count_true = 0
        total_counts = 0
        for data in dataset:
            if data['T'] == self.evaluate(data):
                count_true += data['N']
            total_counts += data['N']
        return 1 - count_true/total_counts
        
# Fitness function:
    def phi_c(self,cs_base,bes):
                     
        # We obtain the cut sets of both trees:
        #c1 = helper.getMCSs(helper.cutsets_from_ft(self),bes)
        c1 = helper.getMCSs(helper.cutsets_from_ft(self,bes))
        c2 = cs_base
               
        if len(c1.shape)==1:
            c1 = c1.reshape((1,c1.shape[0]))
        
        if len(c2.shape)==1:
            c2 = c2.reshape((1,c2.shape[0]))
        
        # Make cut set matrices the same size:
        if c1.shape[1] > c2.shape[1]:
            t1 = np.zeros((c2.shape[0],c1.shape[1]), dtype=int)
            for i in range(0,c2.shape[1]):
                if sum(c2[:,i])>0:
                    t1[:,i] = c2[:,i]
            c2 = t1
        elif c1.shape[1] < c2.shape[1]:
            t1 = np.zeros((c1.shape[0],c2.shape[1]), dtype=int)
            for i in range(0,c1.shape[1]):
                if sum(c1[:,i])>0:
                    t1[:,i] = c1[:,i]
            c1 = t1
            
        # RV-Correlation:
        rv_corr = np.trace(np.matmul(np.matmul(c1.T,c1).T,np.matmul(c2.T,c2)))/np.sqrt(np.trace(np.matmul(np.matmul(c1.T,c1).T,np.matmul(c1.T,c1)))*np.trace(np.matmul(np.matmul(c2.T,c2).T,np.matmul(c2.T,c2))))
        return 1 - rv_corr
    
    def phi_s(self):
        return self.count_connections()
        # return len(self.get_all_bes()) + len(helper.get_gate_statistics(self))

    def phi_r(self,dataset):
        shuffled_dataset = np.copy(dataset)
        np.random.shuffle(shuffled_dataset)
        segment_size = 8
        segment_count = 0
        segment_count_true = 0
        total_counts = 0
        total_counts_true = 0
        for data in shuffled_dataset:
            if (segment_count == segment_size):
                if (segment_count_true == segment_size):
                    total_counts_true += 1
                total_counts += 1
                segment_count = 0
                segment_count_true = 0
            if data['T'] == self.evaluate(data):
                segment_count_true += 1
            segment_count += 1
        if (segment_count == segment_size):
            if (segment_count_true == segment_size):
                total_counts_true += 1
            total_counts += 1
            segment_count = 0
        return 1 - total_counts_true / (total_counts)
    
    def get_unique_list_bes(self):
        bes = []
        for i in self.get_all_bes(): 
            bes.append(str(i))
        return bes

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
# %% Create FT object
def str2ft(ft1):
    
    # Delete the empty spaces:
    ft1 = delete_str_index(ft1,strfind(ft1,' '))
    
    # Example:
    # ft1 = 'OR(BE1,AND(BE3,OR(BE4,AND(BE2,BE5))),BE2,OR(BE1,BE2))'
    _, b, c = helper.get_ft_structure_paths(ft1)

    be = []
    for i in c:
        be = be + i[3]
    used = set()
    be = [x for x in be if x not in used and (used.add(x) or True)]

    # Create the basic events:
    be_s = []
    for i in be:
        be_s.append('BE' + str(i))

    basic_events = {}
    for be in be_s:
        basic_events[be] = BE(be)
    # FT(flattened_representation=ft1_flat, top_event=None, basic_events=basic_events)
    ft = FaultTree.create_from_flat_representation([c[0][0], be_s], basic_events)

    # Delete all the elements below the root event (i.e., top event):
    gate = ft.get_all_gates()[0]
    selected_children = gate.children
    while len(selected_children) > 0:
        gate.children.remove(selected_children[0])

    # Create all the gates:
    dict_gates = {}
    for i in c:
        if i[1][0:i[1].index('_')] == 'AND':
            dict_gates[i[1]] = AND([])
        elif i[1][0:i[1].index('_')] == 'OR':
            dict_gates[i[1]] = OR([])

    # Update the information of the root node:
    dict_gates[c[0][1]] = ft.get_all_gates()[0]

    for i in c:
        # Add the children gates below the parent gate:
        for j in i[2]:
            dict_gates[i[1]].children.append(dict_gates[j])

        # Add BEs under the parent:
        for j in i[3]:
            dict_gates[i[1]].children.append(BE('BE' + str(j)))

    return ft

#%%
def trim_ft(ft):
    if len(ft.get_all_gates()[0].children)>0:
        new_ft = deepcopy(ft)
        master_cont = 1
        while master_cont > 0:
            master_cont = 0
            #%%
            '''
            Task # 1: Remove unnecessary OR gates
            This happens for example when to consecutive OR gates are connected.
            '''  
            idx1 = 1
            while idx1 == 1:
                all_gates_ft = new_ft.get_all_gates()
                len_k = []
                for OR_gate in all_gates_ft:    
                    # Check if the current gate is OR:
                    if type(OR_gate).__name__ == 'OR':
                        # Check if any of the children in the gate is type OR:
                        in_OR_gate = OR_gate.children
                        for i in in_OR_gate:
                            if type(i).__name__ == 'OR':
                                # We move the BEs children of the child OR gate to the parent OR gate:
                                in_OR_gate_child = i.children
                                len_k.append(len(in_OR_gate_child))
                                for j in in_OR_gate_child:    
                                    parent_of_be = new_ft.get_parent_of_be(j)
                                    new_parent   = OR_gate
                                    parent_of_be.children.remove(j)
                                    new_parent.children.append(j)
                                    
                if sum(len_k) == 0:
                    idx1 = 0
                else:
                    idx1 = 1
                    master_cont += 1
                                        
            #%%
            '''
            Task # 2: Delete all the emtpy gates
            Meaning all the gates that have zero children
            '''      
            idx1 = 1             
            while idx1 == 1:
                
                all_gates_ft = new_ft.get_all_gates()
                gate_to_detele = []
                len_k = []
                
                for gate in all_gates_ft: 
                    # Check how many children has the current gate:
                    children_in_gate = gate.children
                    if len(children_in_gate) == 0:
                        gate_to_detele.append(gate)
                        
                len_k.append(len(gate_to_detele))
                # Delete the identified gates
                for i in gate_to_detele:
                    parent_of_gate = new_ft.get_parent_of_be(i)
                    parent_of_gate.children.remove(i)

                
                if sum(len_k) == 0:
                    idx1 = 0
                else:
                    idx1 = 1
                    master_cont += 1
                    
                
            #%%
            '''
            Task # 3: Delete gates with a single basic event
            '''     
            idx1 = 1             
            while idx1 == 1:
                all_gates_ft = new_ft.get_all_gates()
                gate_to_detele = []
                len_k = []
                for gate in all_gates_ft: 
                    # Check how many children has the current gate:
                    children_in_gate = gate.children
                    
                    
                    if len(all_gates_ft) == 1 and len(children_in_gate) == 1:
                        # This is the case where there is a single gate with a single element.
                        break
                        
                    
                    if len(children_in_gate) == 1:
                        if type(children_in_gate[0]).__name__ == 'BE':
                            gate_to_detele.append(gate)
                len_k.append(len(gate_to_detele))
                for i in gate_to_detele:
                    new_parent   = new_ft.get_parent_of_be(i) # Identify the new parent gate
                    new_parent.children.append(i.children[0]) # Add the child to the new parent gate
                    new_parent.children.remove(i) # Delete the whole gate that had a single child
                    
                if sum(len_k) == 0:
                    idx1 = 0
                else:
                    idx1 = 1
                    master_cont += 1
            
            #%%
            '''
            Task # 4: Delete gates with a single child gate
            '''     
            idx1 = 1             
            while idx1 == 1:
                all_gates_ft = new_ft.get_all_gates()
                gate_to_detele = []
                len_k = []
                for gate in all_gates_ft: 
                    # Check how many children has the current gate:
                    children_in_gate = gate.children
                    if len(children_in_gate) == 1:
                        if type(children_in_gate[0]).__name__ != 'BE':
                            gate_to_detele.append(gate)
                len_k.append(len(gate_to_detele))
                for i in gate_to_detele:
                    # If the gate is the root note, it won't have parent gate, first check that
                    if new_ft.get_all_gates()[0]==i:
                        # In this case, we have to eliminate the root gate.
                        # We create a new FT with an AND/OR gate in the root, and a single child, which is going to be deleted later on
                        temporary_basic_event = new_ft.basic_events
                        if type(i.children[0]).__name__ == "AND":
                            ft2_flat = ["AND",new_ft.flattened_representation[1]]
                        else:
                            ft2_flat = ["OR",new_ft.flattened_representation[1]]
                        ft_new = FaultTree(flattened_representation=ft2_flat, top_event=None, basic_events=temporary_basic_event)
                        temporary_basic_event = ft_new.get_all_bes()
                        new_parent   = ft_new.get_all_gates()[0]
                        # BEs of the previous gate:
                        BEs_old_gate = i.children[0].children
                        for j in BEs_old_gate:
                            new_parent.children.append(j)
                        for j in temporary_basic_event:    
                            new_parent.children.remove(j) # Delete the temporary children
                        new_ft = deepcopy(ft_new) # We delete the previous tree and remplace it with the new one.
                        break
                    else:
                        # If the children of the gate is an AND gate, then the children should be moved up:
                        if type(i.children[0]).__name__ == 'AND':
                            new_parent   = new_ft.get_parent_of_be(i) # Identify the new parent gate
                            new_parent.children.append(i.children[0]) # Add the child to the new parent gate
                            new_parent.children.remove(i) # Delete the whole gate that had a single child
                        elif type(i.children[0]).__name__ == 'OR':
                        # If the children of the gate is an OR gate, then the children gate is eliminated and the parent inherit the children:
                            #new_parent   = i # Identify the new parent gate
                            #new_parent.children.append(i.children[0].children) # Add the child to the new parent gate
                            #new_parent.children.remove(i.children[0]) # Delete the whole gate that had a single child
                            new_parent   = new_ft.get_parent_of_be(i) # Identify the new parent gate
                            new_parent.children.append(i.children[0]) # Add the child to the new parent gate
                            new_parent.children.remove(i) # Delete the whole gate that had a single child
                        
                if sum(len_k) == 0:
                    idx1 = 0
                else:
                    idx1 = 1
                    master_cont += 1
            #%%
            '''
            Task # 5: Identify critical elements, and detele them from anywhere else in the tree
            '''
            # idx1 = 1             
            # while idx1 == 1:
            #     list_critical_elements = []
            #     len_k = []
            #     if len(new_ft.get_all_gates())>1:
            #         for i in new_ft.get_all_gates()[0].children:
            #             if type(i).__name__ == 'BE':
            #                 list_critical_elements.append(i.name)
                            
                    
            #         for i in new_ft.get_all_gates()[1:]:
            #             child_name = []
            #             child_obj  = []
            #             for j in i.children:
            #                 try:
            #                     child_name.append(j.name)
            #                     child_obj.append(j)
            #                 except: 
            #                     pass
                        
            #             j = ismember(child_name,list_critical_elements)
            #             if any(j):
            #                 len_k.append(len(j))
            #                 cont = 0
            #                 # Delete the BEs that are critical from the current gate
            #                 for k in j:
            #                     if k is not False:
            #                         if type(i).__name__ == "OR":
            #                             new_parent   = i
            #                             new_parent.children.remove(child_obj[cont])
            #                     cont += 1
                            
                            
            #     if sum(len_k) == 0:
            #         idx1 = 0
            #     else:
            #         idx1 = 1
            #         master_cont += 1
            #%% 
            '''
            Task # 6: Delete repeated BEs under the same parent gate
            '''
            idx1 = 1             
            while idx1 == 1:
                len_k = []
                all_gates_ft = new_ft.get_all_gates()
                for i in all_gates_ft:
                    child_name = []
                    child_obj  = []
                    for j in i.children:
                        try:
                            child_name.append(j.name)
                            child_obj.append(j)
                        except: 
                            pass
        
                    to_del = []
                    while len(child_name)>0:
                        idx1 = helper.get_index_positions(child_name,child_name[0])
                        if len(idx1)>1:
                            for p in idx1[1:]:
                                to_del.append(child_obj[p])
                            for p in sorted(idx1, reverse=True):
                                del(child_name[p])
                                del(child_obj[p])
                        else:
                            del(child_name[0])
                            del(child_obj[0])
                            
                            
                    if len(to_del)>0:
                        new_parent   = i
                        for p in to_del:
                            new_parent.children.remove(p)
                            
                
                    len_k.append(len(to_del))
                    
        
                if sum(len_k) == 0:
                    idx1 = 0
                else:
                    idx1 = 1
                    master_cont += 1
            
            master_cont
    else:
        # The FT has no children.
        new_ft = []
    return new_ft

# %%
def save_results(raw_fts,t,path_save_results,dataset,ft_from_MCSs,multi_objective_function):

    str_fts_sorted = []
    str_fts_sorted_trimmed = []
    #fts_sorted_trimmed = []
    for i in raw_fts[0]:
        str_fts_sorted.append(str(i))
        #fts_sorted_trimmed.append(trim_ft(str2ft(str(i))))
        #str_fts_sorted_trimmed.append(str(fts_sorted_trimmed[-1]))
        
    # Compute the metrics from the FTs
    #fitnesses_trimmed,_ = fitness.compute_metrics_fts(fts_sorted_trimmed,dataset,ft_from_MCSs,multi_objective_function)

    files_in_folder = []
    for name in os.listdir(path_save_results+ '/' ):
        if os.path.isfile(os.path.join(path_save_results+ '/', name)):
            files_in_folder.append(name)
            
    num_files_dir = int(len(files_in_folder))
    saving_dir = path_save_results + '/' + 'gen_' + str(num_files_dir)
    #mdic = {'metrics_fts': raw_fts[1],'metrics_fts_trimmed':fitnesses_trimmed,'fts':str_fts_sorted,'fts_trimmed':str_fts_sorted_trimmed,'time':t}
    mdic = {'metrics_fts': raw_fts[1],'fts':str_fts_sorted,'time':t}
    savemat(saving_dir+'.mat', mdic)
