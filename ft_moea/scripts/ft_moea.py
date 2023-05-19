#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27

FT-MOEA

Main script to infer Fault Tree models via multi-objective evolutionary algorithms. 

This script was built upon the existing implementation made by Linard et al. [1]. 
Relevant sources used to build this code are enlisted below:
    
[1] A. Linard, D. Bucur, and M. Stoelinga, “Fault trees from data: Efficient 
learning with an evolutionary algorithm,” in International Symposium on Dependable 
Software Engineering: Theories, Tools, and Applications. Springer, 2019, pp. 19–37.]
    (https://link.springer.com/chapter/10.1007/978-3-030-35540-1_2)

* From _pymoo by Julian Blank - Multi-Objective Optimization Framework - the [NSGA-II] algorithm
https://github.com/anyoptimization/pymoo/blob/20abef1ade71915352217400c11ece4c2f35163e/pymoo/algorithms/nsga2.py

* From Multi-dimensional pareto front by Jamie Bull: https://code.activestate.com/recipes/578287-multidimensional-pareto-front/

@author: Lisandro A. Jimenez Roa (l.jimenezroa@utwente.nl)
"""

import random
from copy import deepcopy
import numpy as np
import math
import sys 
from sympy.logic import simplify_logic, bool_map
from scipy.io import loadmat, savemat
import os, os.path
import time
import itertools
import re
from sympy.logic.boolalg import to_dnf, to_cnf


class AND():
    """
        Class for an AND gate
    """
    def __init__(self, children=[]):
        """
            Constructor for AND gate.
                children : the children of the AND gate
        """
        self.children = children
        
        
    def valuate(self, be_values=None):
        if len(self.children)>0:
            return all([child.valuate(be_values) for child in self.children])
        else:
            return float("NaN")


class OR():
    """
        Class for an OR gate
    """
    def __init__(self, children=[]):
        """
            Constructor for OR gate.
                children : the children of the OR gate
        """
        self.children = children
        
        
    def valuate(self, be_values=None):
        if len(self.children)>0:
            return any([child.valuate(be_values) for child in self.children])
        else:
            return float("NaN")
        
        
class VoT():
    """
        Class for an VoT gate
    """
    def __init__(self, children=[]):
        """
            Constructor for OR gate.
                children : the children of the OR gate
        """
        self.children = children
        
        
    def valuate(self, be_values=None):
        if len(self.children)>0:
            return any([child.valuate(be_values) for child in self.children])
        else:
            return float("NaN")
        

class BE():
    """
        Class for the Basic Event
    """
    def __init__(self, name=""):
        """
            Constructor for BE.
                name  : 
        """
        self.name  = name
        
        
    def valuate(self,be_values):
        try:
            return be_values[self.name]
        except KeyError:
            return 0
        








class FT():
    """
        Class for Fault Tree (FT)
    """
    
    def __init__(self, flattened_representation=[], top_event=None, basic_events={}):
        """
            Constructor for FT class. Returns a FT
                flattened_representation : the flattened representation of the Fault Tree
                basic_events : set of basic events
        """
        self.flattened_representation = flattened_representation
        self.top_event = top_event
        self.basic_events = basic_events
        
        def construct_children(gate):
            if gate[0] == "OR":
                children = []
                for child in gate[1]:
                    children.append(construct_children(child))
                return OR(children)
            
            elif gate[0] == "AND":
                children = []
                for child in gate[1]:
                    children.append(construct_children(child))
                return AND(children)
            
            # elif "BE" in gate:
            else:
                return basic_events[gate]
        
        #Process flattened_representation
        self.top_event = construct_children(flattened_representation)



    def __str__(self):
        
        def parse_children(gate):
            if isinstance(gate, AND):
                s = 'AND('
                for child in gate.children:
                    s += parse_children(child) + ','
                return s + ')'
            elif isinstance(gate, OR):
                s = 'OR('
                for child in gate.children:
                    s += parse_children(child) + ','
                return s + ')'
            elif isinstance(gate, BE):
                return gate.name
                
        return parse_children(self.top_event).replace(',)',')')



    def __eq__(self, other):
        return str(self) == str(other)


    def __le__(self, other):
        return len(str(self)) <= len(str(self))


    def __gt__(self, other):
        return len(str(self)) > len(str(self))


    def list(self):
    
        def parse_children(gate):
            if isinstance(gate, AND):
                s = ["AND"]
                children = []
                for child in gate.children:
                    children.append(parse_children(child))
                s.append(children)
                return s
            elif isinstance(gate, OR):
                s = ["OR"]
                children = []
                for child in gate.children:
                    children.append(parse_children(child))
                s.append(children)
                return s
            elif isinstance(gate, BE):
                return gate.name
                
        return parse_children(self.top_event)



    def get_all_gates(self):
    
        def parse_children(gate):
            if isinstance(gate, AND):
                s = [gate]
                for child in gate.children:
                    s.extend(parse_children(child))
                return s
            elif isinstance(gate, OR):
                s = [gate]
                for child in gate.children:
                    s.extend(parse_children(child))
                return s
            else:
                return []

        return parse_children(self.top_event)



    def get_parent_of_be(self,be):
        
        def parse_children(gate,be):
            if isinstance(gate, AND):
                if be in gate.children:
                    s = [gate]
                    return s
                s = []
                for child in gate.children:
                    s.extend(parse_children(child,be))
                return s
            elif isinstance(gate, OR):
                if be in gate.children:
                    s = [gate]
                    return s
                s = []
                for child in gate.children:
                    s.extend(parse_children(child,be))
                return s
            else:
                return []

        return parse_children(self.top_event,be)[0]



    def get_all_bes(self):
    
        def parse_children(gate):
            if isinstance(gate, AND):
                s = []
                for child in gate.children:
                    s.extend(parse_children(child))
                return s
            elif isinstance(gate, OR):
                s = []
                for child in gate.children:
                    s.extend(parse_children(child))
                return s
            elif isinstance(gate, BE):
                return [gate]
                
        return parse_children(self.top_event)




    
    def valuate(self,be_values):
        """
            Valuates a FT given BE values
                return the value of the Top Event
        """
        return self.top_event.valuate(be_values)




    def connect_be(ft,p_connect_be):
        """
            Comment
        """
        if np.random.rand() < p_connect_be:
            if not list(set(list(ft.basic_events.values())) - set(ft.get_all_bes())):
                return ft
            new_ft         = deepcopy(ft)
            be_to_connect = random.choice(list(set(list(new_ft.basic_events.values())) - set(new_ft.get_all_bes())))
            gate          = random.choice(new_ft.get_all_gates())
            gate.children.append(be_to_connect)
            return new_ft
        else:
            return ft


    def disconnect_be(ft,p_disconnect_be):
        """
            Comment
        """
        if np.random.rand() < p_disconnect_be:
            new_ft       = deepcopy(ft)
            be_to_del    = random.choice(new_ft.get_all_bes())
            parent_of_be = new_ft.get_parent_of_be(be_to_del)
            parent_of_be.children.remove(be_to_del)
            return new_ft
        else:
            return ft


    def swap_be(ft,p_swap_be):
        """
            Comment
        """
        if np.random.rand() < p_swap_be:
            new_ft       = deepcopy(ft)
            be_to_swap   = random.choice(new_ft.get_all_bes())
            parent_of_be = new_ft.get_parent_of_be(be_to_swap)
            try:
                new_parent   = random.choice(list(set(new_ft.get_all_gates()) - set([parent_of_be])))
            except IndexError:
                return ft
            parent_of_be.children.remove(be_to_swap)
            new_parent.children.append(be_to_swap)
            return new_ft
        else:
            return ft




    def create_gate(ft,p_create_gate):
        """
            Comment
        """
        if np.random.rand() < p_create_gate:
        
            new_ft = deepcopy(ft)
            
            try: 
                gate = random.choice(new_ft.get_all_gates())
                selected_children = random.sample(gate.children, random.randint(1, len(gate.children)-1))
            except ValueError:
                return new_ft
                
            for child in selected_children:
                gate.children.remove(child)
            
            if random.choice([True, False]):
                gate.children.append(AND(selected_children))
            else:
                gate.children.append(OR(selected_children))

            return new_ft
        else:
            return ft    




    def change_gate(ft,p_change_gate):
        """
            Comment
        """
        if np.random.rand() < p_change_gate:
            new_ft = deepcopy(ft)
            gate_to_change = random.choice(new_ft.get_all_gates())
            children       = gate_to_change.children
            if gate_to_change == new_ft.top_event:
                if isinstance(gate_to_change,AND):
                    new_ft.top_event = OR(children)
                else:
                    new_ft.top_event = AND(children)
            else:
                parent_of_gate = new_ft.get_parent_of_be(gate_to_change)
                parent_of_gate.children.remove(gate_to_change)
                if isinstance(gate_to_change,AND):
                    parent_of_gate.children.append(OR(children))
                else:
                    parent_of_gate.children.append(AND(children))
            return new_ft
        else:
            return ft



    def delete_gate(ft,p_delete_gate):
        """
            Comment
        """
        if np.random.rand() < p_delete_gate:
            new_ft         = deepcopy(ft)
            try:
                gate_to_del    = random.choice(list(set(new_ft.get_all_gates()) - set([new_ft.top_event])))
            except IndexError:
                return ft
            children       = gate_to_del.children
            parent_of_gate = new_ft.get_parent_of_be(gate_to_del)
            parent_of_gate.children.remove(gate_to_del)
            parent_of_gate.children.extend(children)
            return new_ft
        else:
            return ft



    def cross_over(ft1,ft2,p_cross_over):
        """
            Comment
        """
        if np.random.rand() < p_cross_over:
            new_ft1      = deepcopy(ft1)
            new_ft2      = deepcopy(ft2)
            gate_of_1    = random.choice(new_ft1.get_all_gates())
            gate_of_2    = random.choice(new_ft2.get_all_gates())
            if gate_of_1 == new_ft1.top_event:
                new_ft1.top_event = gate_of_2
            else:
                parent_of_g1 = new_ft1.get_parent_of_be(gate_of_1)
                parent_of_g1.children.remove(gate_of_1)
                parent_of_g1.children.append(gate_of_2)
            if gate_of_2 == new_ft2.top_event:
                new_ft2.top_event = gate_of_1
            else:
                parent_of_g2 = new_ft2.get_parent_of_be(gate_of_2)
                parent_of_g2.children.remove(gate_of_2)
                parent_of_g2.children.append(gate_of_1)
            return new_ft1
        else:
            return ft1



    def generate_initial_population(be_s):
        """
            Comment
        """    
        basic_events = {}
        for be in be_s:
            basic_events[be] = BE(be)
            
        ft1_flat = ["OR",be_s]
        ft1       = FT(flattened_representation=ft1_flat, top_event=None, basic_events=basic_events)
        
        ft2_flat = ["AND",be_s]
        ft2       = FT(flattened_representation=ft2_flat, top_event=None, basic_events=basic_events)
        
        return [ft1,ft2]



    def fitness(self,dataset):
        count_true = 0
        for data in dataset:
            if data['T'] == self.valuate(data):
                count_true += 1
        return count_true/len(dataset)
        

    def phi_d(self,dataset):
        count_true = 0
        total_counts = 0
        for data in dataset:
            if data['T'] == self.valuate(data):
                count_true += data['N']
            total_counts += data['N']
        return count_true/total_counts
        
# Fitness function:
    def phi_c(self,c2):
        
        #if not isinstance(ft, str):
        #    ft = cutsets2ft(ft)
        
        ft1 = str(self)
        #ft2 = ft
        
        # We obtain the cut sets of both trees:
        _,_,c1 = cutsets_from_ft_string(ft1)
        #_,_,c2 = cutsets_from_ft_string(ft2)
        
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
        c_cos_coef = np.trace(np.matmul(np.matmul(c1.T,c1).T,np.matmul(c2.T,c2)))/np.sqrt(np.trace(np.matmul(np.matmul(c1.T,c1).T,np.matmul(c1.T,c1)))*np.trace(np.matmul(np.matmul(c2.T,c2).T,np.matmul(c2.T,c2))))
        return c_cos_coef
    
    def phi_s(self):
        return len(self.get_all_bes()) + len(get_gate_statistics(self))

def read_dataset_counts(table):
    """
    My method docstring
    """
    # Identify the number of leaf nodes.
    be_s = []
    for i in range(table.shape[1]-2):
        be_s.append('BE'+str(i+1))
    
    # Generates the dataset:
    dataset = []
    for i, be in enumerate(table):
        be_values = {}
        for j in range(be.shape[0]-2):
            be_values[be_s[j]] = bool(be[j])
        be_values['T'] = bool(be[-2])
        be_values['N'] = be[-1]
        dataset.append(be_values)
    return be_s, dataset


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


"""
    Returns the conjunctive normal form and disjunctive normal form of an input FT
"""
def simplify_ft(ft):
    boolean_expression = str(ft)
    boolean_expression = boolean_expression.replace('OR()','').replace('AND()','').replace('OR','Or').replace('AND','And').replace(',)',')')
    # return simplify_logic(expr=boolean_expression,form='cnf'), simplify_logic(expr=boolean_expression,form='dnf')
    return simplify_logic(expr=boolean_expression,form='cnf')


#%% 2-D Correlation 
def mean2(x):
    y = np.sum(x) / np.size(x);
    return y

def corr2(a,b):
    a = a - mean2(a)
    b = b - mean2(b)

    r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum());
    return r
# %% Trimming FT function:
    
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
                        ft_new = FT(flattened_representation=ft2_flat, top_event=None, basic_events=temporary_basic_event)
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
                        idx1 = get_index_positions(child_name,child_name[0])
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
def get_ft_structure_paths(ft1):
    
    op = [m.start() for m in re.finditer('\(', ft1)]
    cp = [m.start() for m in re.finditer('\)', ft1)]
    cm = np.array([m.start() for m in re.finditer(',', ft1)])
    
    m = []
    for i in op:
        m.append(np.array(cp)-i)
    m = np.array(m)
    
    k = []
    cp_ , op_ = deepcopy(cp) , deepcopy(op)
    while len(m)>0:
        idx = np.where(m == np.min(m[m>0]))
        k.append( [ op_[idx[0][0]] , cp_[idx[1][0]] ] )
        m = np.delete(np.delete(m, idx[0][0], 0), idx[1][0], 1)
        cp_.pop(idx[1][0]), op_.pop(idx[0][0])
    del(cp_,op_)
    
    k = np.array(k)
    m = []
    for i in range(0,len(k)):
        m.append(np.logical_and( list(k[i,0]<k[:,0]) , list(k[i,1]>k[:,1]) ))
    m.reverse()    
    
    
    ft_str = {}
    cont_and,cont_or=0,0
    for i in k:
        if any([any(i[0]>cm) , any(i[0]>op)]):
            # Looks for the closest comma:
            if len(cm[i[0]>cm])>0:
                a = [max(cm[i[0]>cm])]
            else:
                a = []
            # Looks for the closest parenthesis:
            b = [max(np.array(op)[i[0]>op])]
            c = max(a+b)
        else:
            c = -1
        if ft1[c+1:i[0]] == 'AND':
            cont_and += 1
            ft_str[ft1[c+1:i[1]+1]] = 'AND_'+str(cont_and)
        elif ft1[c+1:i[0]] == 'OR':
            cont_or += 1
            ft_str[ft1[c+1:i[1]+1]] = 'OR_'+str(cont_or)
    
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
                ft_str[j.replace(i,ft_str[i])] = ft_str[j]
                del ft_str[old_j]
                
            if len(keys_all)>0:
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
                    be_s.append(int(j[2:]))
        
        # Look for the children gates:
        gates = []
        if len(temp1)>0:
            for j in temp1:
                if j.find("BE") == -1:
                # There is at least one gate:
                    gates.append(j) 
    
        table_structure.append([ i[0][0:i[0].index('_')] , i[0] , gates , be_s])

    for i in table_structure:
        
        if i[1][0:i[1].index('_')] == 'AND':
            i[1] = "AND_" + str( - int(i[1][i[1].index('_')+1:]) + cont_and + 1)
        elif i[1][0:i[1].index('_')] == 'OR':
            i[1] = "OR_" + str( - int(i[1][i[1].index('_')+1:]) + cont_or + 1)
        
        if len(i[2])>0:
            new_var = []
            for j in i[2]:
                if j[0:j.index('_')] == 'AND':
                    new_var.append( "AND_" + str( - int(j[j.index('_')+1:]) + cont_and + 1) )
                elif j[0:j.index('_')] == 'OR':
                    new_var.append( "OR_" + str( - int(j[j.index('_')+1:]) + cont_or + 1) )
            i[2] = new_var
    
    
    list_table_str = []
    for i in table_structure:
        for j in i[2]:
            list_table_str.append([i[1] , j ])
            
        for j in i[3]:
            list_table_str.append([i[1] , "BE" + str(j) ])
    
    paths = []
    for j in list_table_str:
        # Check that the first node is a root node:
        if any(np.array(list_table_str, dtype ='object')[:,1] == j[0]) == False:
            master_var = [j]
            cont = 1
            while cont == 1:
                cont = 0
                new_master_var = []
                for k in master_var:
                        indx = [i for i, x in enumerate(np.array(list_table_str, dtype ='object')[:,0].tolist()) if x == k[-1]]
                        if len(indx)>0:
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
                if  j[0:j.index('_')] == 'AND':
                    s2 = s2 + 'AND' + '-'
                elif j[0:j.index('_')] == 'OR':
                    s2 = s2 + 'OR' + '-'
            else:
                s2 = s2 + j + '-'
            
        path_str.append(s2[1:-1])
        path_str_num.append(s[1:-1])
        
    return path_str,path_str_num, table_structure
# %%
def del_list_indexes(l, id_to_del):
    somelist = [i for j, i in enumerate(l) if j not in id_to_del]
    return somelist
# %%
def find_int(l, target):
    sub_list = [x for x in range(len(l)) if target in l[x]][0]
    return sub_list, l[sub_list].index(target)
# %%
def get_ft_feat(ft):
    
    _,feat_cutsets,_ = cutsets_from_ft_string(str(ft))
        
    # Get the number of OR and And gates
    gates_ = ft.get_all_gates()
    cont_or = 0
    cont_and = 0
    for i in gates_:
        if type(i).__name__ == 'OR':
            cont_or += 1 
        elif type(i).__name__ == 'AND':
            cont_and += 1 

    _,_, table_structure = get_ft_structure_paths(str(ft))
    
    
    all_bes = [ i.name for i in ft.get_all_bes()]
    used = set()
    all_bes = [x for x in all_bes if x not in used and (used.add(x) or True)]
    idx_bes = []
    for i in all_bes:
        idx_bes.append(int(i[2:]))

    #Variables | #UniqueBEs | #TotalBEs | #BEs | #OR | #AND | Cutset length | Cutset order | Tree depth | #Branches | #Repeated events | #Superfluous variables
    features = [max(idx_bes), len(idx_bes) , len(ft.get_all_bes()), cont_or, cont_and, len(feat_cutsets), sum(np.array(feat_cutsets)), get_depth(ft) , len(table_structure[0][2]) , abs(len(ft.get_all_bes())-len(idx_bes)) , max(idx_bes) - len(idx_bes)]
    
    return features
#%%
def compare(a, b):
    """
    Compare two base strings, disregarding whitespace
    """
    return re.sub("\s*", "", a) == re.sub("\s*", "", b)
# %% Create FT object
def str2ft(ft1):
    # Example:
    #ft1 = 'OR(BE1,AND(BE3,OR(BE4,AND(BE2,BE5))),BE2,OR(BE1,BE2))'
    _,b,c= get_ft_structure_paths(ft1)
    
    be = []
    for i in c:
        be = be + i[3]
    used = set()
    be = [x for x in be if x not in used and (used.add(x) or True)]
    
    # Create the basic events:
    be_s = []   
    for i in be:
        be_s.append('BE'+str(i))
    
    basic_events = {}
    for be in be_s:
        basic_events[be] = BE(be)
    #FT(flattened_representation=ft1_flat, top_event=None, basic_events=basic_events)
    ft = FT(flattened_representation=[c[0][0],be_s], top_event=None, basic_events=basic_events)
    
    # Delete all the elements below the root event (i.e., top event):
    gate = ft.get_all_gates()[0]
    selected_children = gate.children
    while len(selected_children)>0:
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
            dict_gates[i[1]].children.append(BE('BE'+str(j)))

    return ft
#%%
def flatten(lis):
     flat_list = [item for sublist in lis for item in sublist]
     return flat_list
#%%
def table_to_input_ft(data):
    """
    My method docstring
    """
    # Identify the number of leaf nodes.
    be_s = []
    if len(data.shape)>1:
        for i in range(data.shape[1]):
            be_s.append('BE'+str(i+1))
    else:
        be_s.append('BE1')
    
    if len(data.shape)>1:
        # Generates the dataset:
        dataset = []
        for i, be in enumerate(data):
            be_values = {}
            for j in range(be.shape[0]):
                be_values[be_s[j]] = bool(be[j])
            #be_values['T'] = bool(be[-2])
            #be_values['N'] = be[-1]
            dataset.append(be_values)
    else:
        # Generates the dataset:
        dataset = []
        for i, be in enumerate(data):
            be_values = {}
            for j in range(1):
                be_values[be_s[j]] = bool(be)
            #be_values['T'] = bool(be[-2])
            #be_values['N'] = be[-1]
            dataset.append(be_values)
    return be_s, dataset
# %%
def get_BEs_per_gate(ft):
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
        
        N_BE = []
        for child in gate.children:
            if isinstance(child, BE):
                N_BE.append(child.name)
        
        result.append([gate_name, N_BE])
    
    return result

def get_ID_All_BEs_per_ft(a):
    mylist = []
    for i in a:
        mylist.extend(i[1])
    used = set()
    unique = [x for x in mylist if x not in used and (used.add(x) or True)]
    
    k = []
    for i in unique:
        k.append(int(i[2:]))
    return k

def parse_ft_string_depth(string):
    depth = 0
    max_depth = 0
    for character in string:
        if character == '(':
            depth += 1
            max_depth = max(max_depth, depth)
        elif character == ')':
            depth -= 1
    return max_depth

def get_depth(ft):
    return parse_ft_string_depth(str(ft))

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
# %%
def save_results(raw_fts,t,path_save_results,dataset,ft_from_MCSs,multi_objective_function):

    str_fts_sorted = []
    str_fts_sorted_trimmed = []
    fts_sorted_trimmed = []
    for i in raw_fts[0]:
        str_fts_sorted.append(str(i))
        fts_sorted_trimmed.append(trim_ft(i))
        str_fts_sorted_trimmed.append(str(fts_sorted_trimmed[-1]))
        
    # Compute the metrics from the FTs
    fitnesses_trimmed,_ = compute_metrics_fts(fts_sorted_trimmed,dataset,ft_from_MCSs,multi_objective_function)

    files_in_folder = []
    for name in os.listdir(path_save_results+ '/' ):
        if os.path.isfile(os.path.join(path_save_results+ '/', name)):
            files_in_folder.append(name)
            
    num_files_dir = int(len(files_in_folder))
    saving_dir = path_save_results + '/' + 'gen_' + str(num_files_dir)
    mdic = {'metrics_fts': raw_fts[1],'metrics_fts_trimmed':fitnesses_trimmed,'fts':str_fts_sorted,'fts_trimmed':str_fts_sorted_trimmed,'time':t}
    savemat(saving_dir+'.mat', mdic)

#%% Compute metrics
def compute_metrics_fts(initial_population,dataset,ft_from_MCSs,multi_objective_function):
    
    #show us initial population
    fitnesses = []
    fitness_dict = {}

    #Compute accuracy based on the MCS:
    if multi_objective_function[0] != 0:
        acc_mcs = []
        for ft in initial_population:
            acc_mcs.append([1.0-ft.phi_c(ft_from_MCSs)])
        acc_mcs = np.array(acc_mcs)
    else:
        acc_mcs = np.ones((len(initial_population),2))*-1
        
    #Compute the size:
    if multi_objective_function[1] != 0 or multi_objective_function[1] == 0:
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
            acc_data.append([1.0-ft.phi_d(dataset)])
        acc_data = np.array(acc_data)
    else:
        acc_data = np.ones((len(initial_population),2))*-1
   
    fitnesses = np.column_stack((acc_mcs[:,0],tree_size[:,0],acc_data[:,0]))
    for ft,fi1 in zip(initial_population,fitnesses):
        fitness_dict[str(ft)]  = fi1
    
    return fitnesses, fitness_dict

# %% Fitness function
def cost_function(initial_population,dataset,population_size,ft_from_MCSs,multi_objective_function):

    # Compute the metrics from the FTs
    fitnesses,fitness_dict = compute_metrics_fts(initial_population,dataset,ft_from_MCSs,multi_objective_function)

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
    inputPoints = adjust_array(list(c))
    pt,pt_idx = nsgaII(inputPoints,population_size)
    f=[]
    for i in pt:
        if len(i)>1:
            nds = [fitnesses[i,0] + fitnesses[i,2]][0].argsort()
            f.append(np.array(i)[nds].tolist())
        else:
            f.append(i)
    # Update the order leaving the best accuracy always on top
    pt = f
    pt_idx = list(itertools.chain.from_iterable(f))
    raw_fts = [np.flip(np.array(initial_population)[pt_idx],0).tolist(),np.flip(fitnesses[pt_idx],0),fitness_dict]
    # ------------------------------------------------

    return raw_fts
# %% Execute the main program:
def nsgaII(data,pt_size):
    # Identify all the fronts:
    all_fronts = []
    data_copy = data[:]
    idx_data = list(range(0,len(data)))
    while data:
        _, _, idx = simple_cull(inputPoints = data.copy()) # Identify the pareto front
        c = [e for i, e in enumerate(idx_data) if i in idx]
        c.reverse() # In this way, we give more priority to those FTs with higher accuracy.
        all_fronts.append(c)
        #all_fronts.append(idx[:])
        idx_data = [j for i, j in enumerate(idx_data) if i not in idx]
        data     = [j for i, j in enumerate(data) if i not in idx]
        
    size_per_front = []
    for i in all_fronts:
        size_per_front.append(len(i))
    
    cumsum_size_per_front = np.cumsum(np.array(size_per_front))
    
    final_pt_idx = []# These are the fronts that will pass to the next generation.
    pt_out =[]
    cont = 0
    for i in cumsum_size_per_front<=pt_size:
        if i:
            final_pt_idx.append(cont)
            pt_out.append(all_fronts[cont])
        cont += 1
    
    # Check whether we need to include part of the next front:
    fronts_left_out = [i for i, x in enumerate(cumsum_size_per_front>=pt_size) if x]
    if cumsum_size_per_front[fronts_left_out[0]] > pt_size:
        # Crowding Distance Sorting.
        cws_data = []
        cws_dx = []
        for i in all_fronts[fronts_left_out[0]]: 
            cws_data.append(data_copy[i])
            cws_dx.append(i)
        a = calc_crowding_distance(np.array(cws_data))
        a  = np.flip(a.argsort())
        a  = a[: - cumsum_size_per_front[fronts_left_out[0]] + pt_size ]
        b = []
        for i in a:
            b.append(cws_dx[i])
        pt_out.append(b)
        
    pt_out_index = np.array([])
    for i in pt_out:
        pt_out_index = np.concatenate((pt_out_index,i))
        
    return pt_out, pt_out_index.astype(int)
#%% Pareto front
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

def adjust_array(a): 
    b = []
    for i in a:
        b.append([float(j) for j in i])
    return b
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

#%%
def gen_initial_population_from_cutsets(ft_from_MCSs):
    initial_population = []
    _,_,c = cutsets_from_ft_string(ft_from_MCSs)
    for i in c.tolist():
        fts = 'AND('
        for j in range(0,len(i)):
            if i[j] == 1:
                fts = fts + 'BE' + str(j+1) + ','
        fts = fts[0:-1] + ')' 
        initial_population.append(str2ft(fts))
      
    return initial_population
#%% 
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

#%%
def check_empty_objects(ft):
    if str(ft).find("()") != -1:
        idx = True
    else:
        idx = False
    return idx
#%% Genetic operators:
def genetic_operators(population):
    new_population = []
    go = ['disconnect_be','connect_be','swap_be','create_gate','delete_gate','cross_over','change_gate']
    for ft in population:
        try:
            if ft not in new_population:
                new_population.append(ft)
            
            for i in go:
                if i != 'cross_over':
                    new_ft = eval('FT.'+i+'(ft,'+str(1.0)+')')
                else:
                    new_ft = eval('FT.'+i+'(ft,random.choice(population),'+str(1.0)+')')
                if check_empty_objects(new_ft) == False:
                    if new_ft not in new_population:
                        new_population.append(new_ft)
        except:
            continue
    return new_population
#%%
def cutsets2ft(cut_sets):
    fts = 'OR('
    for i in cut_sets.tolist():
        j = 'AND('
        cont = 0
        for k in i:
            cont += 1
            if k == 1:
                j = j + 'BE' + str(cont) + ','
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
            all_comb_failure_temp = all_comb_failure[:,np.where(ci == True)[0]]
            to_delete = np.where(sum(all_comb_failure_temp.T)  == all_comb_failure_temp.shape[1])
            all_comb_failure = np.delete(all_comb_failure,to_delete,0)
            if len(to_delete)>0:
                cut_sets.append(ci) 
        cut_sets = np.array(cut_sets).astype(int)
    else:
        cut_sets = all_comb_failure.reshape(1,len(all_comb_failure))
    
    return cut_sets
#%%
def check_empty_gates(population):
    c = 0
    idx = []
    for ft in population:
        c += 1
        if str(ft).find("()") != -1:
            idx.append(c)
    idx = (np.array(idx)-1).tolist()
    return idx
#%%
def fts2dcnf(ft_stirng):

    fts_cnf = []
    
    _,_, tab = get_ft_structure_paths(ft_stirng)
    
    cond = {}
    for i in tab:
        str_ = ''
        
        # What is the type of gate?
        if i[0] == 'AND':
            g = '&'
        elif i[0] == 'OR': 
            g = '|'
        
        # Add the basic events in that gate:
        for j in i[-1]:
            str_ = str_ +'BE' + str(j) + ' ' + g + ' '
        
        # Add the other gates within the parent gate.
        for j in i[-2]:
            str_ = str_ + ' ' + j + ' ' + g
        
        cond[i[1]] = str_[0:-2]
    
    tab.reverse()
    
    for i in tab:
        if len(i[2])>0:        
            # Parent string:
            ps = cond[i[1]]
            ps2 = re.split('[ ]',ps)
            # Children string:
            for j in i[2]:
                cs = cond[j]
                #idxs = [i for i, x in enumerate(ps2) if x == j]
                idxs = [k for k, x in enumerate(ps2) if compare(x,j)]
                for k in idxs:
                     ps2[k] = ps2[k].replace(ps2[k], '('+ cs +')')
            ps = ' '.join(ps2)
            cond[i[1]] = ps
                
    output_string = cond[tab[-1][1]]
    cutsets_dnf  = re.split('\|',str(to_dnf(output_string)))
    cutsets_cnf  = re.split('\&',str(to_cnf(output_string)))
    
    # Build the string in the Disjunctive Normal Form:
    fts_dnf = 'OR('
    for i in cutsets_dnf:
        j = re.split('[&()BE]', i)
        be = []
        for k in j:
            try:
                be.append(int(k))
            except:
                continue
        
        aux = 'AND('
        for i in be:
            aux = aux + 'BE'+str(i)+','
        aux = aux[0:-1] + ')'
        fts_dnf = fts_dnf + aux + ','
    
    fts_dnf = fts_dnf[0:-1] + ')'
    
    # Build the string in the Conjunctive Normal Form:
    fts_cnf = 'AND('
    for i in cutsets_cnf:
        j = re.split('|[&()BE]', i)
        be = []
        for k in j:
            try:
                be.append(int(k))
            except:
                continue
        aux = 'OR('
        for i in be:
            aux = aux + 'BE'+str(i)+','
        aux = aux[0:-1] + ')'
        fts_cnf = fts_cnf + aux + ','
    
    fts_cnf = fts_cnf[0:-1] + ')'
    
    return fts_dnf, fts_cnf

#%%
def cutsets_from_ft_string(ft_string):
    _,_, tab = get_ft_structure_paths(ft_string)
    
    cond = {}
    for i in tab:
        str_ = ''
        
        # What is the type of gate?
        if i[0] == 'AND':
            g = '&'
        elif i[0] == 'OR': 
            g = '|'
        
        # Add the basic events in that gate:
        for j in i[-1]:
            str_ = str_ +'BE' + str(j) + ' ' + g + ' '
        
        # Add the other gates within the parent gate.
        for j in i[-2]:
            str_ = str_ + ' ' + j + ' ' + g
        
        cond[i[1]] = str_[0:-2]
    
    tab.reverse()
    
    for i in tab:
        if len(i[2])>0:        
            # Parent string:
            ps = cond[i[1]]
            ps2 = re.split('[ ]',ps)
            # Children string:
            for j in i[2]:
                cs = cond[j]
                idxs = [k for k, x in enumerate(ps2) if compare(x,j)]
                for k in idxs:
                      ps2[k] = ps2[k].replace(ps2[k], '('+ cs +')')
            ps = ' '.join(ps2)
            cond[i[1]] = ps
                
    output_string = cond[tab[-1][1]]
    cutsets  = to_dnf(output_string)
    
    # Compute the length of each cutset
    a = re.split('\|', str(cutsets))
    b = []
    for i in a:
        b.append( len(re.split('\&',i)) )
    
    
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
        a[np.array(i)-1] = 1
        cutset_matrix.append(a.astype(int))
    
    if len(cutset_matrix) == 1:
        cutset_matrix = cutset_matrix[0]
    
    cutset_matrix = getMCSs(np.array(cutset_matrix))
    cutset_length

    return cs2str(cutset_matrix), np.sum(cutset_matrix,axis=1), cutset_matrix

#%%
def cs2str(A):

    all_IE = []
    for i in A:
        j = (np.where(np.array(i) == 1)[0]+1).astype(int).tolist()
        IE = '('
        for k in j:
            IE  += 'BE' + str(k) + ' & '
        IE = IE[0:-2] + ')'
        all_IE.append(IE)   
    str_ = ''
    for i in all_IE:
        str_ +=  i + ' | '
    str_ = str_[0:-2]
    # str_cutsets  = to_dnf(str_)
    
    return str_

#%% FTMOEA
class FTMOEA():
    
    """
        Class for learning a fault tree from data using a genetic algortihm
    """        

    def __new__(self, population_size=100, MCSs=[],parentFT='', ng=100, p_connect_be = 0.1, p_disconnect_be = 0.1, p_swap_be = 0.2, p_create_gate = 0.2, p_change_gate = 0.2, p_delete_gate = 0.2, p_cross_over = 0.2, convergence_criterion=10,multi_objective_function = [1,1,1], file_dataset=None, selection_strategy='elitist',dataset = [], be_s = [],path_save_results = ''):
        """
            Input parameters:
                * Population_size (ps): the minimal population size maintained throughout generations after applying the genetic operations
                * Minimal Cut Sets (MCSs): this is optional, only applicable on noise-free and complete failure data sets. Based on this metric the error based on Minimal Cut Sets (\phi_c) is computed.
                * Parent Fault Tree (parentFT): this is optional, entry a Fault Tree string to be used as parent Fault Tree. If empty, it is going to be initialized with two-parent Fault Trees, one with all BEs connected to an AND gate, and the other with all BEs connected to an OR gate. 
                * Max. number of generations (ng): Terminates the optimization process if the number of generations exceeds $ng$ and none of the other convergence criteria is met. 
                * Max. generations with unchanged best candidate (ug): if after $uc$ number of generations the best individual (i.e., the FT with the smallest size, and smallest error(s) within the best Pareto set) remains unchanged, then we assume the process has converged and is therefore terminated.
                * multi_objective_function: It is a list with three terms that switches on (1) and off (0) the metrics to be optimized.
                    - The first term corresponds to the error based on Minimal Cut Sets (\phi_c)
                    - The second term corresponds to the size of the Fault Tree (\phi_s)
                    - The third term corresponds to the error based on the failure data set (\phi_d)
                * Failure data set (dataset): corresponds to the failure data set.
                * Basic events (be_s): associated with the basic events considered in the inference process.
                * Path to save results externally (path_save_results): when the path of a folder is provided, the results per generation are saved in a .mat file.
            Output:
                * Inferred Fault Trees (fts): this list contains all the objects of inferred Fault Trees of the last generation, sorted from the worst to the best Fault Tree.
                * Convergence time per generation (t): computational time per generation.
                * Metrics: is an array of three columns, corresponding respectively to \phi_c, \phi_s, and \phi_d. And the rows are respectively associated with the inferred Fault Trees contained in (fts).

        """ 
        
        print('... FT-MOEA initialized ...')
        
        if len(MCSs)==0:
            multi_objective_function[0] = 0
        
        if path_save_results != '':
            # Check if the saving folder exist, otherwise create it.
            path_save_results = upfold(path_save_results)

        selection_strategy = 'elitist'
        multi_objective_function = list(-1*np.array(multi_objective_function)) # Minimize the size of the FT
        if MCSs != []:
            ft_from_MCSs = MCSs#cutsets2ft(MCSs)
        t = [time.time()]
        dict_iterations = []

        def index_of(a,list):
            for i in range(0,len(list)):
                if list[i] == a:
                    return i
            return -1


        def sort_by_values(list1, values):
            sorted_list = []
            while(len(sorted_list)!=len(list1)):
                if index_of(min(values),values) in list1:
                    sorted_list.append(index_of(min(values),values))
                values[index_of(min(values),values)] = math.inf
            return sorted_list
        
        
        # ------------------------------------
        # Create the parent fault trees:
        if parentFT == '':
            # Linard's (2019)
            initial_population = FT.generate_initial_population(be_s)
        else:
            # Choosen parent fault tree
            initial_population = [str2ft(parentFT)]

        # ------------------------------------
        while len(initial_population) < population_size:
            initial_population = genetic_operators(initial_population)

        t.append(time.time())
        raw_fts = cost_function(initial_population,dataset,population_size,ft_from_MCSs,multi_objective_function)
        
        if path_save_results != '':
            save_results(raw_fts,t[-1]-t[0],path_save_results,dataset,ft_from_MCSs,multi_objective_function)

        print('Gen. \t Fitness Pop. \t Fitness best \t best individual')
        print('0','\t    phi_d=',"{:.4f}".format(np.mean(raw_fts[1][:,0])),' - phi_c=',"{:.4f}".format(np.mean(raw_fts[1][:,2])),'\t    phi_d=',"{:.4f}".format(raw_fts[1][-1,0]),' - phi_c=',"{:.4f}".format(raw_fts[1][-1,2]),'/ phi_s=',"{:.2f}".format(raw_fts[1][-1,1]),'\t',raw_fts[0][-1])
        
        dict_iterations.append([str(raw_fts[0][-1])] + np.mean(raw_fts[1],axis=0).tolist() + raw_fts[1][-1].tolist()  )
        population = raw_fts[0]
        conv       = 0
        
        for i in range(1,ng):
            
            t.append(time.time())
            st = 0
            
            while st == 0:
                new_population = []
                st = 1
                while len(new_population) < population_size:
                    new_population = genetic_operators(population)
                raw_fts = cost_function(new_population,dataset,population_size,ft_from_MCSs,multi_objective_function)

            if path_save_results != '':
                #Saving dataset
                save_results(raw_fts,t[-1]-t[0],path_save_results,dataset,ft_from_MCSs,multi_objective_function)
            
            dict_iterations.append([str(raw_fts[0][-1])] + np.mean(raw_fts[1],axis=0).tolist() + raw_fts[1][-1].tolist()  )
            print(str(i),'\t    phi_d=',"{:.4f}".format(np.mean(raw_fts[1][:,0])),' - phi_c=',"{:.4f}".format(np.mean(raw_fts[1][:,2])),'\t    phi_d=',"{:.4f}".format(raw_fts[1][-1,0]),' - phi_c=',"{:.4f}".format(raw_fts[1][-1,2]),'/ phi_s=',"{:.2f}".format(raw_fts[1][-1,1]),'\t',raw_fts[0][-1])

            # ------------------------------
            # Convergence criteria:
            if multi_objective_function == [-1,-1,-1]:
                
                if dict_iterations[-2][4] == dict_iterations[-1][4] and dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][6] == dict_iterations[-2][6]:
                    conv += 1
                else:
                    conv = 0
    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1: #or ( dict_iterations[-1][4] == 1.0 and dict_iterations[-1][6] == 1.0 ):
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]
            
            elif multi_objective_function == [-1,0,-1]:
                
                if dict_iterations[-2][4] == dict_iterations[-1][4] and dict_iterations[-1][6] == dict_iterations[-2][6]:
                    conv += 1
                else:
                    conv = 0
    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1 or ( dict_iterations[-1][4] == 0.0 and dict_iterations[-1][6] == 0.0 ):
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]
                
            elif multi_objective_function == [0,0,-1]: 
            
                if dict_iterations[-1][6] == dict_iterations[-2][6]:
                    conv += 1
                else:
                    conv = 0
                    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1 or np.max(raw_fts[1][:,2]) == 0.0:
                    if np.max(raw_fts[1][:,2]) == 0.0:
                        print('Max. Accuracy achieved.')
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]
            
            elif multi_objective_function == [0,-1,0]: 
            
                if dict_iterations[-1][5] == dict_iterations[-2][5]:
                    conv += 1
                else:
                    conv = 0
                    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1:
                    if np.max(raw_fts[1][:,2]) == 0.0:
                        print('Max. Accuracy achieved.')
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]
            
            
            elif multi_objective_function == [0,-1,-1]: 
            
                if dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][6] == dict_iterations[-2][6]:
                    conv += 1
                else:
                    conv = 0
                    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1:
                    if np.max(raw_fts[1][:,2]) == 0.0:
                        print('Max. Accuracy achieved.')
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]
            
            elif multi_objective_function == [-1,-1,0]: 

                if dict_iterations[-1][5] == dict_iterations[-2][5] and dict_iterations[-1][4] == dict_iterations[-2][4]:
                    conv += 1
                else:
                    conv = 0
                    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1:
                    if np.max(raw_fts[1][:,2]) == 0.0:
                    #if np.max(trimmed_fts[1][:,2]) == 1.0:
                        print('Max. Accuracy achieved.')
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]


            elif multi_objective_function == [-1,0,0]: 
                
                if dict_iterations[-1][4] == dict_iterations[-2][4]:
                    conv += 1
                else:
                    conv = 0
                    
                #return best FT is convergence or best fitness reached
                if conv >= convergence_criterion-1 or np.max(raw_fts[1][:,0]) == 0.0:
                    if np.max(raw_fts[1][:,0]) == 0.0:
                        print('Max. Accuracy achieved.')
                    print('... FT-MOEA finalized ...')
                    return raw_fts[0], t, raw_fts[1]


        
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
                            for k in range(0,j):
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
                    subpop = random.sample(list(popul),tournament_size)
                    for f in subpop:
                        if best ==None:
                            best = f
                            continue
                        if fitness_dict[str(f)] > fitness_dict[str(best)]:
                            best = f
                    
                    chosen.append(best)
                
                return chosen
                
            
            def random_select(popul, population_size, fitness_dict):
                return random.sample(list(popul),population_size)
                
                
            def elitism(popul, population_size, fitness_dict):
                return popul[-population_size:]
                

            
            #SELECTION
            sortedPeople,fitness_dict = raw_fts[0],raw_fts[2]
            
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
