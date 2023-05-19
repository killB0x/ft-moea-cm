import sympy
import logging
from sortedcontainers import SortedSet

from ft_learn.modules.modules import Modules
from ft_learn.ft.mcs import CutSet, MinCutSets
from ft_learn.logic.boolean_logic import mcs_to_sympy_formula


def create_from_mcss(mcss, bes, try_and=False):
    """
    Initialize modules from MCSs.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}.
    :param try_and: Whether to try to find modules for a possible AND-gate when the search for modules under an OR-gate was unsuccessful.
    :return: Tuple (Modules, whether the gate over all modules is an OR-gate)
    """
    # Start by trying to find modules with respect to OR-gate
    modules, found_or = find_modules_or(mcss, bes)
    if found_or:
        return modules, True
    if try_and:
        # Try to find modules with respect to AND-gate
        modules, found_and = find_modules_and(mcss, bes)
        return modules, False
    return modules, True


def find_modules_or(mcss, bes):
    """
    Find modules under a possible OR-gate.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}
    :return: Tuple (Modules, whether multiple modules were found)
    """
    modules = Modules()
    all_bes = CutSet(bes.keys())
    for mcs in mcss:
        modules.add_and_merge(mcs)

        # Abort if module contains all BEs
        if len(modules) == 1 and all_bes.issubset(next(iter(modules))):
            return modules, False
    return modules, True


def find_modules_and(mcss, bes):
    """
    Find modules under a possible AND-gate.
    Note that these modules do not yield disjoint cut sets.
    Also note that the performance of this operation can be exponential in the number of BEs as conversion from DNF to CNF takes place.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}
    :return: Tuple (Modules, whether multiple modules were found)
    """
    # Convert to Boolean formula
    formula = mcs_to_sympy_formula(mcss, bes)
    # Convert to Conjunctive Normal Form (CNF)
    formula = sympy.to_cnf(formula, simplify=True)
    logging.debug("Simplified CNF formula from sympy: {}".format(formula))

    # Find modules in clauses of CNF formula
    all_bes = CutSet(bes.keys())
    be_indices = SortedSet({be: index for index, be in bes.items()})
    modules = Modules()

    # Iterate over formula
    if formula.func == sympy.Or or formula.func == sympy.Symbol:
        clauses = [formula]
    else:
        assert formula.func == sympy.And
        clauses = formula.args
    # Iterate over all clauses
    for clause in clauses:
        if clause.func == sympy.Symbol:
            bes_expr = [clause]
        else:
            assert clause.func == sympy.Or
            bes_expr = clause.args
        # Iterate over all BEs in a clause
        bes_mod = []
        for be_expr in bes_expr:
            assert be_expr.func == sympy.Symbol
            bes_mod.append(be_indices[be_expr.name])

        # Add new module
        modules.add_and_merge(CutSet(bes_mod))

        # Abort if module contains all BEs
        if len(modules) == 1 and all_bes.issubset(next(iter(modules))):
            return modules, False
    return modules, True


def find_pseudo_modules_from_symmetry_under_or(mcss, bes, symmetry):
    """
    Tries to split the BEs into two parts (under an OR-gate) such that only shared BEs (symmetric to themselves)
    are occurring in both modules.
    The MCSs is split according to the two modules as well.
    If no split could be found, the symmetric module and symmetric MCS will be empty.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}
    :param symmetry: Symmetry.
    :return: Tuple (Original module, original MCS, symmetric module, symmetric MCS).
    """
    original_module = set()
    symmetric_module = set()
    original_mcss = MinCutSets()
    remaining_mcss = mcss.copy()  # deepcopy() is not necessary as the cut sets are not changed in the following

    # Initialize with first BE
    remaining_bes = SortedSet({SortedSet(bes.keys()).pop()})
    # Iteratively add BEs which are in the same MCS
    while remaining_bes:
        be = remaining_bes.pop()
        # BE must be contained in original module
        original_module.add(be)
        sym_be = symmetry.apply_be(be)
        # Add symmetric BE to symmetric module
        symmetric_module.add(sym_be)
        if sym_be == be:
            # BE must be shared -> ignore
            continue
        if sym_be in original_module:
            # Abort because original BE and symmetric BE cannot be in same module if OR gate is on top
            return CutSet(bes.keys()), mcss, CutSet([]), MinCutSets()
        else:
            # Include all cut sets containing the original BE in the original MCS
            for cut_set in list(remaining_mcss):
                if be in cut_set:
                    original_mcss.add(cut_set)
                    remaining_mcss.remove(cut_set)
                    remaining_bes = remaining_bes.union(cut_set.set.difference(original_module))

    # Check that both modules contain unique elements
    if original_module.issubset(symmetric_module) or symmetric_module.issubset(original_module):
        return CutSet(bes.keys()), mcss, CutSet([]), MinCutSets()

    # Set symmetric counterparts
    symmetric_mcss = MinCutSets(remaining_mcss)
    return CutSet(original_module), original_mcss, CutSet(symmetric_module), symmetric_mcss


def split_mcss_from_symmetry(mcss, bes, symmetry):
    """
    Tries to split the MCSs into two parts (under an OR-gate) with respect to the symmetry.
    Overlaps are allowed for all BEs.
    If no split could be found, the symmetric module and symmetric MCS will be empty.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}
    :param symmetry: Symmetry.
    :return: Tuple (number of BEs in both modules, original module, original MCS, symmetric module, symmetric MCS).
    """

    def try_adding(cut_set, sym_cut_set, mcss_all, module1, mcss1, module2, mcss2, current_best_split):
        """
        Try to add cut set and symmetric cut to modules.
        :param cut_set: Cut set.
        :param sym_cut_set: Symmetric cut set.
        :param mcss_all: All MCS.
        :param module1: First module.
        :param mcss1: MCS of first module.
        :param module2: Second module.
        :param mcss2: MCS of second module.
        :param current_best_split: Number of elements in both modules corresponding to best split so far (less is better).
        :return: New tuple after insertion (feasible, module1, MCS1, module2, MCS2), None if both modules aleady cannot be better than current_best_spit.
            Return value feasible indicates whether the symmetric cut is in principle feasible. Non-feasibility directly indicates that no split exists.
        """
        new_mcss1 = MinCutSets(mcss1.copy())
        assert cut_set not in new_mcss1
        new_mcss1.add(cut_set)
        new_module1 = module1.union(cut_set)
        if current_best_split <= len(new_module1) + len(module2):
            # Cannot beat best split
            return True, None, None, None, None

        # Add symmetric cut set
        if sym_cut_set not in mcss_all:
            # Symmetric cut set does exist
            logging.debug("Symmetric MCS does not exits".format(sym_cut_set))
            return False, None, None, None, None
        elif sym_cut_set in new_mcss1:
            # Original MCS and symmetric MCS belong to the same module
            # -> not allowed under OR-gate
            logging.debug("Original MCS {} and symmetric MCS {} cannot belong to the same module".format(cut_set, sym_cut_set))
            return False, None, None, None, None
        else:
            # Symmetric cut set is feasible
            assert sym_cut_set not in new_mcss1
            new_mcss2 = MinCutSets(mcss2.copy())
            new_mcss2.add(sym_cut_set)
            new_module2 = module2.union(sym_cut_set)
            if current_best_split <= len(new_module1) + len(new_module2):
                # Cannot beat best split
                return True, None, None, None, None
            # Found possibly better split
            return True, new_module1, new_mcss1, new_module2, new_mcss2,

    def recursive_split(mcss, mcss_list, current_index_mcss, original_module, original_mcss, symmetric_module, symmetric_mcss, current_best_split, symmetry):
        """
        Perform recursive split of MCSs.
        :param mcss: Total MCS.
        :param mcss_list: MCS as list (to iterate over).
        :param current_index_mcss: Current index in the mcss_list to consider next.
        :param original_module: Original module.
        :param symmetric_module: Symmetric module.
        :param original_mcss: MCS for original module.
        :param symmetric_mcss: MCS for symmetric module.
        :param current_best_split: Number of elements in both modules corresponding to best split so far (less is better).
        :param symmetry: Symmetry.
        :return: Tuple corresponding to best split (feasible, number of elements in both modules, module1, MCS1, module2, MCS2), None if no split could be found
        """
        # Remember old best split
        old_best_split = current_best_split

        if len(original_module) + len(symmetric_module) >= current_best_split:
            return True, None, None, None, None, None

        # Find next cut set candidate
        search_next = True
        cut_set = None
        while search_next:
            if current_index_mcss >= len(mcss_list):
                # Considered all cut sets
                return True, len(original_module) + len(symmetric_module), original_module, original_mcss, symmetric_module, symmetric_mcss
            # Get next cut set
            cut_set = mcss_list[current_index_mcss]
            current_index_mcss += 1
            # Check if cut set was already considered
            search_next = cut_set in original_mcss or cut_set in symmetric_mcss

        assert cut_set not in original_mcss
        # Compute symmetric cut set
        sym_cut_set = symmetry.apply_cut_set(cut_set)

        # First try: add cut set to original MCS
        feasible, original_module1, original_mcss1, symmetric_module1, symmetric_mcss1 = try_adding(cut_set, sym_cut_set, mcss, original_module, original_mcss, symmetric_module,
                                                                                                    symmetric_mcss, current_best_split)
        if not feasible:
            return False, None, None, None, None, None

        best_split1 = old_best_split + 1  # Init with worse value
        if original_mcss1:
            feasible, best_split1, original_module1, original_mcss1, symmetric_module1, symmetric_mcss1 = recursive_split(mcss, mcss_list, current_index_mcss, original_module1,
                                                                                                                          original_mcss1, symmetric_module1, symmetric_mcss1,
                                                                                                                          current_best_split, symmetry)
            if not feasible:
                return False, None, None, None, None, None
            current_best_split = best_split1 if best_split1 < current_best_split else current_best_split

        # Second try: add cut set to symmetric MCS
        feasible, symmetric_module2, symmetric_mcss2, original_module2, original_mcss2 = try_adding(cut_set, sym_cut_set, mcss, symmetric_module, symmetric_mcss, original_module,
                                                                                                    original_mcss, current_best_split)
        if not feasible:
            return False, None, None, None, None, None

        best_split2 = old_best_split + 1  # Init with worse value
        if original_mcss2:
            feasible, best_split2, original_module2, original_mcss2, symmetric_module2, symmetric_mcss2 = recursive_split(mcss, mcss_list, current_index_mcss, original_module2,
                                                                                                                          original_mcss2, symmetric_module2, symmetric_mcss2,
                                                                                                                          current_best_split, symmetry)
            if not feasible:
                return False, None, None, None, None, None
            current_best_split = best_split2 if best_split2 < current_best_split else current_best_split

        # Return best split
        if old_best_split <= current_best_split:
            # Old split was best
            return True, old_best_split, original_module, original_mcss, symmetric_module, symmetric_mcss
        elif best_split1 <= best_split2:
            # First split was best
            return True, best_split1, original_module1, original_mcss1, symmetric_module1, symmetric_mcss1
        else:
            assert best_split2 < best_split1
            # Second split was best
            return True, best_split2, original_module2, original_mcss2, symmetric_module2, symmetric_mcss2

    # Init
    mcss_list = [cut_set for cut_set in mcss]
    current_best_split = len(bes) * 2 + 1  # Init with worst possible split + 1
    return recursive_split(mcss, mcss_list, 0, CutSet([]), MinCutSets(), CutSet([]), MinCutSets(), current_best_split, symmetry)
