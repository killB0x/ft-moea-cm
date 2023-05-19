import numpy as np
from sortedcontainers import SortedSet
import logging

from ft_learn.ft.mcs import MinCutSets, CutSet


class Symmetry:
    """
    Single symmetry as a mapping from original BE to symmetric BE.
    Missing BEs in the mapping are assumed to map to themselves.
    """

    def __init__(self, symmetry, bes):
        """
        Constructor.
        :param symmetry: Mapping from BE names to names of symmetric BEs.
        :param bes: Dictionary of BEs {index: name}.
        """
        self.sym = symmetry
        # Create mapping from indices as well
        self.index_sym = dict()
        for index, name in bes.items():
            if name in self.sym:
                sym_name = self.sym[name]
                # Get index corresponding to sym_name
                sym_index = -1
                for i, n in bes.items():
                    if sym_name == n:
                        sym_index = i
                        break
                assert sym_index >= 0
                self.index_sym[index] = sym_index

    def apply_ft(self, fault_tree):
        """
        Apply symmetry on fault tree resulting in symmetric fault tree.
        :param fault_tree: Original fault tree.
        :return: Symmetric fault tree.
        """
        sym_ft = fault_tree.copy()
        for be in sym_ft.get_all_bes():
            if be.name in self.sym:
                be.name = self.sym[be.name]
        return sym_ft

    def apply_be(self, index):
        """
        Apply symmetry on index of BE.
        :param index: Index of BE.
        :return: Index of symmetric BE.
        """
        return self.index_sym.get(index, index)

    def apply_cut_set(self, cut_set):
        """
        Apply symmetry on cut set resulting in symmetric cut set.
        :param cut_set: Original cut set.
        :return: Symmetric cut set.
        """
        return CutSet(self.apply_be(i) for i in cut_set.set)

    def apply_mcss(self, mcss):
        """
        Apply symmetry on each cut set of the MCS yielding the symmetric MCS.
        :param mcss: Original MCS.
        :return: MCS after applying the symmetry.
        """
        sym_mcss = MinCutSets()
        for cut_set in mcss:
            sym_mcss.add(self.apply_cut_set(cut_set))
        return sym_mcss

    def is_valid_symmetry(self, mcss, print_cex=False):
        """
        Check that the symmetry is valid with respect to the given MCS.
        :param mcss: Minimal cut sets (MCS).
        :param print_cex: If True, a counterexample for the symmetry is printed.
        :return: True iff the symmetry is valid.
        """
        original_mcss = mcss.copy()  # deepcopy() is not necessary as the cut sets are not changed in the following
        # Apply the symmetry to each cut set and check that the result corresponds to an existing cut set
        # Performing the check on the fly to enable faster abort
        for cut_set in mcss:
            sym_cut_set = self.apply_cut_set(cut_set)
            if sym_cut_set not in original_mcss:
                if print_cex:
                    logging.error("Symmetric cut set {} is invalid (original cut set {})".format(sym_cut_set, cut_set))
                return False
            else:
                original_mcss.remove(sym_cut_set)

        # Check that no further cut sets are contained
        return len(original_mcss) == 0

    def is_subset(self, other):
        """
        Check whether all symmetries are included in other as well.
        :param other: Other symmetry.
        :return: True iff self is a subset of other.
        """
        for index, sym in self.index_sym.items():
            if index not in other.index_sym:
                return False
            if other.index_sym[index] != sym:
                return False
        return True

    def __len__(self):
        return len(self.sym)

    def __str__(self):
        return ", ".join(["{} -> {}".format(o, s) for o, s in self.sym.items()])


class Symmetries(list):
    """
    List of symmetries.
    """

    def is_valid_symmetries(self, mcss, print_cex):
        """
        Check that all symmetries are valid with respect to the given MCS.
        :param mcss: Minimal cut sets (MCS).
        :param print_cex: If True, a counterexample for the symmetry is printed.
        :return: True iff all symmetries are valid.
        """
        for symmetry in self:
            if not symmetry.is_valid_symmetry(mcss, print_cex):
                return False
        return True

    def __str__(self):
        return "\n".join(["\t" + str(sym) for sym in self])


def get_symmetries_from_file(file, bes):
    """
    Get symmetries from given file.
    The first line contains the original order of the BEs.
    Each subsequent line contains one permutation.
    BE names are separated by a whitespace.
    :param file: Text file containing all permutations encoding the symmetries.
    :param bes: Dict of all BEs {index: name}.
    :return: Symmetries.
    """
    symmetries = Symmetries()
    with open(file, 'r') as f:
        lines = f.readlines()
        names = lines[0].split()
        assert sorted(bes.values()) == sorted(names)
        for line in lines[1:]:
            sym_bes = line.split()
            symmetry = dict(zip(names, sym_bes))
            symmetries.append(Symmetry(symmetry, bes))
    return symmetries


def generate_singleton_symmetries(mcss, bes, be_mcs_occurrences):
    """
    Generate all symmetries between two BEs.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}.
    :param be_mcs_occurrences: Number of times a BE occurs in a cut set.
    :return: All symmetries between BEs.
    """
    symmetries = Symmetries()
    remaining_bes = SortedSet(bes.values())
    while len(remaining_bes) > 0:
        original = remaining_bes.pop()
        covered_bes = []
        for cand_name in remaining_bes:
            cand_sym = Symmetry({original: cand_name, cand_name: original}, bes)
            # Check whether both BEs occur equally often in the cut sets.
            if be_mcs_occurrences[original] != be_mcs_occurrences[cand_name]:
                # BEs cannot be symmetric
                continue
            if cand_sym.is_valid_symmetry(mcss):
                # Found a symmetry
                symmetries.append(cand_sym)
                # Add symmetries to other BEs which are already symmetric to the original
                for cand in covered_bes:
                    # Symmetry is by construction valid
                    symmetries.append(Symmetry({cand_name: cand, cand: cand_name}, bes))
                covered_bes.append(cand_name)
        remaining_bes = remaining_bes.difference(covered_bes)
    return symmetries


def generate_all_symmetries(mcss, bes, be_mcs_occurrences):
    """
    Generate all valid symmetries with respect to MCS.
    Note that the approach uses brute-force and just tries out all permutations.
    This operation is therefore very expensive.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}.
    :param be_mcs_occurrences: Number of times a BE occurs in a cut set.
    :return: All symmetries.
    """

    def generate_candidate(symmetries, sym_mapping, avail_bes, be_mcs_occurrences):
        """
        Recursive function to generate and check candidate symmetries
        :param symmetries: Symmetries found so far. This will contain all symmetries in the end.
        :param sym_mapping: Partial mapping for symmetry.
        :param avail_bes: BEs not yet considered for the symmetry mapping.
        :param be_mcs_occurrences: Number of times a BE occurs in a cut set.
        """
        if len(avail_bes) == 0:
            if not sym_mapping:
                # Ignore empty symmetry
                return
            # Check symmetry
            symmetry = Symmetry(sym_mapping, bes)
            if symmetry.is_valid_symmetry(mcss):
                # Check whether other symmetry already captures the new symmetry
                for sym in symmetries:
                    if symmetry.is_subset(sym):
                        return
                symmetries.append(symmetry)
            return

        original = next(iter(avail_bes))
        # Consider all possible mapping for original
        # Consider self-symmetry last
        iter_list_bes = list(avail_bes)[1:] + [original]
        for cand in iter_list_bes:
            remaining = avail_bes.copy()
            extended_sym_map = sym_mapping.copy()
            if original == cand:
                # Identity mapping
                remaining.remove(original)
            else:
                # Check whether both BEs occur equally often in the cut sets.
                if be_mcs_occurrences[original] != be_mcs_occurrences[cand]:
                    # BEs cannot be symmetric
                    continue
                else:
                    # Add mapping
                    extended_sym_map[original] = cand
                    extended_sym_map[cand] = original
                    remaining.remove(original)
                    remaining.remove(cand)
            # Apply recursively
            generate_candidate(symmetries, extended_sym_map, remaining, be_mcs_occurrences)

    symmetries = Symmetries()
    # Generate all candidates
    generate_candidate(symmetries, dict(), SortedSet(bes.values()), be_mcs_occurrences)
    return symmetries


def find_symmetry_between_modules(module1, module2, mcss, bes, be_mcs_occurrences):
    """
    Find a valid symmetry between two modules with respect to MCS.
    Note that the approach uses brute-force and just tries out all mappings from module1 to module2.
    This operation can therefore be very expensive if no symmetry is present.
    :param module1: First module.
    :param module2: Second module.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dict of all BEs {index: name}.
    :param be_mcs_occurrences: Number of times a BE occurs in a cut set.
    :return: A symmetry or None if no valid symmetry exists.
    """

    def generate_candidate(sym_mapping, avail1, avail2, be_mcs_occurrences):
        """
        Recursive function to generate and check candidate symmetries
        :param sym_mapping: Partial mapping for symmetry.
        :param avail1: BEs of first module not yet considered for the symmetry mapping.
        :param avail2: BEs of second module not yet considered for the symmetry mapping.
        :param be_mcs_occurrences: Number of times a BE occurs in a cut set.
        :return: A symmetry or None if no valid symmetry exists.
        """
        if len(avail1) == 0:
            if not sym_mapping:
                # Ignore empty symmetry
                return None
            # Check symmetry
            symmetry = Symmetry(sym_mapping, bes)
            if symmetry.is_valid_symmetry(mcss):
                return symmetry
            return None

        original = next(iter(avail1))
        # Consider all possible mapping for original
        for cand in avail2:
            extended_sym_map = sym_mapping.copy()
            # Check whether both BEs occur equally often in the cut sets.
            if be_mcs_occurrences[original] != be_mcs_occurrences[cand]:
                # BEs cannot be symmetric
                continue
            else:
                remain1 = avail1.copy()
                remain2 = avail2.copy()
                # Add mapping
                extended_sym_map[original] = cand
                extended_sym_map[cand] = original
                remain1.remove(original)
                remain2.remove(cand)
            # Apply recursively
            symmetry = generate_candidate(extended_sym_map, remain1, remain2, be_mcs_occurrences)
            if symmetry:
                return symmetry

    if len(module1) == len(module2):
        # Symmetries can only occur if both modules have the same number of BEs
        bes1 = SortedSet(be for index, be in bes.items() if index in module1)
        bes2 = SortedSet(be for index, be in bes.items() if index in module2)
        return generate_candidate(dict(), bes1, bes2, be_mcs_occurrences)
    else:
        return None
