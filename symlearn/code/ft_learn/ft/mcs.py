import logging
import numpy as np


class CutSet:
    """
    Cut set.
    """

    def __init__(self, entries):
        if type(entries) is frozenset:
            self.set = entries
        else:
            self.set = frozenset(entries)

    def __len__(self):
        return len(self.set)

    def __contains__(self, item):
        return item in self.set

    def __hash__(self):
        return hash(self.set)

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __eq__(self, other):
        if isinstance(other, CutSet):
            return self.set == other.set
        else:
            return False

    def issubset(self, other):
        return self.set.issubset(other.set)

    def isdisjoint(self, other):
        return self.set.isdisjoint(other.set)

    def union(self, other):
        return CutSet(self.set.union(other.set))

    def intersection(self, other):
        return CutSet(self.set.intersection(other.set))

    def difference(self, other):
        return CutSet(self.set.difference(other.set))

    def to_string(self, bes=None):
        """
        Convert cut set to string.
        Use BE names if they are given, otherwise use indices.
        :param bes: (Optional) Dictionary of BEs {index: name}.
        :return: String representation of cut set.
        """
        if bes:
            return "{" + ", ".join([bes[i] for i in self.set]) + "}"
        else:
            return "{" + ", ".join([str(i) for i in self.set]) + "}"

    def __str__(self):
        return self.to_string(bes=None)


class MinCutSets(set):
    """
    Minimal cut sets as sets of CutSet.
    """

    def compute_from_cut_sets(self, matrix):
        """
        Initialize MCSs from failure data.
        :param matrix: Matrix where each row corresponds to a system failure and columns correspond to BEs
        """
        # Two possible approaches for computing all minimal cut sets:
        # 1. directly on matrix
        # 2. via set operations
        # TODO: choose one approach
        approach_matrix = True
        if approach_matrix:
            # Approach 1: directly on matrix
            if len(matrix.shape) == 1:
                # Only one MCS
                self.add(CutSet(np.flatnonzero(matrix.reshape(1, len(matrix)))))
            else:
                while matrix.shape[0] > 0:
                    # Cut set with minimal order is in front
                    min_comb = sum(matrix.T)
                    min_cut_index = np.where(min_comb == min_comb.min())[0][0]
                    # Get corresponding cut set
                    cs = matrix[min_cut_index, :]
                    # Get the values of the BEs in the cut set for the other rows
                    other_values_cs = matrix[:, np.where(cs)[0]]
                    # Delete all rows which are inferred by this cut set (have also 1's at the same positions)
                    to_delete = np.where(sum(other_values_cs.T) == other_values_cs.shape[1])
                    assert len(to_delete) > 0
                    matrix = np.delete(matrix, to_delete, 0)
                    # Add cut set
                    self.add(CutSet(np.flatnonzero(cs)))

        else:
            # Approach 2: using set operations
            for cut_set in matrix:
                # Get indices of non-zero entries
                cs = CutSet(np.flatnonzero(cut_set))

                if len(cs) == 1:
                    # Singleton cut sets are minimal
                    self.add(cs)
                else:
                    # Check whether cut set is already implied by existing (minimal) cut set
                    contained = False
                    for mcs in self:
                        if mcs.issubset(cs):
                            contained = True
                            logging.debug("Cut set {} is already contained in {}".format(cs, mcs))
                            break
                    if not contained:
                        self.add(cs)

    def get_always_failed(self):
        """
        Get the BEs which are required to be failed in all MCS.
        :return: List of BE indices
        """
        first = True
        for mcs in self:
            if first:
                required = mcs
                first = False
            required = required.intersection(mcs)
            if not required:
                break
        return required

    def without_bes(self, bes):
        """
        Remove BEs from MCSs.
        :param bes: List of BEs to remove.
        :return: MinCutSets without given BEs.
        """
        mcss = MinCutSets()
        for mcs in self:
            new_mcs = mcs.difference(bes)
            if new_mcs:
                mcss.add(new_mcs)
        return mcss

    def get_matrix(self, be_indices):
        """
        Return MCS as matrix for given BEs.
        :param be_indices: List of BE indices to use for the columns.
        :return: Matrix where each row corresponds to a MCS.
        """
        data = []
        for mcs in self:
            row = []
            # Check that indices are captured by given BEs.
            for index in mcs.set:
                assert index in be_indices
            # Create row for cut set
            for be in be_indices:
                row.append(1 if be in mcs else 0)
            data.append(row)
        return np.asarray(data)

    def get_be_occurrences(self, bes):
        """
        Return the number of a times a BE occurs in an MCS.
        This can be used to cancel symmetry computations early on.
        :param bes: Dict of all BEs {index: name}.
        :return: Dict {BE name: number of occurrences}.
        """
        # Initialize
        be_occurrences = dict()
        for index, name in bes.items():
            be_occurrences[name] = 0
        # Count occurrences in MCS
        for cutset in self:
            for index in cutset.set:
                be_occurrences[bes[index]] += 1
        return be_occurrences

    def __str__(self):
        return "{" + ", ".join([str(mcs) for mcs in self]) + "}"
