from sortedcontainers import SortedSet

from ft_learn.ft.mcs import MinCutSets


class Modules(SortedSet):
    """
    Modules (independent sub-parts of the failure data).
    """

    def __init__(self, seq=()):
        """
        Constructor.
        :param seq: Iterator of elements to insert into set.
        """
        super().__init__(seq)

    def modularize_mcss(self, mcss, gate_or):
        """
        Split MCSs according to modules.
        :param mcss: Minimal cut sets
        :param gate_or: Whether the top gate is an OR-gate.
        :return: Dict: Module -> MCSs
        """
        mod_mcss = dict()
        for mod in self:
            mod_mcss[mod] = MinCutSets()

        # Add each MCS to corresponding module
        for mcs in mcss:
            if gate_or:
                # For OR-gate: add MCS to containing module
                for mod in self:
                    if mcs.issubset(mod):
                        mod_mcss[mod].add(mcs)
            else:
                # For AND-gate: split MCS according to modules
                for mod in self:
                    new_mcs = mcs.intersection(mod)
                    mod_mcss[mod].add(new_mcs)
        return mod_mcss

    def add_and_merge(self, module):
        """
        Add new module and possibly merge modules into one in case of overlap.
        :param module: New module to add.
        """
        # Find all existing candidates which overlap with new module
        candidates = []
        for mod in self:
            if not module.isdisjoint(mod):
                candidates.append(mod)

        new_module = module
        # Remove candidates and build new module by union
        for cand in candidates:
            self.remove(cand)
            new_module = new_module.union(cand)
        self.add(new_module)

    def to_string(self, bes=None):
        """
        Convert module to string.
        Use BE names if they are given, otherwise use indices.
        :param bes: (Optional) Dictionary of BEs {index: name}.
        :return: String representation of modules.
        """
        return ", ".join([mod.to_string(bes) for mod in self])

    def __str__(self):
        return self.to_string(bes=None)
