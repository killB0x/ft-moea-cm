import sympy

from ft_learn.ft.ft_elements import BE, AND, OR
from ft_learn.ft.mcs import CutSet, MinCutSets


def fault_tree_to_formula_string(fault_tree):
    """
    Return Boolean representation of the fault tree.
    :param fault_tree: Fault tree.
    :return: String encoding the Boolean formula of the fault tree.
    """

    def to_boolean_string(elem):
        if isinstance(elem, BE):
            return elem.name
        else:
            children = [to_boolean_string(child) for child in elem.children]
            if isinstance(elem, AND):
                return "(" + " & ".join(children) + ")"
            elif isinstance(elem, OR):
                return "(" + " | ".join(children) + ")"
            else:
                assert False

    return to_boolean_string(fault_tree.top_event)


def fault_tree_to_sympy_formula(fault_tree):
    """
    Return Boolean representation of the fault tree in sympy format.
    :param fault_tree: Fault tree.
    :return: Boolean formula of the fault tree.
    """

    def to_sympy_formula(elem, elements):
        if str(elem) in elements:
            return elements[str(elem)]
        if isinstance(elem, BE):
            symbol_be = sympy.Symbol(elem.name)
            elements[str(elem)] = symbol_be
            return symbol_be
        else:
            children = [to_sympy_formula(child, elements) for child in elem.children]
            if isinstance(elem, AND):
                return sympy.And(*children, simplify=False)
            elif isinstance(elem, OR):
                return sympy.Or(*children, simplify=False)
            else:
                assert False

    return to_sympy_formula(fault_tree.top_event, dict())


def sympy_formula_to_mcs(formula, bes):
    """
    Get MCS from Boolean formula using sympy.
    :param formula: Boolean formula.
    :param bes: List of all BEs.
    :return: Minimal cut sets (MCS) corresponding to given formula.
    """
    be_indices = {bes[i]: i for i in range(len(bes))}
    # Convert to DNF to allow easy extraction of cut sets
    expression = sympy.to_dnf(formula)

    mcss = MinCutSets()

    if expression.func == sympy.And or expression.func == sympy.Symbol:
        cut_sets_expr = [expression]
    else:
        assert expression.func == sympy.Or
        cut_sets_expr = expression.args
    # Iterate over all cut set formulas
    for cut_set_expr in cut_sets_expr:
        if cut_set_expr.func == sympy.Symbol:
            bes_expr = [cut_set_expr]
        else:
            assert cut_set_expr.func == sympy.And
            bes_expr = cut_set_expr.args
        # Iterate over all BEs in a cut set
        bes_cut_set = []
        for be_expr in bes_expr:
            assert be_expr.func == sympy.Symbol
            bes_cut_set.append(be_indices[be_expr.name])
        mcss.add(CutSet(bes_cut_set))
    return mcss


def mcs_to_sympy_formula(mcss, basic_events):
    # Create variables
    variables = dict()
    for index, be in basic_events.items():
        variables[index] = sympy.Symbol(be)

    # Construct Boolean formula from MCS
    formulas_cutsets = []
    for mcs in mcss:
        # Create formula for each cut set
        vars_cut_set = [variables[be] for be in mcs.set]
        formula_cutset = sympy.And(*vars_cut_set, simplify=False)
        formulas_cutsets.append(formula_cutset)
    return sympy.Or(*formulas_cutsets, simplify=False)
