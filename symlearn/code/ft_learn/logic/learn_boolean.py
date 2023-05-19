import sympy.logic
import pyeda.inter
import logging

from ft_learn.ft.fault_tree import FaultTree
from ft_learn.ft.ft_elements import BE, AND, OR
import ft_learn.logic.boolean_logic as boolean_logic


def fault_tree_from_boolean_prefix(expression_string, basic_events, be_parser):
    """
    Create fault tree from Boolean expression in prefix notation.
    :param expression_string: String representation of Boolean expression in prefix notation.
    :param basic_events: Dictionary of BEs {index: name}.
    :param be_parser: Function to parse BEs.
    :return: Fault tree representing the given expression.
    """

    def element_from_boolean_prefix(formula, be_parser):
        """
        Create element(s) from string of formula in prefix notation.
        :param formula: Rest of formula to parse.
        :param be_parser: Function to parse BEs.
        :return: Next FT element, remaining part of formula without parsed FT element.
        """
        # Handle AND/OR
        if formula.startswith('And('):
            operator = AND
            rest = formula[4:]
        elif formula.startswith('Or('):
            operator = OR
            rest = formula[3:]
        else:
            return be_parser(formula, bes)

        # Parse first child
        child, rest = element_from_boolean_prefix(rest, be_parser)
        children = [child]
        # Parse remaining children
        while rest.startswith(", "):
            child, rest = element_from_boolean_prefix(rest[2:], be_parser)
            children.append(child)
        assert rest.startswith(")")
        # Create gate
        return operator(children), rest[1:]

    # Create BEs
    bes = dict()
    for be in basic_events.values():
        bes[be] = BE(be)
    # Parse formula
    top_event, remaining = element_from_boolean_prefix(expression_string, be_parser)
    assert remaining == ""
    return FaultTree(top_event)


def fault_tree_from_sympy_expression(expression, basic_events):
    """
    Create fault tree from Boolean expression (in sympy format).
    :param expression: Boolean expression.
    :param basic_events: Dictionary of BEs {index: name}.
    :return: Fault tree representing the given expression.
    """

    def parse_be_sympy(formula, bes):
        """
        Parse BE from sympy format.
        Format: Symbol('ABC')
        :param formula: Formula.
        :param bes: Dictionary of BEs.
        :return: BE, rest of formula.
        """
        assert formula.startswith('Symbol(')
        pos = formula.find(")")
        assert pos > 0
        be_name, rest = formula[8:pos - 1], formula[pos + 1:]
        assert be_name in bes
        return bes[be_name], rest

    return fault_tree_from_boolean_prefix(sympy.srepr(expression), basic_events, parse_be_sympy)


def fault_tree_from_pyeda_expression(expression, basic_events):
    """
    Create fault tree from Boolean expression (in pyeda format).
    :param expression: Boolean expression.
    :param basic_events: Dictionary of BEs {index: name}.
    :return: Fault tree representing the given expression.
    """

    def parse_be_pyeda(formula, bes):
        """
        Parse BE from pyeda format.
        Format: ABC
        :param formula: Formula.
        :param bes: Dictionary of BEs.
        :return: BE, rest of formula.
        """
        # Find end of BE name
        pos1 = formula.find(", ")
        pos2 = formula.find(")")
        assert pos2 != 0
        if 0 <= pos1 < pos2:
            pos = pos1
        else:
            pos = pos2
        be_name, rest = formula[:pos], formula[pos:]
        assert "," not in be_name and "(" not in be_name and ")" not in be_name
        assert be_name in bes
        return bes[be_name], rest

    return fault_tree_from_boolean_prefix(str(expression), basic_events, parse_be_pyeda)


def learn_ft_boolean_sympy_tt(mcss, bes):
    """
    Learn fault tree from truth table via Boolean formulas in sympy.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dictionary of BEs {index: name}.
    :return: Fault tree representing MCS.
    """
    # TODO: support for underspecified failure tables
    # formula = sympy.logic.SOPform(be_list, cutset_matrix.tolist())
    assert False


def learn_ft_boolean_sympy_mcss(mcss, bes):
    """
    Learn fault tree from MCS via Boolean formulas in sympy.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dictionary of BEs {index: name}.
    :return: Fault tree representing MCS.
    """
    formula = boolean_logic.mcs_to_sympy_formula(mcss, bes)
    # Simplify formula via sympy (will probably not do much as formula is already in SOP)
    formula = sympy.logic.boolalg.simplify_logic(formula, force=True)
    logging.debug("Simplified formula from sympy: {}".format(formula))
    return fault_tree_from_sympy_expression(formula, bes)


def learn_ft_espresso_mcss(mcss, bes):
    """
    Learn fault tree from MCS via Boolean formulas and Espresso minimization.
    :param mcss: Minimal cut sets (MCS).
    :param bes: Dictionary of BEs {index: name}.
    :return: Fault tree representing MCS.
    """
    # Create variables
    variables = dict()
    for index, be in bes.items():
        variables[index] = pyeda.inter.exprvar(be)

    # Construct Boolean formula from MCS
    formulas_cutsets = []
    for mcs in mcss:
        # Create formula for each cut set
        vars_cut_set = [variables[be] for be in mcs.set]
        formula_cutset = pyeda.inter.And(*vars_cut_set, simplify=False)
        formulas_cutsets.append(formula_cutset)
    formula = pyeda.inter.Or(*formulas_cutsets, simplify=False)

    # Simplify formula via Espresso heuristic logic minimizer
    simpler_formula = pyeda.inter.espresso_exprs(formula.to_dnf())[0]
    logging.debug("Simplified Boolean formula from espresso: {}".format(simpler_formula))
    # Create fault tree from formula
    return fault_tree_from_pyeda_expression(simpler_formula, bes)
