from copy import deepcopy
import random
import logging
import string

from ft_learn.ft.ft_elements import AND, OR
import ft_learn.helper as helper
import multiprocessing
from multiprocessing import Pool

def generate_name():
    """
    Generate random name of length 7.
    :return: Random name.
    """
    return ''.join(random.choices(string.ascii_lowercase, k=7))


def create_be(ft, all_bes, p_create_be, deterministic=False):
    """
    Add a new BE (which is not yet present) to the fault tree.
    :param ft: Fault tree.
    :param all_bes: List of all BEs which can possibly occur in the fault tree.
    :param p_create_be: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_create_be:
        return ft

    new_ft = ft.copy()
    # Possible BEs to add
    possible_bes = list(set(all_bes) - set(new_ft.get_all_bes(sort=False)))
    if not possible_bes:
        return ft
    if deterministic:
        possible_bes.sort(key=lambda x: x.name)
    be_to_connect = deepcopy(random.choice(possible_bes))
    assert len(be_to_connect.parents) == 0
    gate = random.choice(new_ft.get_all_gates(sort=deterministic))
    gate.add_child(be_to_connect)
    assert len(be_to_connect.parents) == 1
    return new_ft


def connect_be(ft, p_connect_be, deterministic=False):
    """
    Make a new connection between a BE and gate.
    :param ft: Fault tree.
    :param p_connect_be: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_connect_be:
        return ft

    new_ft = ft.copy()
    be_to_connect = random.choice(new_ft.get_all_bes(sort=deterministic))
    assert len(be_to_connect.parents) >= 1
    possible_gates = list(set(new_ft.get_all_gates(sort=False)) - set(be_to_connect.parents))
    if not possible_gates:
        return ft
    if deterministic:
        possible_gates.sort(key=lambda x: str(x))
    new_parent = random.choice(possible_gates)
    new_parent.add_child(be_to_connect)
    return new_ft


def disconnect_be(ft, p_disconnect_be, deterministic=False):
    """
    Remove a connection between a random BE and a gate. Removing the last connection removes the BE from the fault tree.
    :param ft: Fault tree.
    :param p_disconnect_be: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_disconnect_be:
        return ft

    new_ft = ft.copy()
    be_to_del = random.choice(new_ft.get_all_bes(sort=deterministic))
    # Choose parent gate from which to remove the BE
    assert len(be_to_del.parents) >= 1
    possible_parents = be_to_del.parents
    if deterministic:
        possible_parents.sort(key=lambda x: str(x))
    parent = random.choice(possible_parents)
    parent.remove_child(be_to_del)
    return new_ft


def delete_be(ft, p_delete_be, deterministic=False):
    """
    Completely remove a BE from the fault tree.
    :param ft: Fault tree.
    :param p_delete_be: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_delete_be:
        return ft

    new_ft = ft.copy()
    be_to_del = random.choice(new_ft.get_all_bes(sort=deterministic))
    assert len(be_to_del.parents) >= 1
    # Remove from all parents
    for parent in list(be_to_del.parents):
        parent.remove_child(be_to_del)
    return new_ft


def move_be(ft, p_move_be, deterministic=False):
    """
    Move BE to other parent gate.
    :param ft: Fault tree.
    :param p_move_be: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_move_be:
        return ft

    new_ft = ft.copy()
    be_to_swap = random.choice(new_ft.get_all_bes(sort=deterministic))
    if len(be_to_swap.parents) > 1:
        # Moving only works if the BE has a single parent
        return ft
    parent_of_be = be_to_swap.parents[0]
    possible_gates = list(set(new_ft.get_all_gates(sort=False)) - {parent_of_be})
    if not possible_gates:
        return ft
    if deterministic:
        possible_gates.sort(key=lambda x: str(x))
    new_parent = random.choice(possible_gates)
    parent_of_be.remove_child(be_to_swap)
    new_parent.add_child(be_to_swap)
    return new_ft


def create_gate(ft, p_create_gate, deterministic=False):
    """
    Insert new gate (of type AND or OR).
    :param ft: Fault tree.
    :param p_create_gate: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_create_gate:
        return ft

    new_ft = ft.copy()
    gate = random.choice(new_ft.get_all_gates(sort=deterministic))
    if len(gate.children) == 1:
        # Creating new gate gives no new fault tree
        return ft
    no_children_to_select = random.randint(1, len(gate.children) - 1)  # Do not select all children
    selected_children = random.sample(gate.children, k=no_children_to_select)
    # Remove children from old gate
    for child in selected_children:
        gate.remove_child(child)
    # Add new gate as input of old gate
    if random.choice([True, False]):
        new_gate = AND(selected_children)
    else:
        new_gate = OR(selected_children)

    if new_gate in gate.children:
        # Gate already exists
        return ft
    else:
        gate.add_child(new_gate)
        return new_ft


def change_gate_type(ft, p_change_gate_type, deterministic=False):
    """
    Change type of gate (AND->OR or OR->AND).
    :param ft: Fault tree.
    :param p_change_gate_type: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_change_gate_type:
        return ft
    new_ft = ft.copy()
    gate_to_change = random.choice(new_ft.get_all_gates(sort=deterministic))
    # Remember parents of gate for later
    parents_of_gate = [parent for parent in gate_to_change.parents]
    is_top_event = gate_to_change == new_ft.top_event
    # Remove old gate
    new_ft.remove_gate(gate_to_change)
    # Create new gate (with the same children and name)
    if isinstance(gate_to_change, AND):
        new_gate = OR(gate_to_change.children, name=gate_to_change.name)
    else:
        new_gate = AND(gate_to_change.children, name=gate_to_change.name)
    if is_top_event:
        new_ft.top_event = new_gate
    else:
        for parent in parents_of_gate:
            parent.add_child(new_gate)
    return new_ft


def delete_gate(ft, p_delete_gate, deterministic=False):
    """
    Delete intermediate gate and add the children of the gate to the parent of the gate.
    :param ft: Fault tree.
    :param p_delete_gate: Probability of performing this action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified fault tree.
    """
    if random.random() >= p_delete_gate:
        return ft
    new_ft = ft.copy()
    possible_gates = list(set(new_ft.get_all_gates(sort=False)) - {new_ft.top_event})  # Do not delete top gate
    if not possible_gates:
        return ft
    if deterministic:
        possible_gates.sort(key=lambda x: str(x))
    gate_to_del = random.choice(possible_gates)
    assert len(gate_to_del.parents) >= 1
    # Remember parents of gate for later
    parents_of_gate = [parent for parent in gate_to_del.parents]
    new_ft.remove_gate(gate_to_del)
    # Add children of gate to parent gate
    for parent in parents_of_gate:
        for child in gate_to_del.children:
            if child not in parent.children:
                parent.add_child(child)
    return new_ft


def cross_over_one_side(ft, old_gate, new_gate, children):
    """
    Perform cross over action by removing old gate and inserting new gate.
    :param ft: Fault tree.
    :param old_gate: Old gate to remove.
    :param new_gate: New gate to use instead.
    :param children: Children of old gate.
    :return: Modified fault tree.
    """

    # Check whether new gate already exists
    duplicate = new_gate in children
    for parent in old_gate.parents:
        for child in parent.children:
            if new_gate.to_string(True) == child.to_string(True):
                duplicate = True
                break
    if duplicate:
        # Generate random name for duplicate child to help distinction
        new_gate.name = generate_name()

    # Modify fault tree
    if old_gate == ft.top_event:
        # Set new top event
        ft.remove_gate(old_gate)
        ft.top_event = new_gate
    else:
        assert len(old_gate.parents) >= 1
        # Remember parents of gate for later
        parents_of_gate = [parent for parent in old_gate.parents]
        ft.remove_gate(old_gate)

        for parent_1 in parents_of_gate:
            parent_1.add_child(new_gate)

    return ft.copy(allow_duplicates=True)


def cross_over(ft1, ft2, p_cross_over, deterministic):
    """
    Perform cross-over swap between two fault tree.
    The cross-over can introduce duplicate children (but with different names) into a fault tree. This is allowed behaviour.
    :param ft1: First fault tree.
    :param ft2: Second fault tree.
    :param p_cross_over: Probability of performing the action.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: Modified first and second fault tree.
    """
    if random.random() >= p_cross_over:
        return ft1, ft2

    new_ft1 = ft1.copy()
    new_ft2 = ft2.copy()
    gate_1 = random.choice(new_ft1.get_all_gates(sort=deterministic))
    gate_2 = random.choice(new_ft2.get_all_gates(sort=deterministic))
    # Abort if both gates are the same
    if str(gate_1) == str(gate_2):
        return ft1, ft2

    # Create copies of children
    children_1 = deepcopy(gate_1.children)
    children_2 = deepcopy(gate_2.children)

    # Create new gates
    new_gate_1 = type(gate_2)(children_2)
    new_gate_2 = type(gate_1)(children_1)

    # Modify fault trees
    new_ft1 = cross_over_one_side(new_ft1, gate_1, new_gate_1, children_1)
    new_ft2 = cross_over_one_side(new_ft2, gate_2, new_gate_2, children_2)
    return new_ft1, new_ft2


class GenOpConfig:
    """
    Configuration for genetic operations.
    The following options exists:
    - p_create_be: Probability of adding a new BE.
    - p_connect_be: Probability of connecting a BE with a gate.
    - p_disconnect_be: Probability of removing a connection.
    - p_move_be: Probability of moving a BE to a new gate.
    - p_create_gate: Probability of creating a new gate.
    - p_change_gate_type: Probability of changing a gate type.
    - p_delete_gate: Probability of removing a gate.
    - p_cross_over: Probability of performing a cross-over between two fault tree.
    """

    def __init__(self, default_prob):
        """
        Constructor.
        :param default_prob: The default value for all probabilities.
        """
        self.p_create_be = default_prob
        self.p_connect_be = default_prob
        self.p_disconnect_be = default_prob
        self.p_delete_be = default_prob
        self.p_move_be = default_prob
        self.p_create_gate = default_prob
        self.p_change_gate_type = default_prob
        self.p_delete_gate = default_prob
        self.p_cross_over = default_prob


def apply_genetic_operators(population, all_bes, prob_config, deterministic=False):
    """
    Apply genetic operators on fault-tree population.
    :param population: Population of fault tree.
    :param all_bes: List of all BEs.
    :param prob_config: Configuration of probabilities for genetic operations.
    :param deterministic: Whether deterministic results should be ensured (useful for debugging).
    :return: New population of fault trees.
    """
    results = []
    for ft in population:
        result = operate_on_ft(ft, population, all_bes, prob_config, deterministic)
        results.append(result)

    new_population = []
    total_assert_errors = 0
    for result in results:
        new_population.extend(result[0])
        total_assert_errors += result[1]
    new_population = list(set(new_population))
    if total_assert_errors > 0:
        logging.debug("Encountered {} assertion errors".format(total_assert_errors))
    return new_population

def operate_on_ft(ft, population, all_bes, prob_config, deterministic=False):
    actions = [create_be, connect_be, disconnect_be, delete_be, move_be, create_gate, change_gate_type, delete_gate, cross_over]
    new_population = [ft]
    assert_errors = 0

    for action in actions:
        prob_action = getattr(prob_config, "p_{}".format(action.__name__))
        try:
            if action == cross_over:
                new_ft, new_ft_2 = action(ft, random.choice(population), prob_action, deterministic)
                if not helper.check_empty_objects(new_ft_2):
                    if new_ft_2 not in new_population:
                        new_ft_2.add_to_history(action.__name__)
                        new_population.append(new_ft_2)
            elif action == create_be:
                new_ft = action(ft, all_bes, prob_action, deterministic)
            else:
                new_ft = action(ft, prob_action, deterministic)
        except AssertionError:
            assert_errors += 1
            new_ft = None

        if new_ft and not helper.check_empty_objects(new_ft):
                # print(new_ft.operation_history)
                new_ft.add_to_history(action.__name__)
                new_population.append(new_ft)

    return new_population, assert_errors

def apply_genetic_operators_multithreaded(population, all_bes, prob_config, deterministic=False):
    with Pool() as p:
        results = p.starmap(operate_on_ft, [(ft, population, all_bes, prob_config, deterministic) for ft in population])

    new_population = []
    total_assert_errors = 0
    for result in results:
        new_population.extend(result[0])
        total_assert_errors += result[1]
    new_population = list(set(new_population))
    if total_assert_errors > 0:
        logging.debug("Encountered {} assertion errors".format(total_assert_errors))
    return new_population