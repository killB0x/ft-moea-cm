import argparse
import logging
from enum import Enum
from timeit import default_timer as timer
import os.path
from pathlib import Path
import datetime
import pandas as pd
from scipy.io import loadmat

from ft_learn.ft.fault_tree import FaultTree
from ft_learn.ft.ft_elements import BE, AND, OR
from ft_learn.ft.mcs import MinCutSets, CutSet
from ft_learn.modules.modules import Modules
from ft_learn.modules.symmetries import Symmetries, get_symmetries_from_file, generate_all_symmetries, generate_singleton_symmetries, find_symmetry_between_modules
from ft_learn.modules.modules_finder import create_from_mcss, find_pseudo_modules_from_symmetry_under_or, split_mcss_from_symmetry
from ft_learn.moea.ft_moea import perform_genetic_ftmoea
from ft_learn.moea.genetic_operators import GenOpConfig
import ft_learn.logic.learn_boolean as learn_boolean
from ft_learn.results import Results
import ft_learn.helper as helper
from ft_learn.dataset import reduce_dataset


class LearnApproach(Enum):
    """
    Different approaches for learning fault trees.
    """
    FTMOEA = "ftmoea"
    SYMPY = "sympy"
    ESPRESSO = "espresso"

    def __str__(self):
        return self.value


class Config:
    """
    Configuration for complete toolchain.
    """

    def __init__(self):
        """
        Constructor.
        """
        # Approach to learn single FT
        self.learn_approach = LearnApproach.FTMOEA
        # Finding and exploiting modules
        self.use_modules = True
        # Finding and exploiting symmetries
        self.use_symmetries = True
        # Recursive learning
        self.use_recursion = True
        # Debugging (enables reproducible random generation)
        self.debug = False
        # Conversion to CNF for searching modules under AND-gates
        self.convert_cnf = False

        # Settings for FT-MOEA
        # Probability to perform each genetic operation
        self.probs_config = GenOpConfig(default_prob=1.0)
        # Population size
        self.population = 400
        # Maximum number of generations
        self.max_generations = 100
        # Stop after not changing in generations
        self.unchanged_generations = 20
        # Multi-objective function: [phi_c (Accuracy MCS), phi_s (Size), phi_d (Accuracy data)]
        self.obj_functions = [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        # Selection strategy
        self.selection_strategy = 'elitist'
        # Saving folder
        self.saving_folder = ""
        # Initial parent FT:
        self.ft_as_input = ''
        self.seg_size = 4


def log_debug(string, recurse_level):
    logger.debug("\t" * recurse_level + string)


def log_info(string, recurse_level):
    logger.info("\t" * recurse_level + string)


def learn_new_fault_tree(mcss, bes, all_bes, config, results, dataset_evaluation=None, recurse_level=0):
    """
    Learn fault tree from MCSs using one of several approaches (FT-MOEA, via Boolean formulas, etc.)
    :param mcss: Minimal cut sets.
    :param bes: Dictionary of BEs {index: name}.
    :param all_bes: All BEs which are globally present.
    :param config: Configuration.
    :param results: Data structure to track results.
    :param dataset_evaluation: Failure data for later evaluation.
    :param recurse_level: Level of nested recursion. Used for pretty printing.
    :return: Fault tree corresponding to MCSs.
    """

    # Learn fault tree
    # --------------------
    start_learn = timer()
    if config.learn_approach == LearnApproach.FTMOEA:
        # Use multi-objective genetic algorithm (FT-MOEA)
        log_debug("Perform FT-MOEA for module {}".format(CutSet(bes.keys()).to_string(bes)), recurse_level)

        # All BEs considered in the problem:
        bes_all = []
        be_all_indices = sorted(all_bes.keys())
        for i in be_all_indices:
            bes_all.append(all_bes[i])

        # Use only BEs contained in this module
        bes_module = []
        bes_module_indices = sorted(bes.keys())
        for i in bes_module_indices:
            bes_module.append(bes[i])

        # Get MCS matrix from MCSs (With ALL BEs)
        mcs_matrix_ftmoea = mcss.get_matrix(be_all_indices)

        dataset_evaluation_ftmoea = reduce_dataset(dataset_evaluation, bes_module)

        log_debug(str(mcs_matrix_ftmoea), recurse_level)
        log_debug(str(bes_module), recurse_level)

        # config.ft_as_input = helper.cutsets2ft(mcs_matrix_ftmoea, bes_all)
        config.ft_as_input = ""

        # Store BEs (globally and per module)
        bes_ftmoea = {'all': bes_all, 'module': bes_module}
        
        fts, _, _ = perform_genetic_ftmoea(dataset=dataset_evaluation_ftmoea, MCSs=mcs_matrix_ftmoea, bes=bes_ftmoea, population_size=config.population,
                                           generations=config.max_generations, convergence_criterion=config.unchanged_generations, multi_objective_function=config.obj_functions,
                                           config_gen_op=config.probs_config, selection_strategy=config.selection_strategy, debugging=config.debug,
                                           path_save_results=config.saving_folder,
                                           ft_as_input=config.ft_as_input, seg_size=config.seg_size, dataset_name=os.path.splitext(os.path.basename(args.file)),
                                           use_multithreading = config.use_multithreading,
                                           use_caching = config.use_caching)
        ft = fts[-1]
    elif config.learn_approach == LearnApproach.SYMPY:
        log_debug("Learn FT via sympy for module {}".format(CutSet(bes.keys()).to_string(bes)), recurse_level)
        # cutset_matrix = cutset_matrix[:, list(mod_bes.keys())]
        # cutset_matrix = cutset_matrix[~np.all(cutset_matrix == 0, axis=1)]
        # ft = learn_boolean.learn_ft_boolean_sympy(cutset_matrix, mod_bes)
        ft = learn_boolean.learn_ft_boolean_sympy_mcss(mcss, bes)
    elif config.learn_approach == LearnApproach.ESPRESSO:
        log_debug("Learn FT via espresso for module {}".format(CutSet(bes.keys()).to_string(bes)), recurse_level)
        ft = learn_boolean.learn_ft_espresso_mcss(mcss, bes)
    else:
        assert False

    results.time_learn_single += timer() - start_learn
    results.fts_single.append((ft, config.learn_approach))
    log_debug("Learned FT: {}".format(ft), recurse_level)
    return ft


def learn_single_module(mcss, bes, all_bes, config=None, results=None, dataset_evaluation=None, symmetries=None, recurse_level=0):
    """
    Learn fault trees for single module (recursively).
    :param mcss: Minimal cut sets.
    :param bes: Dictionary of BEs {index: name}.
    :param all_bes: All BEs which are globally present.
    :param config: Configuration.
    :param results: Data structure to track results.
    :param dataset_evaluation: Failure data for later evaluation.
    :param symmetries: Symmetries.
    :param recurse_level: Level of nested recursion. Used for pretty printing.
    :return: Fault tree corresponding to MCSs.
    """
    # Automatically obtain all symmetries for a module
    if config.use_symmetries:
        start_symmetries = timer()
        be_mcs_occurences = mcss.get_be_occurrences(bes)
        if len(bes) <= 50:  # Avoid magic number
            log_info("Computing all symmetries for module {}...".format(CutSet(bes.keys()).to_string(bes)), recurse_level)
            symmetries = generate_all_symmetries(mcss, bes, be_mcs_occurences)
        else:
            # Fall back to more efficient computation by only considering some symmetries
            logger.info("Computing all BE symmetries for module {}...".format(CutSet(bes.keys()).to_string(bes)))
            symmetries = generate_singleton_symmetries(mcss, bes, be_mcs_occurences)
        if symmetries:
            log_info("Found symmetries:\n" + str(symmetries), recurse_level)
        else:
            log_info("Found no symmetries", recurse_level)
        results.time_symmetries += timer() - start_symmetries
    else:
        symmetries = None

    successful_split = False

    # Try to find some pseudo modules from the symmetries
    start_modules = timer()
    best_symmetry = None
    if symmetries:
        for symmetry in symmetries:
            orig_module, orig_mcss, sym_module, sym_mcss = find_pseudo_modules_from_symmetry_under_or(mcss, bes, symmetry)
            if sym_module:
                # Split into original and symmetric part was possible
                log_info("Use symmetry with shared BEs. Original module: {}, symmetric module: {}".format(orig_module.to_string(bes), sym_module.to_string(bes)), recurse_level)
                successful_split = True
                best_symmetry = symmetry
                # Remember symmetry used
                results.symmetries.append(best_symmetry)
                break

    if not successful_split and symmetries:
        # Try to just split MCS
        best_split = len(bes) * 2 + 1
        orig_module, orig_mcss, sym_module, sym_mcss = None, None, None, None
        for symmetry in symmetries:
            log_debug("Try to split MCS with symmetry: {}".format(symmetry), recurse_level)
            feasible, split_number, cand_orig_module, cand_orig_mcss, cand_sym_module, cand_sym_mcss = split_mcss_from_symmetry(mcss, bes, symmetry)
            if feasible and split_number < best_split:
                # Found better split
                orig_module, orig_mcss, sym_module, sym_mcss = cand_orig_module, cand_orig_mcss, cand_sym_module, cand_sym_mcss
                best_split = split_number
                best_symmetry = symmetry
                log_debug("Found better split of MCS.\n{0}Original MCS: {1}\n{0}Symmetric MCS:{2}".format("\t" * recurse_level, orig_mcss, sym_mcss), recurse_level)
        if best_split < len(bes) * 2 + 1:
            # Split into original and symmetric part was possible
            log_info("Split MCS according to symmetry.\n{0}Original MCS: {1}\n{0}Symmetric MCS:{2}".format("\t" * recurse_level, orig_mcss, sym_mcss), recurse_level)
            successful_split = True
            # Remember symmetry used
            results.symmetries.append(best_symmetry)
    results.time_dependent_modules += timer() - start_modules

    if successful_split:
        # Keep only BEs of module
        orig_bes = {i: be for i, be in bes.items() if i in orig_module}

        if config.use_recursion:
            # Recursively learn fault tree for original part
            # TODO: use subsets for dataset_evaluation
            orig_ft = recursive_learning(orig_mcss, orig_bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, symmetries=None,
                                         recurse_level=recurse_level + 1)  # Symmetries are not helpful for sub-part anymore
        else:
            orig_ft = learn_new_fault_tree(orig_mcss, orig_bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, recurse_level=recurse_level)

        log_debug("Learned original FT: {}".format(orig_ft), recurse_level)
        sym_ft = best_symmetry.apply_ft(orig_ft)
        log_debug("Symmetric FT: {}".format(sym_ft), recurse_level)
        return FaultTree(OR([orig_ft.top_event, sym_ft.top_event]))
    else:
        # Learn new fault tree
        log_info("No split into modules possible.", recurse_level)
        return learn_new_fault_tree(mcss, bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, recurse_level=recurse_level)


def recursive_learning(mcss, bes, all_bes, config=None, results=None, dataset_evaluation=None, symmetries=None, recurse_level=0):
    """
    Perform recursive fault tree learning from MCSs.
    :param mcss: Minimal cut sets.
    :param bes: Dictionary of BEs {index: name}.
    :param all_bes: All BEs which are globally present.
    :param config: Configuration.
    :param results: Data structure to track results.
    :param dataset_evaluation: Failure data for later evaluation.
    :param symmetries: Symmetries.
    :param recurse_level: Level of nested recursion. Used for pretty printing.
    :return: Fault tree corresponding to MCSs.
    """
    log_debug("Recursive learning of FT from MCSs: {}".format(mcss), recurse_level)
    # Obtain modules from data
    # --------------------
    start_modules = timer()
    if config.use_modules:
        mods, module_gate_is_or = create_from_mcss(mcss, bes, try_and=config.convert_cnf)
        log_debug("{} Module(s) under {}-gate: {}".format(len(mods), "OR" if module_gate_is_or else "AND", mods.to_string(bes)), recurse_level)
    else:
        mods = Modules()
        mods.add(CutSet(bes.keys()))
        module_gate_is_or = True

    # Obtain subsets of MCSs according to modules
    # --------------------
    modularized_mcss = mods.modularize_mcss(mcss, module_gate_is_or)
    log_debug("MCSs per module:\n" + "\n".join(["\t" * recurse_level + "\t{}:\t{}".format(mod.to_string(bes), mcs) for mod, mcs in modularized_mcss.items()]), recurse_level)
    results.time_independent_modules += timer() - start_modules

    # Obtain symmetries according to modules
    start_symmetries = timer()
    modules = list(modularized_mcss.keys())
    if len(modules) > 1:
        # Search for symmetries between modules
        be_mcs_occurences = mcss.get_be_occurrences(bes)
        if config.use_symmetries:
            symmetries = Symmetries()
            for i in range(len(modules)):
                for j in range(i + 1, len(modules)):
                    module1 = modules[i]
                    module2 = modules[j]
                    assert module1 != module2
                    symmetry = find_symmetry_between_modules(module1, module2, mcss, bes, be_mcs_occurences)
                    if symmetry:
                        logging.info("Found symmetry between module {} and module {}:\n{}".format(module1.to_string(bes), module2.to_string(bes), symmetry))
                        symmetries.append(symmetry)
    results.time_symmetries += timer() - start_symmetries

    # Perform learning for each set of MCSs independently
    fts = dict()  # Stores correspondence {MCS -> FT}
    covered_mcss = MinCutSets()  # Already covered MCS
    for module, mod_mcss in modularized_mcss.items():
        results.modules.append(module)

        # Check whether symmetric module is already covered
        ft = None
        if symmetries:
            for symmetry in symmetries:
                symmetric_mcss = symmetry.apply_mcss(mod_mcss)
                if str(symmetric_mcss) in fts:
                    # Create symmetric fault tree
                    orig_ft = fts[str(symmetric_mcss)]
                    log_info("Create symmetric fault tree for existing fault tree {}.".format(orig_ft), recurse_level)
                    ft = symmetry.apply_ft(orig_ft)
                    # Remember symmetry used
                    results.symmetries.append(symmetry)
                    break

        # No symmetry present -> learn new fault tree
        if not ft:
            # Keep only BEs of module
            mod_bes = {i: be for i, be in bes.items() if i in module}

            # Obtain BEs which always have to fail
            always_failed_bes = mod_mcss.get_always_failed()

            if always_failed_bes:
                # Create FT by AND over always failed BEs
                log_debug("Always failed BEs: {}".format(always_failed_bes.to_string(mod_bes)), recurse_level)
                if len(always_failed_bes.set) > 1:
                    all_failed_bes = [BE(mod_bes[be]) for be in always_failed_bes.set]
                    ft = FaultTree(AND(all_failed_bes))
                else:
                    single_be = next(iter(always_failed_bes.set))
                    ft = FaultTree(BE(mod_bes[single_be]))
                # Get sub-sets without always failed BEs
                sub_mcss = mod_mcss.without_bes(always_failed_bes)
                sub_bes = {i: be for i, be in mod_bes.items() if be not in always_failed_bes}
                if sub_mcss:
                    # TODO: use subsets for dataset_evaluation
                    if config.use_recursion:
                        # Recursively learn sub-FT
                        sub_ft = recursive_learning(sub_mcss, sub_bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, symmetries=symmetries,
                                                    recurse_level=recurse_level + 1)
                    else:
                        # Directly learn new fault tree
                        sub_ft = learn_new_fault_tree(sub_mcss, sub_bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, recurse_level=recurse_level)
                    ft = FaultTree(AND([ft.top_event, sub_ft.top_event]))

            else:
                # No split of always-failed BEs possible -> learn fault tree for single module
                ft = learn_single_module(mod_mcss, mod_bes, all_bes, config, results=results, dataset_evaluation=dataset_evaluation, symmetries=None,
                                         recurse_level=recurse_level)  # Forget symmetries as they are not computed within a module

        log_info("FT for module {}: {}".format(module.to_string(bes), ft), recurse_level)

        # Remember correspondence between MCS and FT
        fts[str(mod_mcss)] = ft

        # Update the already covered MCS (for later assertion check)
        if module_gate_is_or:
            covered_mcss.update(mod_mcss)
        else:
            if not covered_mcss:
                covered_mcss = mod_mcss
            else:
                newly_covered = MinCutSets()
                for part1 in covered_mcss:
                    for part2 in mod_mcss:
                        newly_covered.add(part1.union(part2))
                covered_mcss = newly_covered

    # Compose fault trees for modules
    # --------------------
    if len(fts) > 1:
        if module_gate_is_or:
            ft = FaultTree(OR([ft.top_event for ft in fts.values()]))
        else:
            ft = FaultTree(AND([ft.top_event for ft in fts.values()]))
    else:
        assert len(fts) == 1
        ft = next(iter(fts.values()))
    assert covered_mcss == mcss
    return ft


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Learn fault tree from data.')

    parser.add_argument('file', metavar='F', type=str, help='Input csv file containing all possible states of BEs and the resulting system status.')
    parser.add_argument('--learn-approach', '-a', type=LearnApproach, choices=list(LearnApproach), default=LearnApproach.FTMOEA)
    parser.add_argument('--symmetries', '-s', type=str, help='Optional file containing symmetries for the BEs.')
    parser.add_argument('--disable-symmetries', help='Disable use of symmetries', action='store_true')
    parser.add_argument('--disable-modules', help='Disable use of modules', action='store_true')
    parser.add_argument('--disable-recursion', help='Disable recursive calls during learning', action='store_true')
    parser.add_argument('--disable-multithreading', help='Disable use of multiple CPU cores', action='store_true')
    parser.add_argument('--disable-caching', help='Disable caching of fault tree metrics', action='store_true')
    parser.add_argument('--result-dir', '-r', type=str, help='Directory to write the results into as matlab files.')
    parser.add_argument('--verbose', '-v', help='Enable verbose output', action='store_true')
    parser.add_argument('--debug', '-d', help='Enable debugging', action='store_true')
    parser.add_argument('--segment-size', '-sz',type=int, help='Set segment size for random segmentation', default='4')
    parser.add_argument('--metric-config','-mc', type=str, help='Set metric configuration', default="11100000000000000000000")

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)

    # Init configuration
    config = Config()
    config.debug = args.debug
    config.learn_approach = args.learn_approach
    config.use_symmetries = not args.disable_symmetries
    config.use_modules = not args.disable_modules
    config.use_recursion = not args.disable_recursion
    config.use_multithreading = not args.disable_multithreading
    config.use_caching = not args.disable_caching
    config.seg_size = args.segment_size
    config.obj_functions = [int(char) for char in args.metric_config]
    print(config.seg_size, config.obj_functions)

    # Set directory for results
    filename, file_extension = os.path.splitext(os.path.basename(args.file))
    if args.result_dir:
        config.saving_folder = os.path.join(args.result_dir, filename)
        Path(config.saving_folder).mkdir(parents=True, exist_ok=True)

    # Start timing
    start = timer()

    # Parse failure data
    # --------------------
    logger.debug("Parsing from {}".format(args.file))
    if file_extension == ".mat":
        # Read data from matlab file
        data_dict = loadmat(args.file)
        data_array = data_dict['D']
        no_bes = len(data_array[0]) - 2  # Last two entries are columns 'T' and 'N'
        df = pd.DataFrame(data_array, columns=['BE' + str(n) for n in range(1, no_bes + 1)] + ['T', 'N'])

        if config.debug:
            # Read MCS from matlab for later validity check
            mcs_array = data_dict['mcs']
            df_mcs_matlab = pd.DataFrame(mcs_array, columns=['BE' + str(n) for n in range(1, no_bes + 1)])
    elif file_extension == ".csv":
        # Read from CSV file
        df = pd.read_csv(args.file, sep=',', header=0)
    else:
        logger.error("Can not parse file extension {}".format(file_extension))
        exit(1)
    end_parse = timer()
    time_parse = end_parse - start

    # Prepare data
    # --------------------
    # Drop last row containing totals
    if df.iloc[-1, 0] == 'Total':
        df.drop(df.tail(1).index, inplace=True)
    if 'Total' in df.index:
        df.drop(labels='Total', axis=0, inplace=True)

    # Drop column 'Cut set(s)' if present
    if 'Cut set' in df.columns:
        df.drop(labels='Cut set', axis=1, inplace=True)
    elif 'Cut sets' in df.columns:
        df.drop(labels='Cut sets', axis=1, inplace=True)

    # Shuffle data
    # --------------------
    if not args.debug:
        df = df.sample(frac=1, axis=0).reset_index(drop=True)  # Shuffle rows
        df = df.sample(frac=1, axis=1).reset_index(drop=True)  # Shuffle columns

    # Set order of cut sets if not present
    if 'Total' not in df.columns:
        df = df.assign(Total=df.sum(axis=1))
    # Sort by row sum as given in column 'Total'.
    df = df.sort_values(by='Total', ascending=True)
    # Drop column containing totals
    df.drop(labels='Total', axis=1, inplace=True)

    # Add top-level column 'T' if not present.
    # In this case we assume that all rows correspond to system failure
    if 'T' not in df.columns:
        df['T'] = 1

    # Add column 'N' (number of observations) if not present
    # In this case we assume that each observation occurs exactly once
    if 'N' not in df.columns:
        df['N'] = 1

    # Create the dataset for later evaluation
    dataset_evaluation = df.to_dict(orient="records")

    # Only keep rows corresponding to system failure
    df = df[df['T'] > 0]

    # Drop columns 'T' and 'N'
    df.drop(labels=['T', 'N'], axis=1, inplace=True)

    # Create mapping from index to BE names
    bes = dict()
    index = 0
    for col in df.columns:
        bes[index] = col
        index += 1
    logger.info("BEs: {}".format(bes))

    # Convert to numpy matrix
    cutset_matrix = df.to_numpy()
    logger.debug("Cut set matrix:\n{}".format(cutset_matrix))

    # Extract MCS from data
    # --------------------
    mcss = MinCutSets()
    mcss.compute_from_cut_sets(cutset_matrix)
    logger.info("MCSs: {}".format(mcss))

    # Check validity of MCS and data if input was matlab file
    if config.debug and file_extension == ".mat":
        mcss_matlab = MinCutSets()
        mcss_matlab.compute_from_cut_sets(df_mcs_matlab.to_numpy())
        while mcss_matlab:
            cut_set = mcss_matlab.pop()
            if cut_set not in mcss:
                logger.error("MCS {} is not contained in data set".format(cut_set))
                exit(2)
        if mcss_matlab:
            logger.error("MCSs {} are not contained in data set".format(mcss_matlab))
            exit(2)

    # Obtain symmetries from file
    # --------------------
    symmetries = None
    if args.symmetries:
        # Parse symmetries from file
        logger.debug("Parsing symmetries from {}".format(args.symmetries))
        symmetries = get_symmetries_from_file(args.symmetries, bes)
        logger.info("Symmetries from file:\n{}".format(symmetries))
        # Check validity of symmetries
        if not symmetries.is_valid_symmetries(mcss, print_cex=True):
            logger.error("Symmetries are not valid!")
            exit(3)
        if not config.use_symmetries:
            logger.info("Symmetries were disabled and will not be used.")

    end_prepare = timer()

    # Learn fault trees recursively
    # --------------------
    results = Results(filename + file_extension, config, bes)
    results.no_mcss = len(mcss)
    results.no_data = len(dataset_evaluation)
    results.time_parse = time_parse
    results.time_prepare = end_prepare - end_parse

    start_learn = timer()
    results.ft = recursive_learning(mcss, bes, bes.copy(), config, results=results, dataset_evaluation=dataset_evaluation, symmetries=symmetries, recurse_level=1)
    end_learn = timer()
    results.time_learn_complete = end_learn - start_learn

    # Simplify fault tree
    try:
        results.ft_simplified = results.ft.copy()
        results.ft_simplified.simplify()
    except:
        pass

    # Evaluate performance on inferred fault tree and print results
    # --------------------
    results.compute_evaluation(mcss, dataset_evaluation)
    end = timer()
    results.time_evaluation = end - end_learn
    results.time_total = end - start
    
    results.print()
    if args.result_dir:
        matlab_file = os.path.join(args.result_dir, 'results-{}-{}.mat'.format(filename, datetime.datetime.now().strftime("%Y%m%d-%H%M%S")))
        results.to_matlab(matlab_file)
