import argparse
import logging
import pandas as pd
import numpy as np

import ft_learn.helper
import ft_learn.ft.fault_tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate data (for learning) from given fault tree.')

    parser.add_argument('input', metavar='fault-tree-file', type=str, help='Input file describing the fault tree.')
    parser.add_argument('output', metavar='output-file', type=str, help='Output file containing the generated sample set.')
    parser.add_argument('--no-samples', '-n', type=int, help='Number of sample points.', default=100000)
    parser.add_argument('--verbose', '-v', help='Enable verbose output', action='store_true')
    parser.add_argument('--debug', '-d', help='Enable debugging', action='store_true')

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)

    no_samples = args.no_samples

    # Parse fault tree
    # --------------------
    with open(args.input, 'r') as f:
        # Will only read first line
        ft_string = f.readline()
    logger.debug("Parsed fault tree: {}".format(ft_string))

    # Failure probability in the BEs (equal for all of them)
    # TODO: let user define in input file
    failure_prob = 0.5

    ft = ft_learn.ft.fault_tree.str2ft(ft_string)
    logger.info("Use fault tree: {}".format(ft))

    # Generating samples for statistically independent dataset
    # --------------------
    # Initialize random number generator
    if args.debug:
        # Use fixed seed for debugging
        rng = np.random.default_rng(42)
    else:
        rng = np.random.default_rng()

    max_be = len(ft.get_all_bes())
    be_names = sorted([str(be) for be in ft.get_all_bes()])

    # Sample data points according to Binomial distribution
    logger.info("Generate {} samples...".format(no_samples))
    data = rng.binomial(1, failure_prob, no_samples).reshape((no_samples, 1))
    for i in range(1, max_be):
        # Going column-wise, because row-wise is orders of magnitude slower
        data = np.hstack((data, np.random.binomial(1, failure_prob, no_samples).reshape((no_samples, 1))))
    assert data.shape[0] == no_samples
    assert data.shape[1] == max_be

    # Create a pandas dataset from numpy array
    df = pd.DataFrame(data=data, columns=be_names)
    # Get evaluations for each BE according to given dataset
    evaluations = df.to_dict(orient="records")
    # Perform evaluation on fault tree
    T = []
    for i in evaluations:
        if ft.evaluate(i):
            T.append(1)
        else:
            T.append(0)
    # Add evaluation results
    df['T'] = T
    logger.info("Generated {} samples".format(no_samples))

    # Keep only unique rows and add the number of observations for each row
    data = df.to_numpy()
    unique, counts = np.unique(data, axis=0, return_counts=True)
    df = pd.DataFrame(data=unique, columns=be_names + ['T'])
    df['N'] = counts
    logger.debug("Final dataset: {}".format(df))
    logger.info("Obtained {} unique cut sets".format(df.shape[0]))

    # Save the dataset
    df.to_csv(args.output, index_label="Cut sets")
    logger.info("Saved dataset to {}".format(args.output))
