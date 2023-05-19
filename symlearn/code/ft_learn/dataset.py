import pandas as pd
import logging


def reduce_dataset(dataset, bes):
    """
    Reduce given dataset to only contain given BEs.
    :param dataset: Dataset.
    :param bes: BEs to use as dict {index: name}
    :return: Reduced dataset.
    """
    # Convert to pandas
    df = pd.DataFrame.from_records(dataset)
    logging.debug("Reduce the following dataset for BEs {}:\n{}".format(bes, df))

    # Drop columns with BEs not present in this module
    be_columns_keep = [be for be in df.columns.values if be in bes]
    be_columns_drop = [be for be in df.columns.values if (be not in bes and be != 'N' and be != 'T')]
    df.drop(labels=be_columns_drop, axis=1, inplace=True)

    # Merge rows with duplicate entries (resulting from removing some BEs)
    # The top event 'T' is given by an AND over all values (represented by the product), i.e. true iff all duplicate rows are true
    # The count 'N' is summed up
    df = df.groupby(be_columns_keep).agg({'T': 'prod', 'N': 'sum'}).reset_index()
    logging.debug("Reduced dataset:\n{}".format(df))

    # Convert to dict again
    return df.to_dict(orient="records")
