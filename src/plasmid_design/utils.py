from tqdm.auto import tqdm
import joblib
import pandas as pd
import yaml
from string import Formatter


def assert_unique(df, *cols):
    """Each argument can be a single column or a list of columns.
    """
    a = df.shape[0]
    for col in cols:
        b = ~df[col].duplicated()
        if a != b.sum():
            counts = df[col].value_counts()
            raise ValueError(
                f'{b.sum()} / {a} entries are unique for column(s) {col}, '
                f'the most frequent duplicate is {counts.index[0]} ({counts.iloc[0]} entries)'
            )
    return df


def load_yaml(filename):
    with open(filename, 'r') as fh:
        return yaml.safe_load(fh)


class ProgressParallel(joblib.Parallel):
    """From SO 37804279
    """
    def __call__(self, *args, **kwargs):
        self.total = len(args)
        with tqdm(total=self.total) as self._pbar:
            return joblib.Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        self._pbar.total = self.n_dispatched_tasks
        # not sure why this doesn't work
        # self._pabar.total = self.total 
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def load_yaml_table(config, verbose=True):
    """Load a table from a YAML description.
    TODO: include drive: option
    """
    if isinstance(config, str):
        config = {'table': config}
    
    load_keys = 'skiprows',
    kwargs = {k: config[k] for k in load_keys if k in config}
    if config['table'].startswith('drive:'):
        from .drive import Drive
        drive = Drive()
        remote = config['table'][len('drive:'):]
        df = drive(remote, **kwargs)
    elif config['table'].startswith('drive_key:'):
        _, sheet_name = config['table'].split(':')
        from .drive import read_csv_from_url
        df = read_csv_from_url(config['drive_key'], sheet_name, **kwargs)
    else:
        df = pd.read_csv(config['table'], low_memory=False, **kwargs)

    if verbose:
        print(f'Loaded {len(df):,} entries from {config["table"]}')
    return filter_yaml_table(df, verbose=verbose, **config)


def filter_yaml_table(df, gate=None, drop_duplicates=None, rename=None, verbose=True, **ignore):
    """Apply some common transformations (e.g., from a YAML description).
    """
    if gate is not None:
        df = df.query(gate)
        if verbose:
            print(f'  Kept {len(df):,} passing gate: {gate}')
    if drop_duplicates is not None:
        df = df.drop_duplicates(drop_duplicates)
        if verbose:
            print(f'  {len(df):,} after dropping duplicates on {drop_duplicates}')
    if rename is not None:
        df = df.rename(columns=rename)
    return df


def format_string_any_field(fmt, fields):
    """Fill in a python format string with invalid field names 
    (i.e., containing any character).
    """
    s = ''
    for prefix, stuffing, _, _ in Formatter().parse(fmt):
        prefix = '' if prefix is None else prefix
        stuffing = '' if stuffing is None else fields[stuffing]
        s += (prefix + stuffing)
    return s
