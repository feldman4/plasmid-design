"""
Generate service file following
https://cloud.google.com/iam/docs/creating-managing-service-account-keys

Place it in the working directory, a parent directory, or your home directory.
"""
from glob import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pygsheets


BASE_URL = 'https://docs.google.com/spreadsheets/d/{key}'
CSV_EXPORT_URL = f'{BASE_URL}/gviz/tq?tqx=out:csv&sheet={{sheet_name}}'


def find_service_file(search='*service*.json'):
    """Searches current path, then parent directories, then user home.
    """
    parent_paths = [os.getcwd()] + [str(x) for x in Path(os.getcwd()).parents]
    # ~/bii/cloud kept for cross-compatibility
    special_paths = [os.environ['HOME'], f'{os.environ["HOME"]}/bii/cloud/']
    paths = parent_paths + special_paths

    matches = []
    for path in paths:
        matches.extend(glob(os.path.join(path, search)))
    
    if len(matches) == 0:
        raise FileNotFoundError(f'no file matching {search}')
    
    return matches[0]


def read_csv_from_url(key, sheet_name, skiprows=0, **kwargs):
    url = CSV_EXPORT_URL.format(key=key, sheet_name=sheet_name)
    try:
        df = pd.read_csv(url, skiprows=skiprows)
    except pd.errors.ParserError as e:
        base_url = BASE_URL.format(key=key)
        extra_info = (
            '\n\nDid you remember to make this sheet public? '
            f'\n{base_url}'
            '\nClick Share > General access > Anyone with the link')
        raise pd.errors.ParserError(str(e) + extra_info)
    return Drive.clean(df, **kwargs)


class Drive():
    def __init__(self):
        if not os.path.exists(SERVICE_FILE):
            SERVICE_FILE = find_service_file()
        
        self.service = pygsheets.authorize(service_file=SERVICE_FILE)
        
    def get_excel(self, name, dropna='all', normalize=True, fix_int=True, 
                  drop_unnamed=True, header=True, skiprows=0):
        """
        """
        if len(name.split('/')) == 2:
            spreadsheet_title, worksheet_title = name.split('/')

        sh = self.service.open(spreadsheet_title)
        ws = sh.worksheet_by_title(worksheet_title)

        # TODO: fix this
        start = f'A{int(skiprows) + 1}'
        df = ws.get_as_df(has_header=header, numerize=True, start=start)
        return self.clean(df, dropna=dropna, normalize=normalize, fix_int=fix_int, 
                     drop_unnamed=drop_unnamed)
    
    @staticmethod
    def clean(df, dropna='all', normalize=True, fix_int=True, drop_unnamed=True, fix_value=True):
        if fix_value:
            df[df == '#VALUE!'] = np.nan
        if dropna:
            df = df.dropna(how=dropna, axis=0)
            df = df.dropna(how=dropna, axis=1)
        if normalize:
            df.columns = [normalize_col_name(x) for x in df.columns]
        if fix_int:
            for c in df.columns:
                try:
                    if (df[c] - df[c].astype(int)).sum() == 0:
                        df[c] = df[c].astype(int)
                except:
                    pass
        if drop_unnamed:
            cols_drop = [c for c in df.columns if str(c).startswith('Unnamed:')]
            df = df.drop(cols_drop, axis=1)
        return df

    def __call__(self, *args, **kwargs):
        return self.get_excel(*args, **kwargs)


def normalize_col_name(s):
    try:
        s = s.replace('# of', 'num')
        s = s.replace('#', 'num')
        s = s.replace('\n', ' ')
        s = s.replace(' / ', '_per_')
        s = s.replace('/', '_per_')
        s = s.replace(' ', '_')
        s = s.replace('-', '_')
        
    except AttributeError: # not a string
        pass
    return s