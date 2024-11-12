#!/usr/bin/env python3

import argparse
import os
import requests


def main(ds: int = 0):
    doi = '10.18419/darus-4225'
    if ds == 0:
        download_dir = 'jet-collision'
    elif ds == 1:
        download_dir = 'jet-collision-ds1'
    elif ds == 2:
        download_dir = 'jet-collision-ds2'
    else:
        raise RuntimeError('Unknown dataset version')

    # Get file list
    r = requests.get('https://darus.uni-stuttgart.de/api/datasets/:persistentId/?persistentId=doi:{}'.format(doi))
    dataset_info = r.json()

    # Sort by directory
    directories = {}
    for f in dataset_info['data']['latestVersion']['files']:
        dirname = f['directoryLabel']
        if dirname not in directories:
            directories[dirname] = []
        directories[dirname].append((f['dataFile']['id'], f['dataFile']['filename']))

    # Download
    os.makedirs(download_dir, exist_ok=True)
    for fileid, filename in directories[download_dir]:
        print(filename)
        with requests.get('https://darus.uni-stuttgart.de/api/access/datafile/{}'.format(fileid), stream=True) as r, \
            open(os.path.join(download_dir, filename), 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                f.write(chunk)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--ds', type=int, choices=range(0, 3), default=0)
    args = parser.parse_args()
    main(args.ds)
