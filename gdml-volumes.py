#!/usr/bin/env python
# Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
# See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
"""
Count the number of logical and physical volumes in a GDML file.
"""

from itertools import count
from collections import deque
import json
import sys
import xml.etree.ElementTree as ET

def ilen(iter):
    # https://github.com/more-itertools/more-itertools/pull/242
    counter = count()
    deque(zip(iter, counter), maxlen=0)
    return next(counter)


def parse_gdml(filename):
    tree = ET.parse(filename)
    structure = next(tree.iter("structure"))

    physical = 0
    for logical, el in enumerate(structure):
        if el.tag in ('volume', 'assembly'):
            physical += ilen(el.iter("volumeref"))
        else:
            raise ValueError(f"Unrecognized structure tag: {el!r}")

    # Account for world volume
    logical += 1
    physical += 1

    return {'filename': filename, 'logical': logical, 'physical': physical}

def main(*args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__, prog="gdml-to-dot")
    parser.add_argument('-o', '--output', default='-')
    parser.add_argument('input', nargs='+')
    ns = parser.parse_args(*args)

    result = []
    for filename in ns.input:
        result.append(parse_gdml(filename))

    if ns.output == '-':
        json.dump(result, sys.stdout, indent=1)
    else:
        with open(ns.output, 'w') as f:
            json.dump(result, f, indent=0)

if __name__ == "__main__":
    main()
