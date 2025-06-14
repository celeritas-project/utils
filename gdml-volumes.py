#!/usr/bin/env python
# Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
# See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
"""
Count the number of logical and physical volumes in a GDML file.
"""

from itertools import count
from collections import namedtuple
import json
import sys
import xml.etree.ElementTree as ET


def get_material(volume_el):
    el = next(volume_el.iter("materialref"))
    return el.attrib["ref"]


def parse_materials(tree):
    materials = next(tree.iter("materials"))

    result = {
        "isotope": count(),
        "element": count(),
        "material": count(),
    }
    for el in materials:
        next(result[el.tag])

    return {k: next(v) for k, v in result.items()}


ChildCount = namedtuple("ChildCount", ["direct", "total", "maxdepth"])


def parse_structure(tree):
    structure = next(tree.iter("structure"))
    child_counts = {}

    # Counters
    border = 0
    skin = 0
    physical = 0

    for logical, el in enumerate(structure):
        if el.tag == "bordersurface":
            border += 1
            continue

        if el.tag == "skinsurface":
            skin += 1
            continue

        if el.tag not in ("volume", "assembly"):
            raise ValueError(f"Unrecognized structure tag: {el!r}")

        indirect = 1
        maxdepth = 0
        direct = count()
        for vrel in el.iter("volumeref"):
            cc = child_counts[vrel.attrib["ref"]]
            indirect += cc.total
            next(direct)
            maxdepth = max(maxdepth, cc.maxdepth)

        cc = ChildCount(direct=next(direct), total=indirect,
                        maxdepth=(maxdepth + 1))
        child_counts[el.attrib["name"]] = cc
        physical += cc.direct

    # Account for world volume
    logical += 1
    physical += 1
    world = tree.findall("./setup/world")[0]

    cc = child_counts[world.attrib["ref"]]

    return {"logical": logical, "physical": physical, "touchable": cc.total,
            "border": border, "skin": skin, "depth": cc.maxdepth}


def parse_gdml(filename):
    result = {"filename": filename}

    tree = ET.parse(filename)
    result.update(parse_materials(tree))
    result.update(parse_structure(tree))
    return result


def main(*args):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__, prog="gdml-volumes")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("input", nargs="+")
    ns = parser.parse_args(*args)

    result = []
    for filename in ns.input:
        try:
            result.append(parse_gdml(filename))
        except Exception as e:
            print(f"While processing {filename}: {e!s}")
            result.append(None)

    if ns.output == "-":
        json.dump(result, sys.stdout, indent=1)
    else:
        with open(ns.output, "w") as f:
            json.dump(result, f, indent=0)

if __name__ == "__main__":
    main()
