#!/bin/env python3
from __future__ import print_function
import os
import subprocess
import re
"""
Very rough interface to isotropy
"""


class IsotropyDisplay:
    """
    Make simple requests to isotropy.
    Minimal implementation capable of making arbitrary requests
    Should eventually be redone such that an isotropy session stays open
    so that everything runs much quicker
    """
    def __init__(self, values=None, shows=None,
                 labels=None, setting=None):
        """
        Args:
        values: list of tuples of the form (key, value)
            key is a string specifying what is being set
                e.g. "basis", "cell", "irrep", "kpoint", "parent"
            value is what this key is being set to
        shows: A list of strings corresponding to data which will be returned
            when display is run.
            Note: some show commands accept additional parameters, for now
            these must be included in the string. Eventually parsing of
            display output should be good enough that they are not needed.
        labels: tuple of two strings (key, notation)
            key corresponds to the object whose notation is being altered
            notation corresponds to the new notation to be used for this object
            for example ("spacegroup", "SCHOENFLIES") will cause returned
            results and entered values to use schoenflies notation
        setting:
            for now this is just a string or list of strings to be passed to
            the setting command can be used to change settings, origin,
            unique axis and/or cell choice. Can also specify if magnetic
            spacegroups are desired.
        """
        iso_location = "/home/john/scripts/isobyu/"
        self.iso_command = ("ISODATA={} ".format(iso_location)
                       + os.path.join(iso_location, 'iso'))
        self.values = values
        self.shows = shows
        if setting:
            if type(setting) == list:
                self.setting = setting
            else:
                self.setting = [setting]
        else:
            self.setting = ["INTERNATIONAL"]
        self.page = "NOBREAK"
        self.screen = 100000  # exploit this too make parsing output easier?

        self.commands = ["SCREEN {}".format(self.screen),
                         "PAGE {}".format(self.page)]
        for s in self.setting:
            self.commands.append("SETTING {}".format(s))
        for value in self.values:
            self.commands.append("VALUE {} {}".format(*value))
        for show in self.shows:
            self.commands.append("SHOW {}".format(show))

    def getDisplayData(self, display, debug=False):
        cmd = (self.iso_command + " <<<'" + "\n"
               + "\n".join(self.commands) + "\n"
               + "DISPLAY {}".format(display) + "\n"
               + "QUIT" + "\n"
               + "\n'")
        output = subprocess.check_output(['bash', '-c', cmd]).decode()
        if debug:
            print(output)
        result_lines = []
        in_result_block = False
        for line in output.split("\n"):
            if in_result_block:
                result_lines.append(line)
            if line[0] == "*":
                in_result_block = not in_result_block
            result = "\n".join(result_lines[:-1])
        return result


def getSymOps(spacegroup, setting=None):
    values = [("parent", spacegroup)]
    shows = ["elements"]
    id = IsotropyDisplay(values, shows, setting=setting)
    res = id.getDisplayData("parent")
    symOps = [(s + ")").replace(",", " ").replace("|", " ").replace("(", "").replace(")", "") for s in res.strip("\n").split("),")]
    return symOps


def getKpoints(spacegroup, setting=None):
    values = [("parent", spacegroup)]
    shows = ["kpoint"]
    id = IsotropyDisplay(values, shows, setting=setting)
    res = id.getDisplayData("kpoint")
    matches = re.findall(r'([A-Z][A-Z]?)\s*\((.*)\)', res)
    kpoints = {lbl: tuple([p for p in loc.split(",")]) for lbl, loc in matches}
    return kpoints


def getIrreps(spacegroup, kpoint=None, setting=None):
    values = [("parent", spacegroup)]
    if kpoint:
        values.append(("kpoint", kpoint))
    shows = ["irrep"]
    id = IsotropyDisplay(values, shows, setting=setting)
    res = id.getDisplayData("irrep")
    return res.split("\n")


def getRepresentation(spacegroup, irrep, element, setting=None):
    values = [("parent", spacegroup), ("irrep", irrep), ("element", element)]
    shows = ["matrix"]
    id = IsotropyDisplay(values, shows, setting=setting)
    res = id.getDisplayData("irrep")
    raw_lines = res.split("\n")
    lines = []
    lines.append(re.match(r'\(.*\)\s*(.*)', raw_lines[0]).groups()[0])
    for l in raw_lines[1:]:
        lines.append(re.match(r'\s*(.*)', l).groups()[0])
    matrix = [[float(i) for i in r.split()] for r in lines if not r == '']
    return matrix


if __name__ == '__main__':
    iso_location = "/home/john/scripts/isobyu/"
    iso_command = ("ISODATA={} ".format(iso_location)
                   + os.path.join(iso_location, 'iso'))
    spacegroup = 221
    elements = getSymOps(spacegroup)
    print(elements)
    print()
    kpoints = getKpoints(spacegroup)
    print(kpoints)
    print()
    kptlist = [k for k in kpoints]
    irreps = getIrreps(spacegroup, kptlist[0])
    print(irreps)
    print()
    for irr in irreps:
        print(irr)
        for element in elements:
            print(getRepresentation(spacegroup, irr, element))
