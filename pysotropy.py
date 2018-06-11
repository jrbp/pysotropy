#!/bin/env python3
from __future__ import print_function
import os
from collections import MutableMapping, MutableSet
from subprocess import PIPE
from sarge import Command, Capture
import re
"""
Very rough interface to isotropy
"""


class Shows(MutableSet):
    def __init__(self, parent, initial_shows):
        self._shows = set()
        self.parent = parent
        if initial_shows:
            for item in initial_shows:
                item = item.upper()
                self.add(item)

    def __contains__(self, item):
        item = item.upper()
        return item in self._shows

    def __iter__(self):
        return iter(self._shows)

    def __len__(self):
        return len(self._shows)

    def add(self, item):
        item = item.upper()
        if item not in self._shows:
            self.parent.sendCommand("SHOW {}".format(item))
            self._shows.add(item)

    def discard(self, item):
        item = item.upper()
        try:
            self._shows.remove(item)
            self.parent.sendCommand("CANCEL SHOW {}".format(item))
        except ValueError:
            pass

    def clearAll(self):
        self.parent.sendCommand("CANCEL SHOW ALL")
        self._shows = set()


class Values(MutableMapping):
    """
    Acts like a dictionary for values set in isotropy,
    when values are set and deleted the appropriate calls
    to the IsotropySession are made
    """
    def __init__(self, parent, initial_values):
        '''Use the object dict'''
        self.parent = parent
        self._vals = dict()
        if initial_values:
            for k, v in initial_values.items():
                k = k.upper()
                self.__setitem__(k, v)

    def __setitem__(self, key, value):
        key = key.upper()
        self.parent.sendCommand("VALUE {} {}".format(key, value))
        self._vals[key] = value

    def __getitem__(self, key):
        key = key.upper()
        return self._vals[key]

    def __delitem__(self, key):
        key = key.upper()
        self.parent.sendCommand("CANCEL VALUE {}".format(key))
        del self._vals[key]

    def __iter__(self):
        return iter(self._vals)

    def __len__(self):
        return len(self._vals)

    def clearAll(self):
        self.parent.sendCommand("CANCEL VALUE ALL")
        self._vals = dict()


class IsotropySession:
    """
    Make simple requests to isotropy.
    isotropy session is kept running in background until closed
    should be used with 'with' statements to ensure isotropy is exited properly
    ex:
    with IsotropySession() as isos:
        do things with isos
    """
    def __init__(self, values=None, shows=None,
                 labels=None, setting=None):
        """
        Args:
        values: dictionary of keys to be set to values
            key is a string specifying what is being set
                e.g. "basis", "cell", "irrep", "kpoint", "parent"
            value is what this key is being set to
        shows: A list of strings corresponding to data which will be returned
            when display is run.
            Note: some show commands accept additional parameters, for now
            these must be included in the string. Eventually parsing of
            display output should be good enough that they are not needed.
        labels: NOT YET IMPLEMENTED
            dictionary where the key corresponds to the object
            whose notation is being altered, the value corresponds
            to the new notation to be used for this object
            for example {"spacegroup": "SCHOENFLIES"} will cause returned
            results and entered values to use schoenflies notation
        setting:
            a string or list of strings to be passed to
            the setting command can be used to change settings, origin,
            unique axis and/or cell choice. Can also specify if magnetic
            spacegroups are desired
            for now the setting options can only be set
            when creating an Isotropy object, not changed later
        """
        iso_location = "/home/john/scripts/isobyu/"  # TODO: don't hard code
        self.iso_process = Command(os.path.join(iso_location, 'iso'),
                                   stdout=Capture(buffer_size=1),
                                   env={"ISODATA": iso_location})
        self.iso_process.run(input=PIPE, async=True)

        # move past initial output
        keep_reading = True
        while keep_reading:
            this_line = self.iso_process.stdout.readline().decode()
            if this_line == 'Use "VALUE IRREP VERSION" to change version\n':
                keep_reading = False

        self.screen = 80  # exploit this too make parsing output easier?
        self.sendCommand("SCREEN {}".format(self.screen))
        self.page = "NOBREAK"
        self.sendCommand("PAGE {}".format(self.page))
        if setting:
            if type(setting) == list:
                self.setting = setting
            else:
                self.setting = [setting]
        else:
            self.setting = ["INTERNATIONAL"]
        for s in self.setting:
            self.sendCommand("SETTING {}".format(s))
        self.values = Values(self, values)
        self.shows = Shows(self, shows)

    def __enter__(self):
        return self

    def __exit__(self, exec_type, exc_value, exc_traceback):
        self.sendCommand("QUIT")

    def sendCommand(self, command):
        self.iso_process.stdin.write(bytes(command + "\n", "ascii"))
        self.iso_process.stdin.flush()

    def getDisplayData(self, display, debug=False):
        self.sendCommand("DISPLAY {}".format(display))
        lines = []
        keep_reading = True
        while keep_reading:
            this_line = self.iso_process.stdout.readline().decode()
            if this_line == '*':
                keep_reading = False
            lines.append(this_line)
        return lines


def getSymOps(spacegroup, setting=None):
    values = {"parent": spacegroup}
    shows = ["elements"]
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData("parent")
        symOps = []
        for line in lines:
            symOps += re.findall(r'\(([A-Za-z0-9]*)\|([0-9,/]*)\)', line)
    return symOps


def getKpoints(spacegroup, setting=None):
    values = {"parent": spacegroup}
    shows = ["kpoint"]
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData("kpoint")
        kpoints = []
        for line in lines:
            kpoints += re.findall(r'([A-Z][A-Z]?)\s*\((.*)\)', line)
        kpt_dict = {label: tuple([p for p in loc.split(",")])
                    for label, loc in kpoints}
    return kpt_dict


def getIrreps(spacegroup, kpoint=None, setting=None):
    values = {"parent": spacegroup}
    if kpoint:
        values["kpoint"] = kpoint
    shows = ["irrep"]
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData("irrep")
        irreps = []
        for line in lines:
            irreps += re.findall(r'^([A-Z].*)\n', line)
    return irreps


def getRepresentations(spacegroup, kpoint_label, setting=None):
    elements = getSymOps(spacegroup, setting)
    irreps = getIrreps(spacegroup, kpoint_label, setting)
    values = {"parent": spacegroup, "kpoint": kpoint_label}
    shows = ["matrix"]
    irrep_dict = {}
    with IsotropySession(values, shows, setting=setting) as isos:
        for irrep in irreps:
            isos.values["irrep"] = irrep
            mat_list = []
            for element in elements:
                elem_str = "{} {}".format(element[0],
                                          element[1].replace(",", " "))
                isos.values["element"] = elem_str
                lines = isos.getDisplayData("irrep")
                temp_lines = []
                temp_lines.append(re.match(r'\(.*\)\s*(.*)',
                                           lines[1]).groups()[0])
                for line in lines[2:-1]:
                    temp_lines.append(re.match(r'\s*(.*)', line).groups()[0])
                matrix = [[float(i) for i in r.split()]
                          for r in temp_lines if not r == '']
                mat_list.append(matrix)
            irrep_dict[irrep] = mat_list
    return irrep_dict


if __name__ == '__main__':
    spacegroup = 221
    print(getSymOps(spacegroup))
    print(getKpoints(spacegroup))
    print(getIrreps(spacegroup))
    print(getRepresentations(spacegroup,
                             list(getKpoints(spacegroup).keys())[0]))
