#!/bin/env python3
import os
import logging
from collections import MutableMapping, MutableSet
from subprocess import PIPE
from sarge import Command, Capture
import re
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
"""
Python interface to isotropy
"""

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

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
        logger.debug("""starting isotropy session in {}
                        using isotropy in: {}""".format(
                            os.getcwd(), iso_location))
        self.iso_process = Command(os.path.join(iso_location, 'iso'),
                                   stdout=Capture(buffer_size=1),
                                   env={"ISODATA": iso_location})
        self.iso_process.run(input=PIPE, async=True)

        # move past initial output
        keep_reading = True
        while keep_reading:
            this_line = self.iso_process.stdout.readline().decode()
            if this_line: # don't log until isotropy responds
                logger.debug("isotropy: {}".format(this_line))
            if this_line == 'Use "VALUE IRREP VERSION" to change version\n':
                keep_reading = False

        self.screen = 999  # exploit this too make parsing output easier?
        self.sendCommand("SCREEN {}".format(self.screen))
        #self.page = "NOBREAK" # still feels the need to periodicly put in labels
        self.page = "999"
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
        logger.debug("sending: {}".format(command))
        self.iso_process.stdin.write(bytes(command + "\n", "ascii"))
        self.iso_process.stdin.flush()

    def getDisplayData(self, display, debug=False):
        self.sendCommand("DISPLAY {}".format(display))
        lines = []
        keep_reading = True
        while keep_reading:
            this_line = self.iso_process.stdout.readline().decode()
            logger.debug("isotropy: {}".format(this_line))
            if this_line in ['*', '']:  # if there is no output '' is returned above
                keep_reading = False
            elif re.match(".*Data base for these coupled subgroups .*", this_line):
                self.iso_process.stdout.readline() # read past Should this
                self.iso_process.stdout.readline() # read past Enter RETURN
                self.sendCommand("")
                self.iso_process.stdout.readline() # read past Adding
                self.iso_process.stdout.readline() # read past Blank
                #self.iso_process.stdout.readline()) # read past Blank
                continue
            #lines.append(this_line.lstrip('*')) #breaks the microdsitortion function atm
            lines.append(this_line)
        return lines


def getSymOps(spacegroup, setting=None):
    values = {'parent': spacegroup}
    shows = ['elements']
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData('parent')
        symOps = []
        for line in lines:
            symOps += re.findall(r'\(([A-Za-z0-9]*)\|([0-9,/]*)\)', line)
    return symOps


def getKpoints(spacegroup, setting=None):
    values = {'parent': spacegroup}
    shows = ['kpoint']
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData('kpoint')
        kpoints = []
        for line in lines:
            kpoints += re.findall(r'([A-Z][A-Z]?)\s*\((.*)\)', line)
        kpt_dict = {label: tuple([p for p in loc.split(',')])
                    for label, loc in kpoints}
    return kpt_dict


def getIrreps(spacegroup, kpoint=None, setting=None):
    values = {'parent': spacegroup}
    if kpoint:
        values['kpoint'] = kpoint
    shows = ['irrep']
    with IsotropySession(values, shows, setting=setting) as isos:
        lines = isos.getDisplayData('irrep')
        irreps = []
        for line in lines:
            irreps += re.findall(r'^([A-Z].*)\n', line)
    return irreps


def getRepresentations(spacegroup, kpoint_label, setting=None):
    elements = getSymOps(spacegroup, setting)
    irreps = getIrreps(spacegroup, kpoint_label, setting)
    values = {'parent': spacegroup, 'kpoint': kpoint_label}
    shows = ['matrix']
    irrep_dict = {}
    with IsotropySession(values, shows, setting=setting) as isos:
        for irrep in irreps:
            isos.values['irrep'] = irrep
            mat_list = []
            for element in elements:
                elem_str = '{} {}'.format(element[0],
                                          element[1].replace(',', ' '))
                isos.values['element'] = elem_str
                lines = isos.getDisplayData('irrep')
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


def getPossiblePhaseTransitions(struct_hs, struct_ls):
    """
    Given two pymatgen structure objects (high symmetry and low symmetry)
    find all irreps describing an order parameter of the phase transition between these structures

    returns a 'phase transition' which is just a dict containing:
            'structure_high_sym'
            'structure_low_sym'
            'parent'
            'subgroup'
            'irrep'
            'direction'
            'basis'
            'origin'
            'secondary_OPs' - other order parameters which do not lower symmetry more than the primary OP
                              a dict of 'irrep', 'direction', 'domain', 'frequency'
    """
    sga_hs = SpacegroupAnalyzer(struct_hs)
    sga_ls = SpacegroupAnalyzer(struct_ls)
    struct_hs_p = sga_hs.get_primitive_standard_structure()
    struct_ls_p = sga_ls.get_primitive_standard_structure()
    sgn_hs = sga_hs.get_space_group_number()
    sgn_ls = sga_ls.get_space_group_number()
    # not limiting results by ratio since the size in isotropy doesn't account for changes of unit cell
    # seems like it might only allow for diagonal basis?? NO! even tutorial has examples with nondiagonal basis
    # ratio = len(struct_ls_p) / len(struct_hs_p)
    # if ratio % 1 != 0:
    #     raise ValueError(("Number of sites in low symmetry structure must be",
    #                       " an integer times the number of sites in high symmetry structure"))
    # ratio = int(ratio)
    # values = {'parent': sgn_hs, 'subgroup': sgn_ls, 'size': ratio}
    values = {'parent': sgn_hs, 'subgroup': sgn_ls}
    shows = ['irrep', 'direction', 'basis', 'origin']
    possible_transitions = []
    with IsotropySession(values, shows) as isos:
        lines = isos.getDisplayData('ISOTROPY')
        for line in lines[1:-1]:
            irrep, direction, basis, origin = line.split()[:4]
            this_transition = {'structure_high_sym': struct_hs_p,
                               'structure_low_sym': struct_ls_p,
                               'parent': sgn_hs,
                               'subgroup': sgn_ls,
                               'irrep': irrep,
                               'direction': direction,
                               'basis': to_array(basis),
                               'origin': to_array(origin)}
            isos.shows.clearAll()
            isos.shows.add('frequency direction')
            isos.values['irrep'] = irrep
            freq_output = isos.getDisplayData('ISOTROPY')
            del isos.values['irrep']
            raw_sec_ops = freq_output[1].split(',')
            sec_ops = []
            for sop in raw_sec_ops:
                freq, s_irrep, raw_s_dir = sop.split()
                freq = int(freq)
                s_dir, s_domain = raw_s_dir.rstrip(')').split('(')
                sec_ops.append({'irrep': s_irrep,
                                'direction': s_dir,
                                'frequency': freq,
                                'domain': s_domain})
            this_transition['secondary_OPs'] = sec_ops
            possible_transitions.append(this_transition)
    return possible_transitions


def getAllowedMicroDistortions(phase_transition):
    """
    Given a 'phase transition' (likely from the getPossiblePhaseTransitions function)
    find allowed displacive modes for all primary and secondary order parameters

    returns a list where each element corresponds to an order parameter
    each element of this list is a tuple where the first element is a dict describing the order parameter
    and the second element is a list of modes
    each mode is a tuple of the form (wyckoff label, [points], [(displacement vectors)])
    the actual distortions can then by projected on to the displacement vectors to
    separate the contributions from each mode
    displacement vectors is a tuple since multidimensional irreps will have more than one at each point
    """
    sga_hs = SpacegroupAnalyzer(phase_transition['structure_high_sym'])
    wyckoffs = ' '.join(set(sga_hs.get_symmetry_dataset()['wyckoffs']))
    values = {'parent': phase_transition['parent'],
              'wyckoff': wyckoffs}
    shows = ['wyckoff', 'microscopic vector']
    dists = []
    with IsotropySession(values, shows) as isos:
        for s_op in phase_transition['secondary_OPs']:
            isos.values['irrep'] = s_op['irrep']
            isos.values['direction'] = s_op['direction']
            raw_dist_out = isos.getDisplayData('DISTORTION')
            distortions = []
            first_dist = True
            for line in raw_dist_out:
                if re.match('.*Wyckoff Point.*', line) or re.match(r'\*\*+', line) or line == '':
                    pass
                # elif line[0] in wyckoffs.split():
                elif re.match('[a-z]', line[0]):
                    if not first_dist:
                        distortions.append(this_dist)
                    wyck_lbl, point= line.split()[:2]
                    proj_vecs = line.split()[2:]
                    this_dist = (wyck_lbl, [to_array(point)], [[to_array(pv.rstrip(',')) for pv in proj_vecs]])
                    first_dist = False
                elif line == '*':
                    distortions.append(this_dist)
                else:
                    point = line.split()[0]
                    proj_vecs = line.split()[1:]
                    this_dist[1].append(to_array(point))
                    this_dist[2].append([to_array(pv.rstrip(',')) for pv in proj_vecs])
            if len(distortions) == 0:
                pass
            else:
                dists.append((s_op, distortions))
    return dists


def to_array(ar_str):
    as_mat = [[mm
               for mm in v.split(',')]
              for v in ar_str.rstrip(')').lstrip('(').split('),(')]
    if len(as_mat) == 1:
        result = as_mat[0]
    else:
        result = as_mat
    return result


if __name__ == '__main__':
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    spacegroup = 221
    logger.info(getSymOps(spacegroup))
    logger.info(getKpoints(spacegroup))
    logger.info(getIrreps(spacegroup))
    logger.info(getRepresentations(spacegroup,
                             list(getKpoints(spacegroup).keys())[0]))
