#!/bin/env python3
import os
import logging
from collections import MutableMapping, MutableSet
from subprocess import PIPE
import re
from fractions import Fraction
import numpy as np
from sarge import Command, Capture
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
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
        # read the '*' that indicates the prompt so they don't build up
        this_line = self.read_iso_line()
        # logger.debug("reading *: {}".format(this_line))
        self.iso_process.stdin.write(bytes(command + "\n", "ascii"))
        self.iso_process.stdin.flush()

    def getDisplayData(self, display, raw=False):
        """
        Args:
        display: what we are asking isotropy to display
            the command sent will be 'DISPLAY {display}'
        raw: if true return a string of the raw output from isotropy
            otherwise the output is automaticly parsed in to a list of dictionaries
        """
        self.sendCommand("DISPLAY {}".format(display))
        lines = []
        keep_reading = True
        while keep_reading:
            this_line = self.read_iso_line()
            if this_line in ['*', '']:  # if there is no output '' is returned above
                keep_reading = False
            elif re.match(".*Data base for these coupled subgroups .*", this_line):
                self.read_iso_line() # read past Should this
                self.read_iso_line() # read past Enter RETURN
                self.sendCommand("")
                self.read_iso_line() # read past Adding
                self.read_iso_line() # read past Blank
                #self.read_iso_line() # read past Blank
            else:
                lines.append(this_line)
        if not raw:
            return self._parse_output(lines)
        return lines

    def read_iso_line(self):
        raw = self.iso_process.stdout.readline().decode()
        this_line = raw.rstrip('\n')
        logger.debug("isotropy: {}".format(this_line))
        return this_line

    def _parse_output(self, lines):
        indexes = detect_column_indexes(lines)
        split_by_ind = [split_line_by_indexes(indexes, line) for line in lines]
        parsed_output = [{key: detect_data_form_and_convert(prop)
                          for key, prop in result.items()}
                         for result in detect_multirows_and_split(split_by_ind)]
        return parsed_output

def detect_data_form_and_convert(prop):
    # if it is a list operate on each element
    if isinstance(prop, list):
        return [detect_data_form_and_convert(p) for p in prop]
    # first split by '|'s (but not '|'s inside paren)
    pipe_split_list = re.split(r'[|]\s*(?![^()]*\))', prop)
    if len(pipe_split_list) > 1:
        return detect_data_form_and_convert(pipe_split_list)
    # first split by commas (but not commas inside paren)
    comma_split_list = re.split(r'[,]\s*(?![^()]*\))', prop)
    if len(comma_split_list) > 1:
        return detect_data_form_and_convert(comma_split_list)
    # remove paren if entirely surrounded in paren with no inner paren
    # return list of paren surrouned bits if there are multiple
    if re.match(r'^\(.*\)$', prop):
        surrounded_by_paren_list = re.findall(r'\((.+?)\)', prop)
        if len(surrounded_by_paren_list) > 1:
            return detect_data_form_and_convert(surrounded_by_paren_list)
        return detect_data_form_and_convert(surrounded_by_paren_list[0])
    # next split by spaces
    space_split_list = prop.split()
    if len(space_split_list) > 1:
        return detect_data_form_and_convert(space_split_list)
    # we leave numbers as strings for ease of giving them
    # back to isotropy (which often wants fractions)
    # they can be converted where needed

    # finally just remove outer whitespace
    return prop.strip()

def detect_multirows_and_split(split_lines):
    result_list = []
    for row in split_lines[1:]:
        if row[0]:
            result = {}
            result_list.append(result)
        for j, prop in enumerate(row):
            if row[0]:
                result[split_lines[0][j]] = prop
            elif prop:
                if isinstance(result[split_lines[0][j]], list):
                    result[split_lines[0][j]].append(prop)
                else:
                    result[split_lines[0][j]] = [result[split_lines[0][j]], prop]
    return result_list

def detect_column_indexes(list_of_lines):
    indexes = [0]
    transitions = [col.count(' ') == len(list_of_lines) for col in zip(*list_of_lines)]
    last = False
    for i, x in enumerate(transitions):
        if (not x
            and last
            and list_of_lines[0][i] != ' '):
            #and all(line[i] != ' ' for line in list_of_lines[1:])):
            # the above commented condition breaks Matricies (where the actual matrix is indented past header)
            # but fixes cases where both the header label has spaces and is longer than the actual data
            # if this indentation only happens for matricies we can apply a fix to that special case
            # example above condition is useful:
            # shows = ['irrep', 'kpoint'], then getDisplayData('irrep')
            indexes.append(i)
        last = x
    return indexes

def split_line_by_indexes(indexes, line):
    tokens = []
    for i1, i2 in zip(indexes[:-1], indexes[1:]): #pairs
        tokens.append(line[i1:i2].rstrip())
    tokens.append(line[indexes[-1]:].rstrip())
    return tokens


def getSymOps(spacegroup, with_matrix=False, lattice_param='1 1 1 90 90 90', setting=None):
    values = {'parent': spacegroup}
    shows = ['elements']
    with IsotropySession(values, shows, setting=setting) as isos:
        if with_matrix:
            isos.values['lattice parameter'] = lattice_param
            isos.shows.add('cartesian')
            symOps = isos.getDisplayData('parent')
        else:
            symOps = isos.getDisplayData('parent')[0]['Elements']
    return symOps

def getKpoints(spacegroup, setting=None):
    values = {'parent': spacegroup}
    shows = ['kpoint']
    with IsotropySession(values, shows, setting=setting) as isos:
        kpoints = isos.getDisplayData('kpoint')
        kpt_dict = {kpt['']: tuple(kpt['k vector']) for kpt in kpoints}
    return kpt_dict

def getIrreps(spacegroup, kpoint=None, setting=None):
    values = {'parent': spacegroup}
    if kpoint:
        values['kpoint'] = kpoint
    shows = ['irrep']
    with IsotropySession(values, shows, setting=setting) as isos:
        results = isos.getDisplayData('irrep')
        irreps = [ir['Irrep (ML)'] for ir in results]
    return irreps

def getRepresentations(spacegroup, kpoint_label, irreps=None, setting=None):
    elements = getSymOps(spacegroup, setting)
    if not irreps:
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
                                          ' '.join(element[1]))
                isos.values['element'] = elem_str
                res = isos.getDisplayData('irrep')[0]
                matrix = res['Matrix'] # unlike previous version leaving elements as strings
                mat_list.append(matrix)
            irrep_dict[irrep] = mat_list
    return irrep_dict

def getDistortion(parent, wyckoffs, irrep, direction=None, cell=None):
    values = {'parent': parent,
              'wyckoff': ' '.join(wyckoffs),
              'irrep': irrep,}
    if direction:
        values['direction'] = direction
    if cell:
        values['cell'] = _matrix_to_iso_string(cell)
    shows = ['wyckoff', 'microscopic vector']
    with IsotropySession(values, shows) as isos:
        dist = isos.getDisplayData('DISTORTION', raw=False)
    # if there is one projected vector for each point we may want to
    # change this so that it is a list of length 1 so data is in same
    # form as cases where there are multiple vectors for each point
    # this currently is not done
    return dist

def _matrix_to_iso_string(mat):
    return ' '.join([','.join([str(i) for i in r]) for r in mat])

def _to_float(number):
    try:
        return float(number)
    except ValueError:
        return float(Fraction(number))

def _list_to_float_array(sl):
    """return float array of (possibly nested) list of strings (or ints/floats) where
    strings can have fractions (i.e. 1/2 will be converted to 0.5)"""
    result = []
    for i in sl:
        if isinstance(i, (list, np.ndarray)):
            result.append(_list_to_float_array(i))
        else:
            result.append(_to_float(i))
    return np.array(result, dtype=float)


def _find_all_equivalent_basis_origin(parent, basis, origin):
    basis = _list_to_float_array(basis)
    origin = _list_to_float_array(origin)
    # we follow a convention that all origin choices are positive with each componenet < 1
    origin = np.array([i % 1 for i in origin])
    symOps = getSymOps(parent, with_matrix=True)
    possible_basis_origins = []
    for symop in symOps:
        rot = _list_to_float_array(symop['Rotation matrix, translation'][:3])
        trans = _list_to_float_array(symop['Rotation matrix, translation'][3])
        new_basis = np.array([np.dot(rot, vec) + trans for vec in basis])
        new_origin = np.dot(rot, origin) + trans
        # we follow a convention that all origin choices are positive with each componenet < 1
        new_origin = np.array([i % 1 for i in new_origin])
        possible_basis_origins.append((new_basis, new_origin))
    # commented line below would remove duplicates, but it might not be worth it
    # no_dupes = list({np.array(bo[0] + bo[1]).tostring(): bo
    #                      for bo in possible_basis_origins}.values())
    return possible_basis_origins


def getPossibleSingleIrrepOPs(parent, subgroup): #, basis=None, origin=None):
    values = {'parent': parent, 'subgroup': subgroup}
    shows = ['irrep', 'direction', 'basis', 'origin']
    with IsotropySession(values, shows) as isos:
        possible_ops = isos.getDisplayData('ISOTROPY')
    return possible_ops

def getPossibleOPs_for_basis(parent, subgroup, basis, origin):
    ops_to_check = getPossibleSingleIrrepOPs(parent, subgroup)
    # TODO: add coupled irreps to ops_to_check
    equivalent_basis = _find_all_equivalent_basis_origin(parent, basis, origin)
    compatible_ops = []
    for op in ops_to_check:
        this_basis = _list_to_float_array(op['Basis Vectors'])
        this_origin = _list_to_float_array(op['Origin']) # assumed to have all 0<x_i<1 for each i
        for b, o in equivalent_basis:
            if (abs(this_basis - b) < 1e-5).all() and (abs(this_origin - o) < 1e-5).all():
                compatible_ops.append(op)
                break
    return compatible_ops



#  def getPossiblePhaseTransitions(struct_hs, struct_ls):
#     """
#     Given two pymatgen structure objects (high symmetry and low symmetry)
#     find all irreps describing an order parameter of the phase transition between these structures
#
#     returns a 'phase transition' which is just a dict containing:
#             'structure_high_sym'
#             'structure_low_sym'
#             'parent'
#             'subgroup'
#             'irrep'
#             'direction'
#             'basis'
#             'origin'
#             'secondary_OPs' - other order parameters which do not lower symmetry more than the primary OP
#                               a dict of 'irrep', 'direction', 'domain', 'frequency'
#     """
#     sga_hs = SpacegroupAnalyzer(struct_hs)
#     sga_ls = SpacegroupAnalyzer(struct_ls)
#     struct_hs_p = sga_hs.get_primitive_standard_structure()
#     struct_ls_p = sga_ls.get_primitive_standard_structure()
#     sgn_hs = sga_hs.get_space_group_number()
#     sgn_ls = sga_ls.get_space_group_number()
#     # not limiting results by ratio since the size in isotropy doesn't account for changes of unit cell
#     # seems like it might only allow for diagonal basis?? NO! even tutorial has examples with nondiagonal basis
#     # ratio = len(struct_ls_p) / len(struct_hs_p)
#     # if ratio % 1 != 0:
#     #     raise ValueError(("Number of sites in low symmetry structure must be",
#     #                       " an integer times the number of sites in high symmetry structure"))
#     # ratio = int(ratio)
#     # values = {'parent': sgn_hs, 'subgroup': sgn_ls, 'size': ratio}
#     values = {'parent': sgn_hs, 'subgroup': sgn_ls}
#     shows = ['irrep', 'direction', 'basis', 'origin']
#     possible_transitions = []
#     with IsotropySession(values, shows) as isos:
#         lines = isos.getDisplayData('ISOTROPY')
#         for line in lines[1:-1]:
#             irrep, direction, basis, origin = line.split()[:4]
#             this_transition = {'structure_high_sym': struct_hs_p,
#                                'structure_low_sym': struct_ls_p,
#                                'parent': sgn_hs,
#                                'subgroup': sgn_ls,
#                                'irrep': irrep,
#                                'direction': direction,
#                                'basis': to_array(basis),
#                                'origin': to_array(origin)}
#             isos.shows.clearAll()
#             isos.shows.add('frequency direction')
#             isos.values['irrep'] = irrep
#             freq_output = isos.getDisplayData('ISOTROPY')
#             del isos.values['irrep']
#             raw_sec_ops = freq_output[1].split(',')
#             sec_ops = []
#             for sop in raw_sec_ops:
#                 freq, s_irrep, raw_s_dir = sop.split()
#                 freq = int(freq)
#                 s_dir, s_domain = raw_s_dir.rstrip(')').split('(')
#                 sec_ops.append({'irrep': s_irrep,
#                                 'direction': s_dir,
#                                 'frequency': freq,
#                                 'domain': s_domain})
#             this_transition['secondary_OPs'] = sec_ops
#             possible_transitions.append(this_transition)
#     return possible_transitions
#
#
#  def getAllowedMicroDistortions(phase_transition):
#      """
#      Given a 'phase transition' (likely from the getPossiblePhaseTransitions function)
#      find allowed displacive modes for all primary and secondary order parameters
#
#      returns a list where each element corresponds to an order parameter
#      each element of this list is a tuple where the first element is a dict describing the order parameter
#      and the second element is a list of modes
#      each mode is a tuple of the form (wyckoff label, [points], [(displacement vectors)])
#      the actual distortions can then by projected on to the displacement vectors to
#      separate the contributions from each mode
#      displacement vectors is a tuple since multidimensional irreps will have more than one at each point
#      """
#      sga_hs = SpacegroupAnalyzer(phase_transition['structure_high_sym'])
#      wyckoffs = ' '.join(set(sga_hs.get_symmetry_dataset()['wyckoffs']))
#      values = {'parent': phase_transition['parent'],
#                'wyckoff': wyckoffs}
#      shows = ['wyckoff', 'microscopic vector']
#      dists = []
#      with IsotropySession(values, shows) as isos:
#          for s_op in phase_transition['secondary_OPs']:
#              isos.values['irrep'] = s_op['irrep']
#              isos.values['direction'] = s_op['direction']
#              raw_dist_out = isos.getDisplayData('DISTORTION')
#              distortions = []
#              first_dist = True
#              for line in raw_dist_out:
#                  if re.match('.*Wyckoff Point.*', line) or re.match(r'\*\*+', line) or line == '':
#                      pass
#                  # elif line[0] in wyckoffs.split():
#                  elif re.match('[a-z]', line[0]):
#                      if not first_dist:
#                          distortions.append(this_dist)
#                      wyck_lbl, point= line.split()[:2]
#                      proj_vecs = line.split()[2:]
#                      this_dist = (wyck_lbl, [to_array(point)], [[to_array(pv.rstrip(',')) for pv in proj_vecs]])
#                      first_dist = False
#                  elif line == '*':
#                      distortions.append(this_dist)
#                  else:
#                      point = line.split()[0]
#                      proj_vecs = line.split()[1:]
#                      this_dist[1].append(to_array(point))
#                      this_dist[2].append([to_array(pv.rstrip(',')) for pv in proj_vecs])
#              if len(distortions) == 0:
#                  pass
#              else:
#                  dists.append((s_op, distortions))
#      return dists
#
#
#def to_array(ar_str):
#    as_mat = [[mm
#               for mm in v.split(',')]
#              for v in ar_str.rstrip(')').lstrip('(').split('),(')]
#    if len(as_mat) == 1:
#        result = as_mat[0]
#    else:
#        result = as_mat
#    return result


if __name__ == '__main__':
    # implement argparse if this main part is ever for more than testing
    import sys
    stream_handler = logging.StreamHandler()
    if len(sys.argv) > 1:
        if sys.argv[1] == 'd':
            stream_handler.setLevel(logging.DEBUG)
    else:
        stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    sg = 221
    logger.info(getSymOps(sg, with_matrix=True))
    #logger.info(getKpoints(sg))
    #logger.info(getIrreps(sg))
    #logger.info(getRepresentations(sg,
    #                               list(getKpoints(sg).keys())[0],
    #                               irreps=['GM5+']))
    #logger.info(getDistortion(sg, 'a b c', 'R4-'))
