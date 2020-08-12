#!/usr/bin/env python3


"""This file should be a temporary measure.
  There are certain changes to pymatgen needed by the code not yet present in pymatgen.
  Instead of insisting that a fork of pymatgen must be used certain class methods are patched here.
  Hopefully these changes make their way upstream and this file can be removed."""

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord import lattice_points_in_supercell
from pymatgen.optimization.linear_assignment import LinearAssignment  # type: ignore
from pymatgen.core.structure import Structure


def _get_lattices(self, target_lattice, s, supercell_size=1, rh_only=False):
    """
    Yields lattices for s with lengths and angles close to the
    lattice of target_s. If supercell_size is specified, the
    returned lattice will have that number of primitive cells
    in it

    Args:
        s, target_s: Structure objects
        rh_only (bool): whether or not to only return lattices with det(basis)>0
    """
    lattices = s.lattice.find_all_mappings(
        target_lattice, ltol=self.ltol, atol=self.angle_tol,
        skip_rotation_matrix=True)
    for l, _, scale_m in lattices:
        is_rh = np.linalg.det(scale_m) > 0
        allowed = is_rh or not rh_only
        if abs(abs(np.linalg.det(scale_m)) - supercell_size) < 0.5 and allowed:
            yield l, scale_m


def _get_supercells(self, struct1, struct2, fu, s1_supercell, rh_only=False):
    """
    Computes all supercells of one structure close to the lattice of the
    other
    if s1_supercell == True, it makes the supercells of struct1, otherwise
    it makes them of s2

    yields: s1, s2, supercell_matrix, average_lattice, supercell_matrix
    """

    def av_lat(l1, l2):
        params = (np.array(l1.parameters) +
                  np.array(l2.parameters)) / 2
        return Lattice.from_parameters(*params)

    def sc_generator(s1, s2):
        s2_fc = np.array(s2.frac_coords)
        if fu == 1:
            cc = np.array(s1.cart_coords)
            for l, sc_m in self._get_lattices(s2.lattice, s1, fu, rh_only=rh_only):
                fc = l.get_fractional_coords(cc)
                fc -= np.floor(fc)
                yield fc, s2_fc, av_lat(l, s2.lattice), sc_m
        else:
            fc_init = np.array(s1.frac_coords)
            for l, sc_m in self._get_lattices(s2.lattice, s1, fu, rh_only=rh_only):
                fc = np.dot(fc_init, np.linalg.inv(sc_m))
                lp = lattice_points_in_supercell(sc_m)
                fc = (fc[:, None, :] + lp[None, :, :]).reshape((-1, 3))
                fc -= np.floor(fc)
                yield fc, s2_fc, av_lat(l, s2.lattice), sc_m

    if s1_supercell:
        for x in sc_generator(struct1, struct2):
            yield x
    else:
        for x in sc_generator(struct2, struct1):
            # reorder generator output so s1 is still first
            yield x[1], x[0], x[2], x[3]


def _strict_match(self, struct1, struct2, fu, s1_supercell=True,
                  use_rms=False, break_on_match=False, rh_only=False):
    """
    Matches struct2 onto struct1 (which should contain all sites in
    struct2).

    Args:
        struct1, struct2 (Structure): structures to be matched
        fu (int): size of supercell to create
        s1_supercell (bool): whether to create the supercell of
            struct1 (vs struct2)
        use_rms (bool): whether to minimize the rms of the matching
        break_on_match (bool): whether to stop search at first
            valid match
        rh_only (bool): whether to only return structures with det(basis)>0
    """
    if fu < 1:
        raise ValueError("fu cannot be less than 1")

    mask, s1_t_inds, s2_t_ind = self._get_mask(struct1, struct2,
                                               fu, s1_supercell)

    if mask.shape[0] > mask.shape[1]:
        raise ValueError('after supercell creation, struct1 must '
                         'have more sites than struct2')

    # check that a valid mapping exists
    if (not self._subset) and mask.shape[1] != mask.shape[0]:
        return None

    if LinearAssignment(mask).min_cost > 0:
        return None

    best_match = None
    # loop over all lattices
    for s1fc, s2fc, avg_l, sc_m in \
            self._get_supercells(struct1, struct2, fu, s1_supercell, rh_only=rh_only):
        # compute fractional tolerance
        normalization = (len(s1fc) / avg_l.volume) ** (1 / 3)
        inv_abc = np.array(avg_l.reciprocal_lattice.abc)
        frac_tol = inv_abc * self.stol / (np.pi * normalization)
        # loop over all translations
        for s1i in s1_t_inds:
            t = s1fc[s1i] - s2fc[s2_t_ind]
            t_s2fc = s2fc + t
            if self._cmp_fstruct(s1fc, t_s2fc, frac_tol, mask):
                inv_lll_abc = np.array(avg_l.get_lll_reduced_lattice().reciprocal_lattice.abc)
                lll_frac_tol = inv_lll_abc * self.stol / (np.pi * normalization)
                dist, t_adj, mapping = self._cart_dists(
                    s1fc, t_s2fc, avg_l, mask, normalization, lll_frac_tol)
                if use_rms:
                    val = np.linalg.norm(dist) / len(dist) ** 0.5
                else:
                    val = max(dist)
                if best_match is None or val < best_match[0]:
                    total_t = t + t_adj
                    total_t -= np.round(total_t)
                    best_match = val, dist, sc_m, total_t, mapping
                    if (break_on_match or val < 1e-5) and val < self.stol:
                        return best_match

    if best_match and best_match[0] < self.stol:
        return best_match


def get_transformation(self, struct1, struct2, rh_only=False):
    """
    Returns the supercell transformation, fractional translation vector,
    and a mapping to transform struct2 to be similar to struct1.

    Args:
        struct1 (Structure): Reference structure
        struct2 (Structure): Structure to transform.
        rh_only (bool): whether to only return structures with det(basis)>0

    Returns:
        supercell (numpy.ndarray(3, 3)): supercell matrix
        vector (numpy.ndarray(3)): fractional translation vector
        mapping (list(int or None)):
            The first len(struct1) items of the mapping vector are the
            indices of struct1's corresponding sites in struct2 (or None
            if there is no corresponding site), and the other items are
            the remaining site indices of struct2.
    """
    if self._primitive_cell:
        raise ValueError("get_transformation cannot be used with the "
                         "primitive cell option")

    struct1, struct2 = self._process_species((struct1, struct2))

    s1, s2, fu, s1_supercell = self._preprocess(struct1, struct2, False)
    ratio = fu if s1_supercell else 1 / fu
    if s1_supercell and fu > 1:
        raise ValueError("Struct1 must be the supercell, "
                         "not the other way around")

    if len(s1) * ratio >= len(s2):
        # s1 is superset
        match = self._strict_match(s1, s2, fu=fu, s1_supercell=False,
                                   use_rms=True, break_on_match=False, rh_only=rh_only)
        if match is None:
            return None
        # invert the mapping, since it needs to be from s1 to s2
        mapping = [list(match[4]).index(i) if i in match[4] else None
                   for i in range(len(s1))]
        return match[2], match[3], mapping
    else:
        # s2 is superset
        match = self._strict_match(s2, s1, fu=fu, s1_supercell=True,
                                   use_rms=True, break_on_match=False, rh_only=rh_only)
        if match is None:
            return None
        # add sites not included in the mapping
        not_included = list(range(len(s2) * fu))
        for i in match[4]:
            not_included.remove(i)
        mapping = list(match[4]) + not_included
        return match[2], -match[3], mapping


def get_s2_like_s1(self, struct1, struct2, include_ignored_species=True, rh_only=False):
    """
    Performs transformations on struct2 to put it in a basis similar to
    struct1 (without changing any of the inter-site distances)

    Args:
        struct1 (Structure): Reference structure
        struct2 (Structure): Structure to transform.
        include_ignored_species (bool): Defaults to True,
            the ignored_species is also transformed to the struct1
            lattice orientation, though obviously there is no direct
            matching to existing sites.

    Returns:
        A structure object similar to struct1, obtained by making a
        supercell, sorting, and translating struct2.
    """
    s1, s2 = self._process_species([struct1, struct2])
    trans = self.get_transformation(s1, s2, rh_only=rh_only)
    if trans is None:
        return None
    sc, t, mapping = trans
    sites = [site for site in s2]
    # Append the ignored sites at the end.
    sites.extend([site for site in struct2 if site not in s2])
    temp = Structure.from_sites(sites)

    temp.make_supercell(sc)
    temp.translate_sites(list(range(len(temp))), t)
    # translate sites to correct unit cell
    for i, j in enumerate(mapping[:len(s1)]):
        if j is not None:
            vec = np.round(struct1[i].frac_coords - temp[j].frac_coords)
            temp.translate_sites(j, vec, to_unit_cell=False)

    sites = [temp.sites[i] for i in mapping if i is not None]

    if include_ignored_species:
        start = int(round(len(temp) / len(struct2) * len(s2)))
        sites.extend(temp.sites[start:])

    return Structure.from_sites(sites)


def PATCH_StructureMatcher():
    # method overrides
    StructureMatcher._get_lattices = _get_lattices
    StructureMatcher._get_supercells = _get_supercells
    StructureMatcher._strict_match = _strict_match
    StructureMatcher.get_transformation = get_transformation
    StructureMatcher.get_s2_like_s1 = get_s2_like_s1
