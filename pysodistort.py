#!/usr/bin/env python
import sys
import logging
from fractions import Fraction
import numpy as np
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
import pysotropy as iso
from sympy import sympify, linsolve, EmptySet
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application)
logger = logging.getLogger("pysotropy")
logger.setLevel(logging.DEBUG)
transformations = (standard_transformations + (implicit_multiplication_application,))


def frac_vec_convert(vec, lat1, lat2):
    """convert from frac coords of lat 1 to frac coords of lat2"""
    cart = lat1.get_cartesian_coords(vec)
    return lat2.get_fractional_coords(cart)

def smallest_disp(s2, s1):
    sign = lambda x: 1 if x > 0 else -1 if x < 0 else 0
    disp  = []
    for d2, d1 in zip(s2, s1):
        d = d2 - d1
        if abs(d) > 0.5:
            d = d - sign(d)
        disp.append(d)
    return disp

def get_sym_info(struct):
    """get spacegroup number and wyckoff set"""
    sga = SpacegroupAnalyzer(struct)
    sgn = sga.get_space_group_number()
    wyckoff = sga.get_symmetry_dataset()['wyckoffs']
    return sgn, wyckoff

def match_structures(s1, s2, scale_lattice=False, rh_only=True):
    """
    Args
        s1: high sym structure
        s2: low sym structure
        scale_lattice (optional): high_sym_superlcell has same lattice vectors as s2 (no strain)
    Returns
        basis: should be the basis that when applied to s1 makes a supercell of the size and orentation of s2
        origin: any additional translation to best match (applied before applying the basis change to match what isotropy does)
        displacements
        high_sym_supercell
"""
    sm = StructureMatcher(attempt_supercell=True, primitive_cell=False)
    basis, origin, mapping = sm.get_transformation(s1, s2, rh_only=rh_only)

    struct_hs_supercell = sm.get_s2_like_s1(s1, s2, rh_only=rh_only)

    # change origin from the supercell basis to the high sym basis
    origin = np.round_(frac_vec_convert(origin,
                                        struct_hs_supercell.lattice,
                                        s2.lattice),
                       decimals=5)
    if scale_lattice:
        struct_hs_supercell = pmg.Structure(s1.lattice,
                                            struct_hs_supercell.species,
                                            [site.frac_coords for site in struct_hs_supercell])
    displacements = []
    for s_hs, s_ls in zip(struct_hs_supercell, s1):
        disp = np.round_(smallest_disp(s_hs.frac_coords, s_ls.frac_coords), decimals=5)
        displacements.append(disp)
    return basis, origin, displacements, struct_hs_supercell

def get_all_distortions(sgn_hs, wyckoff_list, directions, basis, origin):
    directions_dict = {}
    distortions = {}
    with iso.IsotropySession() as isos:
        for direct in directions:
            if iso._kpt_has_params(direct['k vector']):
                logger.warning("low sym kpoint not implemented yet skipping {}".format(direct))
                #TODO: this can be easily implented by setting kvalue correctly using the value by getDirections
                continue
            irrep = direct['Irrep']
            d = "vector,{}".format(','.join(direct['Dir']))
            this_dist = iso.getDistortion(sgn_hs, wyckoff_list,
                                          irrep, cell=basis, origin=origin, direction=d, isos=isos)
            if len(this_dist) > 0:
                distortions[irrep] = this_dist
                directions_dict[irrep] = direct['Dir']
    return distortions, directions_dict

def convert_distortions_basis(distortions, origin,
                              lat1, lat2):
    irreps = {}
    for irrep, wycks in distortions.items():
        irreps[irrep] = []
        for wyck in wycks:
            wyck_sc = {"Wyckoff": wyck["Wyckoff"],
                       "Point": [],
                       "Projected Vectors": []}
    
            # need to check if we have only one site
            if type(wyck["Point"][0]) is not list:
                wyck["Point"] = [wyck["Point"]]
                wyck["Projected Vectors"] = [wyck["Projected Vectors"]]
    
            # also need to check if only 1 proj vector for each site (should we change last step?)
            if type(wyck["Projected Vectors"][0][0]) is not list:
                wyck["Projected Vectors"] = [[pv] for pv in wyck["Projected Vectors"]]
    
            for pt, vcs in zip(wyck["Point"], wyck["Projected Vectors"]):
                pt = np.array([float(Fraction(i)) for i in pt]) + origin
                wyck_sc["Point"].append(list(np.round_(frac_vec_convert(pt, lat1, lat2),
                                                       decimals=5)))
                this_sites_vcs = []
                for vc in vcs:
                    vc = np.array([float(Fraction(i)) for i in vc])
                    vc_sc_basis = list(np.round_(frac_vec_convert(vc, lat1, lat2),
                                                 decimals=5))
                    this_sites_vcs.append(vc_sc_basis)
                wyck_sc["Projected Vectors"].append(this_sites_vcs)
            irreps[irrep].append(wyck_sc)
    return irreps


# TODO: possibly clean this up now that we really only use this for one wyckoff at a time
def get_distortion_dec_struct(wycks, struct_to_match, high_sym_wyckoff, struct_hs):
    coords = []
    species = []
    proj_vecs = []
    wycks_done = []  # dirty maybe wrong fix
    for wyck in wycks:  # CURRENTLY ONLY ONE WYCKOFF IS PASSED AT A TIME, SO THIS IS FINE, BUT UNNECESSARY
        w = wyck["Wyckoff"]
        # PEROVSKITE SPECIFIC, AND I DON'T KNOW IF IT'S THE RIGHT THING TO DO!!!!
        if w in wycks_done:
            # print("SKIPPING duplicate instance of wyckoff {} in irrep {}, this is a perovskite specific thing, and may be wrong even here".format(w, irrep))
            # print("I believe the two sets isotropy returns are equivalent choices, but I'm not certain")
            continue
        wycks_done.append(w)  # dirty maybe wrong fix

        for i, ss in enumerate(high_sym_wyckoff):
            if ss == w:
                sp = struct_hs[i].specie
                break
        for coord, pv in zip(wyck["Point"], wyck["Projected Vectors"]):
            species.append(sp)
            coords.append(coord)
            proj_vecs.append(pv)
    logger.debug("get_distortion_dec_struct:")
    logger.debug("coords:\n{}".format(coords))
    logger.debug("species:\n{}".format(species))
    logger.debug("proj_vecs:\n{}".format(proj_vecs))
    lat = struct_to_match.lattice
    dist_struct = pmg.Structure(lat, species, coords, site_properties={"projvecs": proj_vecs})
    logger.debug("dist_struct:\n{}".format(dist_struct))

    sm_dist = StructureMatcher(ltol = 0.02, primitive_cell=False, allow_subset=True)
    try:
        sc_d, trans_d, mapping = sm_dist.get_transformation(struct_to_match, dist_struct)
    except TypeError:
        logger.warning("couldn't map dist decorated structure to actual structure")
        logger.warning("dist dec struct:\n{}".format(dist_struct))
        logger.warning("struct to match:\n{}".format(struct_to_match))
    dist_struct_matched = dist_struct * sc_d
    dist_struct_matched.translate_sites(list(range(len(dist_struct_matched))), trans_d)
    return dist_struct_matched, mapping

def get_projection_data(displacements, wycks, struct_hs_supercell, high_sym_wyckoff, struct_hs):
    results_by_wyck = {}
    for n, wyck in enumerate(wycks):
        num_proj_vecs = len(wycks[0]["Projected Vectors"][0])
        amplitudes = [0. for i in range(num_proj_vecs)]
        amplitude_as_comps = [0. for i in range(num_proj_vecs)]
        dist_struct_matched, mapping = get_distortion_dec_struct([wyck], struct_hs_supercell, high_sym_wyckoff, struct_hs)
        dist_defs = dist_struct_matched
        full_projvecs = []
        for i, j in enumerate(mapping):
            if j is not None:
                pv = dist_struct_matched[j].properties["projvecs"]
                full_projvecs.append(pv)
            else:
                full_projvecs.append([[0., 0., 0.] for i in range(num_proj_vecs)])
        for i in range(num_proj_vecs):
            sum_cart_squares = 0.
            for pv in full_projvecs:
                pv_cart = struct_hs_supercell.lattice.get_cartesian_coords(pv[i])
                sum_cart_squares += pv_cart.dot(pv_cart)
            norm_factor = sum_cart_squares**(-1/2)
            for disp, pv in zip(displacements, full_projvecs):
                amplitudes[i] += np.dot(disp, pv[i])
                disp_cart = struct_hs_supercell.lattice.get_cartesian_coords(disp)
                pv_cart = struct_hs_supercell.lattice.get_cartesian_coords(pv[i])
                amplitude_as_comps[i] += np.dot(norm_factor * disp_cart, pv_cart)
            amplitude_as = np.sqrt(np.sum([am**2 for am in amplitude_as_comps]))
            amplitude_ap = amplitude_as * np.sqrt(struct_hs.lattice.volume / struct_hs_supercell.lattice.volume)

        results_by_wyck['{}{}'.format(wyck['Wyckoff'], n)] = {
            'amplitude_as': amplitude_as,
            'amplitude_ap': amplitude_ap,
            'amplitudes': amplitudes,
            'dist_defs': dist_defs,
            'full_projvecs': full_projvecs,
            'num_proj_vecs': num_proj_vecs,
            'total_amplitude': np.sqrt(np.sum([a**2
                                               for a in amplitudes]))}
    return results_by_wyck


def get_amps_direction(parent, irrep, irrep_amp, isos=None):
    if isos is None:
        close_after = True
        isos = iso.IsotropySession()
    else:
        close_after = False
    isos.values.update({'parent': parent, 'irrep': irrep})
    isos.shows.update({'direction', 'subgroup'})
    sym_inequiv = isos.getDisplayData('ISOTROPY')
    isos.shows.clearAll()
    inequiv_dir_labels = [s['Dir'] for s in sym_inequiv]
    irrep_domains = {}
    for lbl in inequiv_dir_labels:
        these_domains = iso.getDomains(parent, irrep, lbl, isos=isos)
        these_domains[0]['Dir'] = these_domains[0]['Dir'][-1] # hacky fix
        irrep_domains[lbl] = these_domains
    isos.shows.clearAll()
    isos.values.clearAll()
    if close_after:
        isos.__exit__(None, None, None)

    for lbl, domains in irrep_domains.items():
        for domain in domains:
            syms = set()
            ddir = domain['Dir']
            eqn_set = []
            for pos_d, pos_a in zip(ddir, irrep_amp):
                pos_d_sym = parse_expr(pos_d, transformations=transformations)
                pos_a_sym = sympify(pos_a)
                syms.update(pos_d_sym.free_symbols)
                eqn_set.append(pos_d_sym - pos_a_sym)
            syms = list(syms)
            soln = linsolve(eqn_set, syms)
            if soln == EmptySet():
                continue
            eqns = [0 == eqn.subs([(sy, val)
                                   for sy, val in zip(syms, list(soln)[0])])
                    for eqn in eqn_set]
            if not all(eqns):
                continue
            var_vals = list(zip([str(s) for s in syms], list(soln)[0]))
            logger.debug(var_vals)
            return lbl, ddir, var_vals


def get_mode_decomposition(struct_hs, struct_ls, nonzero_only=False, general_direction=True):
    """
    Args
        struct_hs: high symmetry structure (pymatgen structure)
        struct_ls: low symmetry structure (pymatgen structure)
        nonzero_only (optional): only return modes which have nonzero amplitude TODO: make this a cutoff
    Returns
        dict containing mode_decomposition_data
        {irrep:
               {wyckoff:
                        {'amplitudes': [amp_0, ...], 'total_amplitude: amp',
                         'dist_defs': pymatgen structure with dist as props,
                         'direction': []
                         }
                ...},
          ...}
    """
    # in general need to use value wyckoff xyz if there are free parameters
    # not needed for perovskites here
    sgn_hs, wyckoff_list = get_sym_info(struct_hs)
    if general_direction:
        sgn_ls = None
    else:
        sgn_ls = get_sym_info(struct_ls)[0]

    basis, origin, displacements, struct_hs_supercell = match_structures(struct_ls, struct_hs)

    try:
        directions = iso.getDirections(sgn_hs, basis, origin, subgroup=sgn_ls)
    except iso.IsotropyBasisException:
        # isotropy is picky about some things
        # rh_only option now in match_structures should make this unnecessary
        basis = np.dot(-1 * np.identity(3), basis)
        directions = iso.getDirections(sgn_hs, basis, origin, subgroup=sgn_ls)
        logger.warning("trying with inverted basis {}".format(basis))
    except iso.IsotropySubgroupException:
        # this should also never happen now that we start with geeral directions
        # and sort the rest out later
        logger.warning("Isotropy isn't recognizing the subgroup relation in this basis")
        logger.warning("trying more general projections")
        logger.warning("Perhaps double check this")
        directions = iso.getDirections(sgn_hs, basis, origin)

    logger.info("Basis: \n{}".format(basis))
    logger.info("Origin: {}".format(origin))
    logger.debug("Origin using: {}".format([str(Fraction(i).limit_denominator(10)) for i in origin]))

    this_subgroup_distortions, directions_dict = get_all_distortions(sgn_hs, list(set(wyckoff_list)),
                                                                     directions, basis, origin)

    all_in_sc_basis = convert_distortions_basis(this_subgroup_distortions, origin,
                                                struct_hs.lattice,
                                                struct_hs_supercell.lattice)

    irrep_dist_full_pvecs = {}
    irrep_dist_defs = {}
    irrep_amplitudes = {}

    mode_decomposition_data = {}
    with iso.IsotropySession() as isos:
        for irrep, wycks in all_in_sc_basis.items():
            proj_data_by_wyck = get_projection_data(displacements, wycks,
                                                    struct_hs_supercell, wyckoff_list, struct_hs)
            # TODO: clean this up, also don't find directions for irreps with amp==0
            # should move a lot of the direction stuff to its own function
            for wyck in proj_data_by_wyck.keys():
                this_amp = proj_data_by_wyck[wyck]['amplitudes']
                # temporary work around until cleaned up properly
                if np.sum(np.abs(this_amp)) < 1e-6:
                    proj_data_by_wyck[wyck]['direction'] = ('zero', [0.,0.,0.])
                    continue
                syms = []
                amp_sym = []
                for el in directions_dict[irrep]:
                    el_sym = parse_expr(el, transformations=transformations)
                    amp_sym.append(el_sym)
                    for s in el_sym.free_symbols:
                        if s not in syms:
                            syms.append(s)
                if len(syms) != len(this_amp):
                    logger.warning("WARNING: irrep {} wyck {} has different number of params then amp components".format(irrep, wyck))
                sym_val_pairs = [(sym, val) for sym, val in zip(syms, this_amp)]
                this_amp_conv = [round(float(el_sym.subs(sym_val_pairs)), 4) for el_sym in amp_sym]
                logger.info("{}  {}".format(irrep, wyck))
                dir_lbl, dir_vec, var_vals = get_amps_direction(sgn_hs, irrep, this_amp_conv, isos=isos)
                proj_data_by_wyck[wyck]['direction'] = (dir_lbl, dir_vec)
                # TODO: really seperate components properly and give seperate amplitudes for each freee param
                # currently As and Ap will only be right with one free param
                # otherwise they will be both mixed together and sign will be always positive
                if len(var_vals) == 1:
                    sign = float(var_vals[0][1] / abs(var_vals[0][1]))
                    proj_data_by_wyck[wyck]['amplitude_as'] = sign * proj_data_by_wyck[wyck]['amplitude_as']
                    proj_data_by_wyck[wyck]['amplitude_ap'] = sign * proj_data_by_wyck[wyck]['amplitude_ap']
                proj_data_by_wyck[wyck]['param_vals'] = var_vals
            mode_decomposition_data[irrep] = proj_data_by_wyck
    if nonzero_only:
        nonzero_mode_decomp = {}
        for irrep, wycks in mode_decomposition_data.items():
            tot_amp = 0.
            for wyck, data in wycks.items():
                 tot_amp += np.sum(np.abs(data['amplitudes']))
            if tot_amp > 1e-5:
                nonzero_mode_decomp[irrep] = wycks
        return nonzero_mode_decomp
    return mode_decomposition_data

if __name__ == '__main__':
    stream_handler = logging.StreamHandler()
    if len(sys.argv) > 3:
        if sys.argv[-1] == 'd':
            stream_handler.setLevel(logging.DEBUG)
    else:
        stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    struct_hs = pmg.Structure.from_file(sys.argv[1])
    struct_ls = pmg.Structure.from_file(sys.argv[2])

    irrep_decomposition_data = get_mode_decomposition(struct_hs, struct_ls, nonzero_only=True)

    logger.info("Mode Definitions:")
    for irrep, wycks in irrep_decomposition_data.items():
        logger.info(irrep)
        for wyck, data in wycks.items():
            logger.info("\t{}".format(wyck))
            logger.info("\t\t{}".format(data["direction"]))
            logger.info("\t\t{}".format(data["dist_defs"]))
    logger.info("Mode Amplitudes")
    for irrep, wycks in irrep_decomposition_data.items():
        logger.info(irrep)
        for wyck, data in wycks.items():
            logger.info('\t{}'.format(wyck))
            logger.info('\t\t{}'.format(data["direction"]))
            logger.info('\t\t{}'.format(np.round_(data["amplitudes"], decimals=5)))
            logger.info('\t\tAs: {}'.format(np.round_(data["amplitude_as"], decimals=5)))
            logger.info('\t\tAp: {}'.format(np.round_(data["amplitude_ap"], decimals=5)))
        # logger.info(np.round_(data["Amplitudes_as"], decimals=5))
