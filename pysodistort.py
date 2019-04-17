#!/usr/bin/env python
import sys
from fractions import Fraction
import numpy as np
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
import pysotropy as iso

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

def match_structures(s1, s2, scale_lattice=True):
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
    basis, origin, mapping = sm.get_transformation(s1, s2)

    struct_hs_supercell = sm.get_s2_like_s1(s1, s2)

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
    for direct in directions:
        if iso._kpt_has_params(direct['k vector']):
            print("low sym kpoint not implemented yet skipping")
            print(direct)
            #TODO: this can be easily implented by setting kvalue correctly using the value by getDirections
            continue
        irrep = direct['Irrep']
        d = "vector,{}".format(','.join(direct['Dir']))
        this_dist = iso.getDistortion(sgn_hs, wyckoff_list,
                                      irrep, cell=basis, origin=origin, direction=d)
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
    for wyck in wycks:
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
    # print(coords)
    # print(species)
    # print(proj_vecs)
    lat = struct_to_match.lattice
    dist_struct = pmg.Structure(lat, species, coords, site_properties={"projvecs": proj_vecs})
    # print(dist_struct)

    sm_dist = StructureMatcher(ltol = 0.02, primitive_cell=False, allow_subset=True)
    try:
        sc_d, trans_d, mapping = sm_dist.get_transformation(struct_to_match, dist_struct)
    except TypeError:
        print(dist_struct)
        print()
        print(struct_to_match)
    # print(transformation)
    dist_struct_matched = dist_struct * sc_d
    dist_struct_matched.translate_sites(list(range(len(dist_struct_matched))), trans_d)
    return dist_struct_matched, mapping

def get_projection_data(displacements, wycks, struct_hs_supercell, high_sym_wyckoff, struct_hs):
    results_by_wyck = {}
    for n, wyck in enumerate(wycks):
        num_proj_vecs = len(wycks[0]["Projected Vectors"][0])
        amplitudes = [0. for i in range(num_proj_vecs)]
        amplitudes_as = [0. for i in range(num_proj_vecs)]
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
            # norm_factor = (np.sqrt(np.sum([struct_hs_supercell.lattice.get_cartesian_coords(pv[i]).dot(struct_hs_supercell.lattice.get_cartesian_coords(pv[i]))
            #                                for pv in full_projvecs])))**-1
            for disp, pv in zip(displacements, full_projvecs):
                amplitudes[i] += np.dot(disp, pv[i])

                pv_np = norm_factor * np.array(pv[i], dtype=float)
                proj_disp_cart = struct_hs_supercell.lattice.get_cartesian_coords(np.dot(disp, pv_np) * pv_np)
                # print("displacement:\t{}".format(disp))
                # print("projvec:\t{}".format(pv_np))
                # print("pv_norm:\t{}".format(norm_factor))
                # print("proj_disp_cart:\t{}".format(proj_disp_cart))
                # print("distance:\t{}".format(proj_disp_cart.dot(proj_disp_cart)))
                # print()
                # amplitudes_as[i] += norm_factor * proj_disp_cart.dot(proj_disp_cart)
                amplitudes_as[i] += proj_disp_cart.dot(proj_disp_cart)
            amplitudes_as[i] += np.sqrt(amplitudes_as[i])
        results_by_wyck['{}{}'.format(wyck['Wyckoff'], n)] = {
            'amplitudes_as': amplitudes_as,
            'amplitudes': amplitudes,
            'dist_defs': dist_defs,
            'full_projvecs': full_projvecs,
            'num_proj_vecs': num_proj_vecs,
            'total_amplitude': np.sqrt(np.sum([a**2
                                               for a in amplitudes]))}
    return results_by_wyck

def get_mode_decomposition(struct_hs, struct_ls, nonzero_only=False):
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
    sgn_ls = get_sym_info(struct_ls)[0]

    basis, origin, displacements, struct_hs_supercell = match_structures(struct_ls, struct_hs)

    try:
        directions = iso.getDirections(sgn_hs, basis, origin, subgroup=sgn_ls)
    except iso.IsotropyBasisException:
        # isotropy is picky about some things
        basis = np.dot(-1 * np.identity(3), basis)
        directions = iso.getDirections(sgn_hs, basis, origin, subgroup=sgn_ls)
        print("trying with inverted basis {}".format(basis))
    except iso.IsotropySubgroupException:
        print("Isotropy isn't recognizing the subgroup relation in this basis")
        print("trying more general projections")
        print("Perhaps double check this")
        directions = iso.getDirections(sgn_hs, basis, origin)

    print("Basis: \n{}".format(basis))
    print("Origin detected: {}".format(origin))
    print("Origin using: {}".format([str(Fraction(i).limit_denominator(10)) for i in origin]))

    this_subgroup_distortions, directions_dict = get_all_distortions(sgn_hs, list(set(wyckoff_list)),
                                                                     directions, basis, origin)

    all_in_sc_basis = convert_distortions_basis(this_subgroup_distortions, origin,
                                                struct_hs.lattice,
                                                struct_hs_supercell.lattice)

    irrep_dist_full_pvecs = {}
    irrep_dist_defs = {}
    irrep_amplitudes = {}

    mode_decomposition_data = {}
    for irrep, wycks in all_in_sc_basis.items():
        proj_data_by_wyck = get_projection_data(displacements, wycks,
                                                struct_hs_supercell, wyckoff_list, struct_hs)
        for wyck in proj_data_by_wyck.keys():
            proj_data_by_wyck[wyck]['direction'] = directions_dict[irrep]
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
    struct_hs = pmg.Structure.from_file(sys.argv[1])
    struct_ls = pmg.Structure.from_file(sys.argv[2])

    irrep_decomposition_data = get_mode_decomposition(struct_hs, struct_ls, nonzero_only=False)

    print("Mode Definitions:")
    for irrep, wycks in irrep_decomposition_data.items():
        print(irrep)
        for wyck, data in wycks.items():
            print(wyck)
            print(data["direction"])
            print(data["dist_defs"])
    print()
    print()
    print("Mode Amplitudes")
    for irrep, wycks in irrep_decomposition_data.items():
        print(irrep)
        for wyck, data in wycks.items():
            print(wyck)
            print(data["direction"])
            print(np.round_(data["amplitudes"], decimals=5))
            print(np.round_(data["total_amplitude"], decimals=5))
        # print(np.round_(data["Amplitudes_as"], decimals=5))
