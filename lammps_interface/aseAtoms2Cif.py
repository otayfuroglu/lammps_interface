#
from ase.io.cif import (CIFLoop, format_cell,
                        chemical_formula_header,
                        format_generic_spacegroup_info,
                        atoms_to_loop_data)


def get_cell_info(atoms):
    cell_info = atoms.get_cell_lengths_and_angles()
    return (f'_cell_length_a        {cell_info[0]}\n'
            f'_cell_length_b        {cell_info[1]}\n'
            f'_cell_length_c        {cell_info[2]}\n'
            f'_cell_angle_alpha     {cell_info[3]}\n'
            f'_cell_angle_beta      {cell_info[4]}\n'
            f'_cell_angle_gamma     {cell_info[5]}\n')


def get_space_grpP1():
    return (
        '_symmetry_space_group_name_H-M    "P 1"\n'
        '_symmetry_int_tables_number       1\n')


def aseAtoms2Cif(atoms, cif_format=None, wrap=True,
                 labels=None, loop_keys=None):

    cif_lines = []
    cif_lines.append('data_image0\n')
    #  cif_lines.append(get_cell_info(atoms))
    #  cif_lines.append(get_space_grpP1())
    #  cif_lines.append(chemical_formula_header(atoms))

    rank = atoms.cell.rank
    if rank == 3:
        cif_lines.append(format_cell(atoms.cell))
        cif_lines.append(format_generic_spacegroup_info())
    elif rank != 0:
        raise ValueError('CIF format can only represent systems with '
                         f'0 or 3 lattice vectors.  Got {rank}.')

    loop_keys = {}
    loopdata, coord_headers = atoms_to_loop_data(atoms, wrap,
                                                 labels, loop_keys)

    headers = [
        '_atom_site_type_symbol',
        '_atom_site_label',
        '_atom_site_symmetry_multiplicity',
        *coord_headers,
        '_atom_site_occupancy',
    ]

    headers += ['_' + key for key in loop_keys]

    loop = CIFLoop()
    for header in headers:
        array, fmt = loopdata[header]
        loop.add(header, array, fmt)

    cif_lines.append(loop.tostring())

    # split to list item all lines
    list_cif_lines = []
    for cif_line in cif_lines:
        list_cif_lines += cif_line.split("\n")

    return list_cif_lines
