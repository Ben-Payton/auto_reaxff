##########################################################################################
"""
This file is meant to hold tools that will take an xyz_molecule file and turn it in to a
base input molecule file for lammps
"""
##########################################################################################
import MDAnalysis as mda
import os
import numpy as np



def check_trr(TPR_PATH,TRR_PATH):
    u = mda.Universe(TPR_PATH,TRR_PATH)
    positions = u.atoms.positions * np.array([10,10,10]) # nanometers to angstroms
    velocities = u.atoms.velocities * np.array([10,10,10]) # nanometers to angstroms
    charges = u.atoms.charges
    atom_labels = u.atoms.types
    molecule_ids = [atom.residue.resid for atom in u.atoms]
    masses = [atom.mass for atom in u.atoms]

    data_dict = {
        "positions":positions,
        "velocities":velocities,
        "charges":charges,
        "atom_type":atom_labels,
        "molecule_ids":molecule_ids,
        "masses":masses,
    }
    return data_dict

def make_velocity_section(data_dict):
    string_to_format = "{atom_id} {vx} {vy} {vz}\n"
    final_string = "Velocities\n\n"
    for index, velocity in enumerate(data_dict["velocities"]):
        final_string = final_string + string_to_format.format(atom_id=index+1,vx=velocity[0],vy=velocity[1],vz=velocity[2])
    return final_string

def make_atoms_section(data_dict,swapping_dict):
    string_to_format = "{atom_id} {molecule_id} {atom_type} {charge} {x_pos} {y_pos} {z_pos}\n"
    final_string = "Atoms\n\n"
    for index, postion in enumerate(data_dict["positions"]):
        final_string = final_string + string_to_format.format(atom_id=index+1,
                                                              molecule_id=data_dict["molecule_ids"][index],
                                                              atom_type=swapping_dict[data_dict["atom_type"][index]],
                                                              charge=data_dict["charges"][index],
                                                              x_pos=postion[0],
                                                              y_pos=postion[1],
                                                              z_pos=postion[2],)
    return final_string


def make_mass_dict(data_dict):
    acc = []
    values = []
    
    # First we get all the unique atom types and masses
    for atom_type, atom_mass in zip(data_dict["atom_type"],data_dict["masses"]):
        if not atom_type in acc:
            acc.append(atom_type)
            values.append((atom_type,atom_mass))
    # Sort them by mass
    values = sorted(values,key=lambda x: x[1])

    # add atom labels for lammps based off similar masses full_paired_has format (numeric_atom_type, mass, gromacs_atom_label)
    label_num = 1
    full_paired = []
    for index,value in enumerate(values[:-1]):
        full_paired.append((label_num,value[1],value[0]))
        if value[1] != values[index+1][1]:
            label_num = label_num + 1
    full_paired.append((label_num,values[-1][1],values[-1][0]))

    # Build Dictionaries for writting masses and for subbing atom types

    ## Build dict for writing masses values have form (mass, gromacs label)
    masses_dict = {}
    for group in full_paired:
        if not group[0] in masses_dict.keys():
            masses_dict[group[0]] = [group[1],group[2]]
        else:
            masses_dict[group[0]][1] = masses_dict[group[0]][1] + f", {group[2]}"
    
    # Build dict for swapping atom types
    swap_dict = {}
    for group in full_paired:
        swap_dict[group[2]] = group[0]
            
    return masses_dict, swap_dict




def make_preamble_section(data_dict,masses_dict):
    num_atoms = len(data_dict["charges"])
    num_atom_types = len(np.unique(data_dict["atom_type"]))
    low_x = min(data_dict["positions"],key=lambda x: x[0])[0] - 3
    hi_x = max(data_dict["positions"],key=lambda x: x[0])[0] + 3
    low_y = min(data_dict["positions"],key=lambda x: x[0])[0] - 3
    hi_y = max(data_dict["positions"],key=lambda x: x[0])[0] + 3
    low_z = min(data_dict["positions"],key=lambda x: x[0])[0] - 3
    hi_z = max(data_dict["positions"],key=lambda x: x[0])[0] + 3

    begining_string = f"""LAMMPS input generated using MDAnalysis and code written by Ben Payton. PLEASE CHECK INPUTS ARE CORRECT!!

{num_atoms} atoms

{num_atom_types} atom types

{low_x:.1f} {hi_x:.1f} xlo xhi
{low_y:.1f} {hi_y:.1f} ylo yhi
{low_z:.1f} {hi_z:.1f} zlo zhi

Masses

"""
    masses_format_string = "{atom_type_number} {atom_mass:.6f}       # gromacs atom labels represented: {fomer_atoms}\n"
    for key in masses_dict.keys():
        value = masses_dict[key]
        begining_string = begining_string + masses_format_string.format(atom_type_number=key, atom_mass=value[0],fomer_atoms=value[1])
    return begining_string


def trr_tpp_to_data(TRR_file_path,TPR_file_path,final_file_path):
    data_dict = check_trr(TPR_file_path,TRR_file_path)
    masses_dict, swapping_dict = make_mass_dict(data_dict)
    final_string = ""
    final_string = final_string + make_preamble_section(data_dict,masses_dict) + "\n"
    final_string = final_string + make_atoms_section(data_dict,swapping_dict) + "\n"
    final_string = final_string + make_velocity_section(data_dict)
    with open(final_file_path,"w") as file:
        file.write(final_string)

    
def main():


    # HARD CODED PATHS TO FILES DO NOT FORGET TO CHANGE LATER
    TPR = "/home/ben/01_Projects/work_projects/auto_reaxff/test_files/gromacs_scripts_testing/salt_water-npt-prod.tpr"
    TRR = "/home/ben/01_Projects/work_projects/auto_reaxff/test_files/gromacs_scripts_testing/salt_water-npt-prod.part0001.trr"

    trr_tpp_to_data(TRR,TPR,"./test_files/test_dat.data")

if __name__ == "__main__":
    main()