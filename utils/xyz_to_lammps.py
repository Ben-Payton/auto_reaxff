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
    positions = u.atoms.positions
    velocities = u.atoms.velocities
    charges = u.atoms.charges
    atom_labels = u.atoms.types
    molecule_ids = [atom.residue.resid for atom in u.atoms]

    data_dict = {
        "positions":positions,
        "velocities":velocities,
        "charges":charges,
        "atom_type":atom_labels,
        "molecule_ids":molecule_ids,
    }
    return data_dict

def make_velocity_section(data_dict):
    string_to_format = "{atom_id} {vx} {vy} {vz}\n"
    final_string = "Velocities\n\n"
    for index, velocity in enumerate(data_dict["velocities"]):
        final_string = final_string + string_to_format.format(atom_id=index+1,vx=velocity[0],vy=velocity[1],vz=velocity[2])
    return final_string

def make_atoms_section(data_dict):
    string_to_format = "{atom_id} {molecule_id} {atom_type} {charge} {x_pos} {y_pos} {z_pos}\n"
    final_string = "Atoms\n\n"
    for index, postion in enumerate(data_dict["positions"]):
        final_string = final_string + string_to_format.format(atom_id=index+1,
                                                              molecule_id=data_dict["molecule_ids"][index],
                                                              atom_type=data_dict["atom_type"][index],
                                                              charge=data_dict["charges"][index],
                                                              x_pos=postion[0],
                                                              y_pos=postion[1],
                                                              z_pos=postion[2],)
    return final_string


def main():


    # HARD CODED PATHS TO FILES DO NOT FORGET TO CHANGE LATER
    TPR = "/home/ben/01_Projects/work_projects/auto_reaxff/test_files/gromacs_scripts_testing/salt_water-npt-prod.tpr"
    TRR = "/home/ben/01_Projects/work_projects/auto_reaxff/test_files/gromacs_scripts_testing/salt_water-npt-prod.part0001.trr"

    data_dict = check_trr(TPR,TRR)



if __name__ == "__main__":
    main()