#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob
import subprocess
import random

def write_script_header(cluster, script, event_id, walltime, working_folder):
    if cluster == "nersc":
        script.write(
"""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J photon_%s
#SBATCH -t %s
#SBATCH -L SCRATCH
#SBATCH -C haswell
""" % (event_id, walltime))
    elif cluster == "wsugrid":
        script.write(
"""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J photon_%s
#SBATCH -t %s
""" % (event_id, walltime))
    elif cluster == "guillimin":
        script.write(
"""#!/usr/bin/env bash
#PBS -N photon_%s
#PBS -l nodes=1:ppn=1
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -d %s
""" % (event_id, walltime, working_folder))
    elif cluster == "McGill":
        script.write(
"""#!/usr/bin/env bash
#PBS -N photon_%s
#PBS -l nodes=1:ppn=1:irulan
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -d %s
""" % (event_id, walltime, working_folder))
    else:
        print("Error: unrecoginzed cluster name :", cluster)
        print("Available options: nersc, wsugrid, guillimin, McGill")
        exit(1)


def generate_script(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '10:00:00'
    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir UrQMD_results
for iev in `ls OSCAR_events`
do
    cd osc2u
    ./osc2u.e < ../OSCAR_events/$iev
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh
    mv particle_list.dat ../UrQMD_results/particle_list_`echo $iev | cut -f 2 -d _`
    cd ..
done
""")
    script.close()


def generate_event_folder(cluster_name, working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    mkdir(path.join(event_folder, 'hydro_events'))
    generate_script_iSS(cluster_name, event_folder)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('hydro_photonEmission.e'),
        path.join(path.abspath(event_folder), "hydro_photonEmission.e")),
        shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('ph_rates'),
        path.join(path.abspath(event_folder), "ph_rates")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('parameters.dat'),
        path.join(path.abspath(event_folder), "parameters.dat")), shell=True)


def copy_hydro_events(number_of_cores, input_folder, working_folder):
    events_list = glob('%s/hydro_results*/evolution_all_xyeta.dat'
                       % input_folder)
    for iev in range(len(events_list)):
        event_id = events_list[iev].split('/')[-2].split('_')[-1][0]
        folder_id = iev % number_of_cores
        working_path = path.join(working_folder, 'event_%d' % folder_id,
                                 'hydro_events')
        folder_path = path.join(working_path, events_list[iev].split('/')[-1])
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[iev]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)
        shutil.copy(path.join(input_folder,
                              'music_input_event_%s' % event_id),
                    working_path)


if __name__ == "__main__":
    try:
        from_folder = str(sys.argv[1])
        folder_name = str(sys.argv[2])
        cluster_name = str(sys.argv[3])
        ncore = int(sys.argv[4])
        mode = int(sys.argv[5])
    except IndexError:
        print("Usage:")
        print("  %s input_folder working_folder cluster_name num_of_cores mode"
              % str(sys.argv[0]))
        print("")
        exit(0)

    for icore in range(ncore):
        generate_event_folder_iSS(cluster_name, folder_name, icore)
    copy_hydro_events(ncore, from_folder, folder_name)
