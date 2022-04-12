import mdtraj as md
import os, sys

gmx = '/usr/local/gromacs/bin/gmx'

# snapshots
'''for i in range(5, 140, 5):
    for T in [160, 180, 190]:
        os.system(f'mkdir /home/student/bmim/anneal/{T}.{i}')
        os.system(f'echo 0|{gmx} trjconv -f /home/md_beast/Svetlana/IL_EPR/bmim/6nm/md_nvt3/md_nvt3.part0002.trr -s /home/md_beast/Svetlana/IL_EPR/bmim/6nm/md_nvt3/md_nvt3.tpr -o /home/student/bmim/anneal/{T}.{i}/{i}.gro -b {62500 + i*1000} -e {62500 + i*1000}')
        os.system(f'echo 0|{gmx} trjconv -f /home/md_beast/Svetlana/IL_EPR/bmim/6nm/md_nvt3/md_nvt3.part0002.trr -s /home/md_beast/Svetlana/IL_EPR/bmim/6nm/md_nvt3/md_nvt3.tpr -o /home/student/bmim/anneal/{T}.{i}/{i}.trr -b {62500 + i*1000} -e {62500 + i*1000}')

for i in range(5, 200, 5):
    os.system(f'echo 0 | {gmx} trjconv -f /home/student/bmim/300K_nvt/4md_nvt_0/4md_nvt.xtc -s /home/student/bmim/300K_nvt/4md_nvt_0/4md_nvt.tpr -o /home/student/bmim/anneal/{135 + i}.gro -b {i * 1000} -e {i * 1000}')
    if os.path.exists(f'/home/student/bmim/anneal/{135 + i}.gro'):
        for T in [160, 180, 190]:
            os.system(f'mkdir /home/student/bmim/anneal/{T}.{135 + i}')
            os.system(f'cp /home/student/bmim/anneal/{135 + i}.gro /home/student/bmim/anneal/{T}.{135 + i}/{135 + i}.gro')
            os.system(f'echo 0|{gmx} trjconv -f /home/student/bmim/300K_nvt/4md_nvt_0/4md_nvt.trr -s /home/student/bmim/300K_nvt/4md_nvt_0/4md_nvt.tpr -o /home/student/bmim/anneal/{T}.{135 + i}/{135 + i}.trr -b {i*1000} -e {i*1000}')
        os.system(f'rm /home/student/bmim/anneal/{135 + i}.gro')

for i in range(0, 400, 5):
    os.system(f'echo 0 | {gmx} trjconv -f /home/student/bmim/300K_nvt/4md_nvt_1/4md_nvt.xtc -s /home/student/bmim/300K_nvt/4md_nvt_1/4md_nvt.tpr -o /home/student/bmim/anneal/{190 + i}.gro -b {i * 1000} -e {i * 1000}')
    if os.path.exists(f'/home/student/bmim/anneal/{190 + i}.gro'):
        for T in [160, 180, 190]:
            os.system(f'mkdir /home/student/bmim/anneal/{T}.{190 + i}')
            os.system(f'cp /home/student/bmim/anneal/{190 + i}.gro /home/student/bmim/anneal/{T}.{190 + i}/{190 + i}.gro')
            os.system(f'echo 0|{gmx} trjconv -f /home/student/bmim/300K_nvt/4md_nvt_1/4md_nvt.trr -s /home/student/bmim/300K_nvt/4md_nvt_1/4md_nvt.tpr -o /home/student/bmim/anneal/{T}.{190 + i}/{190 + i}.trr -b {i * 1000} -e {i * 1000}')
        os.system(f'rm /home/student/bmim/anneal/{190 + i}.gro')'''

"""for tmc in range(0, 400, 5):
    if os.path.exists(f'/home/student/bmim/anneal/190.{tmc}/{tmc}.gro'):
        os.system(f'mkdir /home/student/bmim/anneal/220.{tmc}')
        os.system(f'cp /home/student/bmim/anneal/190.{tmc}/{tmc}.gro /home/student/bmim/anneal/210.{tmc}')
        os.system(f'cp /home/student/bmim/anneal/190.{tmc}/{tmc}.trr /home/student/bmim/anneal/210.{tmc}')"""

# md
"""for i in range(120, 400, 5):
    for T in [160, 180]:
        if os.path.exists(f'/home/student/bmim/anneal/{T}.{i}/{i}.gro'):
            os.chdir(f'/home/student/bmim/anneal/{T}.{i}')
            if not os.path.exists(f'anneal.gro'):
                os.system(
                    f'{gmx} grompp -f /home/student/bmim/gmx_inp/5anneal{T}.mdp -c {i}.gro -t {i}.trr -p /home/student/bmim/topol.top -o anneal.tpr')
                os.system(f'{gmx} mdrun -ntmpi 1 -nb gpu -bonded gpu -deffnm anneal')

            if (not os.path.exists(f'eql.gro')) and os.path.exists(f'anneal.gro'):
                os.system(
                    f'{gmx} grompp -f /home/student/bmim/gmx_inp/6md_eql_{T}.mdp -c anneal.gro -t anneal.cpt -p /home/student/bmim/topol.top -o eql.tpr')
                os.system(f'{gmx} mdrun -deffnm eql')

            if (not os.path.exists(f'nvt.gro')) and os.path.exists(f'eql.gro'):
                os.system(
                    f'{gmx} grompp -f /home/student/bmim/gmx_inp/7md_nvt_{T}.mdp -c eql.gro -t eql.cpt -p /home/student/bmim/topol.top -o nvt.tpr')
                os.system(f'{gmx} mdrun -deffnm nvt')"""



