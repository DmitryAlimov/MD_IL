import os, sys


file_stats = os.stat('/home/student/bmim/anneal/190.5/nvt.xtc')


print(file_stats.st_size/1024/1024/1024)


for T in [160, 180, 200, 220]:
    for tmc in range(5, 310, 5):
        if not os.path.exists(f'/home/student/bmim/anneal/{T}.{tmc}'):
            os.system(f'mkdir /home/student/bmim/anneal/{T}.{tmc}')

        print(f'COPYING... tmc = {tmc}ns, T = {T}K')
        files = os.listdir(f'/home/student/bmim/anneal/{T}.{tmc}')
        for file in ['anneal.cpt', 'eql.edr', 'anneal.edr', 'eql_prev.cpt', 'nvt.gro', 'anneal.log', 'nvt.cpt',
                     'mdout.mdp', 'eql.gro', 'nvt.edr', 'anneal.gro', 'eql.trr', 'anneal.trr', 'nvt.xtc', 'eql.tpr',
                     'nvt.trr', 'nvt.log', 'eql.log', 'anneal.tpr', 'anneal.xtc', 'eql.xtc', 'eql.cpt', 'nvt.tpr',
                     'nvt_prev.cpt', 'anneal_prev.cpt']:
            if not os.path.exists(f'/home/student/bmim/anneal/{T}.{tmc}/{file}'):
                print(f'    {file}')
                os.system(
                    f'sshpass -p tcm6zqcy scp -P 22 -r enigma@192.168.1.108:/data/bmim/anneal/{T}.{tmc}/{file} /home/student/bmim/anneal/{T}.{tmc}/')
            elif os.stat(f'/home/student/bmim/anneal/{T}.{tmc}/{file}').st_size < os.stat(
                    f'/home/student/bmim/anneal/160.240/{file}').st_size * 1.05:
                print(f'    {file}')
                os.system(
                    f'sshpass -p tcm6zqcy scp -P 22 -r enigma@192.168.1.108:/data/bmim/anneal/{T}.{tmc}/{file} /home/student/bmim/anneal/{T}.{tmc}/')



