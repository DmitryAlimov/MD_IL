from rdf import average
"""for analysis_type in ['TMP.N-BMI.C06', 'TMP.N-BF.B', 'TMP.N-BF.F', 'TMP.N-BMI.C0N', 'TMP.N-BMI.C0J', 'TMP.N-BMI.COG',
              'TMP.N-BMI.N02', 'TMP.N-BMI.N00', 'TMP.N-BMI.C03_C04', 'TMP.N-BMI', 'TMP.N-BF', 'BF.B-BF.B',
             'BMI.C0N-BMI.C0N', 'BMI.C01-BMI.C01', 'BMI.C0N-BF.B', 'BMI.C0N-BF.F',
             'BMI.N02-BF.B', 'BMI.N02-BF.F', 'BMI.C01-BF.B', 'BMI.C01-BF.F']:
    average(analysis_type)"""

from TMP_L import analysis, average
################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ПРОВЕРИТЬ КУСОК КОДА, ОБОЗНАЧЕННЫЙ ЭТИМ ЗНАКОМ: ЗАМЕНИТЬ ИНДЕКС 5000 НА 100000
#analysis()
#average('DL_SL')
#average('NO_FFT')
#average('NO_DL_van_hoff')
#average('NO_SL_van_hoff')
#average('NO_rotacf')


from BMI_L import analysis, average, BMI_rotacf
#BMI_rotacf('C01_C04')
################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ПРОВЕРИТЬ КУСОК КОДА, ОБОЗНАЧЕННЫЙ ЭТИМ ЗНАКОМ: ЗАМЕНИТЬ ИНДЕКС 5000 НА 100000
#analysis('C01_C04')
#average('C01_C04', 'DL_SL')
#average('C01_C04', 'DL_van_hoff')
#average('C01_C04', 'SL_van_hoff')
#average('C01_C04', 'rotacf')

#BMI_rotacf('C0N_N02')
################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ПРОВЕРИТЬ КУСОК КОДА, ОБОЗНАЧЕННЫЙ ЭТИМ ЗНАКОМ: ЗАМЕНИТЬ ИНДЕКС 5000 НА 100000
#analysis('C0N_N02')
#average('C0N_N02', 'DL_SL')
#average('C0N_N02', 'DL_van_hoff')
#average('C0N_N02', 'SL_van_hoff')
#average('C0N_N02', 'rotacf')

from edr import average
### ПРОВЕРИТЬ ИНДЕКСАЦИЮ, ВОЗМОЖНО НЕ НУЖЕН ИНДЕКС [0]
#average('shear_viscosity')

from clusters import analysis, average
################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ПРОВЕРИТЬ КУСОК КОДА, ОБОЗНАЧЕННЫЙ ЭТИМ ЗНАКОМ: ЗАМЕНИТЬ ИНДЕКС 500 НА 10000
#analysis('Nodes_to_Comms')
#analysis('my_RandtIndex')
#analysis('RandtIndex')
#average('mean')
#average('my_RandtIndex')
#average('RandtIndex')


from BMI_TMP_correlation import analysis, average
analysis('clothest_BMI_chains')
#analysis('clothest_BMI_rings')
#average()







