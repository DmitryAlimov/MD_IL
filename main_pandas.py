import os, sys
import pandas as pd
import numpy as np

tmc_list = [5*(i//5) for i in range(1*5, 400)]
T_list = [[160, 180, 190, 200, 210][i%5] for i in range(1*5, 400)]

def check_data(data_name):
    return [int(os.path.exists(f'/home/student/bmim/analysis/{data_name}/{T_list[i]}.{tmc_list[i]}.npy')) for i in range(len(T_list))]

def data_path(data_name, tmc, T):
    return f'/home/student/bmim/analysis/{data_name}/{T}.{tmc}.npy'

checkDF = pd.DataFrame(
    {
        " |": "|",
        "tmc": tmc_list,
        "T": T_list,
        "|": "|",
        "rdf/TMP.N-BMI.C0N": check_data("rdf/TMP.N-BF.B"),
        "rdf/TMP.N-BMI.C06": check_data("rdf/TMP.N-BMI.C06"),
        "rdf/TMP.N-BMI.C0N/": check_data("rdf/TMP.N-BMI.C0N/"),
        "rdf/TMP.N-BF.B/": check_data("rdf/TMP.N-BF.B/"),
        "rdf/BMI.C01-BF.F/": check_data("rdf/BMI.C01-BF.F/"),
        "rdf/BMI.C01-BF.B/": check_data("rdf/BMI.C01-BF.B/"),
        "rdf/BMI.N02-BF.F/": check_data("rdf/BMI.N02-BF.F/"),
        "rdf/BMI.N02-BF.B/": check_data("rdf/BMI.N02-BF.B/"),
        "rdf/BMI.C0N-BF.F/": check_data("rdf/BMI.C0N-BF.F/"),
        "rdf/BMI.C0N-BF.B/": check_data("rdf/BMI.C0N-BF.B/"),
        "rdf/BMI.C01-BMI.C01/": check_data("rdf/BMI.C01-BMI.C01/"),
        "rdf/BMI.C0N-BMI.C0N/": check_data("rdf/BMI.C0N-BMI.C0N/"),
        "rdf/BF.B-BF.B/": check_data("rdf/BF.B-BF.B/"),
        "rdf/TMP.N-BMI/": check_data("rdf/TMP.N-BMI/"),
        "rdf/TMP.N-BF/": check_data("rdf/TMP.N-BF/"),
        "rdf/TMP.N-BMI.C03_C04/": check_data("rdf/TMP.N-BMI.C03_C04/"),
        "rdf/TMP.N-BMI.N00/": check_data("rdf/TMP.N-BMI.N00/"),
        "rdf/TMP.N-BMI.N02/": check_data("rdf/TMP.N-BMI.N02/"),
        "rdf/TMP.N-BMI.COG/": check_data("rdf/TMP.N-BMI.COG/"),
        "rdf/TMP.N-BMI.C0J/": check_data("rdf/TMP.N-BMI.C0J/"),
        "rdf/TMP.N-BF.F/": check_data("rdf/TMP.N-BF.F/"),
    }
)




print(checkDF)



