import os, sys, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Process, Manager
man = Manager()
run_bool = man.list([0 for _ in range(100)])



from BMI_L import analysis, average, BMI_rotacf
analysis('C01_C04')
