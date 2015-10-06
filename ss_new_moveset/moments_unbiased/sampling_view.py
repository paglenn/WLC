import numpy as np
import matplotlib.pyplot as plt
import os

dirs = os.listdir('.')
dirs = [d for d in dirs if '_0_0' == d[-4:]]
print dirs

for d in dirs:


