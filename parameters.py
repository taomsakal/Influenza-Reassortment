"""
This file contains the rules and tables for what can infect what, along with the contact rates.
"""

import numpy as np
import itertools

viruses = set([tuple(x) for x in itertools.product(range(16),range(9))])

infection_table=np.array(
[[[1,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,1,1,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0]],
[[1,1,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0],
[0,1,1,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0],
[0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0]],
[[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1]],
[[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0]]])