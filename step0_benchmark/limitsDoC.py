#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys

cov=float(sys.argv[1])

low=int(max(8, cov/3))

up=int(2*cov)

print(str(low)+' '+str(up))
