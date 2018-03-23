#!/usr/bin/python

import os
import shutil
import sys
import subprocess


predicted_value = subprocess.Popen(['./stability_predictor.py', '5PTI', 'A', '35', 'Y', 'D'])
print(predicted_value)
