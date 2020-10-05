#!/usr/bin/env bash
rm -r  ~/Work_Space/anaconda3/envs/py_3.7/lib/python3.7/site-packages/GTF/
python setup.py sdist
pip install dist/*.tar.gz
