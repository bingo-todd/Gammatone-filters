#!/bin/sh

dir=$(dirname $0)
python ${dir}/setup.py ${dir}/sdist
pip install ${dir}/dist/*.tar.gz
