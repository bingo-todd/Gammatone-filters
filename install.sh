#!/usr/bin/env bash
python setup.py sdist
pip install dist/*.tar.gz
