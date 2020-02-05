cd ../
python setup.py sdist
pip install dist/GTF-1.0.tar.gz
cd examples
python efficiency.py
