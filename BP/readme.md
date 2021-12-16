# How to install packages and dependencies?

Execute the following commands line by line:

1. conda create -n python373 python=3.7.3
2. conda activate python373
3. pip install h5py (might not be necessary/ actually used)
4. conda install obspy
5. conda install jupyter
6. git clone https://github.com/sergiocallegari/PyDSM
7. cd PyDSM
8. git submodule init
9. git submodule update
10. python setup.py install (this still might return an error! In that case you probably need to install cython via 
`conda install cython`)
11. cd ..
12. pip install geopy
13. jupyter-notebook (starts Jupyter Notebook)

Note 1: Downlaod BP Exercise and unzip in the same folder where you downloaded your PyDSM package.

Note 2: Due to some strange reason execute import cell twice in order to import everything successfully.
