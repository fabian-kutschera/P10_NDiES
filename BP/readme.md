# How to install packages and dependencies?

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
11. pip install geopy

Note: due to some strange reason execute import cell twice in order to import everything successfully.
