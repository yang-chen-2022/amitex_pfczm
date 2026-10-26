tests by LG on marquises
------------------------

ENV MARQUISES: 
  python3.10.4 + pip3 install numpy (+scipy+matplotlib+h5py+pandas)
  intel oneapi2021

31/07/2025  
py4amitex 	- 0753d9 (25/07) : pip3 install . (OK)
amitex_fftp	- 82cc51 (25/07) : .install (OK)
in py4amitex_tests :
	python3 run_tests.py (OK 6 tests)

Seems that PYTHONPATH towards py4amitex is not mandatory : how does it works?
  GB : when installing a package with pip it is stored in ~/.local/lib/python3.X/site-packages and python knows it
       when installing with cmake -DCMAKE_INSTALL_PREFIX=DOSSIER_INSTALL, there PYTHONPATH must be given PYTHONPATH=INSTALL_DIR/lib.
