HOWTO DEFINE A NEW VERSION (i.e. commit) 2DECOMP SUBMODULE?

in 2decomp_fft :
git checkout dev
git pull
git checkout idcommit

checks : compilation + validation tests

AFTER clean:
git add lib_extern/2decomp
git push

HOWTO 'GET' THE NEW VERSION ?

in an existing local repo : 
git checkout lib_extern/2decomp_fft
git pull
git submodule update

for a new repo :
git clone (amitex_fftp)
cd amitex_fftp
git submodule init
git submodule update

=====================================================================
STATUS : 
90b87 (25/07/25) : NOK(test 6/10/2025)
a80e7 (8/7/2025) : OK (submodule update 6/10/2025)


