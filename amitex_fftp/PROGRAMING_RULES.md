PRAOGRAMMING RULES
=====================
(in progress)
examples here : https://www-git-cad.intra.cea.fr/DEC/pleiades/documentation/pgcl-pleiades (CEA users)
                https://github.com/2decomp-fft/2decomp-fft/blob/dev/Contribute.md


RULES for the development of modules in AMITEX
•	Use ‘private’ at the beginning of the module so that everything is private by default
•	Declare ‘public’ only the variables that need to be public.
•	‘public’ variables, if not 'parameter', must be pointers
•	In simu_mod.f90, for every new ‘public‘ pointer a corresponding member (used as target for the pointer) must be added to the object SIMU_AMITEX,
•	In simu_mod.f90 a mechanism allows to move pointer association from one target to another, it must be updated for every new ‘public’ pointer.
If omitting this mechanism of pointer association, the risk is that every simulation will be done with the same variable (which can be an issue if it is modified during the simulation).
