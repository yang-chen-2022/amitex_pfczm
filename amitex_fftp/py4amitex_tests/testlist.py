# -*- coding:utf-8 -*-
# List of all tests

# [(test name, number of processes, type: 'mechanics|diffusion')]
# Default values: type='mechanics', nbprocs=0
TESTS = [
    {"name": "bilayer_2mat_2zone_xmla", "nbprocs": 4},
    {"name": "diffusion_fluximp", "type": "diffusion"},
    {"name": "lin_elas_def_imp_65"},
    {"name": "polyx_elasiso"},
    {"name": "polyxGD_triax2p5_cauchy"},
    {"name": "torsion_beam"},
]
