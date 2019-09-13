This is the code of VSD-MOEA that was used for the paper "VSD-MOEA: A Dominance-Based Multi-Objective Evolutionary Algorithm with Explicit Variable Space Diversity Management" submitted to the IEEE Transactions on Cybernetics.
########Compile
This version only works with linux enviroments, to compile please locate in the dir "Code" and run the following in the terminal:

./make


###############Run
The command to execute VSD-MOEA is the following:
Ejecutable --n POPULATION_SIZE --nfes NUMBER_FUNCTION_EVALUATOINS --nvar NUMBER_VARIABLES --Instance NAME_PROBLEM --Path CURRENT_DIR --Dist_factor INITIAL_DISTANCE_FACTOR --nobj NUMBER_OBJECTIVES --Seed NUMBER_SEED --param_l DISTANCE_PARAMETERS(ONLY WFG PROBLEMS) --param_k 4 POSITION_PARAMETERS(ONLY WFG PROBLEMS)

All the executions carried out can be performed by just specifying the appropiate parameters. For instance, in order to execute VSD-MOEA wich DTLZ1 with the aim of generating the results of the first experiment (see Table 2 of the paper) the following command line must be used: 

./Ejecutable --n 100 --nfes 25000000 --nvar 6 --Instance DTLZ1 --Path . --Dist_factor 0.4 --nobj 2 --Seed 1 

