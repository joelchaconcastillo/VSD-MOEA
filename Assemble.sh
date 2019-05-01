  instance=$1
  nobj=$2
  nvar=$3
  #PATH1=/home/joel.chacon/Chacon/Tesis/MOEA-D_Diversity/Code_CEC09/moead_based_diversity_SBX
  #PATH2=/home/joel.chacon/Chacon/Tesis/jMetalCpp-1.7/SBX-Gen250th-NObjs-eta-2
  #PATH3=/home/joel.chacon/R2MOEA/Results_Main
  #PATH4=/home/joel.chacon/VSD-MOEA/VSD-MOEA
  #PATH4=/home/joel.chacon/Chacon/Tesis/MOEA-D_Diversity/Code_CEC09/moead_based_diversity_SBX
  PATH4=/home/joel.chacon/VSD-MOEA/VSD-MOEA/tmp
  #./distances \"\"$PATH1/POS/*${instance}_*_${nobj}_nvar_*eta-2.dat_2\"\" 11 100 ${nvar}  > moea-d
  #./distances \"\"$PATH2/NSGAII-Obj-${nobj}/VAR_NSGAII_SBX_POLYNOMIAL_EVALUATIONS_25000000_${instance}_*.txt_2\"\" 11 100 ${nvar} | cut -f2 -d' '  > nsga-ii
  #./distances \"\"$PATH3/${instance}_${nobj}_*VAR_2\"\" 11 100 ${nvar} | cut -f2 -d' ' > r2-moea
  
#  ./distances \"\"$PATH4/POS/Scalability_ShortTerm/POS_*${instance}_*${nobj}_nvar_${nvar}_*\"\" 11 100 ${nvar} | cut -f2 -d' ' 
  ./distances \"\"$PATH4/POS/mindist/POS_*${instance}_*var_${nvar}_nobj_${nobj}.dat_*0.4*\"\" 11 100 ${nvar} | cut -f2 -d' '


  #paste moea-d nsga-ii r2-moea vsd-moea | column -s $' ' -t > out
 

