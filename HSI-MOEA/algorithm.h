#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <queue>
#include <iomanip>
#include <cfloat>
#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"

class MOEA
{

public:
	MOEA();
	virtual ~MOEA();

	void init_population();                  // initialize the population

	void load_parameter();

	void evol_population();                                      
	void exec_emo(int run);

	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_pos(char savefilename[1024]);

        void penalize_nearest(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void select_farthest_penalized(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void select_best_candidate(vector<CIndividual *> &survivors, vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void compute_distances(vector<CIndividual *> &candidates, vector<CIndividual *> &survivors);
        void binary_tournament_selection(vector<CIndividual > &population, vector<CIndividual> &child_pop);
        void recombination(vector<CIndividual> &child_pop);
        void reproduction(vector<CIndividual> &population, vector<CIndividual> &child_pop);
        void update_diversity_factor();
        void computing_dominate_information(vector <CIndividual*> &pool);
        void select_first_survivors(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates);
        void update_domianted_information(vector<CIndividual*> &survivors, vector<CIndividual*> &current);
        void update_population(vector<CIndividual*> &survivors, vector<CIndividual> &population);

	void fast_non_dominated_sorting(vector <CIndividual*> &survivors);

	double distance( vector<double> &a, vector<double> &b);

	double distance_improvement( vector<double> &a, vector<double> &b);
	vector <CIndividual> population;
	vector<CIndividual> child_pop;	// memory solutions
	void operator=(const MOEA &moea);

public:
//
//	// algorithm parameters
	int nfes;
//	int     nfes, max_nfes;          //  the number of function evluations

};

MOEA::MOEA()
{

}

MOEA::~MOEA()
{

}
double MOEA::distance( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i])/(vuppBound[i] - vlowBound[i]);
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double MOEA::distance_improvement( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
	double maxd = -INFINITY;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = max(0.0,a[i]-b[i]);
	   dist += factor*factor;
	   maxd = max(maxd, max(b[i]-a[i],0.0));
	}
        if(dist == 0.0) return -maxd; //in case that this indicator is zero, this mean that it is a dominated individual...
   return dist;
   return sqrt(dist);
}
void MOEA::init_population()
{

    for(int i=0; i<pops; i++)
	{
		CIndividual indiv1, indiv2;
		// Randomize and evaluate solution
		indiv1.rnd_init();
		indiv1.obj_eval();
		// Save in the population
		population.push_back(indiv1);

		indiv2.rnd_init();
		indiv2.obj_eval();

		child_pop.push_back(indiv2);
		nfes++;
	}
}
void MOEA::operator=(const MOEA &alg)
{
	//population = alg.population;
}
void MOEA::evol_population()
{
	vector<CIndividual *> penalized, survivors;
//	cout << "generation "<<endl;
	//join the offspring and parent populations
	vector<CIndividual *> candidates;
	for(int i = 0; i < pops; i++)
	{
	  candidates.push_back( &(population[i]));
	  candidates.push_back( &(child_pop[i]));
	}
	computing_dominate_information(candidates); //computing the dominate count of each candidate individual...
	//select the "best" individuals that owns to candidate set and are moved in survivors set...
	select_first_survivors(survivors, candidates);
	//update the diversity-factor-parameter...	
	update_diversity_factor();
	//Pre-computing the neares distances both objective and decision spaces..
	compute_distances(candidates, survivors);
       	while( survivors.size() < pops )
	{
	  penalize_nearest(candidates, penalized);//penalize the nearest individuals.. 
	  if(candidates.empty())	  
	     select_farthest_penalized(survivors, penalized);//in case that all the individuals are penalized pick up the farstest and add it to survirvors
	  else
	    {
	     update_domianted_information(survivors, candidates); //update the rank of each candidate whitout penalized
	     select_best_candidate(survivors, candidates, penalized); // the best candidate is selected considering the improvemente distance, and the rank..
	    }
	}
	fast_non_dominated_sorting(survivors);//rank the survivors individuals..
	//this procedure is necesary since the penalized individuals
	update_population(survivors, population); //update the parent population 
	
	reproduction(population, child_pop); //generate a new population considering the survivors individuals...
}
void MOEA::fast_non_dominated_sorting(vector <CIndividual*> &survivors)
{
   vector< vector < int > > dominate_list(survivors.size()); //in the worst case the number of fronts is the same as the survivors size
   vector< int > dominated_count (survivors.size(), 0), currentfront;
   for(int i = 0; i < survivors.size(); i++)
   {
	   for(int j = 0; j < survivors.size(); j++)
	  {
		if(i==j) continue;
	       if( *(survivors[i]) < *(survivors[j]))
	   	    dominate_list[i].push_back(j);
		else if (*(survivors[j]) < *(survivors[i]))
		   dominated_count[i]++;
 	  }
	if(dominated_count[i] == 0 ) currentfront.push_back(i);// get the first front
   }
   int rank = 0;
   while(!dominate_list[rank].empty())
   {
	vector<int> nextFront;
	for(int i = 0; i < currentfront.size(); i++)
	{
	   survivors[currentfront[i]]->rank = rank;
	   for(int j = 0; j < dominate_list[currentfront[i]].size(); j++)
	   {
		dominated_count[dominate_list[currentfront[i]][j]]--;
		if( dominated_count[dominate_list[currentfront[i]][j]] == 0) nextFront.push_back(dominate_list[currentfront[i]][j]);
		
	   }
	}	
	rank++;
	currentfront = nextFront;
   }
}
void MOEA::update_population(vector<CIndividual*> &survivors, vector<CIndividual> &population)
{
   for(int i = 0; i < population.size(); i++) population[i] = *(survivors[i]);
}
void MOEA::update_domianted_information(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates)
{
     bool firstfrontcurrent = false; 
     for(int i = 0; i < candidates.size(); i++) if(candidates[i]->times_dominated==0) firstfrontcurrent = true; //check if there exists at least one candidate in the lowest current front
     
     if( !firstfrontcurrent) //this indicates that there is not able a current in the lowest front, so the next front is to be considered
	{	
	  
	   for(int i = 0; i < survivors.size(); i++)
	   {
		if(survivors[i]->times_dominated == 0)
		{
		      for(int j = 0; j < survivors[i]->ptr_dominate.size(); j++)
		  	   {
		  		survivors[i]->ptr_dominate[j]->times_dominated--;
		   	   }
		  	   survivors[i]->times_dominated--;
		}
	   }
	}
}
void MOEA::select_first_survivors(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates)
{
	vector<int> BestIndex;
	///Select the best improvement distance candidates....
	vector<bool> grid(candidates.size(), false);
	for(int m = 0; m < nobj; m++)
	{
		int indxmaxim;
		double bestvector = INFINITY;
		for(int i = 0; i <  candidates.size(); i++)
		 {	
		 	  if(grid[i])continue;
		        double s = 0.0;	
		        double maxv = -INFINITY;
		        for(int k = 0; k < nobj; k++)
		        {
		      	   double fi = fabs(candidates[i]->y_obj[k]);
		      	   s += fi;
		      	   double ti = (k==m)?fi:1e5*fi;
			    if(ti > maxv)   maxv=ti;
		        }
		         maxv = maxv + 0.0001*s;
			//maxv = candidates[i]->y_obj[m];
		        if(bestvector > maxv)
		        { indxmaxim = i; bestvector = maxv;}
		  }
		grid[indxmaxim] = true;
		BestIndex.push_back(indxmaxim);
	}
		std::vector<int>::iterator End;
		End = std::unique (BestIndex.begin(), BestIndex.end());
		BestIndex.resize( std::distance(BestIndex.begin(),End));
		sort(BestIndex.begin(), BestIndex.end()); //sort the indexes and remove from candidates

		for(int i = BestIndex.size()-1; i >= 0 ; i--)
		{
			int index = BestIndex[i];
			survivors.push_back( candidates[index]);
			for(int m =0; m < nobj; m++)
			cout << candidates[index]->y_obj[m] << " ";
			cout << endl;
			iter_swap(candidates.begin()+index, candidates.end()-1);
			candidates.pop_back();
		}
//		getchar();
}
//get the rank of each individual...
void MOEA::computing_dominate_information(vector <CIndividual*> &pool)
{
    for(int i = 0; i < pool.size(); i++)
    {
	pool[i]->times_dominated = 0;
	pool[i]->ptr_dominate.clear();
	for(int j = 0; j < pool.size(); j++)
	{
	    if(i == j) continue;
	    if( *(pool[i]) < *(pool[j]) ) //the check if pop[i] dominates pop[j], tht '<' is overloaded
	    {
		pool[i]->ptr_dominate.push_back(pool[j]);
	    }
	    else if( *(pool[j]) < *(pool[i]) )
	   {
		pool[i]->times_dominated++;	
	   }
	}
    }
}
//updates the lowest distance factor of the diversity explicitly promoted
void MOEA::update_diversity_factor()
{
	double ratio = ((double) nfes)/max_nfes;
	lowestDistanceFactor = Initial_lowest_distance_factor - Initial_lowest_distance_factor*(ratio/0.9);
}
void MOEA::reproduction(vector<CIndividual> &population, vector<CIndividual> &child_pop)
{
   //binary tournament selction procedure
   binary_tournament_selection(population, child_pop);
   //recombination of the individuals, through the SBX code (taken from the nsga-ii code), also the evaluation of the population is performed
   recombination(child_pop); 
}
void MOEA::recombination(vector<CIndividual> &child_pop)
{
   vector<CIndividual> child_pop2 = child_pop;
	
   for(int i = 0; i < child_pop.size(); i+=2)
    {
       int indexa = i;//int(rnd_uni(&rnd_uni_init)*pops);
       int indexb = i+1;//int(rnd_uni(&rnd_uni_init)*pops);	
       real_sbx_xoverA( child_pop2[indexa], child_pop2[indexb], child_pop[i], child_pop[i+1]);//the crossover probability and index distribution eta are configured in the global.h file
       realmutation(child_pop[i]); //the index distribution (eta) and  mutation probability are configured in the global.h file
       realmutation(child_pop[i+1]);
       child_pop[i].obj_eval();
       child_pop[i+1].obj_eval();
    }
}
void MOEA::binary_tournament_selection(vector<CIndividual > &population, vector<CIndividual> &child_pop)
{
   for(int i = 0; i < population.size(); i++)
	{
	   int indexa = int(rnd_uni(&rnd_uni_init)*pops);
	   int indexb = int(rnd_uni(&rnd_uni_init)*pops);
	   if(population[indexa].rank < population[indexb].rank)
	      child_pop[i] = population[indexa];
	   else if(population[indexa].rank > population[indexb].rank)
	      child_pop[i] = population[indexb];
	   else 
	   {
	      child_pop[i] = (rnd_uni(&rnd_uni_init) < 0.5  )? population[indexa] : population[indexb];
	   }	
	}
}
void MOEA::compute_distances(vector<CIndividual *> &candidates, vector<CIndividual *> &survivors)
{
	for(int i = 0; i < candidates.size(); i++)
	{
	    candidates[i]->nearest_variable_distance = INFINITY;
	    candidates[i]->neares_objective_distance = INFINITY;
	   for(int j = 0; j < survivors.size(); j++)
	   {
		candidates[i]->nearest_variable_distance = min( candidates[i]->nearest_variable_distance, distance(candidates[i]->x_var, survivors[j]->x_var));
		candidates[i]->neares_objective_distance = min( candidates[i]->neares_objective_distance, distance_improvement(survivors[j]->y_obj, candidates[i]->y_obj));
	   }
	}	

}
void MOEA::select_best_candidate(vector<CIndividual *> &survivors, vector<CIndividual *> &candidates, vector<CIndividual *> &penalized)
{
	int best_index_lastfront = -1;//the index of current with the farthes improvement distance
	double max_improvement = -INFINITY;
	  for(int i = 0 ; i < candidates.size(); i++)
	    {
		   if(candidates[i]->times_dominated != 0) continue;
			if(  max_improvement < candidates[i]->neares_objective_distance  )
			{
				max_improvement = candidates[i]->neares_objective_distance;
				best_index_lastfront= i;
			}
	    }
	 if(best_index_lastfront == -1) return; //this occurs when the first m-survirvors are dominated bewteen them, thus there are not candidates availables to pick, therefore this iteration is skiped, so in the next iteration will be available some candidates...

	//update distances of Current and penalized
	  for(int i = 0 ; i < candidates.size(); i++)
	   {
		if( i != best_index_lastfront) // Avoid to be updated by itself..
	        {
		 candidates[i]->nearest_variable_distance = min( candidates[i]->nearest_variable_distance, distance(candidates[i]->x_var, candidates[best_index_lastfront]->x_var ) );
		 candidates[i]->neares_objective_distance = min( candidates[i]->neares_objective_distance, distance_improvement(candidates[best_index_lastfront]->y_obj, candidates[i]->y_obj));
		}
	   }
	  for(int i = 0 ; i < penalized.size(); i++)
	  {
		penalized[i]->nearest_variable_distance = min( penalized[i]->nearest_variable_distance, distance( penalized[i]->x_var, candidates[best_index_lastfront]->x_var ) )  ;
		penalized[i]->neares_objective_distance  =  min( penalized[i]->neares_objective_distance, distance_improvement(candidates[best_index_lastfront]->y_obj, penalized[i]->y_obj));
	  }
	  survivors.push_back(candidates[best_index_lastfront]);
	  iter_swap(candidates.begin()+best_index_lastfront, candidates.end()-1);
	  candidates.pop_back();
}
void MOEA::select_farthest_penalized(vector<CIndividual *> &survivors, vector<CIndividual *> &penalized)
{
    	double largestDCN = -INFINITY;
	int index_largestDCN=0;
	for(int i = 0; i < (int)penalized.size(); i++) // get the index of penalized with larges DCN
	{
		if(penalized[i]->nearest_variable_distance >  largestDCN )
		{
			index_largestDCN = i;
			largestDCN = penalized[i]->nearest_variable_distance;
		}
	}

	for(int i = 0 ; i < (int)penalized.size(); i++) //update the nearest distance once that the penalized is moved to candidate (thereafter to survivors)
	{
		if( i != index_largestDCN )
		penalized[i]->nearest_variable_distance = min( penalized[i]->nearest_variable_distance, distance( penalized[i]->x_var, penalized[index_largestDCN]->x_var));
	}	

	survivors.push_back(penalized[index_largestDCN]);
	iter_swap(penalized.begin()+index_largestDCN, penalized.end()-1);
	penalized.pop_back();
}
void MOEA::penalize_nearest(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized)
{
   	for(int i = candidates.size()-1; i >=0; i--)
	{	
		if( candidates[i]->nearest_variable_distance < lowestDistanceFactor )
		{
			penalized.push_back(candidates[i]);
			for(int j = 0; j < candidates[i]->ptr_dominate.size(); j++)
			{
				candidates[i]->ptr_dominate[j]->times_dominated--; //decreasing the times in which survivors is dominated, this since penalized individuals are not considered..
			}
			//remove the candidate with index "i"
			iter_swap(candidates.begin()+i, candidates.end()-1);
			candidates.pop_back();
		}
	}
}
void MOEA::exec_emo(int run)
{
       char filename1[5024];
       char filename2[5024];
		seed = run;
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	nfes      = 0;
	init_population(); //Initialize individuals...

	sprintf(filename1,"%s/POS/POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",currentPATH, strTestInstance,run, seed, nobj);
	sprintf(filename2,"%s/POF/POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",currentPATH, strTestInstance,run, seed, nobj);
	while(nfes < max_nfes )
	{
		evol_population();
		nfes += pops;
        save_front(filename2); //save the objective space information
	cout << lowestDistanceFactor <<endl;
	}
	save_pos(filename1); //save the decision variable space information
        save_front(filename2); //save the objective space information
	population.clear();
}
/*
Load the configuration of  each instance
*/
void MOEA::load_parameter()
{
	char filename[1024];
}
void MOEA::save_front(char saveFilename[1024])
{

    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename,fstream::app|fstream::out );
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population[n].y_obj[k]<<"  ";
//	for(int k=0;k<nobj;k++)
//			fout<<child_pop[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void MOEA::save_pos(char saveFilename[1024])
{
    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename, fstream::app|fstream::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].x_var[k] << "  ";
			//fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<child_pop[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}



#endif
