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

	double distance( vector<double> &a, vector<double> &b);
	vector <CIndividual> population;
	vector<CIndividual> child_pop;	// memory solutions
	void operator=(const MOEA &moea);

//public:
//
//	// algorithm parameters
//	int     pops;          //  the population size
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
	   double factor = (a[i]-b[i]);
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
	   maxd = max(maxd, max(B-A,0.0));
	}
        if(dist == 0.0) return -maxd; //in case that this indicator is zero, this mean that it is a dominated individual...
   return sqrt(dist);
}
void MOEA::init_population()
{

    for(int i=0; i<pops; i++)
	{

		CSubproblem indiv1, indiv2;

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
	vector<CIndividual> survivors, penalized;
	//join the offspring and parent populations
	vector<CIndividual> candidates(population.size()+child_pop.size());
	candidates.insert(candidates.end(), population.begin(), population.end());
	candidates.insert(candidates.end(), child_pop.begin(), child_pop.end());

	fast_sort_non_dominated(candidates); //computing the rank of each candidate individual...

	//select the "best" individuals of candidate in survivors...
	select_first_survivors(survivors, candidates);
	//update the diversity-factor-parameter...	
	update_diversity_factor();


 
       //Pre-computing the neares distances both objective and decision spaces..
	
	for(int i = 0; i < current.size(); i++)
	{
	    current[i].nearest_variable_distance = INFINITY;
	    current[i].neares_objective_distance = INFINITY;
	   for(int j = 0; j < survivors.size(); j++)
	   {
		current[i].nearest_variable_distance = min( current[i].nearest_variable_distance, distance(current[i].x_var, survivors[j].x_var));
		current[i].neares_objective_distance = min( current[i].neares_objective_distance, distance_improvement(current[i].y_obj, survivors[j].y_obj));

	   }
	}	
	while( survivors.size() < pops )
	{
	  penalize_nearest(current, penalized);//penalize the nearest individuals.. 
	  if(current.empty())	  
	     select_farthest_penalized(current, penalized);//in case that all the individuals are penalized pick up the farstest

	  fast_update_non_domianted_front(survivors, current, ); //update the rank of each candidate whitout penalized
	  select_best_candidate(survivors, current, penalized); // the best candidate is selected considering the improvemente distance, and the rank..
	}
	fast_sort_non_dominated(survivors); //computing the rank of each candidate individual...
	survivors_reproduction(survivors, ); //generate a new population considering the survivors individuals...
}
void select_best_candidate(vector<CIndividual> &survivors, vector<CIndividual> &current, vector<CIndividual> &penalized, vector<CIndividual> &lowes_front_current)
{
	int best_index_lastfront= -1;//the index of current with 
	 double max_improvement = -INFINITY;
	  for(int i = 0 ; i < lowes_front_current.size(); i++)
		{
			int inner_index = lowes_front_current[i];	
			if(  max_improvement < curret[inner_index].neares_objective_distance  )
			{
				max_improvement =  current[best_index_lastfront].neares_objective_distance;
				best_index_lastfront= inner_index;
			}
		}
	//update distances of Current and penalized
	  for(int i = 0 ; i < current.size(); i++)
	   {
		if( i != best_index_lastfront) // Avoid updated itself..
	        {
		 current[i].nearest_variable_distance = min( current[i].nearest_variable_distance, distance(current[i], current[best_index_lastfront] ) );
		 current[i].neares_objective_distance = min( current[i].neares_objective_distance, distance_improvement(current[best_index_lastfront], current[i]));
		}
	   }
	  for(int i = 0 ; i < penalized.size(); i++)
	  {
		penalized[i].nearest_variable_distance = min( penalized[i].nearest_variable_distance, distance( penalized[i], current[best_index_lastfront] ) )  ;
		penalized[i].neares_objective_distance  =  min( penalized[i].neares_objective_distance, distance_improvement(current[best_index_lastfront], penalized[i]));
	  }
	  survivors.push_back(current[best_index_lastfront]);
	  iter_swap(current.begin()+best_index_lastfront, current.end()-1);
	  current.pop_back();
}
void select_farthest_penalized(vector<CIndividual> &current, vector<CIndividual> &penalized)
{
    	double largestDCN = penalized[0].nearest_variable_distance;
	int index_largestDCN=0;
	for(int i = 0; i < (int)penalized.size(); i++) // get the index of penalized with larges DCN
	{
		if(penalized[i].nearest_variable_distance >  largestDCN )
		{
			index_largestDCN = i;
			largestDCN = penalized[i].nearest_variable_distance;
		}
	}

	for(int i = 0 ; i < (int)penalized.size(); i++) //update the nearest distance once that the penalized is moved to candidate (thereafter to survivors)
	{
		if( i != index_largestDCN )
		penalized[i].nearest_variable_distance = min( penalized[i].nearest_variable_distance, distance( penalized[i].x_var, penalized[index_largestDCN].x_var));
	}	

	for(int i = 0; i < penalized[index_largestDCN].index_dominate.size(); i++) //update the dominate count 
          penalized[index_dominate].index_dominate[i].times_dominated++;

	current.push_back(penalized[index_dominate]);
	iter_swap(penalized.begin()+index_dominate, penalized.end()-1);
	penalized.pop_back();
}
void penalize_nearest(vector<CIndividual> &current, vector<CIndividual> &penalized)
{
   	for(int i = survivors.size()-1; i >=0; i--)
	{	
		if( survivors[i].nearest_variable_distance < lowestDistanceFactor )
		{
			penalized.push_back(survivors[i]);
			for(int j = 0; j < survivors[i].index_dominate.size(); j++)
			{
				survivors[i].index_dominate[j]->times_dominated--; //decreasing the times in which survivors is dominated, this since penalized individuals are not considered..
			}
			//remove the survivor with index "i"
			iter_swap(survivors.begin()+i, survivors.end()-1);
			survivors.pop_back();
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

	sprintf(filename1,"/home/joel.chacon/Current/MyResearchTopics/MOEA-Improvement/HSI-MOEA/POS/POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",strTestInstance,run, seed, nobj);
	sprintf(filename2,"/home/joel.chacon/Current/MyResearchTopics/MOEA-Improvement/HSI-MOEA/POF/POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",strTestInstance,run, seed, nobj);

	while(nfes < max_nfes )
	{
		evol_population();
		nfes += pops;
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

cout << max_nfes<<endl;
max_gen= max_nfes/pops;
prob=0.9;
rate=0.8;

}
void MOEA::save_front(char saveFilename[1024])
{

    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename,fstream::app|fstream::out );
	for(int n=0; n<pops; n++)
	{
	//	for(int k=0;k<nobj;k++)
	//		fout<<best[n].y_obj[k]<<"  ";
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.y_obj[k]<<"  ";
	for(int k=0;k<nobj;k++)
			fout<<child_pop[n].indiv.y_obj[k]<<"  ";
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
			fout<<population[n].indiv.x_var[k] << "  ";
			//fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<best[n].x_var[k]<<"  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<child_pop[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}



#endif
