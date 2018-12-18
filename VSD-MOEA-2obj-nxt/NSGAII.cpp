/**
	Autores:
		 Carlos Segura González
		 Joel Chacón Castillo
	Fecha: 31/10/2016
	Descripción:
**/
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <queue>
#include <algorithm>
#include "Measures.hpp"
#include "Benchmark.hpp"
#include "Individual.hpp"
#include "NSGAII.hpp"
inline ostream & operator << ( ostream &out,vector<Individual> Pool)
{
    for(int i = 0; i < (int)Pool.size(); i++)
    {
        for(int j = 0; j < Pool[i].getObjectives().getNumberObjectives() ; j++)
         out << Pool[i].getObjectives().SpaceObjectives[j].Fitness<<" ";
         out << endl;
    }
    return out;
}
inline ostream & operator << ( ostream &out,vector<double> data)
{
    for(int i = 0; i < (int)data.size(); i++)
    {
        out << data << endl;
    }
    return out;
}
inline vector<Individual> operator+(vector<Individual> PoolA, vector<Individual> PoolB )
{
    vector<Individual> C;
    for(int i=0 ; i < (int)PoolA.size(); i++ )
        C.push_back(PoolA[i]);
    for(int i = 0; i < (int)PoolB.size(); i++)
        C.push_back(PoolB[i]);
    return C;
}
NSGAII::NSGAII(int SizePool, double ProbCross, int NBitsMut, Benchmark &ObjBenchmark, int NumberGenerations, double ParameterDInit ,string ProblemName, int Semilla, string Ruta, string Label, int IdCrossover, int IdMutation, int PeriodReport  )
{
    srand(Semilla);

    this->D = 0;
    this->SizePool = SizePool;
    this->NumberObjectives = ObjBenchmark.getNObjectives();
    this->Dimension = ObjBenchmark.getDimension();
    this->ProbCross = ProbCross;
    this->NBitsMut = NBitsMut;
    this->ObjBenchmark = &ObjBenchmark;
    this->NumberGenerations = NumberGenerations;
    this->CurrentGeneration=0;
    this->IdCrossover = IdCrossover;
    this->IdMutation = IdMutation;
    this->PeriodReport = PeriodReport;
    this->DI = ParameterDInit;
    this->VariableRepresentation = BINARY_ENCODE;
    this->FilenameGenotype = Ruta+"/"+ProblemName+"/Genotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFenotype =  Ruta+"/"+ProblemName+"/Fenotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFitness = Ruta+"/"+ProblemName+"/Fitness_"+Label+"_"+to_string(Semilla)+".txt";
    this->FileNameSummary = Ruta+"/"+ProblemName+"/SummaryFront_"+Label+"_"+to_string(Semilla)+".txt";
    this->ProblemName = ProblemName;
	//Eliminar el archivo en caso de que existan
	ofstream SummaryGen, SummaryFen, SummaryFit;
    	SummaryGen.open (this->FilenameGenotype);
    	SummaryFen.open (this->FilenameFenotype);
    	SummaryFit.open (this->FilenameFitness);
	SummaryGen << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;
	SummaryFen << "Total_Generations\t" << this->NumberGenerations<< "\tDimension\t"<< this->Dimension <<endl;
	SummaryFit << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;
 //this->FileNameSummary = FileNameSummary ;

}
NSGAII::NSGAII(int SizePool, double ProbCross, double ProbMut ,Benchmark &ObjBenchmark, int NumberGenerations, double ParameterDInit , string ProblemName, int Semilla, string Ruta, string Label, int IdCrossover, int IdMutation, int PeriodReport  )
{
    srand(Semilla);
    this->DI = ParameterDInit;
    this->D = 0;
    this->SizePool = SizePool;
    this->NumberObjectives = ObjBenchmark.getNObjectives();
    this->Dimension = ObjBenchmark.getDimension();
    this->ProbCross = ProbCross;
    this->ProbMutation = ProbMut;
    this->ObjBenchmark = &ObjBenchmark;
    this->NumberGenerations = NumberGenerations;
    this->CurrentGeneration=0;
    this->IdCrossover = IdCrossover;
    this->IdMutation = IdMutation;
    this->PeriodReport = PeriodReport;
    this->VariableRepresentation = REAL_ENCODE;
    this->FilenameGenotype = Ruta+"/"+ProblemName+"/Genotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFenotype =  Ruta+"/"+ProblemName+"/Fenotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFitness = Ruta+"/"+ProblemName+"/Fitness_"+Label+"_"+to_string(Semilla)+".txt";
    this->FileNameSummary = Ruta+"/"+ProblemName+"/SummaryFront_"+Label+"_"+to_string(Semilla)+".txt";
    this->ProblemName = ProblemName;
	//Eliminar el archivo en caso de que existan
	ofstream SummaryGen, SummaryFen, SummaryFit;
   if(VariableRepresentation == BINARY_ENCODE)
    	SummaryGen.open (this->FilenameGenotype);
    	SummaryFen.open (this->FilenameFenotype);
    	SummaryFit.open (this->FilenameFitness);
   if(VariableRepresentation == BINARY_ENCODE)
	SummaryGen << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;
	SummaryFen << "Total_Generations\t" << this->NumberGenerations<< "\tDimension\t"<< this->Dimension <<endl;
	SummaryFit << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;
    //IdealVector.resize(this->NumberObjectives, INFINITY);
    //NadirVector.resize(this->NumberObjectives, -INFINITY);

 //this->FileNameSummary = FileNameSummary ;
}
NSGAII::~NSGAII()
{
   FreePool(Pool_Q);
   FreePool(Survivors);
   //FreePool(archive);

//   free(ReferencePoint);
}
void NSGAII::ControlParameterD()
{

	double TElapsed = this->CurrentGeneration;
	double TEnd = this->NumberGenerations;
	double ratio = (sqrt(Dimension))*TElapsed/TEnd;
	//En esta sección se controla
	//el decremento del parámetro D
	this->D = DI - DI * ( TElapsed / (TEnd*0.9)  );

  //this->D = DI - (DF - DI )* ( 2*TElapsed / TEnd  );
  this->D = ( 0 > this->D)? -1.0: this->D ;
//	if( TElapsed < TEnd*0.5)
//	this->D = ratio;else D = 1.0 - ratio;
}
/**
	Se regresa a los individuos más cercanos
	al vector ideal..
	Closest Ideal Vector
**/
void NSGAII::CIV(vector<ptrIndividual> &Set, vector<ptrIndividual> &Extremal)
{
    Extremal.clear();
	  
	//Definir maximizar o minimizar para cada objetivo
	vector<int> KindObjective =  ObjBenchmark->getTypeDuality();
	vector<double> Best(this->NumberObjectives);
	//vector<int> BestIndex(this->NumberObjectives);
	vector<int> BestIndex;
	

	///Select the best improvement distance survivor....
	vector<bool> grid(Set.size(), false);
	for(int m = 0; m < this->NumberObjectives; m++)
	{
	        vector<int> subset;
		int indxmaxim;
		   double bestvector = INFINITY;
		   for(int i = 0; i <  Set.size(); i++)
		   {	
		   	  if(grid[i])continue;
			  double s = 0.0;	
			  double maxv = -INFINITY;
			  for(int k = 0; k < this->NumberObjectives; k++)
			  {
				double fi = fabs(Set[i]->getObjectiveValue(k));
				s += fi;
				double ti = (k==m)?fi:1e5*fi;
				if(ti > maxv)	maxv=ti;
			  }
			   maxv = maxv + 0.0001*s;
			  if(bestvector > maxv)
			  { indxmaxim = i; bestvector = maxv;}
		
		    }
		grid[indxmaxim] = true;
		BestIndex.push_back(indxmaxim);
	}
	
	/**
		Agregar los individuos que corresponden a los índices...
	**/
		std::vector<int>::iterator End;
		End = std::unique (BestIndex.begin(), BestIndex.end());
		BestIndex.resize( std::distance(BestIndex.begin(),End));
		//Ordenar los indices en orden ascendente y no tener problemas al momento de eliminar los elementos del vector original
		sort(BestIndex.begin(), BestIndex.end());

		for(int i = BestIndex.size()-1; i >= 0 ; i--)
		{
			int index = BestIndex[i];
			Extremal.push_back( Set[index]);
//			Es el mismo resultado que "Set.erase(Set.begin()+index);"
			iter_swap(Set.begin()+index, Set.end()-1);
			Set.pop_back();
		}

}
double NSGAII::EuclideanVariableDistance(ptrIndividual &indA, ptrIndividual &indB)
{
	double TElapsed = this->CurrentGeneration;
	double TEnd = this->NumberGenerations;
	  double ratio =  ( 2.0*TElapsed / TEnd  );


	double Sum = 0.0;
	bool flag = true;
	for(int d = 0; d < indA->getDimension(); d++)
	{
		double A = (indA->getVariable(d)) / (indA->getMaximum(d) - indA->getMinimum(d));
		double B = (indB->getVariable(d)) / (indA->getMaximum(d) - indA->getMinimum(d));
		Sum += pow(A-B,2);
	}

	return sqrt(Sum);
	double normA =0.0, normB=0.0, dotAB=0.0;

	for(int d = 0; d < indA->getDimension() ; d++)
	{
//		Sum += pow((indA->getVariable(d)- indB->getVariable(d)) / (indA->getMaximum(d) - indA->getMinimum(d)),2);
		double A = (indA->getVariable(d)) / (indA->getMaximum(d) - indA->getMinimum(d));
		double B = (indB->getVariable(d)) / (indA->getMaximum(d) - indA->getMinimum(d));
		normA += pow(A,2);
		normB += pow(B,2);
		dotAB += A*B;
	}
	normA=sqrt(normA);	
	normB=sqrt(normB);	
	double theta = acos(dotAB/(normA*normB));
   	return theta;
}
double NSGAII::EuclideanObjectiveDistance(ptrIndividual &indA, ptrIndividual &indB)
{
	double Sum = 0;
	for(int d = 0; d < indA->getNumberObjectives(); d++)
	{
		Sum += pow(indA->getObjectiveValue(d)- indB->getObjectiveValue(d),2);
	}
	return sqrt(Sum);
}
double NSGAII::ImprovementDistance(ptrIndividual &Reference, ptrIndividual &X)
{
	double Sum = 0.0;
	double maxd = -INFINITY;
	for(int d = 0; d < Reference->getNumberObjectives(); d++)
	{
		double A = Reference->getObjectiveValue(d);
		double B = X->getObjectiveValue(d);

		Sum += pow( max( 0.0, A-B),2);
		maxd = max(maxd, max(B-A,0.0));
	}
        if(Sum ==0.0) return -maxd; 
	return Sum;
	double normA =0.0, normB=0.0, dotAB=0.0;

	for(int d = 0; d < Reference->getNumberObjectives(); d++)
	{
		double A = Reference->getObjectiveValue(d);
		double B = X->getObjectiveValue(d);
		normA += pow(A,2);
		normB += pow(B,2);
		dotAB += A*B;
	}
	normA=sqrt(normA);	
	normB=sqrt(normB);	
	double theta = acos(dotAB/(normA*normB));
   	return theta;
	return Sum+theta;
}
double NSGAII::ImprovementDistance2(ptrIndividual &Reference, ptrIndividual &X)
{
	double Sum = 0.0;
	double maxd = -INFINITY;
	for(int d = 0; d < Reference->getNumberObjectives(); d++)
	{
		double A = Reference->getObjectiveValue(d);
		double B = X->getObjectiveValue(d);

		Sum += pow( max( 0.0, A-B),2);
		maxd = max(maxd, max(B-A,0.0));
	}
	return Sum;	
}
/**
	Distance Closest Neighbor
  Se busca la distancia mínima entre el conjunto de referencia y el conjunto, los
  individuos con menor distancia se mueven a penalized.
  Nota:  Se puede implementar una algoritmo mas eficiente este es O(n^2)
**/
void NSGAII::DCN(vector<ptrIndividual> &Set, vector<ptrIndividual> &Penalized)
{

	for(int i = Set.size()-1; i >=0; i--)
	{	
		if( Set[i]->getVarDistance() < D)
		{
			Penalized.push_back(Set[i]);
			for(int j = 0; j < Set[i]->WhoDominate.size(); j++)
			{
			////	if(Set[i]->getRank() < Set[i]->WhoDominate[j]->getRank() )
				//if(Set[i]->getObjectiveValue(0) !=  Set[i]->WhoDominate[j]->getObjectiveValue(0) && Set[i]->getObjectiveValue(1) !=  Set[i]->WhoDominate[j]->getObjectiveValue(1))
				Set[i]->WhoDominate[j]->TimesIsDominated--;
			}
			iter_swap(Set.begin()+i, Set.end()-1);
			Set.pop_back();
		}
	}
}




/**
	Objective Distance Closest Negihbor
	Descripcion:
		Obtiene el individuo con mayor DCN en el espacio objetivo, 
		los indices de los individuos considerados pertenecen a Current y al menor 
		rank (frente mas cercano)
**/
void NSGAII::BestODCN( vector<ptrIndividual> &CurrentSet, vector<ptrIndividual> &Survivors, vector<ptrIndividual> &Penalized, vector<int> &FirstFrontCurrent)
{
	 int BestIndex= -1;//Improvement_Current(FirstFrontCurrent, CurrentSet, Survivors);
	 double maximpr = -INFINITY;
	  for(int i = 0 ; i < FirstFrontCurrent.size(); i++)
		{
			int indexf = FirstFrontCurrent[i];	
			if(  maximpr < CurrentSet[indexf]->getObjDistance()  )
			{
				 maximpr =  CurrentSet[indexf]->getObjDistance();
				BestIndex = indexf;
			}
		}
//         if(maximpr == 0.0) BestIndex = FirstFrontCurrent[0];
	
	  
	//update distances of Current and penalized
	  for(int i = 0 ; i < CurrentSet.size(); i++)
	   {
		if( i != BestIndex) // Avoid updated itself..
	        {
		 CurrentSet[i]->setVariableDistance( min( CurrentSet[i]->getVarDistance() , EuclideanVariableDistance(CurrentSet[i], CurrentSet[BestIndex] ) )  );
		 CurrentSet[i]->setObjDistance( min( CurrentSet[i]->getObjDistance() , ImprovementDistance(CurrentSet[BestIndex], CurrentSet[i])));
		}
	   }
	  for(int i = 0 ; i < Penalized.size(); i++)
	  {
		Penalized[i]->setVariableDistance( min( Penalized[i]->getVarDistance() , EuclideanVariableDistance( Penalized[i], CurrentSet[BestIndex] ) )  );
		Penalized[i]->setObjDistance( min( Penalized[i]->getObjDistance() , ImprovementDistance(CurrentSet[BestIndex], Penalized[i])));
	  }
	  Survivors.push_back(CurrentSet[BestIndex]);
	  iter_swap(CurrentSet.begin()+BestIndex,CurrentSet.end()-1);
	  CurrentSet.pop_back();

//	  CurrentSet.erase(CurrentSet.begin()+BestIndex);
}
/**
	Get the best penalized individual
  Se agrega el individuo con mayor distancia de los considerados en el DCN.
**/
void NSGAII::AddCurrentBestIndividualPenalized(vector<ptrIndividual> &Penalized, vector<ptrIndividual> &CurrentSet)
{
	if(Penalized.empty()) return ;
	double BestPenalized = Penalized[0]->getVarDistance();
	int IndexBest=0;
	for(int i = 0; i < (int)Penalized.size(); i++)
	{
		if(Penalized[i]->getVarDistance() >  BestPenalized )
		{
			IndexBest = i;
			BestPenalized = Penalized[i]->getVarDistance();
		}
	}
	for(int i = 0 ; i < (int)Penalized.size(); i++)
	{
		if( i != IndexBest )
		Penalized[i]->setVariableDistance( min( Penalized[i]->getVarDistance() , EuclideanVariableDistance( Penalized[i], Penalized[IndexBest] ) )  );
	//	Penalized[i]->setObjDistance( min( Penalized[i]->getObjDistance() , ImprovementDistance(Penalized[IndexBest], Penalized[i]))); //is needed for the binary tournament..
	}	
	for(int i = 0; i < Penalized[IndexBest]->WhoDominate.size(); i++) Penalized[IndexBest]->WhoDominate[i]->TimesIsDominated++;
	
//	for(int i = 0 ; i < CurrentSet.size(); i++)
//	   {
//		 CurrentSet[i]->setVariableDistance( min( CurrentSet[i]->getVarDistance() , EuclideanVariableDistance(CurrentSet[i], Penalized[IndexBest] ) )  );
//		 CurrentSet[i]->setObjDistance( min( CurrentSet[i]->getObjDistance() , ImprovementDistance(Penalized[IndexBest], CurrentSet[i])));
//	   }
	CurrentSet.push_back(Penalized[IndexBest]);
	//Survivors.push_back(Penalized[IndexBest]);
	//Penalized.erase(Penalized.begin()+IndexBest);
	iter_swap(Penalized.begin()+IndexBest, Penalized.end()-1);
	Penalized.pop_back();
}
/**
	Se implementa una union entre el set1 y el set2, la informacion que se almacena en Result son
	los apuntadores de los individuos y no se realizan copias de individuos.
**/
void NSGAII::Union(vector<ptrIndividual> &Set1, vector<ptrIndividual> &Set2, vector<ptrIndividual> &Result)
{
    Result.resize(Set1.size()+Set2.size(), NULL);
    for(int i = 0; i < (int)Set1.size(); i++)
    {
        Result[i] = Set1[i];
    }
    for(int j = 0; j < (int)Set2.size(); j++)
    {
        Result[j+Set1.size()] = Set2[j];
    }
}
void NSGAII::Init_NSGAII()
{
    InitializePool(Survivors);
    //InitializePool(archive);
//    this->ReferencePoint = new Individual;
//    (* ReferencePoint) = (*(Survivors[0]));
    //Inicialmente Pool_Q es un conjunto vacío...
    while(CriterionStop())
    {

		//cout << "Generation: "<<CurrentGeneration<<endl;

          vector<ptrIndividual> CurrentSet1, CurrentSet;

	  
          //Se realiza la unión entre Survivors y PoolQ
          //se almacena en CurrentSet (solo los apuntadores)
          Union(Survivors, Pool_Q, CurrentSet);
//          Union(CurrentSet1, archive, CurrentSet);
			  
 
          //Se insertan los individuos extremos
          //en survivors (primero se vacía survivors)
          CIV(CurrentSet, Survivors);
          //Control del parámetro
          ControlParameterD();

	  
          vector<ptrIndividual> Penalized;
	 //Compute the improvement distance matrix.....
	   for(int i = 0 ; i < CurrentSet.size(); i++)
	   {
		CurrentSet[i]->setVariableDistance(INFINITY);
		CurrentSet[i]->setObjDistance(INFINITY);
	   	for(int j = 0 ; j < Survivors.size(); j++)
		{
		  CurrentSet[i]->setVariableDistance( min(  CurrentSet[i]->getVarDistance(), EuclideanVariableDistance( CurrentSet[i],  Survivors[j])));
		  CurrentSet[i]->setObjDistance( min(  CurrentSet[i]->getObjDistance(),  ImprovementDistance( Survivors[j], CurrentSet[i])));
		}
	   }
		vector<int>KFrontCurrent;
          	Fast_Non_Dominated_Sort_Current(CurrentSet, Survivors, KFrontCurrent);

	  int Rank = 0;
          while( (int)Survivors.size() < this->SizePool)
          {
		double TElapsed = this->CurrentGeneration;
		double TEnd = this->NumberGenerations;
	  	double ratio =  ( TElapsed / TEnd  );

	//	if(Survivors.size()>SizePool*ratio && Penalized.empty()) D=100;
                DCN(CurrentSet, Penalized);
              if(CurrentSet.empty())
		{
		    AddCurrentBestIndividualPenalized(Penalized, CurrentSet);
		}
	  	    //Penalize the individuals...
		    KFrontCurrent.clear();
	  	    FastComputingKfront(Survivors, CurrentSet, KFrontCurrent, Rank);	
		    BestODCN(CurrentSet, Survivors, Penalized, KFrontCurrent);
          }
	    vector<ptrIndividual> nothing;
            Fast_Non_Dominated_Sort_Current(nothing, Survivors, KFrontCurrent);
	    Recycle_Memory(Penalized, CurrentSet, Pool_Q,  archive);
//            updateArchive(Survivors, Pool_Q, archive);
            New_Pop(Survivors, Pool_Q);

          if( ! (this->CurrentGeneration % PeriodReport) )//|| this->CurrentGeneration < 25000)
          {
              ExportIndividualsFile(this->Survivors, this->FilenameGenotype, this->FilenameFenotype, this->FilenameFitness, this->CurrentGeneration);
          //    ExportIndividualsFile(this->archive, this->FilenameGenotype+"_archive", this->FilenameFenotype+"_archive", this->FilenameFitness+"_archive", this->CurrentGeneration);
             // ExportIndividualsFile(this->Pool_Q, this->FilenameGenotype, this->FilenameFenotype, this->FilenameFitness, this->CurrentGeneration);
              //this->PlotInterfaceRSpaceObjective();
              //this->PlotInterfaceRSpaceVariables();
              //cout << "Generacion...."<<CurrentGeneration<<endl;
		/**
			Para realizar un analisis de los penalizados y del parametro D teorico
		**/
		//cout << CurrentGeneration << " "<<Penalized.size() << " " <<this->D<<endl;
		cout << "Generation: "<<CurrentGeneration<<endl;
          }
    }

          End();
}
void NSGAII::updateArchive(vector<ptrIndividual> &Survivors, vector<ptrIndividual> &Pool_Q, vector<ptrIndividual> &archive)
{
    vector<ptrIndividual> joinall;
        for(int i = 0; i < archive.size(); i++)
    {
	joinall.push_back(Survivors[i]);
	if(!Pool_Q.empty())
	joinall.push_back(Pool_Q[i]);
	joinall.push_back(archive[i]);
    }
 	vector<ptrIndividual> survivors; 
	int k = 0;
        vector<double> mem_improvement_distance(joinall.size(), INFINITY);
        vector<double> mem_improvement_distance2(joinall.size(), INFINITY);
        vector< vector<double> > objind (archive.size(), vector<double>(NumberObjectives));
        vector< vector< double > > xvarind(archive.size(), vector<double>(Dimension));

	vector<bool> grid(joinall.size(), false); //select best extremes...
        vector<bool> selected(joinall.size(), true);
	for(int m = 0; m < this->NumberObjectives; m++)
	{
	        vector<int> subset;
		priority_queue< pair<double, int> > pq;
		   for(int i = 0; i <  joinall.size(); i++)
		   {	
		   	  if(grid[i])continue;
		   	pq.push( make_pair(-joinall[i]->getObjectiveValue(m), i));
		   }
		double last = -pq.top().first;
		subset.push_back(pq.top().second);
		pq.pop();
		while( last+1e-3 > -pq.top().first && !pq.empty())
		{
		      subset.push_back(pq.top().second);
		      pq.pop();
		}
		int indxmaxim = -1;

		double MAXI=-INFINITY;
		for(int i = 0; i < subset.size(); i++)
		{
		    double minnearest = 0.0;
		   for(int j = 0; j < subset.size(); j++)
		   {
			minnearest +=ImprovementDistance2( joinall[subset[j]], joinall[subset[i]] );
		   }
			if( minnearest > MAXI)
			{
			   MAXI = minnearest;
			   indxmaxim = subset[i];
			}
		}
		grid[indxmaxim] = true;
		survivors.push_back(joinall[indxmaxim]);
		selected[indxmaxim] = false;
		for(int i = 0 ; i < NumberObjectives; i++) objind[k][i] = joinall[indxmaxim]->getObjectiveValue(i);
	 	for(int i = 0 ; i < Dimension; i++) xvarind[k][i]= joinall[indxmaxim]->getVariable(i);
	   k++;
	//  break;
	}

   //memory nearest improvement distances...
   for(int i = 0 ; i < joinall.size(); i++)
	   {
		if(!selected[i])continue;
	   	for(int j = 0 ; j < survivors.size(); j++)
		{
		  mem_improvement_distance[i] = ( min(  mem_improvement_distance[i],  ImprovementDistance( survivors[j], joinall[i])));
		  mem_improvement_distance2[i] = ( min(  mem_improvement_distance2[i],  ImprovementDistance2( survivors[j], joinall[i])));
		}
	   }
    while(k < archive.size())
    {
	double maxd = -INFINITY;
	double maxd2 = -INFINITY;
	int index = -1;
	int index2 = -1;
	for(int i = 0 ; i < mem_improvement_distance.size(); i++)
	{
		if(!selected[i])continue;
		if( maxd < mem_improvement_distance[i]  )
		{
		   index = i;
		   maxd = mem_improvement_distance[i];
		}
		if( maxd2 < mem_improvement_distance2[i]  )
		{
		   index2 = i;
		   maxd2 = mem_improvement_distance2[i];
		}
	}
        selected[index] = false;
	survivors.push_back(joinall[index]);
        //update vector distances...
	for(int i = 0 ; i < mem_improvement_distance.size(); i++)
	   {
		if( i != index ) // Avoid updated itself..
	        {
		 mem_improvement_distance[i] = ( min( mem_improvement_distance[i] , ImprovementDistance(joinall[index], joinall[i])));
		}
		if( i != index2 ) // Avoid updated itself..
	        {
		 mem_improvement_distance2[i] = ( min( mem_improvement_distance2[i] , ImprovementDistance2(joinall[index], joinall[i])));
		}
	   }
	
	 for(int i = 0 ; i < NumberObjectives; i++) objind[k][i] = joinall[index]->getObjectiveValue(i);
	 for(int i = 0 ; i < Dimension; i++) xvarind[k][i]= joinall[index]->getVariable(i);
 	 
	k++;	
    }
	for(int j = 0 ; j < archive.size(); j++)
	{
	  for(int i = 0 ; i < NumberObjectives; i++)  archive[j]->setObjective(i, objind[j][i] );
	  for(int i = 0 ; i < Dimension; i++) archive[j]->setFenotype(i, xvarind[j][i] );
	}

	  
   
}
void NSGAII::Recycle_Memory(vector<ptrIndividual> &Penalized, vector<ptrIndividual> &CurrentSet, vector<ptrIndividual> &Pool_Q,  vector<ptrIndividual> &archive)
{
		vector<ptrIndividual> recycling_address;
		for(int i = 0; i < Penalized.size(); i++)  recycling_address.push_back(Penalized[i]);
		for(int j = 0; j < CurrentSet.size(); j++) recycling_address.push_back(CurrentSet[j]);

		int k = 0;
		while( k < Pool_Q.size() ) //recycling memory address available
		{
		    Pool_Q[k] = recycling_address[k];
		    k++;
		}
	//	int index = 0;
	//	while( index < archive.size() ) //recycling memory address available
	//	{
	//	    archive[index] = recycling_address[k];
	//	    k++;
	//	    index++;
	//	}
}
/**
    Agregar información al final del archivo...
**/
void NSGAII::InitializePool(vector<ptrIndividual> &Survivors)
{
    Survivors.resize(this->SizePool);
    for(int i = 0; i < this->SizePool; i++)
    {
        Survivors[i] = new Individual;
        Survivors[i]->InitializeIndividual(this->ObjBenchmark);
    }
}
void NSGAII::FreePool(vector<ptrIndividual> &Pool)
{
  for(int i = 0; i < (int)Pool.size(); i++)
    delete (Pool[i]);
    Pool.clear();
}
void NSGAII::New_Pop(vector <ptrIndividual> &Survivors, vector <ptrIndividual> &Pool_Q)
{
    if( Pool_Q.empty() )
      InitializePool(Pool_Q);

    Selection(Survivors, Pool_Q);
    Recombination(Pool_Q);
	if(this->VariableRepresentation == BINARY_ENCODE)
		BitMutation(Pool_Q ,this->NBitsMut);
	else if(this->VariableRepresentation == REAL_ENCODE)
	{
		switch(this->IdMutation)
		{
		   case _POLYNOMIAL:
			PolynomialMutation(Pool_Q, ProbMutation);
		   break;
		   case _NORMAL:
		   	NormallyDistribuitedMutation(Pool_Q);
		   break;
		   case _RANDOM:
			RandomMutation(Pool_Q);
		   break;
		}
	}
    Evaluation(Pool_Q);
}
void NSGAII::Selection(vector <ptrIndividual> &Survivors, vector <ptrIndividual> &Pool_Q)
{

    for(int i = 0; i < (int)Survivors.size(); i++)
    {
        int Index1 = rand()%Survivors.size();
        int Index2 = rand()%Survivors.size();

        if( Survivors[Index1]->getRank() < Survivors[Index2]->getRank() )
        //if( Survivors[Index1]->ObjectiveDistance > Survivors[Index2]->ObjectiveDistance )
        {
            *(Pool_Q[i]) = *(Survivors[Index1]);
            //Pool_Q.push_back(Survivors[Index1]);
        }
        //else if(Survivors[Index1]->ObjectiveDistance == Survivors[Index2]->ObjectiveDistance)
        else if(Survivors[Index1]->getRank() == Survivors[Index2]->getRank())
        {
				if (getRandom(0,1.0) <= 0.5)
			{
			  *(Pool_Q[i]) = *(Survivors[Index1]);
			    // Pool_Q.push_back(Survivors[Index1]);
			}else
			*(Pool_Q[i]) = *(Survivors[Index2]);
			    // Pool_Q.push_back(Survivors[Index2]);
        }
        else
        {
          *(Pool_Q[i]) = *(Survivors[Index2]);
            // Pool_Q.push_back(Survivors[Index2]);
        }
    }
}
void NSGAII::Recombination(vector <ptrIndividual> &Pool_Q)
{
	if(this->VariableRepresentation == BINARY_ENCODE)
		UniformCross(Pool_Q, Pool_Q, this->ProbCross);
	else if(this->VariableRepresentation == REAL_ENCODE)
	{
		
			double TElapsed = this->CurrentGeneration;
			double TEnd = this->NumberGenerations;
			double ProbVariable  = 1.0-( TElapsed / TEnd  );

			switch(this->IdCrossover)
			{
			   case _SBX:
				SBX(Pool_Q, Pool_Q, this->ProbCross);
			   break;
			   case _SBXH:
				SBXHybrid(Pool_Q, Pool_Q, this->ProbCross);
			   break;
			   case _SBXHD:
				ProbVariable = max(0.5, ProbVariable);
				if(CurrentGeneration%100 == 0)
				cout << ProbVariable << endl;
				SBXHybridDynamic(Pool_Q, Pool_Q, this->ProbCross, ProbVariable);
			   break;
			   case _BLX:
				BLX(Pool_Q, Pool_Q, this->ProbCross);
			   break;
			   case _FUZZY:
				FuzzyCrossOver(Pool_Q, Pool_Q, this->ProbCross);
			   break;
			   case _HBX:
				HBX(Pool_Q, Pool_Q, this->ProbCross);
			   break;
			   case _DE:
         		        DE(Pool_Q, Pool_Q, this->ProbCross); //Simulated Binary Crossover
			   break;

			}
	}
}
void NSGAII::Evaluation(vector<ptrIndividual> &Pool_Q)
{
    for(int i = 0; i < (int)Pool_Q.size(); i++)
    {
        Pool_Q[i]->EvalIndividual();
      	Pool_Q[i]->setVariableDistance(0);
      	Pool_Q[i]->setObjDistance(0);
    }
}
double NSGAII::getRandom(double Min, double Max)
{
    double f = (double)rand() / RAND_MAX;
    return Min + f * (Max - Min);
}
bool NSGAII::CriterionStop()
{
    if(this->CurrentGeneration < this->NumberGenerations)
    {
        this->CurrentGeneration++;
        return true;
    }
    else
    {
        return false;
    }
}
void NSGAII::FastComputingKfront(vector<ptrIndividual> &Survivors, vector<ptrIndividual>  &CurrentSet, vector<int> &KFrontCurrent, int &Rank)
{
	
//decrease dominance count...
	if(CurrentSet.size() == 1) KFrontCurrent.push_back(0);

	while(KFrontCurrent.empty())
	{
           vector<ptrIndividual> tmpFront, tmpSurvivors;
	   for(int i = 0; i < Survivors.size(); i++)
		{
		   if(Survivors[i]->TimesIsDominated == 0)
			{ 
			   tmpFront.push_back(Survivors[i]); 
			}
		}

	   for(int i = 0; i < CurrentSet.size(); i++)
		{
		if(CurrentSet[i]->TimesIsDominated == 0) 
			{ 
			   KFrontCurrent.push_back(i);
			}
		}
		if(KFrontCurrent.empty())
		{
			for(int i = 0; i < tmpFront.size(); i++) 
	      		{
		  	   for(int j = 0; j < tmpFront[i]->WhoDominate.size(); j++)
			     {
			        tmpFront[i]->WhoDominate[j]->TimesIsDominated--;
			     }
			        tmpFront[i]->setRank( Rank );
			     tmpFront[i]->TimesIsDominated--;
	      		}	   

		 	Rank++;
		}
	}
}
void NSGAII::Fast_Non_Dominated_Sort_Current(vector<ptrIndividual> &Current, vector<ptrIndividual> &Survivors, vector<int> &CurrentFirstFront)
{
	CurrentFirstFront.clear();
    vector<ptrIndividual> Pop;
    //Realizar la unión de los dos conjuntos
    Union(Current, Survivors, Pop);

    vector< vector<int> > F(1), S(Pop.size());
    vector<int> n(Pop.size(),0);

    int mindominatedcount = 10*Pop.size();
    for(int i = 0; i < (int)Pop.size(); i++)
    {
	
        Pop[i]->TimesIsDominated=0;
	Pop[i]->WhoDominate.clear();
        for(int j = 0; j < (int)Pop.size(); j++)
        {
            if(j==i) continue;
            if( Pop[i]->getObjectives().Dominate(Pop[j]->getObjectives()) )
            {
		   S[i].push_back(j);
                   Pop[i]->WhoDominate.push_back(Pop[j]);
            }
            else if( Pop[j]->getObjectives().Dominate(Pop[i]->getObjectives()))
            {
                Pop[i]->TimesIsDominated++;
		n[i]++;
            }
        }
        if( n[i] == 0)
        {
            F[0].push_back(i);
        }
    }
	
    int k=0;

    while(F[k].size())
    {
        vector<int> Q;
        for(int i = 0; i  < (int)F[k].size(); i++)
        {
            int indexSi = F[k][i];
            for(int j = 0; j < (int)S[indexSi].size(); j++)
            {
                n[S[indexSi][j]]--;
                if(n[S[indexSi][j]]== 0)
                {
                    Q.push_back(S[indexSi][j]);
                }
            }
        }
      k++;
        F.push_back(Q);
    }
    bool FlagFirstFront = true;
     for(int i = 0; i < (int)F.size(); i++)
    {
        for(int j = 0; j < (int)F[i].size(); j++)
        {
            int indexIndividual = F[i][j];
            if(indexIndividual < (int)Current.size())
            {
             // Pop[indexIndividual]->setRank(i);
              //Current[indexIndividual] = Pop[indexIndividual];
		//cout << Pop.size()<<" "<<Current.size() <<" " <<indexIndividual<<endl;

              if( FlagFirstFront )
                  CurrentFirstFront.push_back(indexIndividual);
            }else
            {
              //Survivors[indexIndividual-Current.size()]->setRank(i);
            }
              Pop[indexIndividual]->setRank(i);
        }

        if(CurrentFirstFront.size())
          FlagFirstFront =false;
    }
  //  return CurrentFirstFront;
}
void NSGAII::PlotInterfaceRSpaceObjective()
{
    string Comand, X = "", Y = "";

    for(int i = 0; i < (int)this->Pool_Q.size(); i++)
    {
            std::stringstream temporal;
            temporal << Pool_Q[i]->getObjectives().SpaceObjectives[0].Fitness << " , ";
            X += temporal.str();
            temporal.str("");
            temporal << Pool_Q[i]->getObjectives().SpaceObjectives[1].Fitness << " , ";
            Y += temporal.str();
            temporal.str("");
    }
    X = X.substr(0, X.size()-2);
    Y = Y.substr(0, Y.size()-2);
    Comand = "echo \"  pdf(file = paste('Objective_"+this->ProblemName+"_"+to_string(this->CurrentGeneration)+"','.pdf' , sep = '' )) ;  plot(x = c("+ X +"), y = c("+Y+"), main=('NSGAII Generation "+to_string(this->CurrentGeneration)+"') ,xlab='f1', ylab='f2' ) ;    \" | R --Silent --no-save 2>/dev/null | tail -n 0";

    //Comand = "echo \"plot(x = c(3.47788,4), y = c(4,5) ) \" | R --Silent --no-save 2>/dev/null | tail -n 0";
   // cout << Comand<<endl;
    system(Comand.c_str());
}
void NSGAII::PlotInterfaceRSpaceVariables()
{
    string Comand, X = "", Y = "";

    for(int i = 0; i < (int)this->Pool_Q.size(); i++)
    {
            std::stringstream temporal;
            temporal << Pool_Q[i]->getFenotype(0) << " , ";
            X += temporal.str();
            temporal.str("");
            temporal << Pool_Q[i]->getFenotype(1) << " , ";
            Y += temporal.str();
            temporal.str("");
    }
    X = X.substr(0, X.size()-2);
    Y = Y.substr(0, Y.size()-2);
    Comand = "echo \"  pdf(file = paste('Desicion_"+this->ProblemName+"_"+to_string(this->CurrentGeneration)+"','.pdf' , sep = '' )) ;  plot(x = c("+ X +"), y = c("+Y+"), main=('NSGAII Generation "+to_string(this->CurrentGeneration)+"') ,xlab='f1', ylab='f2' ) ;    \" | R --Silent --no-save 2>/dev/null | tail -n 0";

    //Comand = "echo \"plot(x = c(3.47788,4), y = c(4,5) ) \" | R --Silent --no-save 2>/dev/null | tail -n 0";
   // cout << Comand<<endl;
    system(Comand.c_str());
}
void NSGAII::End()
{
 ofstream save;
  save.open (FileNameSummary);
   save << this->Survivors.size() << " " << this->NumberObjectives<<endl;
  for(int i = 0; i < (int)this->Survivors.size(); i++)
  {
    for(int m = 0; m < this->NumberObjectives; m++)
      save << this->Survivors[i]->getObjectiveValue(m)<< " ";
    save << endl;
  }
  //save << this->Survivors;
}
