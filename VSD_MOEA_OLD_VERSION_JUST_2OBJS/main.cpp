/**
	Autores:
		 Carlos Segura González
		 Joel Chacón Castillo
	Fecha: 18/10/2016
	Descripción:
		En este fichero se almacena la función central de todo
		el programa.
**/
#include <iostream>
#include <omp.h>
#include <cmath>
#include <string>
#include "Measures.hpp"
#include "NSGAII.hpp"
using namespace std;

string DiccionarioProblemas(int Id);
void PrepareDirectory(string Dest, string Instance);
int strInstancetoId(string Instance);
void SetConfiguration(int argc, char*argv[], double &ProbCross, double &ProbMut, string &Path, int &SizePool, int &TotalGenerations, int &Sed, int &PeriodReport, string &Label, string &Instance, string &Crossover, string &Mutation, double &ParameterDInit);
int Str_to_Id_Crossover(string Crossover);
int Str_to_Id_Mutation(string Mutation);
void PrintHelp();
int main(int argc, char *argv[])
{
	if(argc<2)
         {
	    
	    cout << "Unknown Argument.."<<endl;
	    PrintHelp();
	    exit(0);
	 }

	int Sed = 1;
        double ProbCross= 0.9, ProbMut = 1.0/24.0, ParameterDInit= 0.25*sqrt(24);

        int SizePool = 100, TotalGenerations =250, PeriodReport=10;
	string Path="./", Instance="WFG1", Label ="0000", Crossover="SBX", Mutation="Polynomial";


	SetConfiguration(argc, argv, ProbCross, ProbMut, Path, SizePool, TotalGenerations, Sed, PeriodReport, Label, Instance, Crossover, Mutation, ParameterDInit);

 	PrepareDirectory(Path, Instance);
	int IdCrossover = Str_to_Id_Crossover(Crossover);	
	int IdMutation = Str_to_Id_Mutation(Mutation);
	int i = strInstancetoId(Instance);
		clock_t TimeInit = clock();
		string ProblemName = DiccionarioProblemas(i);
		cout << " Iniciando...."<<ProblemName<<endl;
		Benchmark ObjBenchmark(i);
		ParameterDInit = ParameterDInit*sqrt(ObjBenchmark.getDimension()); 	
		ProbMut = 1.0/ObjBenchmark.getDimension();
		NSGAII ObjNSGAII(SizePool, ProbCross, ProbMut,  ObjBenchmark, TotalGenerations, ParameterDInit, ProblemName, Sed, Path, Label, IdCrossover, IdMutation, PeriodReport);
		ObjNSGAII.Init_NSGAII();
		cout << "Problema "+ProblemName+" terminado.. TiempoTotal: " << (clock()-TimeInit)/ (double)  CLOCKS_PER_SEC << endl  ;
    return 0;
}
void PrintHelp()
{
	cout << "Instructions:"<<endl;
	cout << "--Instance NAMEINSTANCE (WFG1)"<<endl;
	cout << "--Sed SED (299)" <<endl;
	cout << "--Label OPERATOR_SBX, is a postfix by filename"<<endl;
	cout << "--Pcross 0.9, is the Crossover probability" <<endl;
	cout << "--Pmut 0.3, is the Mutation Probability " << endl;
	cout << "--Path ./RESULT, is the directory where will save results"<<endl;
	cout << "--SizePool 100, is the number of individual by generation"<<endl;
	cout << "--Generations 25000, is the number of generations"<<endl;
	cout << "--Crossover SBX, is the Crossover operator to use"<<endl;
	cout << "--Mutation Polynomial ,is the mutation operator to use"<<endl;
	cout << "--PeriodReport 10, each N generations will save informaton in the \"--Path\" files"<<endl;
	cout << "--DI 0.75 , initial valor of diversity D"<<endl;
	cout <<"---------Crossover Operators Availables-------------"<<endl;
	cout << "SBX, BLX, HBX, Fuzzy, SBXHybrid"<<endl;
	cout <<"---------Mutation Operators Availables--------------"<<endl;
	cout << "Polynomial, NormalDistribution, Random"<<endl;
}
void PrepareDirectory(string Dest, string Instance)
{
	string term = "mkdir -p ";
	term += Dest;
	term += "/"+Instance;
	system(term.c_str());
}
int Str_to_Id_Crossover(string Crossover)
{
	if(Crossover == "SBX")
		return _SBX;
	if(Crossover == "BLX")
		return _BLX;
	if(Crossover == "HBX")
		return _HBX;
	if(Crossover == "Fuzzy")
		return _FUZZY;
	if(Crossover == "SBXHybrid")
		return _SBXH;
	if(Crossover == "SBXHybridDynamic")
		return _SBXHD;
	if(Crossover == "DE")
		return _DE;
	else
	{
		cout << "Unknown Crossover..."<<endl;
		exit(0);
	}
}
int Str_to_Id_Mutation(string Mutation)
{
	if(Mutation == "Polynomial")
		return _POLYNOMIAL;
	else  if(Mutation == "NormalDistribution")
		return _NORMAL;
	else if(Mutation == "Random")
		return _RANDOM;
	else
	{
		cout << "Unknown Mutation...."<<endl;
		exit(0);
	}
}
string DiccionarioProblemas(int Id)
{
	switch(Id)
	{
		case Type_WFG1:
			return "WFG1";
		break;
		case Type_WFG2:
			return "WFG2";
		break;
		case Type_WFG3:
			return "WFG3";
		break;
		case Type_WFG4:
			return "WFG4";
		break;
		case Type_WFG5:
			return "WFG5";
		break;
		case Type_WFG6:
			return "WFG6";
		break;
		case Type_WFG7:
			return "WFG7";
		break;
		case Type_WFG8:
			return "WFG8";
		break;
		case Type_WFG9:
			return "WFG9";
		break;
		case ZDT1:
			return "ZDT1";
		break;
		case ZDT2:
			return "ZDT2";
		break;
		case ZDT3:
			return "ZDT3";
		break;
		case ZDT4:
			return "ZDT4";
		break;
		case ZDT5:
			return "ZDT5";
		break;
		case ZDT6:
			return "ZDT6";
		break;
		case DTLZ1:
			return "DTLZ1";
		break;
		case DTLZ2:
			return "DTLZ2";
		break;
		case DTLZ3:
			return "DTLZ3";
		break;
		case DTLZ4:
			return "DTLZ4";
		break;
		case DTLZ5:
			return "DTLZ5";
		break;
		case DTLZ6:
			return "DTLZ6";
		break;
		case DTLZ7:
			return "DTLZ7";
		break;
		case UF1:
			return "UF1";
		break;
		case UF2:
			return "UF2";
		break;
		case UF3:
			return "UF3";
		break;
		case UF4:
			return "UF4";
		break;
		case UF5:
			return "UF5";
		break;
		case UF6:
			return "UF6";
		break;
		case UF7:
			return "UF7";
		break;
	}
}
int strInstancetoId(string Instance)
{
	if(Instance == "WFG1" )
		return Type_WFG1;
	else if(Instance == "WFG2")
		return Type_WFG2;	
	else if(Instance == "WFG3")
		return Type_WFG3;	
	else if(Instance == "WFG4")
		return Type_WFG4;	
	else if(Instance == "WFG5")
		return Type_WFG5;	
	else if(Instance == "WFG6")
		return Type_WFG6;
	else if(Instance == "WFG7")
		return Type_WFG7;		
	else if(Instance == "WFG8")
		return Type_WFG8;	
	else if(Instance == "WFG9")
		return Type_WFG9;	
	else if(Instance == "ZDT1")
	        return ZDT1;
	else if(Instance == "ZDT2")
	        return ZDT2;
	else if(Instance == "ZDT2")
	        return ZDT2;
	else if(Instance == "ZDT3")
	        return ZDT3;
	else if(Instance == "ZDT4")
	        return ZDT4;
	else if(Instance == "ZDT5")
	        return ZDT5;
	else if(Instance == "ZDT6")
	        return ZDT6;
	else if(Instance == "DTLZ1")
	        return DTLZ1;
	else if(Instance == "DTLZ2")
	        return DTLZ2;
	else if(Instance == "DTLZ3")
	        return DTLZ3;
	else if(Instance == "DTLZ4")
	        return DTLZ4;
	else if(Instance == "DTLZ5")
	        return DTLZ5;
	else if(Instance == "DTLZ6")
	        return DTLZ6;
	else if(Instance == "DTLZ7")
	        return DTLZ7;
	else if(Instance == "UF1")
	        return UF1;
	else if(Instance == "UF2")
	        return UF2;
	else if(Instance == "UF3")
	        return UF3;
	else if(Instance == "UF4")
	        return UF4;
	else if(Instance == "UF5")
	        return UF5;
	else if(Instance == "UF6")
	        return UF6;
	else if(Instance == "UF7")
	        return UF7;
	else if(Instance == "UF10")
	        return UF10;
	else
	{
		cout << "Unknown instance... " << endl;
		exit(0);
	}
}
void SetConfiguration(int argc, char*argv[], double &ProbCross, double &ProbMut, string &Path, int &SizePool, int &TotalGenerations, int &Sed, int &PeriodReport, string &Label, string &Instance, string &Crossover, string &Mutation, double &ParameterDInit)
{
	for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--Instance")
			Instance = string(argv[++i]);
		else if(Terminal == "--Sed")
			Sed = atoi(argv[++i]);
		else if(Terminal == "--Label")
			Label = string(argv[++i]);
		else if(Terminal == "--Pcross")
			ProbCross = atof(argv[++i]);
		else if(Terminal == "--PMut")
			ProbMut = atof(argv[++i]);
		else if(Terminal == "--Path")
			Path = string(argv[++i]);
		else if(Terminal =="--SizePool")
			SizePool = atoi(argv[++i]);
		else if(Terminal == "--Generations")
			TotalGenerations = atoi(argv[++i]);
		else if(Terminal == "--Crossover")	
			Crossover = string(argv[++i]);
		else if(Terminal == "--Mutation")
			Mutation = string(argv[++i]);
		else if(Terminal == "--DI")
			ParameterDInit = atof(argv[++i]);
		else if(Terminal =="--PeriodReport")
			PeriodReport = atoi(argv[++i]);
		else if(Terminal == "--help" || Terminal == "--h")
			PrintHelp();
		else
		{
			cout << Terminal<<endl;
			cout << "Unknown Argument...";
			exit(0);
		}
	    }
}
