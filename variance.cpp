#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
int nvar, nss, npop;
using namespace std;
vector<string> getfiles(char *pattern)
{
  //look for a list of files given a pattern...
  char tt[500];
  strcpy(tt, "/bin/ls ");
  strcat(tt, pattern);
  FILE *proc = popen(tt,"r");
  char buf[1024];
   vector<string> files;
  while ( !feof(proc) && fgets(buf,sizeof(buf),proc) )
  {
    string line(buf);
	line.pop_back();
    files.push_back(line);
  }     
   return files;
}
double hypot(vector<double> &A, vector<double> &B)
{
   double dist = 0.0;
   for(int i = 0; i < A.size(); i++)// -- all
   //for(int i = 0; i < 4; i++) //--position
   //for(int i = 4; i < A.size(); i++) //--distance
   {
	dist += (A[i]-B[i])*(A[i]-B[i]);
   }
   return sqrt(dist);
}
int main(int argc, char **argv)
{
  vector <string> ff = getfiles(argv[1]);
  nss= atoi(argv[2]);
  npop = atoi(argv[3]);
  nvar = atoi(argv[4]);
  vector< vector< vector< double > > > dcn_run(ff.size(), vector< vector<double> > (nss, vector<double>( nvar, 0.0)));
  for(int i = 0; i < ff.size(); i++)
  {
	
	   FILE *summary = fopen(ff[i].c_str(), "r");
        for(int k = 0; k < nss; k++)
        {
	   vector< vector< double> > information( npop, vector< double> ( nvar, 0.0)  );
	if(summary==NULL)
	{
//		cout <<"The file "<< ff[i].c_str()<<" cann't be opened.. "<<endl;
		exit(0);
	}
	  for(int l = 0; l < npop; l++)
	  {
	    for(int o = 0; o < nvar; o++)
	    {
		fscanf(summary, "%lf ",&information[l][o]);	
	        information[l][o] /= (2.0*(o+1.0));
	    //    information[l][o] =  (information[l][o]-100.0)/200.0;// (2.0*(o+1.0));
	    }
	  }

	  vector<double> expectation(nvar, 0.0);
	  double mean = 0.0;
	  for(int j = 0; j < npop; j++) //Computing the mean..
	  {
	    for(int d = 0 ; d < nvar; d++)
	    {
	     expectation[d] += (information[j][d])/(npop);
	    }
	  }
	  for(int j = 0; j < npop; j++) ///computing the variance...
	  {
	    for(int d = 0 ; d < nvar; d++)
	    {
	     dcn_run[i][k][d] += ((information[j][d] - expectation[d])*(information[j][d] - expectation[d]))/(npop);
	    }
	  }
	}
	  fclose(summary);
  }
  vector<vector<double> > mean( nss, vector<double> (nvar,0.0));
  for(int i = 0; i < dcn_run.size(); i++)
  {
     for(int j = 0 ;j < dcn_run[i].size(); j++)
	{
	   for(int d =0; d < nvar; d++)
		 mean[j][d] +=dcn_run[i][j][d];	
	}
  }

 for(int i = 0; i < nss; i++)
  {
	printf("%d ",i+1);
	for(int d = 0; d < nvar; d++)	
	
	printf("%lf ", mean[i][d]/ff.size());
	printf("\n");
  }
}
