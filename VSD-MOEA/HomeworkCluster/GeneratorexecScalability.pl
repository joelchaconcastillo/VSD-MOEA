#Con este script se genera la lista de ejecuciones....
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $file = "ExecutionFileDiversity";
my $fout;
open($fout, '>' ,$file);
my $PathAlgorithm =  `cd ..; pwd`;#"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt";
chomp $PathAlgorithm;
#my $Path = "/home/joel.chacon/Chacon/Tesis/MOEA-D_Diversity/Code_CEC09/moead_based_diversity_SBX/";
my $Instance=0;
my $Sed=0;

####Realizar la búsqueda del parámetro D inicial que proporcione mejores resultados
my @Conf =(
"UF1 2",
"UF4 2",
"UF5 2",
"UF10 3",
"UF9 3",
"UF8 3",
"DTLZ4 2",
"DTLZ4 3",
);




###foreach(@DI)
##{
##	foreach(@Conf)
##	{
##		for($Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
##		{
##			print $fout "~$PathAlgorithm/Ejecutable $_ $Sed \n";
##		}
##	}
##	
##}

my @Variables= ("50","100","250", "500", "1000");

foreach my $v (@Variables)
{
        foreach my $configuration (@Conf)
        {
                for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
                {
                      my @configuration2 = split ' ', $configuration;#~ s/ /_/g; 

		     print $fout "~$PathAlgorithm/Ejecutable $configuration2[0] $v $configuration2[1] $Sed $PathAlgorithm \n";
                }
        }

}

