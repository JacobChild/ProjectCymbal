FullUoverTimeOutput
This came from taking UoverTime and writing it to a delimited file. Each row contains all of the data for a time step. Time goes from 0 to 2 with 100 steps. The output data for a time step was arranged in a matrix correlating to every r and theta value. The writedlm command took that matrix and wrote it on one line, so I am unsure how it is organized now.

UoverTime45
This came from outputting all of the data at the 45th (.88888sec) time step. It is in the proper format, and correlates to the outputs shown in some of the plots.