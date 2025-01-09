README

######
Update: 20250109

######
Contact: Elise S. Droste (e.droste@uea.ac.uk)

######
This workflow uses Python packages and code. Use the environment.yml file to create a virtual conda environment containing all the required packages. 

To build the exact identical conda environment, create an environment usingn the spec-file.txt. 

See https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html on how to create an environment from environment.ylm and/or spec-file files. 


######
The scripts in this workflow assume that the raw VINDTA datafile abides by the following filenaming conventions. 

Filenaming conventions

For CRMs: 

Filename: 8888_batch_CRMbottlenumber_repeat_replicate
Station = 8888
Cast = CRM batch number
Niskin = CRMbottlenumber
Repeat = repeat
Replicate = replicate


For junk samples: 

Station = 9999
Cast = analysis date