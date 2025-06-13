Autopicking takes in a job with imported movies, and outputs a locally-refined angle-unified volume and particles for each class that contains > 5000 particles in classification.

It requires a cs_config.yml file in your $HOME directory. An example cs_config.yml is provided. 
It requires linking all the csparc_*.py files into the csparc_automation directory, as following
    cd csparc_automation
    ln -s ../*.py .
It requires cryosparc-tools, requirements updated. 

The least error-prone way to run it is to start with imported movies. In that case, it can be run from anywhere, and will change directory to your project directory, and create its outputs there. If the processing is interrupted, it can be restarted by running the script from the project directory, where autopick.log is present. It will redo the last step. In the project directory, a cs_autopick.yml file can be optionally provided to override defaults settings for model files for classification, segment length and tube diameter. 

Models for classification should correspond to the following 6 classes of microtubules: 11-3, 12-3, 13-3, 14-3, 15-4, and 16-4. They should be named accordinly so that sorting them alphabetically gives that same order. 

If pipeline continuation is required, one can fake an autopick.log file in the project directory, and it will attempt to resume from last step that cryosparc logs as 'finished'. An example autopick.log file is provided. If the logfile is found in the project directory, it will accept the variables from the header of the file - segment_length, tube_diameter, particle padding factor, and apix. The logfile has to be written exactly as provided, including leading spaces. The first number, sequence of job submitted, doesn't matter. The P, W, and J entries matter, as well as the name of the step. 
