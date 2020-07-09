This folder was created on 2020-03-10
The contexts are based on an older version of the code which we have used to start a fresh repo

Notes on how to use files in this repo:
---------------------------------------
The scripts contained in the \matlabFiles directory are examples on how to postprocess the data. They are not to be used to analyze data during a production run. 
To create custom made postprocessing scripts, create a new directory and add this directory to the .gitignore file so that when pushing to the remote repo, only the essential machinery of the code is tracked

Questions:
- How does it scale on megabucky?
- How to save output files directly to the "outputFile" directory?
- How to set the fueling rate independent of the leak rate?
- How to we calculate the wall time of the calculation correctly?
