# TFG_MoA_MLalgorithms
--> Here, I include all the codes I used for building the ML algorithms. I can not include the original datasets because anyone who wants to have access to it, must agree on the competition terms and rules. You can find it on https://www.kaggle.com/c/lish-moa/data.

Requirements: Python and packages that are scpecified in the notebooks. For running the SBM, we recommend to use pypy (otherwise, it is really slow).

To run the code, we suggest the following:
1. Go to the link above and download the data, after singing in. Note that you should agree the terms of the competition. The folder (lish-moa) should be on the same directory as the other folder we provide.
2. Go to the folder 'Data' and run the notebook, several files will be created in the same folder.
3. Choose whether using 'Off-The-Shelf' algorithms or the 'StochasticBlockModel':
   a) 'Off-The-Shelf' algorithms:
   Go to the folder 'Off-The-Shelf'. The notebook 0 will be used for checking the log-loss of the models created in the other notebooks. Note that the last       notebooks contain also some parts of analysis (mainly related to cell viability and gene expression discretisations).
   
   b) 'StochasticBlockModel':
   Go to the folder 'StochasticBlockModel' and, if the objective is to create a new SBM, you can use the MMSBM_CMG.py (which does not consider cell viability data) or MMSBM_CMGC.py (which does consider cell viability data). We have used this code with pypy to run it faster (it can take several hours). We have optimized this code, so it can be difficult to compare it to the equations of our model. That is why we include two notebooks (MMSBM_AllData.ipynb and MMSBM_CMG.ipynb) that contains all the functions to create the SBM in a more comprehensive way, although it would be highly time consuming. 
   On the other hand, in order to check the SBM's performance, we have created two files of functions that will help us to read the parameters and compute the log-loss (function_read_param.py and results_functions.py) which should be included in the environment you use to run the code. Finally, there are two examples (for both cases of cell viability presence) that will use the aformentioned functions to finally get the log-loss and other valuable information and graphs. 
   
   
