## Approaches to treatment effect heterogeneity in the presence of confounding
### Sarah C. Anoke, Sharon-Lise Normand, Corwin M. Zigler

#### Code + Results for simulation study from manuscript

This simulation is the same as the one described in section 4.1 of the manuscript (only scenarios A through D, see [`mainSimStudy-nonLinearY`](https://github.com/sanoke/approachesTEH/tree/master/mainSimStudy-nonLinearY) for D* (the nonlinear simulation scenario).

To run this simulation:

1. **Generate the data** by running `dataGeneration.R`.  
    - The script `helperFcns.R` contains the function `datagen()`, which is used to generate a dataset of size `n` (as described in section 4.1.1 of the manuscript). 
    - The script `dataGeneration.R` uses `datagen()` within a loop to generate `numsim` datasets of size `n`, under all four simulation scenarios (A, B, C, D). Example datasets are stored within the directories `dataFilesA`, `dataFilesB`, `dataFilesC`, `dataFilesD`. 
    - Note that in the manuscript `numsim = 1000`... one GBM or BART simulation completes in a few minutes, but one FS simulation takes up to 12 hours. Thus you may want to select a more practical number of simiulations.
    - The outcome variance is defined using the variable `varOutcome` within the script `dataGeneration.R`. 

2. **Apply each of the three TEH identification methods to the data.**
    - The script `simulationStudy.R` applies one of the three estimation procedures to one dataset. Thus the application of all three estimation procedures to all 4 x `numsim` datasets would require calling this script 3 x 400 = 1200 times, resulting in 1200 result files.
       - The results are saved in the `res` directory.
    - The function that does the TEH estimation work `estimation()` (defined in `helperFcns.R`). 
       - This function has three arguments, `ds` (the dataset to be partitioned), `procedure` (an integer indicating which estimation procedure to use in determining the partition, where 2 is GBM, 3 is BART, and 8 is the facilitating score), and `numGrp` (the number of desired subgroups in the estimated partition).
       - This function returns a `list` with the following elements:
           1. `TE` is a vector of length `numGrp`, with the estimated treatment effect in each subgroup (arranged in ascending order). 
           2. `var` is a vector of length `numGrp`, with the estimated variance of each treatment effect in `TE`. Looking at a particular subgroup, and letting `Y1data` represent the vector of outcomes among observations in the treatment group, and `Y0data` represent the control group, this subgroup-specific variance is calculated as
              ```r
              var(Y1data)/length(Y1data) + var(Y0data)/length(Y0data)
              ```
           3. `n0` is a vector of length `numGrp`, containing the number of control observations in each subgroup.
           4. `n1` is a vector of length `numGrp`, containing the number of treatment observations in each subgroup.
           5. `PSgrp` is a vector of length `nrow(ds)`, with an integer between 1 and `numGrp` indicating which subgroup each observation belongs to. The subgroups are ordered, where subgroup 1 has the smallest treatment effect, 2 has the second smallest treatment effect, etc. 
           6. `D` is the dissimilarity matrix (only used in the case of the facilitating score).
           7. `mmt` is a vector of length `nrow(ds)`, with each individual's estimated statistic that is used in the partition procedure. Note that the facilitating score does not produce this.

    - Functions used to apply the Facilitating Score to the data were written by [Xiaogang Su et al](http://www.jmlr.org/papers/v13/su12a.html), and were used + posted here with permission. We ask that you also [seek permission](http://www.math.utep.edu/people/dmprofile.php?username=xsu) before distributing their code.

3. **Visualize the results.**
    - The results contained in the `res` directory are then visualized using the `resultProcessing.R` script. Examples of the expected figures are `mainSimStudy-annot.png` and `mainSimStudy.png`.
    - Sensitivity calculations are done using the `sensitivity.R` script.
    - Forest plots are generated using `resultProcessing-forest.R`. Note: The white dot represents the mean.
    - Example result files have been provided in the `res` directory. 

4. **Other notes**
    - `shellScript.txt` and `shellScriptFS.txt` are examples of scripts we use to run the simulation on our clusters; they may be helpful for you.

#### Bug Reporting

Report any bugs or suggestions as an [issue](https://github.com/sanoke/approachesTEH/issues).

#### Licensing

Our code is open source licensed under the [MIT License](https://github.com/sanoke/approachesTEH/blob/master/LICENSE).
