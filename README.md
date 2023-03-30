# MILP-read-across
The proposed methodology is a regression technique based on mathematical programming. The algorithm solves Mixed Integer Linear Programming (MILP) problems to split the hyperspace into consistent regions and allocate the samples into regions. Also, it identifies linear regression models for each region. Therefore, predictions for untested samples can be generated according to the group they belong and the use of multivariable linear equations, unique for each group. The grouping of the samples is achieved by selecting one feature of the available dataset as the partition feature and arrange the samples into regions depending on their value of  that specific feature (“1D MILP problem”). If the dataset contains two different types of descriptors, the grouping boundaries can be determined in the two-dimensional space (“2D MILP problem”). 

<b>Requirements:</b> 
<ul>
To run  MILP read-across, you will need: 
<ul>
 <li>Matlab</li> 
 <li>YALMIP (open source)</li> 
 <li>MIP solvers: 
<ul>
  <li>Mosek (commercial or academic license)</li>
 <li>Gurobi (commercial or academic license)</li> 
   </ul>
  </li>
  </ul>
</ul>
<ul>
To reproduce results shown in the publication you will need Matlab R2015b, YALMIP R20190425 and Mosek 9.0. 
</ul>

<b>Case study: </b>
<ul>The proposed method is demonstrated on two datasets. The dataset derived from the publication of Walkey et al. (2014) (Gold ENMs), and the dataset derived from the publication of Xia et al., (2011) (MWCNT). The Gold ENMs dataset consists of 84 culture medium incubated gold anionic and cationic ENMs which are characterized by 40 physicochemical and 129 biological descriptors (protein corona fingerprints). The protein corona fingerprints were filtered by GSVA and only 63 were consisted as statistically significant proteins. (Varsou et al. (2017)). There are also available measurements of there cell association with human A549 cells (in mL/μg(Mg)) as a toxicity index. 
</ul>
<ul>The MWCNT dataset, consists of a 40 nm diameter MultiwalledCarbon Nanotube (MWCNT), coated with hydroxyl derivatives and the surface absorption properties. Five nanodescriptors representing the surface adsorption interactions (lipophilicity, hydro-gen bond acidity and basicity, polarity/polarizability, and lone-pair electrons) were used to train a model for predicting the adsorption coefficient (k) of this particular ENM. The dataset contained adsorption coefficient (k) values for 28 probe compounds with various physicochemical properties. The k values were converted to logarithmic scale.
 </ul>
 
 
<b>Default parameters:</b> 
  <ul>
For Gold NPs in “1D MILP problem”: 
<li>lamda=0.01</li>
<li>beta=0.05</li>
<li>epsilon=0.05</li>
<li>U=10</li>
 </ul>
   <ul>
For MWCNT in “1D MILP problem”: 
<li>lamda=0.01</li>
<li>beta=0.05</li>
<li>epsilon=0.05</li>
<li>U=10</li>
 </ul>
   <ul>
For Gold NPs in “2D MILP problem”: 
<li>lamda=0.02</li>
<li>beta=0.05</li>
<li>epsilon=0.05</li>
<li>U=10</li>
</ul>


<b>Solving "1D problem"</b>
<ul>
To solve the "1D problem" you will need OPLRA_1D_trainingtest.m and OPLRA_1D_test.m files. The algorithm examines the addition of regions by recording the values of the objective function of test set between two iterations with different number of regions. OPLRA_1D_test.m calculates the value of the objective function for the test set to check the addition or not of regions. 
</ul>

<b>Solving "2D problem"</b>
<ul>
OPLRA_2D_independent_pf1_final.m, OPLRA_2D_independent_pf1_pf2_final.m and OPLRA _2D_trainingtest_final.m  are used to solve the "2D problem" according to the approach of the selection of partition features. OPLRA_2D_independent_pf1_final.m is used for the sequentially approach, OPLRA_2D_independent_pf1_pf2_final.m selects the partition features of two groups of features independently and OPLRA_2D_trainingtest_final.m is used for the simultaneously approach.

Each one of the above files contains two functions which are necessary for solving the "2D problem". OPLRA_2D_test_final.m applies the produced model to the test set and records the value of the objective function while OPLRA_2D_final.m is a function used for the addiction of regions for the "2D problem". 
</ul>

# Cite us!
If you find this methodology useful, please cite us! :) 

Dimitra-Danai Varsou, Nikoletta-Maria Koutroumpa, and Haralambos Sarimveis, "Automated Grouping of Nanomaterials and Read-Across Prediction of Their Adverse Effects Based on Mathematical Optimization", Journal of Chemical Information and Modeling 2021 61 (6), 2766-2779. <a href="https://doi.org/10.1021/acs.jcim.1c00199">DOI: 10.1021/acs.jcim.1c00199</a>

# GUI
The trained models using MILP-read-across and the above data sets are made publicly available through a user-friendly web service (<a href="https://vythos.jaqpot.org/">https://vythos.jaqpot.org/</a>).

# License
This application is released under <a href="https://www.gnu.org/licenses/gpl.html"> GNU General Public License v.3</a>.
```html

This program is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
You should have received a copy of the GNU General Public License along with this program.  
If not, see here: http://www.gnu.org/licenses/.
