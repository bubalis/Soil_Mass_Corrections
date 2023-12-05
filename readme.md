<h1> Methods for Equivalent Soil-Mass Corrections and Assessments </h1>

This repository relates to Equivalent Soil Mass (ESM) corrections in Soil Carbon Accounting.

ESM methods are used to correct for erroneous calculations of increased carbon stocks caused by changes in soil bulk density. Since at least 1995, various methods have been proposed. See the references section for more explanation of the why and how.

<h2> Directory Structure</h2>
ESM_and_mass_correction_functions.R : A set of functions for generating equivalent soil-mass estimates of SOC content. This script includes modifed code taken from two previous peer-reviewed publications[1,2].

estimate_ESM_errors.R : An analysis of internal validity of linear and spline-based interpolation using Leave-One-or-More Out validation. This analysis is also used to generate estimates of the error of spline-based interpolation, used later when using ESM as a benchmark.

data_sim.R : Functions for simulating fine-resolution soil horizons, modified from [2].

analysis_functions.R : Functions used to simulate the use of different soil sampling strategies and estimation methods.

analyze.R : analysis and simulation of different methods on all of the available datasets.

figures.R : generate figures summarizing the simulations.

\Source_data : directory containing the source data used for the analysis.


<h1> References</h1>
<a id="1">[1]</a> 
von Haden, Adam C and Yang, Wendy H and DeLucia, Evan H (2020). 
Soils' dirty little secret: Depth-based comparisons can be inadequate for quantifying changes in soil organic carbon and other mineral soil properties.
Global Change Biology. 26:7
<a id="1">[2]</a> 
Fowler, Ames F and Basso, Bruno and Millar, Neville and Brinton, William F (2023).
A simple soil mass correction for a more accurate determination of soil carbon stock changes
Scientific Reports. 13:1











