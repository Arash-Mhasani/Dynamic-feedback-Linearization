# Dynamic-feedback-Linearization
Date: Dec. 2021

Author:
Arash Mohammadhasani <Arash.Mhasani@ualberta.ca>
Applied Nonlinear Control Lab 
Dept ECE
University of Alberta
http://ancl.ece.ualberta.ca/

This Maple code implements the Dynamic Extension Algorithm which contructs a dynamics state feedback which extends the dynamics 
so that its relative degree is well-defined.

This code has been tested using Maple 2021 and Matlab 2020b and extends the work in [1].

1) To implement the dynamic extension algorithm using

pre-defined models (4 sample models are provided): Run "DynamicFeedbackLinearization.mw" with other ".mpl" in the same directory.

new models: insert your model at the end of "Models.mw", and then "File|Export As" "Models.mw" into a Maple input file "Models.mpl". Now run "DynamicFeedbackLinearization.mw" with your selected model.

2) For simulating (in Matlab) the controller designed by dynamic feedback linearization algorithn on the SLS system, Run "DFL_SLS_Simulation.mw". You can later 
change the reference trajectory, model parameters, and controller gains in the "PlotResults.m" script and re-run "PlotResults.m" in Matlab.

Description of Files:
mw files:
DynamicFeedbackLinearization.mw - runs DEA for a list of pre-defined models
Models.mw - Maple worksheet version of list of models used in DynamicFeedbackLinearization.mw 
DFL_SLS_Simulation.mw - generates Matlab Simulation of the Control designed based on DEA algorithm for the Slung Load System.

mpl files:
Moldes.mpl - .mpl version of "Models.mw"
dynExt.mpl - .mpl version of "dynExt.mw". This function is used inside "DynamicFeedbackLinearization.mw"
RelDeg.mpl- .mpl version of the "RelDeg.mw". This function is used inside "dynExt.mpl"
lDiff.mpl- .mpl version of the "lDiff.mw". This function is used inside "dynExt.mpl" and "RelDeg.mpl"
printModel.mpl- .mpl version of the "printModel.mw". This function is used for printing extended model at the last iteration of th DEA alg.

worksheets - folder which contains Maple worksheet (mw) versions of files for future development (these files are not required to run DEA)



References
[1] G. Fischer, “Nonlincon: symbolic analysis and design package for nonlinear control systems,” Master’s thesis, Eindhoven University of Technology, 1994.
