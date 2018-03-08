Program: BSMPT version 1.0

Released by: Philipp Basler and Margarete Mühlleitner

Manual: version 1.0 

BSMPT - Beyond the Standard Model Phase Transitions:
The C++ program package BSMPT calculates for models with extended
Higgs sectors the loop-corrected effective potential at finite temperature
including the daisy resummation of the bosonic masses. The program
computes the vacuum expectation value (VEV) \f$ v \f$ of the potential
as a function of the temperature, and in particular the critical VEV
\f$v_c\f$ at the temperature \f$T_c\f$ where the phase transition takes
place. In addition, the loop-corrected trilinear Higgs self-couplings are
provided. We apply an 'on-shell' renormalization scheme in the sense
that the loop-corrected masses and mixing angles are required to be
equal to their tree-level input values. This allows for efficient
scans in the parameter space of the models. The models implemented so far
are the CP-conserving and the CP-violating 2-Higgs-Doublet Models (2HDM) and the
Next-to-Minimal 2HDM (N2HDM). The program structure is such that the
user can easily implement further models.

The program package can be downloaded at:
https://github.com/phbasler/BSMPT

The documentation is given at: https://phbasler.github.io/BSMPT/

Modifications and corrected bugs are reported in the file 'Changelog.md'.

Sample input and output files are provided in the directory 'example'.

For additional information, comments, complaints or suggestions please e-mail
to:  philipp.basler@kit.edu

---
   

##Installation:


1) The program requires the GSL library which the code assumes to be installed in PATH. 

2) Go to the directory sh and type 'chmod +x autogen.sh' and subsequently
   'chmod +x InstallLibraries.sh'.

3) Then type `./InstallLibraries.sh --lib=PathToLib --CXX=YourC++Compiler --CC=YourCCompiler` in order to install the
   eigen and CMAES libraries in the path 'PathToLib'. This might take a while.

4) Afterwards type `./autogen.sh --lib=PathToLib --CXX=YourC++Compiler` in order to create a makefile in the
   main path, with the libraries installed in 'PathToLib'.  

5) Finally go back to the main directory and type 'make'. 
 
---
  
##How to add a new model (for further details, also see the manual):

1) Go to the file IncludeAllModels.h and rename the variable
   'C_ModelTemplate' to the variable name with which the new model shall
   be selected by the program.

2) Go to IncludeAllModels.cpp and add

``` c++    
	  else if(choice == C_ModelTemplate)
     {
       return std::unique_ptr<Class_Potential_Origin> { new Class_Template };
     }
```


   to the function Fchoose.

3) Adjust the functions in ClassTemplate.cpp as needed for the new model.

4) Have fun!

