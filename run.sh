#!/usr/bin/env bash

module load mathematica/14.00


unset DISPLAY

# Every named list is updated key-by-key; unmentioned keys keep their default.
wolframscript -file fluctuations.wls \
  "lyapunovSolverParams={Nint->20, fbox->20}; \
   otherParams={outputDir->\"./output\", saveData->True, plotBox->9}"
