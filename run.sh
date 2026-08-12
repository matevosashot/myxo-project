#!/usr/bin/env bash

# Every named list is updated key-by-key; unmentioned keys keep their default.
wolframscript -file fluctuations.wls \
  "lyapunovSolverParams={Nint->19, fbox->10}; \
   otherParams={outputDir->\"./output\", saveData->True, plotBox->9}"
