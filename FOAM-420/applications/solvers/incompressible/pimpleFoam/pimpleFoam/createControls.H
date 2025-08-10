#include "cfdTools/general/solutionControl/createControl.H"
#include "cfdTools/general/include/createTimeControls.H"

bool ddtCorr
(
    pimple.dict().lookupOrDefault("ddtPhiCorr", true)
);
