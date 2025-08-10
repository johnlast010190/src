#include "cfdTools/general/solutionControl/createControl.H"
#include "cfdTools/general/include/createTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool ddtCorr
(
    pimple.dict().lookupOrDefault("ddtPhiCorr", true)
);
