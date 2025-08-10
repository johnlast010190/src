#include "cfdTools/general/include/readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", false);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);

Switch detailTime
(
    pimple.dict().lookupOrDefault<Switch>("detailTime", false)
);

Switch solveThermal
(
    pimple.dict().lookupOrDefault<Switch>("solveThermal", true)
);

Switch ddtPhiCorr =
    pimple.dict().lookupOrDefault<Switch>("ddtPhiCorr", true);
