#include "cfdTools/general/include/readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);

Switch solveThermal
(
    pimple.dict().lookupOrDefault<Switch>("solveThermal", true)
);
