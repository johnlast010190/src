#include "cfdTools/general/include/readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", mesh.changing());

checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);

bool fastMode = pimple.dict().lookupOrDefault("fastMode", false);
