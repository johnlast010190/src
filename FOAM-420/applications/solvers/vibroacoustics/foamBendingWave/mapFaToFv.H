label pID = mesh.boundary().findPatchID(faPatch);

loadV.boundaryFieldRef()[pID] = load.primitiveField();
wV.boundaryFieldRef()[pID]    = w.primitiveField();

magSf.boundaryFieldRef()[pID]  = mesh.boundary()[pID].magSf();
crit.boundaryFieldRef()[pID]   = plateD.value() * sqr(runTime.deltaTValue()*mesh.boundary()[pID].magSf())();

//wLapV.boundaryFieldRef()[pID] = wLap.primitiveField();

