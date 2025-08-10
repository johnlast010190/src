label pID = mesh.boundary().findPatchID(faPatch);

loadV.boundaryFieldRef()[pID] = load.primitiveField();
wV.boundaryFieldRef()[pID]    = w.primitiveField();

//wLapV.boundaryFieldRef()[pID]    = wLap.primitiveField();

acc = fac::d2dt2(w);
