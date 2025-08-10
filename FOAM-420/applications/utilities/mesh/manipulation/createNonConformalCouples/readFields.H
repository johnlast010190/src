// Read vol fields
PtrList<volScalarField> vsFlds;
if (fields) ReadFields(mesh, objects, vsFlds);

PtrList<volVectorField> vvFlds;
if (fields) ReadFields(mesh, objects, vvFlds);

PtrList<volSphericalTensorField> vstFlds;
if (fields) ReadFields(mesh, objects, vstFlds);

PtrList<volSymmTensorField> vsymtFlds;
if (fields) ReadFields(mesh, objects, vsymtFlds);

PtrList<volTensorField> vtFlds;
if (fields) ReadFields(mesh, objects, vtFlds);


// Read surface fields
PtrList<surfaceScalarField> ssFlds;
if (fields) ReadFields(mesh, objects, ssFlds);

PtrList<surfaceVectorField> svFlds;
if (fields) ReadFields(mesh, objects, svFlds);

PtrList<surfaceSphericalTensorField> sstFlds;
if (fields) ReadFields(mesh, objects, sstFlds);

PtrList<surfaceSymmTensorField> ssymtFlds;
if (fields) ReadFields(mesh, objects, ssymtFlds);

PtrList<surfaceTensorField> stFlds;
if (fields) ReadFields(mesh, objects, stFlds);


// Read point fields
PtrList<pointScalarField> psFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, psFlds);

PtrList<pointVectorField> pvFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, pvFlds);

PtrList<pointSphericalTensorField> pstFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, pstFlds);

PtrList<pointSymmTensorField> psymtFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, psymtFlds);

PtrList<pointTensorField> ptFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, ptFlds);
