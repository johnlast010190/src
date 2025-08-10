#include "graphField/writeCellGraph.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvMesh.H"
#include "graph/graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeCellGraph
(
    const volScalarField& vsf,
    const word& graphFormat
)
{
    fileName path(vsf.time().path()/"graphs"/vsf.time().timeName());
    mkDir(path);

    graph
    (
        vsf.name(),
        "x",
        vsf.name(),
        vsf.mesh().C().primitiveField().component(vector::X),
        vsf.primitiveField()
    ).write(path/vsf.name(), graphFormat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
