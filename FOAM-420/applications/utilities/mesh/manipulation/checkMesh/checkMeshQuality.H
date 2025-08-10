#include "fvMesh/fvMesh.H"

namespace Foam
{
    class surfaceWriter;

    label checkMeshQuality
    (
        const fvMesh&,
        const dictionary&,
        const wordList&,
        const bool&,
        const bool&,
        const autoPtr<surfaceWriter>&
    );
}
