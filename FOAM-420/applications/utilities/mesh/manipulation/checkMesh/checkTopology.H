#include "primitives/ints/label/label.H"
#include "memory/autoPtr/autoPtr.H"
#include "sampledSetWriters/writer.H"

namespace Foam
{
    class polyMesh;
    class surfaceWriter;

    label checkTopology
    (
        const polyMesh&,
        const bool,
        const bool,
        const autoPtr<surfaceWriter>&,
        const autoPtr<writer<scalar>>&
    );
}
