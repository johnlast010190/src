#ifndef faCFD_H
#define faCFD_H

#include "global/argList/parRun.H"

#include "db/Time/Time.H"
#include "faMesh/faMesh.H"
#include "finiteArea/fac/fac.H"
#include "faMatrices/faMatrices.H"
#include "finiteArea/fam/fam.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "interpolation/volSurfaceMapping/volSurfaceMapping.H"
#include "fields/faPatchFields/basic/calculated/calculatedFaPatchFields.H"
#include "fields/faPatchFields/basic/fixedValue/fixedValueFaPatchFields.H"

#include "include/OSspecific.H"
#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

#endif
