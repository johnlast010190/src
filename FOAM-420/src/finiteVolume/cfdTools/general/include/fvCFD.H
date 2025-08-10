#ifndef fvCFD_H
#define fvCFD_H

#include "global/argList/parRun.H"

#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/fvc/fvc.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvm.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedFluxPressure/fixedFluxPressureFvPatchScalarField.H"
#include "cfdTools/general/constrainHbyA/constrainHbyA.H"
#include "cfdTools/general/constrainPressure/constrainPressure.H"
#include "cfdTools/general/adjustPhi/adjustPhi.H"
#include "cfdTools/general/findRefCell/findRefCell.H"
#include "global/constants/constants.H"

#include "include/OSspecific.H"
#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

#endif
