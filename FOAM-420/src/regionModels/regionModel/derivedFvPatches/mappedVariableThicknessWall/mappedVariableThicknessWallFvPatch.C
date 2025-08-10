/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "derivedFvPatches/mappedVariableThicknessWall/mappedVariableThicknessWallFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "regionModel1D/regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedVariableThicknessWallFvPatch, 0);
    addToRunTimeSelectionTable
    (
        fvPatch,
        mappedVariableThicknessWallFvPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedVariableThicknessWallFvPatch::
makeDeltaCoeffs(scalarField& dc) const
{
    const mappedVariableThicknessWallPolyPatch& pp =
        refCast<const mappedVariableThicknessWallPolyPatch>(patch());

    typedef regionModels::regionModel1D modelType;

    const modelType& region1D =
        patch().boundaryMesh().mesh().time().lookupObject<modelType>
        (
            "thermalBaffleProperties"
        );

    dc = 2.0/(pp.thickness()/region1D.nLayers());
}


// ************************************************************************* //
