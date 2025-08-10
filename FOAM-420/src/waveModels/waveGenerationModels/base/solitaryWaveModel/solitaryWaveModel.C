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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2015 IH-Cantabria

\*---------------------------------------------------------------------------*/

#include "waveGenerationModels/base/solitaryWaveModel/solitaryWaveModel.H"
#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"
#include "fields/Fields/Field/SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(solitaryWaveModel, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::solitaryWaveModel::timeCoeff
(
    const scalar t
) const
{
    // Ramping not applicable to solitary waves
    return 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::solitaryWaveModel::solitaryWaveModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveGenerationModel(dict, mesh, patch, false),
    x_
    (
        patch.faceCentres().component(0)*cos(waveAngle_)
      + patch.faceCentres().component(1)*sin(waveAngle_)
    ),
    x0_(gMin(x_))
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::solitaryWaveModel::~solitaryWaveModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::solitaryWaveModel::readDict
(
    const dictionary& overrideDict
)
{
    if (waveGenerationModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::solitaryWaveModel::info(Ostream& os) const
{
    waveGenerationModel::info(os);

    os << "    x0: " << x0_ << nl;
}


// ************************************************************************* //
