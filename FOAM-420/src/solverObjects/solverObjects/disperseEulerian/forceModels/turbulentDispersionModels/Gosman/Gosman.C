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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Gosman.H"
#include "disperseEulerian/phase/phase.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "disperseEulerian/forceModels/dragModels/dragModel/dragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{
    defineTypeNameAndDebug(Gosman, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Gosman,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decoupledEulerian::Gosman::Gosman
(
    const dictionary& dict,
    const phase& phase,
    const bool registerObject
)
:
    turbulentDispersionModel(dict, phase, registerObject),
    Sct_("Sct", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decoupledEulerian::Gosman::~Gosman()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::decoupledEulerian::Gosman::D() const
{
    const fvMesh& mesh(phase_.alphad().mesh());
    const dragModel&
        drag
        (
            mesh.lookupObject<dragModel>
            (
                "dragModel" + Foam::word(phase_.name())
            )
        );

    const turbulenceModel& turbulence(mesh.lookupObject<turbulenceModel>(turbulenceModel::propertiesName));

    return
        0.75
       *drag.CdRe()
       *phase_.muc()
       *turbulence.nut()/Sct_
       /(sqr(phase_.diam())
       *phase_.rhod());
}


// ************************************************************************* //
