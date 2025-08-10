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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "velocityDampingConstraint.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(velocityDampingConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        velocityDampingConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::addDamping(fvMatrix<vector>& eqn)
{
    // Note: we want to add
    //      deltaU/deltaT
    // where deltaT is a local time scale:
    //  U/(cbrt of volume)
    // Since directly manipulating the diagonal we multiply by volume.

    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    scalarField& diag = eqn.diag();

    label nDamped = 0;

    forAll(U, cellI)
    {
        scalar magU = mag(U[cellI]);
        if (magU > UMax_)
        {
            scalar scale = sqr(Foam::cbrt(vol[cellI]));

            diag[cellI] += scale*(magU-UMax_);

            nDamped++;
        }
    }

    reduce(nDamped, sumOp<label>());

    Info<< type() << " " << name_ << " damped "
        << nDamped << " ("
        << 100*scalar(nDamped)/mesh_.globalData().nTotalCells()
        << "%) of cells" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::velocityDampingConstraint::velocityDampingConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addDamping(eqn);
}


void Foam::fv::velocityDampingConstraint::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::velocityDampingConstraint::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        UMax_ = readScalar(coeffs_.lookup("UMax"));

        if (coeffs_.found("UNames"))
        {
            coeffs_.lookup("UNames") >> fieldNames_;
        }
        else if (coeffs_.found("UName"))
        {
            word UName(coeffs_.lookup("UName"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
            fieldNames_ = wordList(1, "U");
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
