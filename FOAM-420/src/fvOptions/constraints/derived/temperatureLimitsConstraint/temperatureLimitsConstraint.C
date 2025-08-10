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
    (c) 2012-2013 OpenFOAM Foundation

\*----------------------------------------------------------------------------*/

#include "temperatureLimitsConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(temperatureLimitsConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        temperatureLimitsConstraint,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::temperatureLimitsConstraint::temperatureLimitsConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    Tmin_(readScalar(coeffs_.lookup("Tmin"))),
    Tmax_(readScalar(coeffs_.lookup("Tmax")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::temperatureLimitsConstraint::initialise()
{
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);
    fieldNames_.setSize(1, thermo.he().name());
    applied_.setSize(1, false);
    return true;
}


bool Foam::fv::temperatureLimitsConstraint::alwaysApply() const
{
    return true;
}


void Foam::fv::temperatureLimitsConstraint::correct(volScalarField& he)
{
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);

    if (he.name() == thermo.he().name())
    {
        scalarField Tmin(cells_.size(), Tmin_);
        scalarField Tmax(cells_.size(), Tmax_);
        const scalarField pSubSet(thermo.p(), cells_);
        scalarField heMin(thermo.he(pSubSet, Tmin, cells_));
        scalarField heMax(thermo.he(pSubSet, Tmax, cells_));

        scalarField& hec = he.ref();

        forAll(cells_, i)
        {
            label cellI = cells_[i];
            hec[cellI] = max(min(hec[cellI], heMax[i]), heMin[i]);
        }

        // handle boundaries in the case of 'all'
        if (selectionMode_ == smAll)
        {
            volScalarField::Boundary& bf = he.boundaryFieldRef();

            forAll(bf, patchI)
            {
                fvPatchScalarField& hep = bf[patchI];
                if (hep.fixesValue())
                {
                    // not over-riding fixed conditions
                    continue;
                }

                const scalarField& pp = thermo.p().boundaryField()[patchI];

                scalarField Tminp(pp.size(), Tmin_);
                scalarField Tmaxp(pp.size(), Tmax_);

                scalarField heMinp(thermo.he(pp, Tminp, patchI));
                scalarField heMaxp(thermo.he(pp, Tmaxp, patchI));

                forAll(hep, faceI)
                {
                    hep[faceI] =
                        max(min(hep[faceI], heMaxp[faceI]), heMinp[faceI]);
                }
            }
        }
    }
}


void Foam::fv::temperatureLimitsConstraint::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::temperatureLimitsConstraint::read(const dictionary& dict)
{
    Info<< "IN READ" << endl;
    if (option::read(dict))
    {
        coeffs_.readIfPresent("Tmin", Tmin_);
        coeffs_.readIfPresent("Tmax", Tmax_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
