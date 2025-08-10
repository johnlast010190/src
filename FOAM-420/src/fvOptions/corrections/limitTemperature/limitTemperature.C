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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "limitTemperature.H"
#include "fvMesh/fvMesh.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTemperature, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTemperature::limitTemperature
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    TName_(coeffs_.lookupOrDefault<word>("field", "T")),
    Tmin_(readScalar(coeffs_.lookup("min"))),
    Tmax_(readScalar(coeffs_.lookup("max")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitTemperature::initialise()
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());
    }
    else
    {
    fieldNames_.setSize(1, TName_);
    }

    applied_.setSize(1, false);

    return true;
}


bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("min") >> Tmin_;
        coeffs_.lookup("max") >> Tmax_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitTemperature::correct(volScalarField& he)
{

    scalarField Tmin(cells_.size(), Tmin_);
    scalarField Tmax(cells_.size(), Tmax_);


    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);
    const scalarField pSubSet(thermo.p(), cells_);
    scalarField heMin(thermo.he(pSubSet, Tmin, cells_));
    scalarField heMax(thermo.he(pSubSet, Tmax, cells_));

    scalarField& hec = he.primitiveFieldRef();

    forAll(cells_, i)
    {
        label celli = cells_[i];
        hec[celli]= max(min(hec[celli], heMax[i]), heMin[i]);
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& bf = he.boundaryFieldRef();

        forAll(bf, patchi)
        {
        fvPatchScalarField& hep = bf[patchi];

        if (!hep.fixesValue())
        {
            const scalarField& pp = thermo.p().boundaryField()[patchi];

            scalarField Tminp(pp.size(), Tmin_);
            scalarField Tmaxp(pp.size(), Tmax_);

            scalarField heMinp(thermo.he(pp, Tminp, patchi));
            scalarField heMaxp(thermo.he(pp, Tmaxp, patchi));

            forAll(hep, facei)
            {
            hep[facei] =
                max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
            }
        }
        }
    }
    }
    else
    {
    scalarField& Tec = he.primitiveFieldRef();

    forAll(cells_, i)
    {
        label celli = cells_[i];
        Tec[celli]= max(min(Tec[celli], Tmax[i]), Tmin[i]);
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& bf = he.boundaryFieldRef();

        forAll(bf, patchi)
        {
        fvPatchScalarField& Tep = bf[patchi];

        if (!Tep.fixesValue())
        {
            scalarField Tminp(he.size(), Tmin_);
            scalarField Tmaxp(he.size(), Tmax_);

            forAll(Tep, facei)
            {
            Tep[facei] =
                max(min(Tep[facei], Tmaxp[facei]), Tminp[facei]);
            }
        }
        }
    }

    }
}


// ************************************************************************* //
