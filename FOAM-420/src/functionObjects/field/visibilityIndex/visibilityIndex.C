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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "visibilityIndex/visibilityIndex.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(visibilityIndex, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        visibilityIndex,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::visibilityIndex::calc()
{
    if
    (
        foundObject<volScalarField>(fieldName_)
    )
    {
        const volScalarField& rho = lookupObject<volScalarField>(rhoName_);
        const volScalarField& smoke = lookupObject<volScalarField>(fieldName_);

        // constant object viewed across smoke
    dimensionedScalar C("C", dimensionSet(0,0,0,0,0,0,0), C_);

        // mass specific extinction coeff
    dimensionedScalar km("km", dimensionSet(-1,2,0,0,0,0,0), kappam_);

    //Maximum visibility clipped to 30 m (CBSE guidelines)
    dimensionedScalar Lmax("Lmax", dimensionSet(0,1,0,0,0,0,0), Lmax_);

    Info<<"The visibility index C/K has been now computed."<<endl;

        return store
        (
            resultName_,
            min(Lmax, C/(km*max(smoke,SMALL/100)*rho))
        );
    }
    else
    {
        return false;
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::visibilityIndex::visibilityIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "smoke"),
    rhoName_("rho"),
    kappam_(8700),
    C_(3),
    Lmax_(30)
{
    setResultName("visibility", "smoke");
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::visibilityIndex::~visibilityIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::visibilityIndex::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    kappam_ = dict.lookupOrDefault<scalar>("kappam", 8700);
    C_ = dict.lookupOrDefault<scalar>("C", 3);
    Lmax_ = dict.lookupOrDefault<scalar>("Lmax", 30);

    Info<<"Visibility Index Calculation:"<<endl;
    Info<<"---Density field name:                     "<<rhoName_<<endl;
    Info<<"---Mass specific extinction coeff.[m2/kg]: "<<kappam_<<endl;
    Info<<"---C constant:                             "<<C_<<endl;
    Info<<"---Max visibility Length:                  "<<Lmax_<<endl;

    return true;
}


// ************************************************************************* //
