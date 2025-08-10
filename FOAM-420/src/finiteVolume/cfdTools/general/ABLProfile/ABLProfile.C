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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "cfdTools/general/ABLProfile/ABLProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(ABLProfile, 0);
defineRunTimeSelectionTable(ABLProfile, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ABLProfile::ABLProfile(const fvPatch& patch)
:
    patch_(patch),
    Z_()
{}


ABLProfile::ABLProfile(const fvPatch& patch, const dictionary& dict)
:
    patch_(patch),
    Z_(new patchDistanceFunction(patch, dict))
{}

ABLProfile::ABLProfile(const ABLProfile& abl)
:
    patch_(abl.patch_),
    Z_(abl.Z_, false)
{}
// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<ABLProfile> ABLProfile::New
(
    const fvPatch& patch,
    const dictionary& dict
)
{
    word profileType = dict.lookup("profileType");

    const auto ctor = ctorTableLookup("function type", dictionaryConstructorTable_(), profileType);
    return autoPtr<ABLProfile>(ctor(patch, dict));
}

// * * * * * * * * * * * * * * * * * Functions* * * * * * * * * * * * * * * //

void ABLProfile::write(Ostream& os) const
{
    Z_->write(os);
}

// * * * * * * * * * * * * * * * * * Operators* * * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
