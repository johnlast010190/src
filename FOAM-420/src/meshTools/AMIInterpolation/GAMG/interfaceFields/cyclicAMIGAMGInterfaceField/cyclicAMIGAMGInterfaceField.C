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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/GAMG/interfaceFields/cyclicAMIGAMGInterfaceField/cyclicAMIGAMGInterfaceField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicAMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicAMIGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicAMIGAMGInterfaceField::cyclicAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicAMIInterface_(refCast<const cyclicAMIGAMGInterface>(GAMGCp)),
    rank_(0)
{
    const cyclicAMILduInterfaceField& p =
        refCast<const cyclicAMILduInterfaceField>(fineInterface);

    rank_ = p.rank();
}


Foam::cyclicAMIGAMGInterfaceField::cyclicAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, rank),
    cyclicAMIInterface_(refCast<const cyclicAMIGAMGInterface>(GAMGCp)),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::cyclicAMIGAMGInterfaceField::~cyclicAMIGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicAMIGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Get neighbouring field
    scalarField pnf
    (
        cyclicAMIInterface_.nbrPatch().interfaceInternalField(psiInternal)
    );

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    scalarField pf(psiInternal,cyclicAMIInterface_.faceCells());

    if (cyclicAMIInterface_.owner())
    {
        if (cyclicAMIInterface_.AMI().applyLowWeightCorrection())
        {
            const vectorField& nHat = cyclicAMIInterface_.
                AMI().srcNf();
            Field<tensor> tensField(this->size(), I);
            Field<tensor> trans(tensField - 2.0*sqr(nHat));
            Field<scalar> mirrorField(Foam::transform(trans, pf));

            UList<scalar> pif
            (
                mirrorField
            );
            pnf = cyclicAMIInterface_.AMI().interpolateToSource
            (
                pnf,
                pif
            );
        }
        else
        {
            pnf = cyclicAMIInterface_.AMI().interpolateToSource(pnf);
        }
    }
    else
    {
        if (cyclicAMIInterface_.nbrPatch().AMI().applyLowWeightCorrection())
        {
            const vectorField& nHat = cyclicAMIInterface_.nbrPatch().
                AMI().tgtNf();
            Field<tensor> tensField(this->size(), I);
            Field<tensor> trans (tensField - 2.0*sqr(nHat));
            Field<scalar> mirrorField(Foam::transform(trans, pf));

            UList<scalar> pif
            (
                mirrorField
            );
            pnf = cyclicAMIInterface_.nbrPatch().AMI().interpolateToTarget
            (
                pnf,
                pif
            );
        }
        else
        {
            pnf = cyclicAMIInterface_.nbrPatch().AMI().
                interpolateToTarget(pnf);
        }
    }

    this->addToInternalField(result, !add, coeffs, pnf);
}


// ************************************************************************* //
