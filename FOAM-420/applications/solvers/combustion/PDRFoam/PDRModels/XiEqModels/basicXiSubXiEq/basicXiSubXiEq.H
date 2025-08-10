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
    (c) 2011-2012 OpenFOAM Foundation

Class
    Foam::XiEqModels::basicSubGrid

Description
    Basic sub-grid obstacle flame-wrinking enhancement factor model.
    Details supplied by J Puttock 2/7/06.

    <b> Sub-grid flame area generation </b>

    \f$ n = N - \hat{\dwea{\vec{U}}}.n_{s}.\hat{\dwea{\vec{U}}} \f$
    \f$ n_{r} = \sqrt{n} \f$

    where:

        \f$ \hat{\dwea{\vec{U}}} = \dwea{\vec{U}} / \vert \dwea{\vec{U}}
        \vert \f$

        \f$ b = \hat{\dwea{\vec{U}}}.B.\hat{\dwea{\vec{U}}} / n_{r} \f$

    where:

        \f$ B \f$ is the file "BT".

        \f$ N \f$ is the file "N".

        \f$  n_{s} \f$ is the file "ns".

    The flame area enhancement factor \f$ \Xi_{sub} \f$ is expected to
    approach:

    \f[
        \Xi_{{sub}_{eq}} =
            1 + max(2.2 \sqrt{b}, min(0.34 \frac{\vert \dwea{\vec{U}}
            \vert}{{\vec{U}}^{'}}, 1.6)) \times min(\frac{n}{4}, 1)
    \f]


SourceFiles
    basicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef basicSubGrid_H
#define basicSubGrid_H

#include "XiEqModel/XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class basicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class basicSubGrid
:
    public XiEqModel
{
    // Private data

        //- tblock
        volSymmTensorField BT_;

        //- Equilibrium Xi model due to turbulence
        autoPtr<XiEqModel> XiEqModel_;


    // Private Member Functions

        //- Disallow copy construct
        basicSubGrid(const basicSubGrid&);

        //- Disallow default bitwise assignment
        void operator=(const basicSubGrid&);


public:

    //- Runtime type information
    TypeName("basicSubGrid");


    // Constructors

        //- Construct from components
        basicSubGrid
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~basicSubGrid();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
