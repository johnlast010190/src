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

Class
    Foam::TomiyamaSwarm

Description
    Swarm correction of Tomiyama et al.

    Reference:
    \verbatim
        "Drag Coefficients of Bubbles. 2nd Report. Drag Coefficient for a Swarm
         of Bubbles and its Applicability to Transient Flow."
        Tomiyama, A., Kataoka, I., Fukuda, T., and Sakaguchi, T.,
        Nippon Kikai Gakkai Ronbunshu
        Volume 61, Issue 588, 1995, pp. 2810-2817
    \endverbatim

SourceFiles
    TomiyamaSwarm.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaSwarm_H
#define TomiyamaSwarm_H

#include "interfacialModels/swarmCorrections/swarmCorrection/swarmCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{

/*---------------------------------------------------------------------------*\
                         Class TomiyamaSwarm Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaSwarm
:
    public swarmCorrection
{
    // Private data

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;

        //- Constant exponent
        const dimensionedScalar l_;


public:

    //- Runtime type information
    TypeName("Tomiyama");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaSwarm
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~TomiyamaSwarm();


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace swarmCorrections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
