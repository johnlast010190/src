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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::TimeScaleModels::equilibrium

Group
    grpLagrangianIntermediateMPPICTimeScaleSubModels

Description
    Equlibrium model for the time scale over which properties of a dispersed
    phase tend towards the mean value.

    Reference:
    \verbatim
        "A model for collisional exchange in gas/liquid/solid fluidized beds"
        P O'Rourke, P Zhao and D Snider
        Chemical Engineering Science
        Volume 64, Issue 8, Pages 1784-1797, April 2009
    \endverbatim

SourceFiles
    equilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef equilibrium_H
#define equilibrium_H

#include "submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{

/*---------------------------------------------------------------------------*\
                        Class equilibrium Declaration
\*---------------------------------------------------------------------------*/

class equilibrium
:
    public TimeScaleModel
{
public:

    //- Runtime type information
    TypeName("equilibrium");


    //- Constructors

        //- Construct from components
        equilibrium(const dictionary& dict);

        //- Construct a copy
        equilibrium(const equilibrium& hc);

        //- Construct and return a clone
        virtual autoPtr<TimeScaleModel> clone() const
        {
            return autoPtr<TimeScaleModel>
            (
                new equilibrium(*this)
            );
        }


    //- Destructor
    virtual ~equilibrium();


    //- Member Functions

        //- Time scale
        tmp<FieldField<Field, scalar>> oneByTau
        (
            const FieldField<Field, scalar>& alpha,
            const FieldField<Field, scalar>& r32,
            const FieldField<Field, scalar>& uSqr,
            const FieldField<Field, scalar>& f
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TimeScaleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
