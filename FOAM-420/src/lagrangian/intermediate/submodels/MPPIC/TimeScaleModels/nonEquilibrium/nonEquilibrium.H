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
    Foam::TimeScaleModels::nonEquilibrium

Group
    grpLagrangianIntermediateMPPICTimeScaleSubModels

Description
    Non-Equlibrium model for the time scale over which properties of a
    dispersed phase tend towards the mean value.

    Reference:
    \verbatim
        "An improved collision damping time for MP-PIC calculations of dense
        particle flows with applications to polydisperse sedimenting beds and
        colliding particle jets"
        P O'Rourke and D Snider
        Chemical Engineering Science
        Volume 65, Issue 22, Pages 6014-6028, November 2010
    \endverbatim

SourceFiles
    nonEquilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef nonEquilibrium_H
#define nonEquilibrium_H

#include "submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{

/*---------------------------------------------------------------------------*\
                        Class nonEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class nonEquilibrium
:
    public TimeScaleModel
{
public:

    //- Runtime type information
    TypeName("nonEquilibrium");


    //- Constructors

        //- Construct from components
        nonEquilibrium(const dictionary& dict);

        //- Construct a copy
        nonEquilibrium(const nonEquilibrium& hc);

        //- Construct and return a clone
        virtual autoPtr<TimeScaleModel> clone() const
        {
            return autoPtr<TimeScaleModel>
            (
                new nonEquilibrium(*this)
            );
        }


    //- Destructor
    virtual ~nonEquilibrium();


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
