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
    Foam::TimeScaleModels::isotropic

Group
    grpLagrangianIntermediateMPPICTimeScaleSubModels

Description
    Model for the time scale over which the velocity field of a dispersed phase
    tends towards an isotropic distribution.

    Reference:
    \verbatim
        "Inclusion of collisional return-to-isotropy in the MP-PIC method"
        P O'Rourke and D Snider
        Chemical Engineering Science
        Volume 80, Issue 0, Pages 39-54, December 2012
    \endverbatim

SourceFiles
    isotropic.C

\*---------------------------------------------------------------------------*/

#ifndef isotropic_H
#define isotropic_H

#include "submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{

/*---------------------------------------------------------------------------*\
                        Class isotropic Declaration
\*---------------------------------------------------------------------------*/

class isotropic
:
    public TimeScaleModel
{
public:

    //- Runtime type information
    TypeName("isotropic");


    //- Constructors

        //- Construct from components
        isotropic(const dictionary& dict);

        //- Construct a copy
        isotropic(const isotropic& hc);

        //- Construct and return a clone
        virtual autoPtr<TimeScaleModel> clone() const
        {
            return autoPtr<TimeScaleModel>
            (
                new isotropic(*this)
            );
        }


    //- Destructor
    virtual ~isotropic();


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
