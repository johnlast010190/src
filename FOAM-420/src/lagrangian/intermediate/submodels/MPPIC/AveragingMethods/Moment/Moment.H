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
    Foam::AveragingMethods::Moment

Group
    grpLagrangianIntermediateMPPICAveragingMethods

Description
    Moment lagrangian averaging procedure.

    Point values and moments from the cell centroid are summed over
    computational cells. A linear function is generated which has the same
    integrated moment as that of the point data.

    The computed linear function is used to interpolate values within a cell.
    The gradient is calculated from the coefficients of the function, and is
    assumed constant over the cell.

SourceFiles
    Moment.C

\*---------------------------------------------------------------------------*/

#ifndef Moment_H
#define Moment_H

#include "submodels/MPPIC/AveragingMethods/AveragingMethod/AveragingMethod.H"
#include "meshes/pointMesh/pointMesh.H"
#include "meshes/polyMesh/polyMeshTetDecomposition/tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace AveragingMethods
{

/*---------------------------------------------------------------------------*\
                           Class Moment Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Moment
:
    public AveragingMethod<Type>
{
public:

    //- Public typedefs

        //- Gradient type
        typedef typename AveragingMethod<Type>::TypeGrad TypeGrad;


private:

    //- Private data

        //- Data mean
        Field<Type>& data_;

        //- X-data moment
        Field<Type>& dataX_;

        //- Y-data moment
        Field<Type>& dataY_;

        //- Z-data moment
        Field<Type>& dataZ_;

        //- Transform tensor from moment to gradient
        Field<symmTensor> transform_;

        //- Length scale for moment values
        Field<scalar> scale_;


    //- Private member functions

        //- Re-calculate gradient
        virtual void updateGrad();


public:

    //- Runtime type information
    TypeName("moment");


    //- Constructors

        //- Construct from components
        Moment
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        Moment(const Moment<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<AveragingMethod<Type>> clone() const
        {
            return autoPtr<AveragingMethod<Type>>
            (
                new Moment<Type>(*this)
            );
        }


    //- Destructor
    virtual ~Moment();


    //- Member Functions

        //- Add point value to interpolation
        void add
        (
            const point position,
            const tetIndices& tetIs,
            const Type& value
        );

        //- Interpolate
        Type interpolate
        (
            const point position,
            const tetIndices& tetIs
        ) const;

        //- Interpolate gradient
        TypeGrad interpolateGrad
        (
            const point position,
            const tetIndices& tetIs
        ) const;

        //- Return an internal field of the average
        tmp<Field<Type>> primitiveField() const;

        //- Return an internal field of the gradient
        tmp<Field<TypeGrad>> internalFieldGrad() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AveragingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/MPPIC/AveragingMethods/Moment/Moment.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
