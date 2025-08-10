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
    (c) 2014 Hannes Kroeger <hannes@kroegeronline.net>
 *
 */

/*---------------------------------------------------------------------------*/

#ifndef volFieldSampler_H
#define volFieldSampler_H

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"

namespace Foam
{

//- Class used for sampling volFields
template<class Type>
class volFieldSampler
: public List<Field<Type>>
{
    //- Name of this collection of values
    const word name_;

public:

    //- Construct interpolating field to the sampleSets
    volFieldSampler
    (
        const word& interpolationScheme,
        const GeometricField<Type, fvPatchField, volMesh>& field,
        const PtrList<sampledSet>& samplers
    )
        :
        List<Field<Type>>(samplers.size()),
        name_(field.name())
    {
        autoPtr<interpolation<Type>> interpolator
        (
            interpolation<Type>::New(interpolationScheme, field)
        );

        forAll(samplers, setI)
        {
            Field<Type>& values = this->operator[](setI);
            const sampledSet& samples = samplers[setI];

            values.setSize(samples.size());
            forAll(samples, sampleI)
            {
                const point& samplePt = samples[sampleI];
                label cellI = samples.cells()[sampleI];
                label faceI = samples.faces()[sampleI];

                if (cellI == -1 && faceI == -1)
                {
                    // Special condition for illegal sampling points
                    values[sampleI] = pTraits<Type>::max;
                }
                else
                {
                    values[sampleI] = interpolator().interpolate
                                      (
                                          samplePt,
                                          cellI,
                                          faceI
                                      );
                }
            }
        }
    }

    volFieldSampler
    (
        const List<Field<Type>>& values,
        const word& name
    )
        :
        List<Field<Type>>(values),
        name_(name)
    {}

    inline const word& name() const
    {
        return name_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
