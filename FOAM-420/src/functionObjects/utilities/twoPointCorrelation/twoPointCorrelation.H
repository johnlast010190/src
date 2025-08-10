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

#ifndef functionObjects_twoPointCorrelation_H
#define functionObjects_twoPointCorrelation_H

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "sampledSet/cloud/cloudSet.H"
#include "sampledSet/sampledSets/sampledSets.H"
#include "meshSearch/meshSearch.H"
#include "coordinate/systems/coordinateSystem.H"
#include "twoPointCorrelation/volFieldSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//#include <boost/concept_check.hpp>

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/


class twoPointCorrelation
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    //- Mesh search engine
    autoPtr<meshSearch> searchEngine_;

    point p0_;
    vector directionSpan_;
    label np_;

    vector homogeneousTranslationUnit_;
    label nph_;

    autoPtr<coordinateSystem> csys_;

    PtrList<sampledSet> lines_;

    // Merging structures
    autoPtr<scalarField> x_;
    PtrList<coordSet> masterSampledSets_;
    labelListList indexSets_;

    autoPtr<tensorField> correlationCoeffs_;
    scalar totalTime_;

    //- Output file header information
    virtual void writeFileHeader(Ostream& os) const;

    void createInterpolators();

    void combineSampledSets
    (
        PtrList<coordSet>& masterSampledSets,
        labelListList& indexSets
    );

    void resetAveraging();

    template<class T>
    void combineSampledValues
    (
        const volFieldSampler<T>& sampledField,
        const labelListList& indexSets,
        autoPtr<volFieldSampler<T>>& masterField
    );

    template<class Type>
    autoPtr<volFieldSampler<Type>>
    sample
    (
        PtrList<sampledSet>& sets,
        const GeometricField<Type, fvPatchField, volMesh>& field,
        const labelListList& indexSets
    );

public:
    //- Runtime type information
    TypeName("twoPointCorrelation");

    // Constructors
    twoPointCorrelation
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~twoPointCorrelation();

    //- Read and set the function object if its data have changed
    virtual bool read(const dictionary&);

    virtual bool execute();

    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
