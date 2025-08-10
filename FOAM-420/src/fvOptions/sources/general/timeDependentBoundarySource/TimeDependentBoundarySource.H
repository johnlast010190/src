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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2021 Esi Ltd

Class
    Foam::fv::TimeDependentBoundarySource

Group
    grpFvOptionsSources

Description
    Time dependent boundary flux source. The specified values need to be used
    by an appropriate boundary condition, otherwise they have no effect and
    a warning will be given.

    Example of the source specification:

    \verbatim
    <Type>TimeDependentBoundarySourceCoeffs
    {
        surfaceMode      absolute; // specific
        fluxSource
        {
            h           table (( 0 0.005 ) (3.14159 0.01) (13.02 0.002));
            epsilon     constant ( 1e5 );
        }
        patches ( <patchNames> );
    }
    \endverbatim

    Notes

    Valid options for the \c surfaceMode entry include:
    - absolute: values are total flux across the patches specified (in \<equation units\> * m^3)
    - specific: values are flux density (in \<equation units\> * m)

    All patches will be used if unspecified.

SeeAlso
    Foam::fvOption

SourceFiles
    TimeDependentBoundarySource.C

\*---------------------------------------------------------------------------*/

#ifndef TimeDependentBoundarySource_H
#define TimeDependentBoundarySource_H

#include "primitives/Tuple2/Tuple2.H"
#include "interpolations/interpolationTable/interpolationTable.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declaration of classes

template<class Type>
class TimeDependentBoundarySource;


// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeDependentBoundarySource<Type>&
);


/*---------------------------------------------------------------------------*\
                     Class TimeDependentBoundarySource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeDependentBoundarySource
:
    public option
{
public:

    // Public data

        //- Enumeration for volume types
        enum surfaceModeType
        {
            smAbsolute,
            smSpecific
        };

        //- Word list of volume mode type names
        static const wordList surfaceModeTypeNames_;


protected:

    // Protected data

        //- Volume mode
        surfaceModeType surfaceMode_;

        //- Surface normalisation
        scalar areaNorm_;

        //- The time series being used for each equation
        HashTable<autoPtr<Function1<Type>>> timeSeriesList_;


    // Protected functions

        //- Helper function to convert from a word to a surfaceModeType
        surfaceModeType wordToSurfaceModeType(const word& stName) const;

        //- Helper function to convert from a surfaceModeType to a word
        word surfaceModeTypeToWord(const surfaceModeType& stType) const;

        //- Set the field time series
        void setTimeSeries(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("TimeDependentBoundarySource");


    // Constructors

        //- Construct from components
        TimeDependentBoundarySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        // Access

            //- Return const access to the volume mode
            inline const surfaceModeType& surfaceMode() const;


        // Edit

            //- Return access to the surface mode
            inline surfaceModeType& surfaceMode();

        // Evaluation

            //- Add flux contribution
            virtual void addBoundarySource
            (
                const word& fieldName,
                const label patchID,
                const Field<Type>& pf,
                Field<Type>& f,
                Field<Type>& df
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TimeDependentBoundarySource.C"
    #include "TimeDependentBoundarySourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeDependentBoundarySourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
