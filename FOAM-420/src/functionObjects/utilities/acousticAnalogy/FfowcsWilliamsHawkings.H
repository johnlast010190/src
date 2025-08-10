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

Class
    Foam::functionObjects::FfowcsWilliamsHawkings

Group
    grpUtilityFunctionObjects

Description
    Implementation of the Ffowcs-Williams-Hawkings acoustic analogy for
    moving or stationary walls in incompressible flow as
    described in:
    \verbatim
        F. Farassat (2007).
        Derivation of Formulations of 1 and 1A of Farassat.
        NASA Technical Memorandum 214853.
    \endverbatim

    Example of function object specification:
    \verbatim
    FfowcsWilliamsHawkings1
    {
        type                    FfowcsWilliamsHawkings;
        functionObjectLibs      ("libutilityFunctionObjects.so");
        log                     true;

        patches                 ( cylinder1 cylinder2 );
        rhoRef                  1.225;
        cRef                    343;

        observers
        {
            patch1
            {
                type patch;
                name outlet;
                writeFile false;
                writeAsField true;
            }
            points1
            {
                type set;
                positions (
                            ( 1.5 0 0)
                            ( 1.5 1 0)
                            ( 1.5 2 0)
                            ( 1.5 3 0)
                          );
            }
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description | Required | Default value
        type         | type name: FfowcsWilliamsHawkings | yes |
        log          | write acoustic data to standard output | no | no
        patches      | patches included in the acoustic calculation | yes |
        rhoRef       | reference density | yes |
        cRef         | reference speed of sound | yes |
        pName        | pressure field name | no | p
        UName        | velocity field name | no | U
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    FfowcsWilliamsHawkings.C

\*---------------------------------------------------------------------------*/

#ifndef FfowcsWilliamsHawkings_H
#define FfowcsWilliamsHawkings_H
#include <map>
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcDdt.H"
//#include "acousticAnalogy/soundObserver.H"
#include "acousticAnalogy/observer/observer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class FfowcsWilliamsHawkings Declaration
\*---------------------------------------------------------------------------*/

class FfowcsWilliamsHawkings
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data
		word modelType_;
		bool writeComponent_;
        //- Initialised flag
        bool initialised_;
	    word modelName_;
	    label tAcc_; //temporal accuracy

		scalar dtAco_;
        //- Patches to sample
        labelHashSet patches_;
        //- Name of pressure field (optional)
        word pName_;

        //- Name of velocity field (optional)
        word UName_;

        //- Reference density needed for incompressible cases
        scalar rhoRef_;

        //- Reference speed of sound
        scalar cRef_;
		scalar dt_;
		scalar dto_;
		bool hasOfield_; //whether the old field values exists
		bool hasOOfield_; //whether the old-old field values exists

		std::map<label,vectorField> Up_;
		std::map<label,vectorField> Up0_;
		std::map<label,vectorField> Up00_;

		std::map<label,scalarField> pp_;
		std::map<label,scalarField> pp0_;
		std::map<label,scalarField> pp00_;

        //- List of observers
        PtrList<observer> observerSets_;
        //List<SoundObserver> observers_;

        //- Prime Pressure ( observer sets -> observers -> time steps )
        List<List<scalarList>> pPrime_;
		List<List<scalarList>> pThickness_;

        //- Start simulation time
        scalar startSamplingTime_;

        //- End simulation time
        scalar endSamplingTime_;

    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(Ostream& os) const;
		void setFields();
		void updateFields();
		void calUdotPdot
		(
			std::map<label,vectorField>& Udot,
			std::map<label,scalarField>& pDot
		);

		void updateState();

        //- Helper function to write acoustic data
        virtual void writeFfowcsWilliamsHawkings();

        //- Initialise the fields
        virtual void initialise();

        //- Return pressure field
        tmp<volScalarField> p() const;

        //- Return pressure field time derivative
        tmp<volScalarField> dpdt() const;

        //- Return velocity field
        tmp<volVectorField> U() const;

        //- Return velocity field time derivative
        tmp<volVectorField> dUdt() const;

        //- Return surface normal vector
        tmp<surfaceVectorField> n() const;

        tmp<scalarField> getSurfaceNormalVelocityGradient(const label& patchI)
        const;

        //- Disallow default bitwise copy construct
        FfowcsWilliamsHawkings(const FfowcsWilliamsHawkings&);

        //- Disallow default bitwise assignment
        void operator=(const FfowcsWilliamsHawkings&);


public:

    //- Runtime type information
    TypeName("FfowcsWilliamsHawkings");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        FfowcsWilliamsHawkings
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word &modelName="FWH"
        );

    //- Destructor
    virtual ~FfowcsWilliamsHawkings();


    // Member Functions

        //- Read the FfowcsWilliamsHawkings data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the acoustic data
        virtual bool write();

        //- Calculate the acoustic data
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
