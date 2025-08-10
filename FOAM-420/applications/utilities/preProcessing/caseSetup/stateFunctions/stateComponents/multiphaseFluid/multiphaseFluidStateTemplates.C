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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
Xfer<HashTable<Type, stateFunction::interfacePair, stateFunction::interfacePair::hash>>
stateFunctions::multiphaseFluidState::assembleBinaryData
(
    word property,
    const Type& defaultValue,
    const List<dictionary>& matprops
)
{
    label nPhases = phases_.size();
    label nBinaryPairs = (nPhases*nPhases - nPhases)/2;

    HashTable<Type, stateFunction::interfacePair, stateFunction::interfacePair::hash>
        binaryData(nBinaryPairs);


    for (label n = 0; n < nPhases-1; n++)
    {
        for (label m = n+1; m < nPhases; m++)
        {
            //check for reciprocal data definitions, primary phase wins
            //if phases define incompatible reciprocals
            Type cDat(defaultValue);

            //do not insert entry at all if not available
            bool found = false;

            if (matprops[n].found("binaryPhaseData"))
            {
                if (matprops[n].subDict("binaryPhaseData").found(phases_[m]))
                {
                    const dictionary& bpdDict
                    (
                        matprops[n].subDict("binaryPhaseData")
                        .subDict(phases_[m])
                    );

                    if (bpdDict.found(property))
                    {
                        cDat = bpdDict.lookupOrDefault<Type>
                        (
                            property,
                            defaultValue
                        );

                        found = true;
                    }
                }
            }
            if (matprops[m].found("binaryPhaseData") && !found)
            {
                if (matprops[m].subDict("binaryPhaseData").found(phases_[n]))
                {

                    const dictionary& bpdDict
                    (
                        matprops[m].subDict("binaryPhaseData")
                        .subDict(phases_[n])
                    );

                    if (bpdDict.found(property))
                    {
                        cDat = bpdDict.lookupOrDefault<Type>
                        (
                            property,
                            defaultValue
                        );
                    }
                }
            }

            binaryData.insert(interfacePair(phases_[n], phases_[m]), cDat);
        }
    }

    return binaryData.xfer();
}

} // End namespace Foam

// ************************************************************************* //
