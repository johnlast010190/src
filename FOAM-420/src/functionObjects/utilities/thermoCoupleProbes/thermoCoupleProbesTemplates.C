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
    (c) 2016 OpenCFD.Ltd

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::thermoCoupleProbes::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    if (Pstream::master())
    {
        if (newFileFormat_)
        {
            columnatedFileWriter& fileWriter = *fileWriters_[vField.name()];
            fileWriter.writeTime();
            forAll(*this, probeI)
            {
                fileWriter.writeDelimited(Ttc_[probeI]);
            }
            fileWriter.endLine();
        }
        else
        {
            unsigned int w = IOstream::defaultPrecision() + 7;
            OFstream& probeStream = *probeFilePtrs_[vField.name()];
            probeStream
                << setw(w)
                << vField.time().timeOutputValue();

            forAll(*this, probeI)
            {
                probeStream << ' ' << setw(w) << Ttc_[probeI];
            }
            probeStream << endl;
        }
    }
}


// ************************************************************************* //
