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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/
#include "layersProperties/layersProperties.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layersProperties::calculateTotalResistance() const
{
    rtot_ = 0.0;

    if (layerNames_.size()>0)
    {
        forAll(layerNames_, li)
        {
            // sum layer and contact resistances
            rtot_ += (t_[li]/max(VSMALL, kappa_[li]) + rcontact_[li]);
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layersProperties::layersProperties
(
    const dictionary& dict
)
:
    layerNames_(),
    t_(),
    kappa_(),
    rcontact_(),
    rtot_()
{
    read(dict);

    calculateTotalResistance();
}

Foam::layersProperties::layersProperties
(
    const layersProperties& lp
)
:
    layerNames_(lp.layerNames_),
    t_(lp.t_),
    kappa_(lp.kappa_),
    rcontact_(lp.rcontact_),
    rtot_(lp.rtot_)
{
    // copy constructor
}

// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

bool Foam::layersProperties::read(const dictionary& dict)
{

    if (dict.found("layers"))
    {
        const dictionary& layerDicts(dict.subDict("layers"));

        layerNames_.setSize(layerDicts.size());

        t_.setSize(layerDicts.size());
        kappa_.setSize(layerDicts.size());
        rcontact_.setSize(layerDicts.size());

        label nLayers = 0;

        forAllConstIter(dictionary, layerDicts, iter)
        {

            if (!iter().isDict())
            {
                continue;
            }

            const word& key = iter().keyword();
            const dictionary& ldict = iter().dict();

            layerNames_[nLayers] = key;

            ldict.lookup("thickness") >> t_[nLayers];
            ldict.lookup("kappa") >> kappa_[nLayers];
            ldict.lookup("rcontact") >> rcontact_[nLayers];

            nLayers++;
        }
    }

    return true;

}

void Foam::layersProperties::write(Ostream& os) const
{
    os.beginBlock("layers");
    forAll(t_, li)
    {
        os.beginBlock(layerNames_[li]);
        os.writeEntry("thickness", t_[li]);
        os.writeEntry("kappa", kappa_[li]);
        os.writeEntry("rcontact", rcontact_[li]);
        os.endBlock();
    }
    os.endBlock();

}


// ************************************************************************* //
