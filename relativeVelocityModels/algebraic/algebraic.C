/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "algebraic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(algebraic, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, algebraic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::algebraic::algebraic
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    residualRe_("residualRe", dimless, dict),
    turbulenceCorrect_(dict.getOrDefault("turbulenceCorrect", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::algebraic::~algebraic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::algebraic::correct()
{
    volScalarField betac(alphac_*rhoc_);
    volScalarField betad(alphad_*rhod_);

    // Calculate the relative velocity of the continuous phase w.r.t the mean
    volVectorField Ucm(betad*Udm_/betac);

    // particle Reynolds number
    volScalarField Re_p(alphac_ * rhoc_ * mixture_.dd() * mag(Udm_ - Ucm) / (mixture_.nucModel().nu() * rhoc_) + residualRe_);

    volScalarField fd(1. + 0.15 * pow(Re_p, 0.687));

    forAll(Re_p, i)
    {
        if (Re_p[i] >= 1000)
        {
            fd[i] = 0.0183 * Re_p[i];
        }
    }

    Info << "max fd: " << max(fd) << endl;
    Info << "min fd: " << min(fd) << endl;

    const meshObjects::gravity& g = meshObjects::gravity::New(mixture_.U().db().time());

    volVectorField a(g - mixture_.U()*fvc::div(mixture_.U()) - fvc::ddt(mixture_.U()));

    volVectorField Udc((rhod_-rho()) * sqr(mixture_.dd()) / (mixture_.nucModel().nu() * rhoc_) / 18. / fd * a);

    const volScalarField &nut = mixture_.U().mesh().lookupObject<volScalarField>("nut");

    if(turbulenceCorrect_)
    {
        Udc -= nut * (fvc::grad(alphad_)/(alphad_+1e-8) - fvc::grad(alphac_)/(alphac_+1e-8));
    }
    
    //volScalarField tau_p(rhod_* mixture_.dd() * mixture_.dd() / 18 / (mixture_.nucModel().nu() * rhoc_));

    //volScalarField K(3 / 4 / mixture_.dd() * rhoc_ * Cd * pow(alphac_, -1.65) * (1 - alphac_) * mag(Udm_ - Ucm));

    Udm_ = betac / rho() * Udc;
    //Udm_ = Udc - alphac_ *Ucm;

    Info << "max Udm: " <<max(Udm_) << endl;
    Info << "min Udm: " <<min(Udm_) << endl;
}


// ************************************************************************* //
