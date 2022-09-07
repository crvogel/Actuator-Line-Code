/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::fv::flowBlockingSource

Description
    A large force is added in the areas where a solid boundary is supposed to be located

SourceFiles
    flowBlockingSource.C

\*---------------------------------------------------------------------------*/

#ifndef flowBlockingSource_H
#define flowBlockingSource_H

#include "fvOption.H"
#include "fvCFD.H" // provides vector type, e.g. for rotatePoint



//#include "fvOption.H"
//#include "fvMesh.H"
//#include "fvMatrix.H"
//#include "addToRunTimeSelectionTable.H"
//#include "HashPtrTable.H"
//#include "OFstream.H"
//#include "meshSearch.H"
//#include "interpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class flowBlockingSource Declaration
\*---------------------------------------------------------------------------*/

class flowBlockingSource
:
    public option
{

protected:

        //- Runtime pointer.
        const Time& runTime_;


       //- Current simulation time.
       word time;
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        flowBlockingSource(const flowBlockingSource&);

        //- Disallow default bitwise assignment
        void operator=(const flowBlockingSource&);
    
    //- Radius of the surface at a specific point in axis
    scalar radiusAtS(scalar s,List<List<scalar> > surfaceNonDimProfile);
    
    //- scalar field of volume fraction f
    volScalarField f;

    //- vector field of body Force 
    volVectorField bodyForce;

    //- scalar field of volume fraction f
    volScalarField sp;
    
    //- averageVelocity_
    vector averageVelocity_;
    

    //- Calculation of the (independent of time) sp field that takes
    //zero value outside the surface, large value inside, and 
    //a volume averaged value in the cells that have both solid and fluid
    
    void computeSpFieldCircle();    

    void computeSpFieldSquare();    

    void spSquareCylinder();

    void spCube(vector sentre,scalar sideX,scalar sideY, scalar sideZ);

    void spSphere(vector centre, scalar radius);

    void spEllipse(vector axisStart, vector axisEnd, scalar radius);
    
    void spAxisymmetric(vector axisStart,vector axisEnd,List<List<scalar> > surfaceNonDimProfile);

    // - CharacteristicLength
    scalar characteristicLength_;    

    // - Characteristic Velocity
    scalar characteristicVelocity_;    
    
    // - List of names of geometry components
    wordList geometryNames_;

    // - Max edge
    scalar maxEdge_;            

    // - Min edge
    scalar minEdge_;            
    
public:

    //- Runtime type information
    TypeName("flowBlockingSource");

    // Constructors

        //- Construct from components
        flowBlockingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~flowBlockingSource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
