/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::fv::actuationLineSource

Description
    Actuation disk source

    \heading Source usage

    Example usage:
    \verbatim
    actuationLineSourceCoeffs
    {
        fieldNames      (U);        // names of fields to apply source
    }
    \endverbatim


SourceFiles
    actuationLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuationLineSource_H
#define actuationLineSource_H

#include "fvOption.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "fvCFD.H" // provides vector type, e.g. for rotatePoint
#include "addToRunTimeSelectionTable.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "meshSearch.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuationLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuationLineSource
:
    public option
{

protected:

    // Protected data

        // Constants

            //- Degrees to radians conversion factor.
            const scalar degRad;

            //- Revolutions per minute to radians per second conversion factor.
            const scalar rpmRadSec;

            //- Runtime pointer.
            const Time& runTime_;

        //- Current time step size.
        scalar dt;

        //- Current simulation time.
        word time;

        // Velocity sampling

            //- Instantiate octree mesh searching for polyMesh
            meshSearch ms;

            //- Dictionary to define interpolation method for velocity
            dictionary interpolationDict;

        //- Print debug info
        Switch debug_;

        //- Body force field applied to fluid by turbine.
        volVectorField bodyForce;

        //- Mark cells where sampling occurs or focres are applied
        volScalarField cellFinder;

        // Turbine Array Level Data (all variables start with lower case letter)
        // .

            //- Restart flag for rotational blade
            DynamicList<int> restart_flag;

            //- List of names of turbines in array.
            List<word> turbineName;

            //- Number of turbines in array.
            int numTurbines;

            //- List of names of turbine types in array.
            DynamicList<word> turbineType;

            //- List of locations of bases of turbines in array relative to
            //- origin (m).
            DynamicList<vector> baseLocation;

            //- List of number of actuator line points on blades of turbines
            //- in array (set to uniform for now--here for future upgrades).
            DynamicList<int> numBladePoints;

            //- List of description of actuator line point distribution types
            //- for each turbine.
            DynamicList<word> pointDistType;

            //- List of the sampling method of each type of turbine
            DynamicList<word> samplingMethod; 

            //- Sampling distance = this * smearing radius. (Default = 1.1)
            DynamicList<scalar> sampleDistScalar;

            //- Sampling distance with a constant value.
            DynamicList<scalar> sampleDistParam;

            //- List of number of sample points for each blade point of turbines
            //- in array (set to uniform for now--here for future upgrades).
            DynamicList<int> numSamplePoints;

            //- Switch for bladeForceCorrector.
            DynamicList<int> bladeForceCorrector_;

            //- List of body force normalization parameters for each turbine.
            DynamicList<scalar> epsilon;
            DynamicList<scalar> smearRadius; //(m)

            //- List of factors to scale sphere of action radius for each turbine.
            DynamicList<scalar> sphereRadiusScalar;

            //- Smearing model for the blade forces ('Global' = use the 
            //- fixed values that are provided globally, otherwise we can 
            // use a 'LocalChord' based or 'Circulation' based model)
            DynamicList<word> smearingModel;

            //- Fraction of the local chord that we smear over (epsilon =
            // localChordFraction * local Chord)
            DynamicList<scalar> localChordFraction;

            //- Power for the function of epsilon
            DynamicList<scalar> epsilonFunctionPow;

            //- List of tip/root loss correction type for each turbine.
            DynamicList<word> tipRootLossCorrType;
            DynamicList<scalar> c1Faxi;
            DynamicList<scalar> c2Faxi;
            DynamicList<scalar> c3Faxi;
            DynamicList<scalar> c1Ftan;
            DynamicList<scalar> c2Ftan;
            DynamicList<scalar> c3Ftan;

            //- Tip speed ratio for use in Shen tip loss correction
            DynamicList<scalar> tipSpeedRatioShen;

            //- List of dynamic gust model for each type of turbine
            DynamicList<word> dynamicGustModel; 

            //- List of the airfoil polar interpolation method of each type of turbine
            DynamicList<word> airfoilInterpolationMethod; 

            //- Rotor rotation direction as viewed from upwind.
            // Options are "cw" for clockwise and "ccw" for counter-clockwise.
            DynamicList<word> rotationDir;

            //- Write every "outputInterval" time steps or seconds.
            // Options are "timeStep" or "runTime".  "runTime" writes out as
            // closely to every "outputInterval" seconds as possible, but
            // doesn't adjust the time step to write on the exact interval.
            word outputControl;

            //- The inteveral over which to write out turbine data.
            scalar outputInterval;

            //- Last time when output written.
            scalar lastOutputTime;

            //- Last time step when output written.
            label outputIndex;

            //- Initial or fixed rotor speed (rpm).
            // A positive value means clockwise rotation for a clockwise
            // rotating turbine (see rotationDir above) or counter-clockwise
            // rotation for a counter-clockwise rotating turbine.
            DynamicList<scalar> rotSpeed;

            //- Initial blade 1 azimuth angle (degrees)
            // (looking from upwind to downwind, a positive azimuth angle makes
            // a clockwise movement if this is a clockwise rotating turbine (see
            // rotationDir above) or a counterclockwise movement if this is a
            // counter-clockwise rotating turbine).
            DynamicList<scalar> azimuth;

            //- Initial or fixed nacelle yaw angle.
            // Direction that the turbine is pointed in cardinal directions
            // (i.e. 0 = north, 90 = east, 180 = south, 270 = west) (degrees).
            // This is converted to radians in the more standard mathematical
            // convention of 0 degrees on the + x axis and positive degrees in
            // the counter-clockwise direction.
            DynamicList<scalar> nacYaw;

            //- Specify the fluid density (kg/m^3).
            // This turbine model is to be used with an incompressible solver,
            // so density divides out of the momentum equations.  Therefore,
            // turbine forces are given to the solver as force/density. To get
            // actual forces, torques, and power written to file, provide a
            // density by which to multiply.
            DynamicList<scalar> fluidDensity;

            //- Incoming Induced velocity.
            DynamicList<vector> rotorUref;

            //- Incoming Turbulent velocity.
            DynamicList<scalar> rotorUturb;

            //- pitchZero
            DynamicList<scalar> pitchZero;

            //- Number of distinct turbines in array.
            int numTurbinesDistinct;

            //- List of distinct names of turbine types in array.
            DynamicList<word> turbineTypeDistinct;

            //- ID label given to each distinct type of turbine in the array.
            DynamicList<label> turbineTypeID;

        // Turbine Level Data (all variables start with a capital letter).
            //*** THE FOLLOWING VARIABLES MATCH FAST INPUT FILE ***

            //- Number of blades;
            DynamicList<int> NumBl;

            //- Distance from rotor apex to blade tip (m).
            DynamicList<scalar> TipRad;

            //- Distance from rotor apex to blade root (m).
            DynamicList<scalar> HubRad;

            //- Distance from teeter pin to rotor apex (m).
            DynamicList<scalar> UndSling;

            //- Distance from nacelle yaw axis to teeter pin or rotor apex (m).
            DynamicList<scalar> OverHang;

            //- Height of tower top above ground (m).
            DynamicList<scalar> TowerHt;

            //- Vertical distance from tower-top to rotor shaft centerline (m).
            DynamicList<scalar> Twr2Shft;

            //- Shaft tilt-up angle (degrees).
            DynamicList<scalar> ShftTilt;

            //- Coning angle of blades (degrees) (one for each blade).
            DynamicList<List<scalar> > PreCone;

            //*** END OF FAST INPUT FILE VARIABLES ***

            //- Pitch rate of blades.
            DynamicList<scalar> PitchRate;

            //- Pitch amplitude of blades.
            DynamicList<scalar> PitchAmplitude;

            //- Pitch controller type.
            DynamicList<word> PitchControllerType;

            //- List of airfoils that compose turbine blade;
            DynamicList<List<word> > AirfoilType;

            // Lists of blade data for each turbine type.

                //- Overall blade data array.
                DynamicList<List<List<scalar> > > BladeData;

                //- Overall rotor design array.
                DynamicList<List<List<scalar> > > rotorDesign;

                //- Station along blade in which Cl and Cd table is given (m).
                DynamicList<DynamicList<scalar> > BladeStation;

                //- Reynolds number at this station.
                DynamicList<DynamicList<scalar> > BladeP1;

                //- Tubulence intensity at this station.
                DynamicList<DynamicList<scalar> > BladeP2;

                //- Airfoil type ID at this station.
                DynamicList<DynamicList<label> > BladeAirfoilTypeID;

                //- Station along blade in which design information is given (m).
                DynamicList<DynamicList<scalar> > DesignStation;

                //- Chord at the design station (m).
                DynamicList<DynamicList<scalar> > DesignChord;

                //- Twist at the design station (degrees).
                DynamicList<DynamicList<scalar> > DesignTwist;

                //- Twist at the design station (degrees).
                DynamicList<DynamicList<scalar> > DesignTC;

            //- List of locations of the intersection of the tower axis and the shaft
            //- centerline relative to the origin (m).
            DynamicList<vector> towerShaftIntersect;

            //- List of locations of the rotor apex relative to the origin (m).
            DynamicList<vector> rotorApex;

        // Airfoil Level Data (all variables start with a lower-case letter).

            //- Number of distinct airfoils being used in turbines in array.
            int numAirfoilsDistinct;

            //- List of distinct type of airfoils amongst all turbines in array.
            DynamicList<word> airfoilTypesDistinct;

            //- Overall airfoil data array.
            DynamicList<List<List<scalar> > > airfoilData;

            //- Angle-of-attack.
            DynamicList<DynamicList<scalar> > airfoilAlpha;

            //- Lift Coefficient.
            DynamicList<DynamicList<scalar> > airfoilCl;

            //- Drag Coefficient.
            DynamicList<DynamicList<scalar> > airfoilCd;

        //- List of list of labels or cells within sphere of action.
        DynamicList<List<DynamicList<label> > > allSphereCellsI;

        //- Total actuator line points in array.
        int totBladePoints;

        //- An indicator of shaft direction.
        // The convention is that when viewed from upwind, the rotor turns
        // clockwise for positive rotation angles, regardless of if it is an
        // upwind or downwind turbine. uvShaft is found by subtracting the rotor
        // apex location from the tower shaft intersection point. This vector
        // switches direciton depending on if the turbine is upwind or downwind,
        // so this uvShaftDir multiplier makes the vector consistent no matter
        // what kind of turbine.
        DynamicList<scalar> uvShaftDir;

        //- Unit vector pointing along the rotor shaft (axis of blade rotation).
        DynamicList<vector> uvShaft;

        //- Unit vector pointing along the tower (axis of yaw).
        DynamicList<vector> uvTower;

        //- Actuator element width.
        DynamicList<DynamicList<scalar> > db;

        //- Actuator line point locations with respect to origin.
        DynamicList<List<List<vector> > > bladePoints;

        //- Blade radius away from rotor apex.  Must take into account coning.
        DynamicList<List<List<scalar> > > bladeRadius;

        //- Chord length at each actuator point.
        DynamicList<List<List<scalar> > > Chord;

        //- Twist angle at each actuator point.
        DynamicList<List<List<scalar> > > Twist;

        //- TC ratio at each actuator point.
        DynamicList<List<List<scalar> > > TCratio;

        // function

            //- Change in azimuth each time step.
            DynamicList<scalar> deltaAzimuth;

        // function

            //- Current pitch.
            DynamicList<scalar> pitch;

        // function

            //- deltaNacYaw
            DynamicList<scalar> deltaNacYaw;

            //- Blade-aligned vectors.
            // Three vectors for each blade of each turbine that define the local
            // blade-aligned coordinate system. Vector 0 is along the blade pointed
            // from root to tip, vector 1 is in the tangential direction (direction
            // of blade rotation) where positive points in the direction opposite
            // rotation if the rotor turns clockwise as viewed from upstream, and
            // points orthogonal to vector 0 and 1 and points toward downstream (but
            // vector 2 is not perfectly aligned with downstream due to rotor coning
            // and nacelle tilt).
            DynamicList<List<List<vector> > > bladeAlignedVectors;

            //- Chord-aligned vectors.
            // [turbine][blade][point][3 base-vectors] unit vectors
            DynamicList<List<List<List<vector> > > > chordAlignedVectors;

        // function

            //- Actuator line point locations with respect to origin.
            // [turbines][blades][AL-points][3 sample-points] coords
            DynamicList<List<List<List<vector> > > > samplePoints;

        // function

            //- samplePointsI
            // [turbines][blades][AL-points][3 sample-points]
            DynamicList<List<List<List<label> > > > samplePointsI;

        // function

            //- samplePointsU
            // [turbines][blades][AL-points][3 sample-points] velocity vectors
            DynamicList<List<List<List<vector> > > > samplePointsU;

            //- samplePointsUChord
            // [turbines][blades][AL-points][3 sample-points] velocity vectors in chord frame
            DynamicList<List<List<List<vector> > > > samplePointsUChord;

        // function

            //- Sample distance, rs, at each actuator point.
            DynamicList<List<List<scalar> > > rs;

            //- q at each actuator point.
            DynamicList<List<List<scalar> > > qC;

            //- Relative velocity at each actuator point.
            DynamicList<List<List<scalar> > > VmagC;

            //- Angle of attack at each actuator point.
            DynamicList<List<List<scalar> > > alphaC;

            //- Circulation at each actuator point.
            DynamicList<List<List<scalar> > > GammaC;

            //- Gaussian width parameter at each actuator point. 
            DynamicList<List<List<scalar> > > epsilonC;

            //- Gaussian smearing radius at each actuator point. 
            DynamicList<List<List<scalar> > > smearingRadiusC;

        // function

            //- Thrust/density on turbine.
            DynamicList<scalar> thrust;

            //- Total torque/density on turbine.
            DynamicList<scalar> torqueRotor;

            //- Power/density on turbine.
            DynamicList<scalar> powerRotor;

        // function

            //- Coefficient of lift at each actuator point.
            DynamicList<List<List<scalar> > > Cl;

            //- Coefficient of drag at each actuator point.
            DynamicList<List<List<scalar> > > Cd;

            //- Lift/density at each actuator point.
            DynamicList<List<List<scalar> > > lift;

            //- Drag/density at each actuator point.
            DynamicList<List<List<scalar> > > drag;

        // function

            //- Blade force at each actuator point.
            DynamicList<List<List<vector> > > bladeForce;

            //- Blade force corrector at each actuator point.
            DynamicList<List<List<scalar> > > bladeForceCorrector;

            //- Axial force/density at each actuator point
            DynamicList<List<List<scalar> > > axialForce;

            //- Tangential force/density at each actuator point.
            DynamicList<List<List<scalar> > > tangentialForce;

            //- Normal force/density at each actuator point.
            DynamicList<List<List<scalar> > > normalForce;

            //- Chordwise force/density at each actuator point.
            DynamicList<List<List<scalar> > > chordwiseForce;


        // Output files

        //- List of setup parameters.
        OFstream* setupFile_;

            // blade


                //- List of output files for rs.
                OFstream* radiusCFile_;

                //- List of output files for rs.
                OFstream* rsFile_;

                //- List of output files for angle of attack.
                OFstream* alphaCFile_;

                //- List of output files for relative velocity.
                OFstream* VmagCFile_;

                //- List of output files for q.
                OFstream* qCFile_;

                //- List of output files for epsilon.
                OFstream* epsilonCFile_;

                //- List of output files for smearing radius.
                OFstream* smearingRadiusCFile_;

                //- List of output files for circulation strength.
                OFstream* gammaCFile_;

                //- List of output files for coefficient of lift.
                OFstream* ClFile_;

                //- List of output files for coefficient of drag.
                OFstream* CdFile_;

                //- List of output files for lift/density.
                OFstream* liftFile_;

                //- List of output files for drag/density.
                OFstream* dragFile_;

                //- List of output files for axial force/density.
                OFstream* axialForceFile_;

                //- List of output files for tangential/density.
                OFstream* tangentialForceFile_;

                //- List of output files for normal force/density.
                OFstream* normalForceFile_;

                //- List of output files for chordiwse force/density.
                OFstream* chordwiseForceFile_;

                //- List of output files for chordiwse force/density.
                OFstream* bladeForceCorrectorFile_;

            // turbine

                //- List of output files for total aerodynamic torque/density.
                OFstream* torqueRotorFile_;

                //- List of output files for total thrust/density.
                OFstream* thrustFile_;

                //- List of output files for total power/density.
                OFstream* powerRotorFile_;

//              //- List of output files for probe index.
//              OFstream* probeCellIndexFile_;

//              //- List of output files for samplePoints locations.
                OFstream* samplePointsFile_;

//              //- List of output files for samplePoints velocities.
//              OFstream* samplePointsUFile_;


private:

    // Private Member Functions

        //- Print variables for debugging.
        void printDebug();

        //- Calculate the blade rotation rate
        //  (fixed rate based on user-specified initial rotor speed).
        void computeRotSpeed();

        //- Calculate the blade pitch.
        void computePitch();

        //- Rotate the blades.
        void rotateBlades();

        //- Calculate the nacelle yaw position
        //  (fixed based on user-specified initial yaw direction).
        void computeNacYaw();

        //- Yaw the nacelle.
        void yawNacelle();

        //- Compute chord-aligned vectors.
        void computeChordAlignedVectors();

        //- Compute sample points.
        void computeSamplePoints();

        //- Compute blade force.
        void computeBladeForce();

        //- Compute body force.
        void computeBodyForce();

        //- Compute body force.
        void computeBladeForceCorrector();

        //- Rotates a point about a rotation axis and rotation point by the
        // specified angle in radians.
        vector
        rotatePoint
        (
            vector point,
            vector rotationPoint,
            vector axis,
            scalar angle
        );

        //- Perform interpolation.
        scalar
        interpolate
        (
            scalar xNew,
            DynamicList<scalar>& xOld,
            DynamicList<scalar>& yOld
        );
        label
        interpolate
        (
            scalar xNew,
            DynamicList<scalar>& xOld,
            DynamicList<label>& yOld
        );
    List<label>    
        extrapolate
        (
            scalar xNew,
            DynamicList<scalar>& xOld,
            DynamicList<label>& yOld
        );


        //- Open turbine data output files.
        void openOutputFiles();

        //- Write turbine information to file.
        void printOutputFiles();

        //- Disallow default bitwise copy construct
        actuationLineSource(const actuationLineSource&);

        //- Disallow default bitwise assignment
        void operator=(const actuationLineSource&);


public:

    //- Runtime type information
    TypeName("actuationLineSource");


    // Constructors

        //- Construct from components
        actuationLineSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~actuationLineSource()
    {}


    // Member Functions

        // Public Functions

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


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
