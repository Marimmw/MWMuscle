rm -rf *
cmake ..
make


clear
make -j32
LIBGL_ALWAYS_SOFTWARE=1 ./MWmuscleWrapper

FIEPATH: \\wsl.localhost\Ubuntu\home\mariu\projects\MWmuscleWrapper

Git push
1. ~/projects/MWmuscleWrapper
2. git checkout dev
3. git status
4. git add -A    ODER git add .
5. git commit -m "update.........."
6. git push origin dev

0. git reset

UMITTELBAR
-> 

WEITER MIT:
- "exact" hand model
- Setup first hand model
    -> add parametrizatiable finger to hand
    -> ellipsoids as thick as bone in the middle part
- Moment Arm conpuation (per muscle or per joint)
- complex initial guess path through tori(+bodies)
- evtl scalable torus


- universal Joints- ball joint
- distance funktion zu körper -> variante zu mesh-mittelpunkt und meshoberflächefür nähesten body


Own Cpp Programm

- Moment arms
- Objects
	- Bodies
		# more meshes/warapper srufaces
		# add synthetic (anatomic) origin and frame
	- Joints
		- Spherical
		- (Translation)
		# revolute
		- (saddle)
		- live movment q
	- Surfaces
	- Muscles
		# reference bodies
		+ origins on muscles
		- Torus, ellisopid, cylinder, (capsule)
		- obj files or similar loading
- Penetration tests
- dynamic local, local muscle Point tracker
- constraint equations of each mesh/shapes
+ solver
+ System handling
# live adding bodies 3D
# Import/Export System
# Import/Export Data
# Results/Solving page
# Body/Joints/muscles overview topdownlist

+ System for parralellization
+ interpolating results for smoother movement
+ Show coordaxis


Further ToDO:
	- see derivatives of anybody studies



ERKENNTNISSE:

Via-Points:
	- its hard for the solver to get the point (viaPointTol=0.02, 30points -> 0.103m pointdistance)
		into the via point?!
Initial-Guess:
	- if muscle point is inside mesh, it might fuck up the solver
Muscle_Points:
	- the optimal points need to be found

Local Gammas
	- good -> muscle points should move less from step x to step x+1 -> stay where they are



push in lokal ellipsCoord:
-y -> yes
x -> yes
-z -> yes




________________________________________________________________________

PROGRAM FLOW DOCUMENTATION: MUSCLE SIMULATION FRAMEWORK
=======================================================

1. INITIALIZATION PHASE
-------------------------------------------------------
    1.1. Application Setup

    1.2. Scene Construction [Function: setupScene]
        - Mesh Instantiation:
            - Create geometric objects (Ellipsoids, Tori, Joints).
            - Define physical properties (Dimensions A/B/C, Colors, Names).
            - Set initial Global Position and Orientation matrices.
            - Set Flags: 'bIsViaPoint' (soft constraint) or 'bIsJointMesh'.
        - Muscle Instantiation [Class: SSMuscle]:
            - Define topology: Node count, Origin Mesh/Offset, Insertion Mesh/Offset.
            - Assign interaction lists (Obstacles, Via-Points).
        - Muscle Node Generation [Function: createMusclePoints]:
            - Calculate initial linear path between global Origin and Insertion.
            - Initialize 'MNodes' (Muscle Nodes) with global positions.
            - Assign fixed parents to start/end nodes; set internal nodes as free.

    1.3. Solver System Initialization [Class: CasadiSystem]
        - Instantiate Wrapper for IPOPT NLP solver.
        - Symbolic Setup [Function: setupCasadi...]:
            - Define Decision Variables (x): Node coordinates (x,y,z) and contact forces (eta).
            - Define Parameters (p): Mesh poses (Pos + RotMatrix), Muscle Anchors.
            - Define Constraints (g):
                * Obstacles: Signorini conditions (Complementarity: distance * force = 0).
                * Via-Points: (Soft)Min distance aggregation or hard tolerance constraints.
                * Physics: Euler-Lagrange equilibrium equations.
            - Define Objective (f): Minimize muscle length (or potential energy).
        - Compilation of NLP function.

    1.4. Data Container Preparation
        - Resize multidimensional vectors for result storage (Positions, Rotations, Colors, Guesses).
        - Memory Reservation [Function: mus->initializeSimulationMuscle]:
            - Reserve memory for internal node history vectors.


2. SIMULATION LOOP (TIME STEPPING)
-------------------------------------------------------
    Iterate 't' from 0 to TotalSteps:

    2.1. Kinematic Scene Update [Function: updateSceneMovement]
        - Calculation of time progress ratio (0.0 -> 1.0).
        - Mesh Manipulation:
            - Forward Kinematics (FK) based on Scene ID.
            	-> Update of 'PositionGlobal' and 'OrientationGlobal' for all meshes.
        		-> Saving for visualization.
        - Discretization [Function: m->discretizeMesh]:
            - Generation of surface point clouds for distance heuristics.

    2.2. Muscle Anchor & Guess Prediction
        - Anchor Update:
            - Recalculation of Global Origin/Insertion points based on parent Mesh movements.
        - Initial Guessing:
            - Prediction of positions for internal nodes (Linear interpolation or previous step velocity).
            #- Heuristic Application (Optional): "Snapping" of nodes near Via-Points to exact mesh locations.

    2.3. Non-Linear Optimization [Function: sys->solveStepX]
        - Input Assembly:
            - Flattening of Initial Guess (x0) and Parameters (p).
            - Setting of Bounds for Variables (lbx, ubx) and Constraints (lbg, ubg).
        - IPOPT Execution:
            - NLP Solving: Finding optimal Node positions satisfying constraints (while minimizing length.)
        - Result Extraction:
            - Update of 'MNodes' with optimized Global Positions.
            - Storage of Contact Forces (eta) (e.g. for warm-starting the next step).

    2.4. Post-Optimization Logic
        - Parent Association [Function: updateMusclePointsParents]:
            - Nearest Neighbor Search: Calculate distance to all meshes for each node.
            - State Update: Assign parent mesh and color based on proximity.
        - Validation:
            - Collision Check [Function: checkCollision]: Verification of node penetration (Logging warnings).
            - Via-Point Check [Function: getViaPointNodeInfo]: Logging of distance deviations.

    2.5. Data Archiving


3. VISUALIZATION & OUTPUT
-------------------------------------------------------
    3.1. Data Aggregation
        - Retrieval of full history vectors (Global Points, Guesses, Colors) from all Muscle instances.
        - Compilation of Mesh trajectory history.

    3.2. Performance Logging
        - Output of total simulation duration.
        - Printing of final Muscle Length profiles to console.

    3.3. Rendering [Class: VTKSimViewerSimple]
        - Initialization of VTK pipeline.
        - Loading of aggregated history data.
        - Execution: Launching of Qt Event Loop.
        - Display: Interactive 3D playback with slider control.
        - Export (Optional): Rendering frames to image files.





________________________________________________________________________
1. iter (Iteration):
   Die Nummer des Rechenschritts.
   ACHTUNG: Das "r" hinter der Zahl (z.B. 14r bis 47r) steht für "Restoration Phase".
   Bedeutung: Der Solver hat sich "verrannt". Die normalen Schritte haben die
   Nebenbedingungen so stark verletzt, dass er in einen Notfallmodus schaltet.
   Er versucht nun krampfhaft, nur noch eine physikalisch gültige Lösung zu
   finden, ignoriert dabei aber weitgehend das eigentliche Minimierungsziel.
   Dass er bis zum Ende im "r"-Modus bleibt, ist ein schlechtes Zeichen.

2. objective (Zielfunktion):
   Wert: 1.0000000e+00
   Bedeutung: Dieser Wert ändert sich über 47 Schritte nicht.
   Das heißt: Der Solver ist so sehr damit beschäftigt, die geometrischen
   Zwänge zu erfüllen, dass er die eigentliche Zielfunktion (z.B. Muskel-
   länge oder Energie) überhaupt nicht verbessern kann.

3. inf_pr (Primal Infeasibility) -> WICHTIGSTER WERT!
   Bedeutung: Dies ist der "Fehler" in den Nebenbedingungen.
   Wie weit werden die Regeln (Abstand > 0, Radius < Tol) verletzt?
   - Start: 2.04e+01 (Sehr großer Fehler am Anfang).
   - Ende:  2.87e-02 (ca. 0.028).
   Problem: Der Wert muss nahe 0 sein (z.B. 1e-08). Er bleibt aber bei
   0.028 stecken. Das bedeutet, irgendeine Bedingung wird permanent um
   ca. 2.8 cm (wenn Einheiten Meter sind) verletzt.

4. lg(mu):
   Der Logarithmus des Barriere-Parameters (für Ungleichungen).
   Er sollte gegen negativ Unendlich gehen. -5.2 ist okay, aber der Abbruch
   erfolgt zu früh.

5. ||d|| (Schrittweite):
   Wie stark ändern sich die Variablen?
   Am Ende (8.86e-05) fast gar nicht mehr. Der Solver "steckt fest".

6. "14r":
	Der Solver hat sich "verrannt". Die reguläre Newton-Methode hat einen 
	Schritt gemacht, der die Bedingungen so stark verletzt hat, dass er 
	nicht mehr weiterkommt. Er schaltet in den "Notfallmodus" 
	(Restoration Phase), um irgendeinen Punkt zu finden, der die Constraints 
	besser erfüllt, selbst wenn das Zielfunktions-Minimum ignoriert wird