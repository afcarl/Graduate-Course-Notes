# Graduate-Course-Notes
I scanned a few of my graduate course notes into pdf files.

<h1> Mathematical Fluid Mechanics </h1>
<ul>
<li> What is a Fluid?
<li> Geometry, R^3, Cartesian Tensors
<ul>
<li> 2nd Order Tensors and nth Order Tensors
<li> Quotient Rule, Contraction, Determinants, Isotropic
<li> Tensor Calculus: Gradient, Divergence, Curl
<li> Integrals: Stoke's Theorem, Green's Theorem
</ul>
<li> Kinematics: Velocity, Accelaration
<li> Density
<li> Material Derivative
<li> Particle Paths, Streamlines, Streaklines
<li> Reynolds Transport Theorem
<li> Conservation of Mass
<li> Conservation of Momentum, Euler's Equation, Hydrostatic Equation
<li> Principle of Local Stress Equilibrium
<li> Motion = Translation + Deformation + Rotation
<li> Cauchy Representation Theorem
<li> Cayley-Hamilton Theorem
<li> Compressible Navier-Stokes Equation
<li> Scaling of the Navier-Stokes Equation
<li> Derivation of the Boundary Layer Equation from Navier Stokes
<li> Stream Function in 2D flow
<li> Vorticity Equation, Vortex Stretching
<li> Bernouilli's Equation
<li> Velocity Potential
<li> Complex Variables To Solve 2D Steady Irrotational Inviscid Flows
<li> Circulation Results in Lift
<li> van Karman Vortex
<li> Very Viscous Flows - Stokes Flow
</ul>


<h1> Methods of Applied Mathematics </h1>
<ul>
<li> Eigenvalues and Eigenvectors of Continuous Systems
<ul>
<li> Head Conduction (uniform temperature)
<li> Shallow Water Waves
<li> Rotating FLow
</ul>
<li> Differential Operators
<li> Regular Sturm-Liouville Problems (Separated Boundary Conditions)
<li> Heat Conduction in a Nonuniform Medium
<li> Special Functions
<ul>
<li> Eigenfunction Expansions
<li> Chebyshev's Equation
<li> Singular Points of Differential Equations
<li> Euler-Cauchy Equation
<li> Weierstrauss Convergence Criteria
<li> Hypergeometric Series
<li> Confluent Hypergeometric Eqn (Kummer's Eqn)
<li> Back to Chebyshev's Equation
<li> Bessel Functions, Gamma Function, Applications to Membranes
<li> Modified Bessel Functions
</ul
<li> Green's Functions
<ul>
<li> Motivation: Solving Nonhomogeneous Differential Equations with Boundary Conditions
<li> Variation of Parameters to Derive Green's Functions
<li> Various Examples, Including Non-Separated Boundary Conditions
<li> Modified Green's Functions
<li> Differential Operators, Adjoint
</ul>
<li> Calculus of Varitions
<ul>
<li> Motivating Examples Including Catenary
<li> First Variation of Functional, Euler-Lagrange Equation and Solution
<li> Various Examples, Including Minimal Surface of Revolution
<li> Natural Boundary Conditions
<li> Transition Conditions
<li> Functions of More Than One Variable
<li> Find the Minimum of the Functional
<li> Calculus of Variations With Constraints
<li> Vibrational Problem of Buckling
<li> Variable Endpoints: Transversality Condition
<li> Finite Constraints
<li> Differential Equations as Constraints
</ul>
<li> Rayleigh-Ritz Method
<ul>
<li> Galerkin Method
<li> Approximate Methods Using Gram-Schmidt Orthogonolization of Basis Functions
</ul>
<li> Hamiltonian Systems
</ul>


<h1>Numerical Solution of Ordinary Differential Equations</h1>
<ul>
<li> Initial Value Problems
<ul>
<li> Pendulum
<li> Predator-Prey
<li> Biochemical Kinetics
<li> Diffusion Problem
</ul>
<li> Regularity Result
<li> Boundary Value Problems
<ul>
<li> A Diffusion Problem
<li> Singular Perturbation Problem
<li> Blasius Problem (Boundary Layer)
</ul>
<li> Differential Algebraic Problems
<ul>
<li> Mechanics
<li> Hamilton's Principle
<li> Pendulum example
</ul>
<li> Initial Value Problems
<ul>
<li> "Stability"
<li> Test Equation, Full Equation
<li> Variable Coefficients, Non-Homogeneous
<li> Nonlinear Case
<li> Hamiltonian Systems
</ul>
<li> Numerical Methods for IVPs
<ul>
<li> Euler's Method
<li> Local Truncation Error, Consistent, Convergence, O-Stability
<li> Local Error
<li> Absolute Stability
<li> Spring Equation: Overdamped and Underdamped
</ul>
<li> Stiffness and Implicit Methods
<ul>
<li> Backward Euler
<li> Newton's Method
<li> Trapezoidal Method
<li> Single-Step Methods
</ul>
<li> Runge-Kutta Methods
<ul>
<li> Derivation
<li> General S-Stage
<li> Order of Accuracy
<li> Special Cases and Classical RK4 Method
<li> Region of Absolute Stability and Error Control
<li> Implicity Runge-Kutta Methods
<li> Diagonally Implicit
</ul>
<li> Linear Multistep Methods
<ul>
<li> Adams Family
<li> Backwards Difference Formulas
<li> Order of Accuracy
<li> Root Condition, O-Condition, Absolute Stability
<li> Implementation
<ul>
<li> Predictor-Corrector Methods
<li> Error Estimates
<li> Variable Step Size
</ul>
</ul>
<li> Boundary Value Problems
<ul>
<li> Green's Functions
<li> Problem Stability For Linear BVPs
<li> Stiff BVPs
<li> Shooting Methods for BVPs
<li> Multiple Shooting
</ul>
<li> Finite Difference Methods for BVPs
<ul>
<li> Midpoint Methods
<li> Nonlinear BVPs
<li> Consistency, O-Stability, Convergence
<li> Higher Order Methods: Collocation Methods, Richardson Extrapolation, Continuation Methods
</ul>
</ul>



<h1> Numerical Solution of Partial Differential Equations </h1>
<ul>
<li> Classification of PDEs
<ul>
<li> hyperbolic, parabolic, elliptic
<li> linar/nonlinear
<li> classification via characteristics
<li> first order systems
</ul>
<li> Finite Difference Methods for PDEs
<ul>
<li> Simple example involving heat equation: introduce the mesh; Taylor series to approximate derivatives
<li> Method of lines
<li> Solve using forward difference
<li> What are some questions one might ask: accuracy, stability, cost
</ul>
<li> Consistency, Order of Accuracy
<ul>
<li> Difference Operators
<li> Error grid function
<li> Local Truncation Error
<li> Definition of consistency
</ul>
<li> Computational Cost
<li> Neumann Problem for the Heat Equation
<ul>
<li> Neumann Boundary Conditions
<li> Introduce "ghost lines"
<li> Integral Conservation for Neumann Problem
</ul>
<li> Generating Discrete Approximations
<ul>
<li> Taylor Series Approach
<li> Interpolation Approach
<li> Finite Volume Approach
</ul>
<li> Convergence, Consistency and Stability
<ul>
<li> Lax Theorem
<li> Fourier Stability Analysis
<li> Fourier Mode Analysis
<li> Stability Analysis for Initial BVP
</ul>
<li> Finite Difference Methods for Parabolic Methods
<ul>
<li> Crank Nicholson Method
<li> Alternating Direction Implicit (ADI)
<li> Heat Equation in Multi-Spatial Dimensions including Polar Coordinates
<li> Nonlinear Heat Equation
</ul>
<li> Hyperbolic Partial Differential Equations
<ul>
<li> Linear Advection Equation
<li> Behavior of Discontinuities for Linear Equations
<li> Lax-Friedrich, Lax-Wendroff
<li> Courant-Friedrich-Lewy (CFL) Condition
<li> Non-constant coefficients
<li> Linear systems: Upwind methods, Boundary Conditions
</ul>
<li> Hyperbolic Conservation Laws
<ul>
<li> Scalar Conservation Laws
<li> Zero limit of "viscous" solution
<li> Weak Solutions of the Integral Form
<li> Jump Notation
<li> Burger's Equation
<li> Riemann Problems
</ul>
<li> Finite Volume Formulation of Conservation Laws
<ul>
<li> Conservative Finite Volume Scheme
<li> Lax-Wendroff Theorem
<li> Godunov Methods
<li> High Resolution Methods - Flux Limiters
<li> Total Variation Diminishing
<li> High Resolution Methods - Slope Limiters
</ul>
<li> Systems of Conservation Laws
<ul>
<li> Finite Volume Method
<li> Godunov Method
<li> Multiple Spatial Dimensions
<li> Directional Splitting
</ul>
<li> Elliptic Equations
<ul>
<li> Properties of Laplace's Equation
<li> The Numerical Problem: Solvability of the Linear System
<li> Convergence of Finite Difference Approximation
<li> Direct Methods: Direct Factorization, Block-Tridiagonal Solvers
<li> Iterative Methods: Jacobi, Gauss-Seiel, Successive Over Relaxation (SOR)
<li> Analysis of Residual Correction Schemes
<li> Multigrid Methods
<li> Two Grid Algorithm
</ul>
</ul>



<h1> Partial Differential Equations II: Special Topics </h1>
<ul>
<li> Wave Equation
<ul>
<li> Linear Elasticity
<ul>
<li> Isotropic Case
<li> Constant Coefficient and Isotropic Case
</ul>
<li> Incompressible Case
<li> Electromagnetic Eqns (Constant Coefficient Case)
<li> Weak Solutions
<ul>
<li> Derivation via Multiplication by Test Function and Integration
<li> How to deal with Initial Conditions
<li> Use Energy Estimates to Prove Convergence
<li> Prove Existence of Weak Solution
<li> Definition of Weak Convergence
</ul>
<li> Propagation of Disturbaces (Finite Propagation)
<li> Lax-Milgram Theorem
<li> Helmholtz Equation
<li> Complex Lax-Milgram Theorem
</ul>
<li> Nonlinear Parabolic Equations
<ul>
<li> Derivation of Eikonal Equation
<li> Definition of Weak Solution
<li> Definition of Viscosity Solution
</ul>
<li> Semigroup Theory
<ul>
<li> Definition and Elementary Properties of Semigroup
<li> Differentiatial Properties of Semigroups
<li> Definition of Infinitesimal Generators
<li> Hill-Yoshida Theorem
</ul>
</ul>



<h1> Perturbation Methods </h1>
<ul>
<li> Examples: Nonlinear Oscillators (Quadratic Equations)
<ul>
<li> Comparison of Exact Solution and Expansion for simple case
<li> Example using Gauge (scaling) functions
<li> Example with singularity
</ul>
<li> Background Theory
<ul>
<li> Order Notation - "Big Oh", "Little Oh"
<li> Gauge functions (scale functions, basis functions)
<li> Transcendentally small terms
</ul>
<li> Asymptotic Expansion of Functions
<ul>
<li> Uniform Asymptotic Expansions
<li> Differentiation
<li> Convergent Series vs Asymptotic Series
</ul>
<li> Regular ODEs
<ul>
<li> Projectile Motion
<li> Interior Boundary Layer
<li> Boundary Layer on Both Sides
<li> Lubrication Theory / Slider Bearing
<li> Mass-Spring Damper System (using 3-2 Van Dyke Matching)
<li> Modified Van Dyke Principle
<li> Far Field / Switchback (using Van Dyke Matching)
<li> Classical Problem From Kevorkian and Cole (BVP)
</ul>
<li> Weakly Nonlinear Oscillators
<ul>
<li> Duffing Equation, Linear Spring With Damping, van Der Pol Eqn, Rayleigh's Equation
<li> The Naive Expansion Fails
<li> Naive Expansion - Rescaled
<li> Three approaches: Renormalization, Strained Coordinates, Multiple Scales
<li> Duffing Equation Using Poincarre-Lighthill Method
</ul>
<li> Multiple Scales
<ul>
<li> Introduce Two Scales - Calculate Derivatives
<li> Example Involving Alternative Forms Of The Homogeneous Solution
<li> Various Forms of Slow/Fast Multiple Scales
<li> Multiple Scales and Boundary Layer Problems
<il> General Weakly Nonlinear Oscillator
</ul>
<li> Phase Plane / Limit Cycles
<ul>
<li> Rayleigh's Equation using Multiple Scales Approach
<li> Rayleigh's Equation using WKB(J) Approximation
</ul>
<li> Asymptotic Integration
<ul>
<li> Methods: Small/Large Parameters, Stationary Phase, Laplace Method
<li> Integration by Parts
<li> Laplace Method
<li> Special Cases
<li> Fourier Integral / Method of Stationary Phase
<li> Bessel Function of Order n
</ul>
</ul>














