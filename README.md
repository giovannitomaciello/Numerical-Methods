

<div align="center">
    <h1>Numerical Methods for Molecular Simulation</h1>
    
![MATLAB](https://img.shields.io/badge/MATLAB-e86e05?style=for-the-badge&logo=Octave&logoColor=white)
![POLIMI](https://img.shields.io/badge/POLIMI-050065?style=for-the-badge)
    
</div> 

<h4 align="right">Andrea Somma, Giovanni Tomaciello, Paolo Zagardo</h4> 

## Group Assignment 1 
<div>
    <h3>Mathematical simulators of a double elastic pendulum and a frozen Argon crystal in MATLAB</h3> <h5 align="right">March 2024</h5>
</div>


### Implemented numerical methods:
-  [Forward Euler method](./@INT/euleroavanti.m)
    - Explicit,
    - 1st order,
    - Instable solution - Amplification,
    - Not reversible.
-  [Backward Euler method](./@INT/euleroindietro.m)
    - Implicit,
    - 1st order,
    - Stable solution - Damping,
    - Not reversible.
-  [Crank–Nicolson method](./@INT/crankNick.m)
    - Implicit,
    - 2nd order in _t_, permitting long timesteps,
    - Stable solution, but out of phase,
    - Reversible.
-  [Velocity Verlet method](./@INT/velVerlet.m)
    - Explicit,
    - 2nd order in _t_, permitting long timesteps,
    - Reversible,
    - Symplectic.

### Double elastic pendulum, [Matlab code](./pendolo.m)
<div align ="center">
  
  <img src="https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/8d5d75d5-94d7-4aba-9df0-2483e19cfc79" width="700" height="400">
</div>

#### Features:
- two or more masses;
- forces: gravity, elastic force;
- no friction contribution;
- relies on the above numerical methods.

#### Graphs (velocity Verlet method)
- Energy conservation
<div align ="center">
  
<img src="https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/aec5d3d8-2fb7-435f-a9c3-d831f43b617a" width="500" height="400">
</div>

- Pendulum motion from 0.00 to 1.00
<div align ="center">
  
<img src="https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/ef6c880b-a0ba-4f0d-9665-f291e62ecd0e" width="1200" height="500">
</div>





### Frozen Argon Crystal, [Matlab code](./argonCry.m)

<div align ="center">
  
  ![image](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/bc101e96-b38f-4e75-af52-911f671965af)
</div>

This code models the interaction of seven argon atoms in a plane, where six of them are arranged symmetrically around a centre atom.
The chosen mathematical model is the following Hamiltonian:
<div align ="center">
  
  ![image](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/de418490-2fa3-4a15-aca5-c1c941ec0bf5)
</div>

Where V<sub>i,j</sub>(r) represents the Lennard–Jones potential:
<div align ="center">
  
  ![image](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/7693fc49-a3dc-4b0c-a5b8-861af47f32f2)
</div>

Additionally, temperature was computed using the following formula:
<div align ="center">
  
  ![image](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/e43fa3c6-c8b4-45e9-af35-62d05631a1d8)
</div>

#### Results
- Velocity Verlet:
  <div align ="center">
  
  ![argon E e T vv](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/aba258e9-311a-4281-98c7-d832459a354f)
  </div>
- Eulero Forward:
  <div align ="center">
  
  ![argon E e T fe](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/ff39f45b-cd7f-4bd7-b412-111f4ad85bb0)
  </div>
- Eulero Backward:
  <div align ="center">
    
  ![argon E e T be](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/6bb6d436-ad0d-413b-becf-163793a6277b)
  </div>
- Crank-Nicolson:
  <div align ="center">
  
  ![argon E e T cn](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/03eeed65-335c-4702-b670-9cf1606f90f4)
  </div>
