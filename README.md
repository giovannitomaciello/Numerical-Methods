

<div align="center">
    <h1>Numerical Methods for Molecular Simulation</h1>
    
![MATLAB](https://img.shields.io/badge/MATLAB-e86e05?style=for-the-badge&logo=Octave&logoColor=white)
![POLIMI](https://img.shields.io/badge/POLIMI-050065?style=for-the-badge)
    
</div> 

<h4 align="right">Andrea Somma, Giovanni Tomaciello, Paolo Zagardo</h4> 

<!--- ## Group Assignment 1 
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
#
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
  
<img src="https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/d57d93fb-ddba-4bd2-9a1d-cca091ed73b3" width="500" height="400">

<img src="https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/c02f6b72-4fa7-4b92-b54a-72385c8769c6" width="500" height="400">

</div>

#
- Pendulum motion from 0.00 to 10s
<div align ="center">

https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/66b9a125-7caa-47aa-b47a-42f81def5554
  
https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/ed220e49-8909-4c92-9d47-1402c7d5abe9

</div>




#
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

<div align ="center">

https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/9a921526-6490-465e-83b5-1ac2ed745873

</div>

- Velocity Verlet:
  <div align ="center">
  
  ![vv40](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/741e25e3-a56d-450e-a0ca-6e139587749c)
 ![vv80](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/77fc90c8-72e7-4422-8ce5-b825b0453fe4)

  </div>
- Eulero Forward:
  <div align ="center">
  
  ![feul10](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/debfde79-0479-4572-a912-bc6162cacc27)
 ![feul20](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/b8e2056f-1f63-40e0-bd10-a3476ffd9a07)

  </div>
- Eulero Backward:
  <div align ="center">
  
    ![beul20](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/1045ee55-6b6a-44f0-9e3e-31948211db65)

  </div>
- Crank-Nicolson:
  <div align ="center">
  
    ![nick20](https://github.com/giovannitomaciello/Numerical-Methods/assets/120776791/1a4d9169-1772-441c-9f56-90cc630d9f28)

  </div>

  --->
