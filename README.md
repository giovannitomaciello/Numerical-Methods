# Numerical-Methods

<div align="center">

![MATLAB](https://img.shields.io/badge/MATLAB-e86e05?style=for-the-badge&logo=Octave&logoColor=white)
![POLIMI](https://img.shields.io/badge/POLIMI-050065?style=for-the-badge)
</div>

<!-- PROJECT LOGO -->
<br />
<div align="center">  
   Trova una immagine 
      <br />
    <h1 align="center">Group Assignment 1</h1>
    <h3 align="center">Mathematical simulators of a double elastic pendulum and an Argon crystal in MATLAB</h3>
    <h4>Team: Andrea Somma, Giovanni Tomaciello, Paolo Zagardo</h4>
</div>

## Implemented numerical methods:
- Crank–Nicolson method [Matlab code](./@INT/crankNick.m)
  - Piccola descrizione vantaggi e svantaggi 
- Velocity Verlet method [Matlab code](./@INT/velVerlet.m)
  - Piccola descrizione vantaggi e svantaggi 
- Forward Euler method [Matlab code](./@INT/euleroavanti.m)
  - Piccola descrizione vantaggi e svantaggi 
- Backward Euler method [Matlab code](./@INT/euleroindietro.m)
  - Piccola descrizione vantaggi e svantaggi 

## Double elastic pendulum 
[Matlab code](./pendolo.m)
<div align ="center">

  ![image](https://github.com/giovannitomaciello/Numerical-Methods/assets/162450790/961cf984-801a-4cd4-93c3-0330fdb019f9)
</div>
  
- two or more masses;
- gravity effect;
- no friction contribution;
- relies on the above numerical methods.

### Results
Grafici e immagini

## Frozen Argon Crystal

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

### Results
Grafici e immagini