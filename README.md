# Variable Delay Line Algorithm for Simulating the Doppler effect created by a moving object
Descreiption of the Project
--------------------------------------------------------------------------------------------
This project aims to simulate a variable delay line algorithm and demonstrate its use in simulating the Doppler effect caused by a moving sound source. The moving object is modeled as following a circular path, and the resulting sound effect is simulated at the center of the circle.

The mathematical details of the delay line algorithm can be found on "Damiano, S., Bondi, L., Guntoro, A. et al. A framework for the acoustic simulation of passing vehicles using variable length delay lines. J AUDIO SPEECH MUSIC PROC. 2024, 49 (2024)".

Two microphones placed at the center of the circle receive the sound signal, and the Time Difference of Arrival (TDOA) between the signals is estimated. This TDOA is used to align the received signals at the two microphones. The estimation of the TDOA is performed using a correlation-based method.

For a more detailed explanation of the estimation methods used, please refer to:
"Spagnolini, Umberto. Statistical Signal Processing in Engineering. John Wiley & Sons, 2018".
