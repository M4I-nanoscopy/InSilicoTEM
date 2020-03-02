# InSilicoTEM_2020

Simulator of a Transmission Electron Microscope (TEM) optimized for imaging of biological specimen.

Implemented in MATLAB using the DIPimage Toolbox for MATLAB.

Adapted by Yue Zhang and Ruben Tammaro from the code of Vulovic et al., 2013 (https://doi.org/10.1016/j.jsb.2013.05.008).

Changes to the original code include:
  -

1) Addition of an ideal phase plate to the microscope's optics
    - Phase plate adds phase shift to the electron wave which provides higher phase contrast when the scattered waves recombine on the image plane
    
2) Addition of a Point of Entry script: 'RunTEMsim.m'
    - Allows for the iterative generation of multiple micrographs in a single run with either varying or fixed parameters
    - Provides output to command line on the status of the simulation
    
3) Optimization of memory usage for parallel processing of multiple micrographs
    - Producing a single 4K micrograph would use up to 210 GB of RAM before optimization
    - After optimization it only uses up to 60 GB
    - Up to 3 simulations at the same time can therefore now be run on a single machine with 256 GB RAM without outofmemory errors

How to use InSilicoTEM simulator: testing updates
-

!! BEFORE USING THE SIMULATOR, MAKE SURE THAT A WORKING VERSION OF DIPimage (http://www.diplib.org/) FOR MATLAB IS INSTALLED ON THE COMPUTER !! 

1) Add the PDB file to be imaged to the '/PDBs' folder
2) Start MATLAB and open 'RunTEMsim.m'
3) Here some simulation parameters can be set:
    - Number of micrographs to generate
    - Phase plate usage
    - Defocus range
    - Correction factor [multiple values can be specified so that for each micrograph nr. as many micrographs as correction factor values will be produced]
    - Motion blur (radiation damage) [same as for correction factor]
    - Size of the micrographs (pixels)
    - Pixel size
    - Minimum particle distance
    - Directory where to save the micrographs
4) Other parameters (such as acceleration voltage, specimen thickness, detector type, etc.) can be set within the 'src/TEMsim.m' function
5) Run 'RunTEMsim.m' with the chosen parameters
6) Find the generated micrographs in the '/Micrographs' folder

How the InSilicoTEM simulator works:
-
Here only a brief description of the simulator is presented.

For a the full theoretical treatment please refer to OUR PAPER (when we have the DOI we will update it) as well as Vulovic et al.,2013 (https://doi.org/10.1016/j.jsb.2013.05.008).

The functioning of InSilicoTEM can be broken down into 4 steps:

1) Interaction Potential
    - First, a specimen (PDB file) is loaded into matlab
    - The positions of every atom are recorded
    - The electrostatic potential of the specimen is calculated
    - A volumetric potential map of the specimen is constructed
    
2) Electron Wave Propagation and Specimen Interaction (Multislice algorithm)
    - The 3D potential map of the specimen is divided in a number of slices of fixed thickness
    - An incident planar electron wave is generated 
    - The wave is propagated through the specimen one slice at a time
    - On every slice, the electron wave is diffracted according to the specimen's interaction potential
    - This goes on until all slices have been done and an exit wave is formed
    
3) TEM Optics
    - The exit wave from (2) is passed through the microscopes contrast transfer function (CTF)
    - The CTF describes how the optics of the TEM process the incoming electron wave
      - Here is where the shift in phase provided by the Phase Plate is introduced, if desired
    - After these steps we obtain the electron intensity on the image plane, where the detector is
    
4) Detector Response
    - The signal obtained in (3) is captured by the detector
    - Detection of the whole signal in real microscopes is not possible, therefore to produce a realistic micrograph the detector is modelled according to:
      - Conversion factor
      - Detective Quantum Efficiency (DQE)
      - Modulation Transfer Function (MTF)
      - Noise Transfer Function (NTF)
      - Poisson noise
      

