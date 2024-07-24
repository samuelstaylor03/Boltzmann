**Description:** Generates a velocity.inp file for Varga Group's TDDFT code. Uses the Boltzmann distribution to generate velocities for individual atoms given an initial temperature and seed.

**Program Input:** The program expects three different inputs. The latter two can be given as command line arguments or runtime input after running the program with no arguments:
        
        dft.inp file: Must be in the same directory that you are in when you run the program. Needed in order to access the number of atoms, species, and positions of atoms in a molecule. 
        
        Temperature (in Kelvin): Provide as a command line argument or as runtime input. Any real number

        Seed: Provide as command line argument or as runtime input. Initial value used to generate random numbers for the Boltzmann distribution.

**Running the Program:** The program must be compiled in order to create an executable. There are two methods of doing this: either via the Makefile or manually compiling and then linking each file. See the .txt files included in the 'instruction' directory.

**Acknowledgements:** Kalman Varga, Cody Covington
