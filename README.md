# Boltzmann
**Description:** Generates a velocity.inp file for Varga Group's TDDFT code. Uses the Boltzmann distribution to generate velocities for individual atoms given an initial temperature and seed.

**Program Input:** The program expects 4 different inputs. The former two must be in the form of input files. The latter two can be given as command line arguments or runtime input after running the program with no arguments:
        
        control_boltzmann.inp file: Must be in your current directory when you run the program. Input parameters for the program such as where to access the dft.inp file, where to output the velocity.inp file, and if an extra line for a proton's velocity should be printed.

        dft.inp file: Input the path to this file in the control_boltzmann.inp file. The path should be relative to your current directory.. Needed in order to access the number of atoms, species, and positions of atoms in a molecule. 
        
        Temperature (in Kelvin): Provide as a command line argument or as runtime input. Any real number.

        Seed: Provide as command line argument or as runtime input. Initial value used to generate random numbers for the Boltzmann distribution. Any Integer.

**Running the Program:** The program must be compiled in order to create an executable. There are two methods of doing this: either via the Makefile or manually compiling and then linking each file. See the .txt files included in the 'instruction' directory. Run with the following command:

        ./<program_executable> <temperature_in_K> <seed>

**Single File:** The directory single_file contains one source file that can be compiled individually and performs the same calculation as the program. It is simply all the modules moved to one singular main file. It is much less organized, but performs the same calculation and is easier to compile (see single_file/compile_instructions.txt).

**Acknowledgements:** Kalman Varga, Cody Covington
