# ceaToFoam
Python script that joins the thermochemical data outputs from NASA Chemical Equilibrium with Applications and OpenFOAM. Automatically generates an engine geometry based on user specified performance parameters and surrounding firing conditions.

To download the inital files need to run the ceaToFoam program download the following Python files:

          chart_cea.py
          CEAtoFOAM.py
          foamCaseGeneration.py 



The three files listed above are responsible for parsing data from the CEAWrap package developed by civilwargeeky. The ceaToFoam.py file is required to run foamCaseGeneration.py since it calls directly from CEAWrap to gather data from the rocket module inputs.

To download CEAWrap, use the terminal and enter:
          
          pip install CEA-Wrap==1.7.4


Also verify that the following packages are installed in the latest version of Python:

          Matplotlib
          Numpy
          Cantera

Install the latest version of each of these packages:

          pip install matplotlib

          pip install numpy

          py -m pip install cantera



To run the ceaToFoam program:

          1. Download the repository and extract the folder "ceaToFoam".

          2. Ensure all previous packages have been downloaded.

          3. Run the program in an IDE like VSCode and follow the prompts in the console.

The only inputs required for this program are:
          1. Chamber Pressure (Pc)
          2. O/F Ratio
          3. Propellant Combination (Fuel and Oxidizer) 
          4. Altitude operating conditions (nozzle expansion depends on this)
          5. Design Thrust

Other values that are needed for the geometry of the engine contour itself are up to the user and recommended in the console based on the information available in Modern Engineering for the Design of Liquid Propellant Rocket Engines by Huzel and Huang.

Final output of the program will be a folder with a name of the user's choice that can be placed into OpenFOAM and compiled in the Linux environment.
