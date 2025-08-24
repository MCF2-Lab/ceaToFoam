import numpy as np


def extract_rows_in_range(input_file, output_file, start_row, end_row):
    #if you want to run a new propellant combination, you must change the inputs of within the CEA_Code.py file to
    #the specifications of your engine, this file will handle the output from my_output.out and parse it into a format
    #which is applicable to the OpenFOAM simulation
    """
    Extracts rows within a specified range from the input file and writes them to the output file.

    :param input_file: Path to the input text file.
    :param output_file: Path to the output text file.
    :param start_row: Starting row number (1-based index).
    :param end_row: Ending row number (1-based index).
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
    
    # Extract the specified range of rows
    extracted_lines = lines[start_row - 1:end_row]
    
    with open(output_file, 'w') as outfile:
        outfile.writelines(extracted_lines)

import numpy as np

def text_file_to_array(input_file):
    """
    this function converts the text file into a 2D array
    Removes asterisks from the text column.

    :param input_file: Path to the input text file.
    :return: A 2D numpy array with 9 rows and 4 columns.
    """
    array = []
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.split()
            # Remove asterisks from the first column
            text = parts[0].replace('*', '')
            numbers = list(map(float, parts[1:]))
            array.append([text] + numbers)
    
    return np.array(array)

# Example usage
if __name__ == "__main__":
    input_file = "input.txt"  # Replace with your input file path
    output_file = "output.txt"  # Replace with your desired output file path
    start_row = 187  # Replace with your desired start row
    end_row = 195  # Replace with your desired end row

    extract_rows_in_range('my_output.out', 'myoutput2.txt', start_row, end_row)

    input_file = "myoutput2.txt"  # Replace with your input file path

    result_array = text_file_to_array(input_file)
    

species_from_cea = result_array[:, 0] #outputs an array for the species names which make up the largest propellant compositions in the reaction

#ALL OF THE FOLLOWING ARE THE MOLE FRACTIONS OF EACH SPECIES AT THE DIFFERENT NODES OF THE ENGINE AS OUTPUTTED BY CEA
chamber_mole_fractions = result_array[:, 1].astype(float) #mole fractions of species in the combustion chamber

throat_mole_fractions = result_array[:, 2].astype(float) #mole fractions of species in the throat section

nozzle_exit_mole_fractions = result_array[:, 3].astype(float) #mole fractions of species in the nozzle section



print("Species from CEA:", species_from_cea)

molecular_weights = {
    'CO': 28.01, #g/mol 
    'CO2': 44.01, #g/mol 
    'COOH': 46.03, #g/mol
    'H': 1.008, #g/mol 
    'HCO': 28.01, #g/mol 
    'H2': 2.016, #g/mol
    'H2O': 18.015, #g/mol 
    'H2O2': 34.02, #g/mol
    'O' : 15.999, #g/mol 
    'OH': 17.01, #g/mol
    'O2': 32.00, #g/mol
}

print("Chamber Mole Fractions:", chamber_mole_fractions)
print("Throat Mole Fractions:", throat_mole_fractions)
print("Nozzle Exit Mole Fractions:", nozzle_exit_mole_fractions)


molecular_weight_array = np.array([molecular_weights [species] for species in species_from_cea])

print("Molecular Weights: ", molecular_weight_array)
