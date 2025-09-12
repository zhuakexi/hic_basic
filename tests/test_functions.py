# tests/test_functions.py
import pandas as pd
import os
from pathlib import Path

def simple_copy(input_path, output_path):
    """
    Simple function that copies content from input to output file.
    
    Parameters:
        input_path (str): Path to input file
        output_path (str): Path to output file
        
    Returns:
        str: Success message
    """
    # Read input file
    with open(input_path, 'r') as f:
        content = f.read()
    
    # Write to output file
    with open(output_path, 'w') as f:
        f.write(content)
    
    return f"Copied {input_path} to {output_path}"

def process_multiple_inputs(input1, input2, output_path):
    """
    Function that processes two inputs and creates an output.
    
    Parameters:
        input1 (str): First input value
        input2 (str): Second input value
        output_path (str): Path to output file
        
    Returns:
        str: Success message with inputs concatenated
    """
    result = f"{input1}_{input2}"
    
    with open(output_path, 'w') as f:
        f.write(result)
    
    return result

def process_multiple_outputs(input_path, output1_path, output2_path):
    """
    Function that processes an input and creates two outputs.
    
    Parameters:
        input_path (str): Path to input file
        output1_path (str): Path to first output file
        output2_path (str): Path to second output file
        
    Returns:
        tuple: Paths of the two output files
    """
    with open(input_path, 'r') as f:
        content = f.read()
    
    mid_point = len(content) // 2
    part1 = content[:mid_point]
    part2 = content[mid_point:]
    
    with open(output1_path, 'w') as f:
        f.write(part1)
    
    with open(output2_path, 'w') as f:
        f.write(part2)
    
    return (output1_path, output2_path)

def return_series(input_path, output_path):
    """
    Function that returns a pandas Series.
    
    Parameters:
        input_path (str): Path to input file
        output_path (str): Path to output file
        
    Returns:
        pd.Series: A series with file information
    """
    # Read input file
    with open(input_path, 'r') as f:
        content = f.read()
    
    # Write to output file
    with open(output_path, 'w') as f:
        f.write(content)
    
    # Return a Series
    return pd.Series({
        'input_file': input_path,
        'output_file': output_path,
        'content_length': len(content)
    })

def failing_function(input_path, output_path):
    """
    Function that intentionally fails for testing error handling.
    
    Parameters:
        input_path (str): Path to input file
        output_path (str): Path to output file
        
    Raises:
        ValueError: Always raises an error for testing
    """
    raise ValueError("This function always fails for testing purposes")