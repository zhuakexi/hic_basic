import os
import concurrent.futures
import importlib
import traceback

from functools import wraps
from tqdm import tqdm
import pandas as pd

def task_wrapper(module_name, func_name, *args, **kwargs):
    """
    Wrapper function to import and execute a function from a module.
    
    Parameters:
        module_name (str): Name of the module to import from
        func_name (str): Name of the function to execute
        *args: Positional arguments to pass to the function
        **kwargs: Keyword arguments to pass to the function
        
    Returns:
        Result of the executed function
    """
    module = importlib.import_module(module_name)
    func = getattr(module, func_name)
    #print(f"Executing {func_name} with args: {args}")
    return func(*args, **kwargs)

def mt(
    module_name: str
):
    """
    A decorator factory to create parallelized function decorators.

    This factory returns a decorator that parallelizes the function it decorates,
    providing features like file skipping, path templating, progress bars, and result handling.

    Args:
        module_name (str): The name of the module where the function to be decorated is located.

    Returns:
        A decorator function that can be applied to a function.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(input_data: pd.DataFrame, force=False, nproc=None, outcol=None, concat=False, *args, **kwargs):
            """
            Parallelized function wrapper with input/output pattern support.
            
            The function uses a pattern-based system for input and output file handling:
            - Input patterns: input_pattern, input_pattern1, input_pattern2, etc.
            - Output patterns: output_pattern, output_pattern1, output_pattern2, etc.
            
            Patterns are formatted using {sample_name} placeholder which is replaced with
            the DataFrame index value for each row.
            
            Parameters:
                input_data (pd.DataFrame): Input DataFrame containing data to process
                force (bool): Whether to overwrite existing output files (default: False)
                nproc (int): Number of parallel processes (default: all CPUs)
                outcol (str): Column name to store output file paths in the input DataFrame
                concat (bool): If True, concatenate Series results into a DataFrame
                *args: Additional positional arguments passed to the decorated function
                **kwargs: Additional keyword arguments including input/output patterns
                
            Pattern Parameters (passed as kwargs):
                input_pattern: Template for input file paths
                input_pattern1, input_pattern2, ...: Additional input patterns
                output_pattern: Template for output file paths  
                output_pattern1, output_pattern2, ...: Additional output patterns
                input_col: Column name containing input paths (default: "input_path")
                
            Returns:
                If concat=True and results are pandas Series: DataFrame with concatenated results
                Otherwise: Dictionary of results keyed by DataFrame index
            """
            njobs = nproc or 4

            # Extract main input/output patterns
            input_patterns = []
            output_patterns = []

            # Extract input_pattern, input_pattern1, input_pattern2...
            idx_num = 0
            while True:
                key = "input_pattern" if idx_num == 0 else f"input_pattern{idx_num}"
                if key in kwargs:
                    input_patterns.append(kwargs.pop(key))
                    idx_num += 1
                else:
                    break

            # Extract output_pattern, output_pattern1, output_pattern2...
            idx_num = 0
            while True:
                key = "output_pattern" if idx_num == 0 else f"output_pattern{idx_num}"
                # print(key, key in kwargs)
                if key in kwargs:
                    output_patterns.append(kwargs.pop(key))
                    idx_num += 1
                else:
                    break
            # print(kwargs)
            # input_col check — only one allowed
            input_col = kwargs.pop("input_col", "input_path")
            if isinstance(input_col, (list, tuple)):
                raise ValueError("Only one input_col is allowed.")

            # check output
            if not output_patterns and outcol is None:
                assert concat, "If no output_pattern and outcol are provided, concat must be True"

            # Store output paths and results
            output_paths = {}
            results = {}

            # Generate tasks with row indices
            tasks = []
            for idx, row in input_data.iterrows():
                sample_name = idx  # Assuming the index is the sample identifier

                # Resolve all input paths
                resolved_inputs = []
                if input_patterns:
                    for pat in input_patterns:
                        resolved_inputs.append(str(pat).format(sample_name=sample_name))
                else:
                    resolved_inputs.append(row[input_col])

                # Resolve all output paths
                resolved_outputs = []
                if output_patterns:
                    for pat in output_patterns:
                        resolved_outputs.append(str(pat).format(sample_name=sample_name))
                else:
                    # resolved_outputs.append(None)
                    pass

                # Skip check — if *any* output exists and force=False, skip
                if any(resolved_outputs) and not force:
                    if all(op is None or os.path.exists(op) for op in resolved_outputs):
                        continue

                tasks.append((idx, resolved_inputs, resolved_outputs))
            # create directories for all output paths
            for _, _, resolved_outputs in tasks:
                for op in resolved_outputs:
                    if op:
                        os.makedirs(os.path.dirname(op), exist_ok=True)

            # Execute tasks in parallel
            with concurrent.futures.ProcessPoolExecutor(max_workers=njobs) as executor:
                futures = {}
                for idx, resolved_inputs, resolved_outputs in tasks:
                    # Arguments order: inputs..., outputs..., *args, **kwargs
                    task_args = resolved_inputs + resolved_outputs + list(args)
                    futures[
                        executor.submit(
                            task_wrapper,
                            module_name,
                            func.__name__[3:],  # Remove 'mt_' prefix
                            *task_args,
                            **kwargs
                        )
                    ] = idx

                # Track progress and collect results
                with tqdm(total=len(futures), desc="Processing") as pbar:
                    for future in concurrent.futures.as_completed(futures):
                        try:
                            result = future.result()
                            results[futures[future]] = result
                        except Exception as e:
                            error_message = f"[ERROR] Task failed: {e}\n{traceback.format_exc()}"
                            print(error_message)
                        finally:
                            pbar.update(1)

            # Map first output path back to DataFrame (outcol mechanism unchanged)
            if outcol:
                input_data[outcol] = pd.NA
                for idx, _, resolved_outputs in tasks:
                    input_data.at[idx, outcol] = resolved_outputs[0] if resolved_outputs else None

            # Concatenate results if required
            if concat and results:
                if all(isinstance(r, pd.Series) for r in results.values()):
                    results_df = pd.concat(results, axis=1).T
                    return results_df
                else:
                    raise ValueError("All results must be pandas Series when concat=True")

            return results

        return wrapper
    return decorator


# # In your decorator:
# futures.append(executor.submit(task_wrapper, 
#                                "my_functions", 
#                                "mt_count_chrom_phased", 
#                                input_path, output_path))