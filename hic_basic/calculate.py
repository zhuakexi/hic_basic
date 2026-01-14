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
        def wrapper(input_data: pd.DataFrame, force=False, nproc=None, concat=False, new_cols=None, *args, **kwargs):
            """
            Parallelized function wrapper with input/output pattern support.
            
            The function uses a pattern-based system for input and output file handling:
            - Input sources: input_cols, input_pattern, input_pattern1, input_pattern2, etc.
            - Output sources: output_cols, output_pattern, output_pattern1, output_pattern2, etc.
            
            Parameters are processed in order: input_cols first, then input_patterns.
            Patterns are formatted using {sample_name} placeholder which is replaced with
            the DataFrame index value for each row.
            
            Parameters:
                input_data (pd.DataFrame): Input DataFrame containing data to process
                force (bool): Whether to overwrite existing output files (default: False)
                nproc (int): Number of parallel processes (default: all CPUs)
                concat (bool): If True, concatenate Series results into a DataFrame
                new_cols (list): Column names for pattern-based output files (non-concat mode only)
                *args: Additional positional arguments passed to the decorated function
                **kwargs: Additional keyword arguments including input/output patterns
                
            Pattern Parameters (passed as kwargs):
                input_cols: Column name(s) containing input paths (single or list)
                output_cols: Column name(s) containing output paths (single or list)
                input_pattern: Template for input file paths
                input_pattern1, input_pattern2, ...: Additional input patterns
                output_pattern: Template for output file paths  
                output_pattern1, output_pattern2, ...: Additional output patterns
                
            Returns:
                If concat=True and results are pandas Series: DataFrame with concatenated results
                Otherwise: DataFrame with pattern-based output paths and execution status
            """
            njobs = nproc or 4

            # Extract input and output parameters in order
            input_sources = []
            output_sources = []
            
            # Extract input_cols first
            if "input_cols" in kwargs:
                input_cols = kwargs.pop("input_cols")
                if isinstance(input_cols, str):
                    input_sources.append(("col", [input_cols]))
                else:
                    input_sources.append(("col", input_cols))
            
            # Extract input_pattern, input_pattern1, input_pattern2...
            idx_num = 0
            while True:
                key = "input_pattern" if idx_num == 0 else f"input_pattern{idx_num}"
                if key in kwargs:
                    input_sources.append(("pattern", kwargs.pop(key)))
                    idx_num += 1
                else:
                    break

            if concat:
                assert new_cols is None, (
                    "new_cols parameter is not allowed when concat=True. This parameter is "
                    "only relevant for file-based workflows (concat=False)."
                )
                
                # Inform user that file existence checks are skipped in concat mode
                print("Note: In concat=True mode, all tasks will be executed regardless of "
                      "any existing output files. File existence checks are disabled as "
                      "concat mode is designed for computational results rather than file I/O.")

            # Extract output_cols next
            output_cols_list = []
            if "output_cols" in kwargs:
                output_cols = kwargs.pop("output_cols")
                if isinstance(output_cols, str):
                    output_cols_list = [output_cols]
                    output_sources.append(("col", output_cols_list))
                else:
                    output_cols_list = output_cols
                    output_sources.append(("col", output_cols_list))
            
            # Extract output_pattern, output_pattern1, output_pattern2...
            output_patterns = []
            idx_num = 0
            while True:
                key = "output_pattern" if idx_num == 0 else f"output_pattern{idx_num}"
                if key in kwargs:
                    pattern = kwargs.pop(key)
                    output_patterns.append(pattern)
                    output_sources.append(("pattern", pattern))
                    idx_num += 1
                else:
                    break
            if concat:
                # Prevent output pattern usage in concat mode
                assert len(output_patterns) == 0, (
                    "Output patterns (output_pattern, output_pattern1, etc.) are not allowed "
                    "when concat=True. Concat mode is designed for computational results, not "
                    "file generation. Use concat=False for file-based workflows."
                )

            # Check if we have any output sources
            if not output_sources and not concat:
                raise ValueError("Either output patterns/columns must be provided or concat must be True")
            
            # For non-concat mode, validate new_cols parameter
            if not concat and output_patterns:
                assert isinstance(new_cols, (list, type(None))), "new_cols must be a list or None"
                if new_cols is None:
                    raise ValueError("new_cols must be provided when using output patterns in non-concat mode")
                if len(new_cols) != len(output_patterns):
                    raise ValueError(f"new_cols length ({len(new_cols)}) must match number of output patterns ({len(output_patterns)})")

            # Store task information and results
            tasks_info = {}  # Will store task info by index
            execution_results = {}  # Will store execution results by index

            # Generate tasks with row indices
            tasks = []
            for idx, row in input_data.iterrows():
                sample_name = idx  # Assuming the index is the sample identifier

                # Resolve all input paths in order
                resolved_inputs = []
                for source_type, source_value in input_sources:
                    if source_type == "col":
                        for col in source_value:
                            resolved_inputs.append(row[col])
                    else:  # pattern
                        resolved_inputs.append(str(source_value).format(sample_name=sample_name))

                # Resolve all output paths in order
                resolved_outputs = []
                pattern_outputs = []  # Store only pattern-based outputs
                for source_type, source_value in output_sources:
                    if source_type == "col":
                        for col in source_value:
                            resolved_outputs.append(row[col])
                    else:  # pattern
                        output_path = str(source_value).format(sample_name=sample_name)
                        resolved_outputs.append(output_path)
                        pattern_outputs.append(output_path)

                # Determine if task should be skipped
                skip_task = False
                # Only perform skip check if concat=False AND outputs exist AND force=False
                if not concat and resolved_outputs and not force:
                    if all(os.path.exists(op) for op in resolved_outputs):
                        skip_task = True
                # When concat=True, skip_task remains False for all tasks

                # Store task information
                tasks_info[idx] = {
                    'resolved_inputs': resolved_inputs,
                    'resolved_outputs': resolved_outputs,
                    'pattern_outputs': pattern_outputs,
                    'skip_task': skip_task
                }
                
                if not skip_task:
                    tasks.append((idx, resolved_inputs, resolved_outputs))
            
            # Create directories for all output paths
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
                        idx = futures[future]
                        try:
                            result = future.result()
                            execution_results[idx] = {
                                'success': True,
                                'result': result,
                                'error': None
                            }
                        except Exception as e:
                            error_message = f"[ERROR] Task failed: {e}\n{traceback.format_exc()}"
                            print(error_message)
                            execution_results[idx] = {
                                'success': False,
                                'result': None,
                                'error': str(e)
                            }
                        finally:
                            pbar.update(1)

            # For skipped tasks, mark as success but don't store result
            for idx, info in tasks_info.items():
                if info['skip_task'] and idx not in execution_results:
                    execution_results[idx] = {
                        'success': True,
                        'result': None,
                        'error': None
                    }

            # Concatenate results if required
            if concat and execution_results:
                if all(isinstance(execution_results[idx]['result'], pd.Series) for idx in execution_results):
                    results_df = pd.concat([execution_results[idx]['result'] for idx in execution_results], axis=1).T
                    results_df.index = list(execution_results.keys())
                    return results_df
                else:
                    raise ValueError("All results must be pandas Series when concat=True")
            
            # For non-concat mode, create a DataFrame with pattern-based outputs
            if not concat:
                # Create result DataFrame with the same index as input_data
                result_df = pd.DataFrame(index=input_data.index)
                
                # Add pattern-based output columns
                if new_cols and output_patterns:
                    for col_name in new_cols:
                        result_df[col_name] = pd.NA
                
                # Check output files and populate result DataFrame
                for idx, info in tasks_info.items():
                    exec_result = execution_results.get(idx, {'success': False, 'result': None, 'error': 'Not executed'})
                    
                    # Check if all pattern outputs were created
                    outputs_created = True
                    if exec_result['success']:
                        for output_path in info['pattern_outputs']:
                            if not os.path.exists(output_path):
                                outputs_created = False
                                break
                    
                    # Populate result DataFrame
                    if outputs_created and exec_result['success']:
                        for i, col_name in enumerate(new_cols):
                            if i < len(info['pattern_outputs']):
                                result_df.at[idx, col_name] = info['pattern_outputs'][i]
                    else:
                        # Mark all outputs as NA if any check fails
                        for col_name in new_cols:
                            result_df.at[idx, col_name] = pd.NA
                
                # Append result_df to input_data
                input_data_with_results = pd.concat([input_data, result_df], axis=1, join='outer')
                return input_data_with_results

            return execution_results

        return wrapper
    return decorator
# # In your decorator:
# futures.append(executor.submit(task_wrapper, 
#                                "my_functions", 
#                                "mt_count_chrom_phased", 
#                                input_path, output_path))