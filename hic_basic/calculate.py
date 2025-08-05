import os
import concurrent.futures
from functools import wraps
from tqdm import tqdm
import pandas as pd

def mt(
    force: bool = False,
    nproc: int = None,
    outcol: str = None,
    concat: bool = False
):
    """
    A decorator to parallelize functions with file skipping, path templating,
    progress bars, and result handling.

    Parameters:
        force (bool): Whether to overwrite existing output files.
        nproc (int): Number of parallel processes (default: all CPUs).
        outcol (str): Column name to store output file paths in the input DataFrame.
        concat (bool): If True, concatenate Series results into a DataFrame.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(input_data: pd.DataFrame, *args, **kwargs):
            nproc = 16
            nproc = nproc or 16  # Default to 16 processes if not specified
            # Extract dynamic input/output patterns from **kwargs
            input_pattern = kwargs.pop("input_pattern", None)
            output_pattern = kwargs.pop("output_pattern", None)
            input_col = kwargs.pop("input_col", "input_path")

            # check output
            if output_pattern is None and outcol is None:
                assert concat, "If no output_pattern and outcol are provided, concat must be True"

            # Store output paths and results
            output_paths = []
            results = []

            # Generate tasks with row indices
            tasks = []
            for idx, row in input_data.iterrows():
                sample_name = idx  # Assuming the index is the sample identifier
                input_path = input_pattern.format(sample_name=sample_name) if input_pattern else row[input_col]
                output_path = output_pattern.format(sample_name=sample_name) if not output_pattern is None else None

                # Skip existing files if not forced
                if output_path is not None:
                    if not force and os.path.exists(output_path):
                        continue

                tasks.append((idx, input_path, output_path))

            # Execute tasks in parallel
            with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
                futures = []
                for idx, input_path, output_path in tasks:
                    futures.append(executor.submit(func, input_path, output_path, *args, **kwargs))

                # Track progress and collect results
                with tqdm(total=len(futures), desc="Processing") as pbar:
                    for future in concurrent.futures.as_completed(futures):
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as e:
                            print(f"[ERROR] Task failed: {e}")
                        finally:
                            pbar.update(1)

            # Map results/output paths back to the original DataFrame
            for idx, input_path, output_path in tasks:
                output_paths.append((idx, output_path))

            # Sort by index to maintain order
            output_paths.sort(key=lambda x: x[0])
            output_paths_only = [path for _, path in output_paths]
            results_only = [res for res in results]

            # Update input DataFrame with output paths
            if outcol:
                input_data[outcol] = pd.NA
                for idx, output_path in output_paths:
                    input_data.at[idx, outcol] = output_path

            # Concatenate results if required
            if concat and results_only:
                if all(isinstance(r, pd.Series) for r in results_only):
                    results_df = pd.concat(results_only, axis=1).T
                    return results_df
                else:
                    raise ValueError("All results must be pandas Series when concat=True")

            return results_only

        return wrapper
    return decorator

# import importlib

# def task_wrapper(module_name, func_name, *args):
#     module = importlib.import_module(module_name)
#     func = getattr(module, func_name)
#     return func(*args)

# # In your decorator:
# futures.append(executor.submit(task_wrapper, 
#                                "my_functions", 
#                                "mt_count_chrom_phased", 
#                                input_path, output_path))