import concurrent
from functools import wraps

from tqdm import tqdm

### --- single-pc --- ###
# def mt(n_workers=4):
#     """
#     e.g: mt(4)(func)([],c=3)
#     """
#     def decorator(func):
#         def wrapper(iterable, *args, **kwargs):
#             with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
#                 futures = [executor.submit(func, item, *args, **kwargs) for item in iterable]
#                 results = []
#                 for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
#                     results.append(future.result())
#                 return results
#         return wrapper
#     return decorator

def mt(func):
    @wraps(func)
    def wrapper(iterable, *args, num_workers=4, **kwargs):
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(func, item, *args, **kwargs) for item in iterable]
            results = []
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
                results.append(future.result())
            return results
    return wrapper