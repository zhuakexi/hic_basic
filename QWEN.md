# Code Style Guide

This part specifies the coding standards for Qwen3 codebase. All contributions **must** adhere to these guidelines to ensure consistency, readability, and maintainability.

---

## 1. Docstring Standards (Google Style)

All functions, classes, and modules **must** have docstrings following the **Google style** format. All documentation must be in **English**.

### 1.1 Basic Structure

```python
def function_name(param1: type, param2: type) -> return_type:
    """One-line summary of the function.

    Extended description of the function's behavior, purpose, and usage.

    Args:
        param1 (type): Description of parameter 1.
        param2 (type): Description of parameter 2.

    Returns:
        return_type: Description of return value.

    Raises:
        ValueError: When invalid input is provided.

    Example:
        >>> function_name(5, 10)
        15
    """
    pass
```

### 1.2 Key Rules

- Always use triple quotes `"""` for docstrings
- Always include:
  - One-line summary
  - Extended description
  - `Args` section (with types)
  - `Returns` section
  - `Raises` section (if applicable)
  - `Example` section (with doctest)
- Never omit type hints in docstrings
- Always end docstring with a period
- Never use `#` comments inside docstrings

## 2. Code Organization

### 2.1 File Structure

- Group related functions by logical functionality
- Label each block with `### --- Block Title --- ###`
- Add two blank lines above and below each block title.

```python
### --- Data Processing Functions --- ###

def clean_data(raw_data: list) -> list:
    ...

def transform_data(data: list) -> list:
    ...


### --- Mathematical Operations --- ###


def calculate_mean(values: list) -> float:
    ...

def calculate_variance(values: list) -> float:
    ...
```

### 2.2 Function Internal Structure

- Group related logic within a function
- Separate logic blocks with one blank line
- Label each block with `# --- Block Title --- #`

```python
def process_data(data: list) -> dict:
    """Process raw data into structured format."""
    
    # --- Validate Input Data ---
    if not data:
        raise ValueError("Input data cannot be empty")
    
    # --- Filter Valid Entries ---
    valid_entries = [item for item in data if item['value'] > 0]
    
    # --- Calculate Statistics ---
    total = sum(item['value'] for item in valid_entries)
    count = len(valid_entries)
    
    # --- Generate Report ---
    return {
        'total': total,
        'count': count,
        'average': total / count if count else 0
    }
```
## 3. Inline Comments

### 3.1 When to Add Comments

- **Complex algorithms**  
  (e.g., mathematical formulas, state machines)

- **Non-obvious logic**  
  (e.g., edge cases, optimizations)

- **Critical calculations**  
  (e.g., financial formulas, physics models)

### 3.2 Comment Style

- Always explain **why** the code exists, not just **what** it does
- Include examples when possible
- Keep comments concise (2-3 lines max)

```python
# Calculate BMI using formula: weight(kg) / (height(m) ** 2)
# Example: 70kg person with 1.75m height -> 70 / (1.75**2) = 22.86
bmi = weight / (height ** 2)

# Handle edge case: height must be > 0 (prevents division by zero)
if height <= 0:
    raise ValueError("Height must be greater than 0")
```

## 4. Commit Guidelines

All commits **must** follow a consistent format to ensure clear change tracking and review. Each commit message should:

- List modified functions **one per line**
- Specify modification type (**bug fix**, **feature**, **refactor** or **documentation**) in the first line
- Each function should be listed only once, with all changes related to it summarized in a single, concise sentence.

```bash
bug fix - calculate_discount: Fix edge case for discount_rate > 100
```