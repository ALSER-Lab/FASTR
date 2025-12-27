import re
import numpy as np


def get_scaling_equation(scaling_method: str, custom_formula=None, phred_alphabet_max=41) -> str:
    """Returns the mathematical equation string for each scaling method"""
    if scaling_method == 'custom' and custom_formula:
        return custom_formula
    elif scaling_method == 'none':
        return "x"
    elif scaling_method == 'log':
        return f"1 + 62*(ln(x)/ln({phred_alphabet_max + 1}))"
    elif scaling_method == 'log_reverse':
        return f"62-62*(ln({phred_alphabet_max + 1}-x)/ln({phred_alphabet_max + 1}))"
    elif scaling_method == 'linear':
        return "1+62*(x-40)/53"
    else:
        return "x"


def create_phred_quality_map(phred_offset=33, phred_alphabet_max=41):
    """Create mapping from ASCII quality characters to numeric quality scores"""
    phred_map = np.zeros(128, dtype=np.uint8)
    
    for ascii_val in range(128):
        quality_score = ascii_val - phred_offset
        quality_score = max(0, min(quality_score, (phred_alphabet_max-1)))
        phred_map[ascii_val] = quality_score
    
    return phred_map


def parse_custom_formula(formula: str, quality_scores: np.ndarray) -> np.ndarray:
    """
    Parse and evaluate a custom mathematical formula.
    Supports: +, -, *, /, **, (), ln(), log(), log10(), exp(), sqrt(), abs(), min(), max()
    Variable 'x' represents quality scores
    
    Example formulas:
      "1 + 62 * (x - 40) / 53" - Linear
      "ln(x) * 10" - Logarithmic variant
      "x ** 2 / 100" - Power scaling
      "1 + 62 * (ln(x - 39) / ln(54))" - Log
    """
    # Remove f(x)= prefix if present
    formula = re.sub(r'^\s*f\s*\(\s*x\s*\)\s*=\s*', '', formula.strip())
    
    # Create safe namespace with numpy functions
    safe_dict = {
        'x': quality_scores,
        'ln': np.log,
        'log': np.log,
        'log10': np.log10,
        'exp': np.exp,
        'sqrt': np.sqrt,
        'abs': np.abs,
        'min': np.minimum,
        'max': np.maximum,
        'np': np,
        '__builtins__': {}
    }
    
    try:
        # Replace common math notation
        formula = formula.replace('^', '**')
        
        # Evaluate formula
        with np.errstate(divide='ignore', invalid='ignore'):
            result = eval(formula, safe_dict)
        
        # Handle scalar results (broadcast to array)
        if np.isscalar(result):
            result = np.full_like(quality_scores, result, dtype=np.float32)
        
        return result.astype(np.float32)
    
    except Exception as e:
        print(f"ERROR: Failed to parse custom formula '{formula}'")
        print(f"Error details: {e}")
        print("\nSupported operations: +, -, *, /, ** (power), (), ln(), log(), log10(), exp(), sqrt(), abs()")
        print("Variable 'x' represents quality scores")
        print("\nExample formulas:")
        print("  '1 + 62 * (x - 40) / 53'")
        print("  'ln(x - 39) / ln(54) * 62 + 1'")
        print("  'x ** 2 / 100'")
        raise


def apply_clamping_to_equation(arr):
    """Clamp values to valid range [1, 63]"""
    return np.maximum(np.minimum(arr, 63), 1)


def apply_quality_to_bases(base_values, quality_scores, base_map, scaling_method='none',
                           custom_formula=None, phred_alphabet_max=41):
    """
    Apply quality scaling to base values.
    Transforms base encoding based on quality scores using various scaling methods.
    """
    # Straight up one-hot encoding
    if scaling_method == 'one_hot':
        return base_values
    
    # Calculate scale factors based on method
    if scaling_method == 'custom':
        if custom_formula is None:
            print("ERROR: custom scaling method requires --custom_formula argument")
            return base_values
        
        # Parse and evaluate custom formula
        scale_factors = apply_clamping_to_equation(parse_custom_formula(custom_formula, quality_scores))

    elif scaling_method == 'log':
        # Logarithmic scaling - makes larger quality scores closer in difference
        # f(x) = 1 + 62 * (log(x) / log(phred_alphabet_max+1))
        with np.errstate(divide='ignore', invalid='ignore'):
            LOG_DICT = np.clip(
                1 + 62 * (np.log(np.arange(0, phred_alphabet_max + 1, dtype=np.float32)) / 
                         np.log(phred_alphabet_max+1)), 
                1, 63
            ).astype(np.uint8)
            scale_factors = LOG_DICT[quality_scores]

    elif scaling_method == 'log_reverse':
        # Reverse logarithmic scaling - high quality -> small input to log
        with np.errstate(divide='ignore', invalid='ignore'):
            x = np.arange(0, phred_alphabet_max + 1, dtype=np.float32)
            LOG_REV_DICT = np.clip(
                62 - 62 * (np.log((phred_alphabet_max + 1) - x) / 
                          np.log(phred_alphabet_max + 1)), 
                1, 63
            ).astype(np.uint8)
            scale_factors = LOG_REV_DICT[quality_scores]

    elif scaling_method == 'linear':
        # Linear scaling
        scale_factors = 1 + 62 * (quality_scores - 40) / 53
        
    elif scaling_method == 'adaptive':
        # Linear scaling with min/max as bounds
        min_q = quality_scores.min()
        max_q = quality_scores.max()
        if max_q == min_q:
            return np.full_like(base_values, 32, dtype=np.uint8)
        scale_factors = 1 + 62 * (quality_scores - min_q) / (max_q - min_q)
        
    else:
        return base_values
    
    # Clip and convert to uint8
    scale_factors = np.clip(scale_factors, 0, 62).astype(np.uint8)
    
    # Lookup table approach for fast base scaling
    base_min_lookup = np.zeros(256, dtype=np.uint8)

    # N gets special handling w/ only 3 values (0-2)
    # For N we'll scale the quality to 0-2 range instead of 0-62
    n_mask = base_values == base_map[ord('N')]
    base_min_lookup[base_map[ord('N')]] = 0     # N: 0-2
    base_min_lookup[base_map[ord('A')]] = 3     # A: 3-65
    base_min_lookup[base_map[ord('G')]] = 66    # G: 66-128
    base_min_lookup[base_map[ord('C')]] = 129   # C: 129-191
    base_min_lookup[base_map[ord('T')]] = 192   # T: 192-254
    
    result = base_min_lookup[base_values] + scale_factors
    
    # We apply a mask to n that scales it down (ONLY TO N)
    if np.any(n_mask):
        n_scale = (scale_factors[n_mask] * 2 // 62).astype(np.uint8)  # Map 0-62 to 0-2
        result[n_mask] = base_min_lookup[base_map[ord('N')]] + n_scale
    
    return result.astype(np.uint8)


def validate_and_adjust_formula(formula: str, phred_alphabet_max: int) -> str:
    """
    Validate that a custom formula produces values in 1-63 range.
    Tests formula across quality score range and reports any clamping.
    """
    # Test the formula across the quality range
    test_qualities = np.arange(1, phred_alphabet_max + 1, dtype=np.float32)
    
    try:
        output = parse_custom_formula(formula, test_qualities)
        
        # Check the range AFTER clipping
        min_val = np.min(output)
        max_val = np.max(output)
        
        print(f"\nFormula validation for '{formula}':")
        print(f"  Input range: 1-{phred_alphabet_max}")
        print(f"  Output range: {min_val:.2f} to {max_val:.2f}")
        
        # Re-evaluate without clipping to see original range
        cleaned = re.sub(r'^\s*f\s*\(\s*x\s*\)\s*=\s*', '', formula.strip()).replace('^', '**')
        safe_dict = {
            'x': test_qualities,
            'ln': np.log, 'log': np.log, 'log10': np.log10, 'exp': np.exp,
            'sqrt': np.sqrt, 'abs': np.abs, 'min': np.minimum, 'max': np.maximum,
            'np': np, '__builtins__': {}
        }
        
        with np.errstate(divide='ignore', invalid='ignore'):
            raw_output = eval(cleaned, safe_dict)
        
        if np.isscalar(raw_output):
            raw_output = np.full_like(test_qualities, raw_output, dtype=np.float32)
        else:
            raw_output = raw_output.astype(np.float32)
        
        # Handle infinities and NaNs
        raw_output_clean = np.nan_to_num(raw_output, nan=0.0, posinf=100, neginf=-100)
        
        raw_min = np.min(raw_output_clean)
        raw_max = np.max(raw_output_clean)
        
        # Calculate how many values are clipped
        clamped_low = np.sum(raw_output_clean < 1)
        clamped_high = np.sum(raw_output_clean > 63)
        total_values = len(test_qualities)
        
        if raw_min < 1 or raw_max > 63:
            print(f"  WARNING: Formula produces values outside 1-63 range before clipping!")
            print(f"  Raw output range: {raw_min:.2f} to {raw_max:.2f}")
            print(f"  Values clamped: {clamped_low + clamped_high}/{total_values} " +
                  f"({100 * (clamped_low + clamped_high) / total_values:.1f}%)")
            
            if clamped_low > 0:
                print(f"    - {clamped_low} values below 1 (will be clamped to 1)")
            if clamped_high > 0:
                print(f"    - {clamped_high} values above 63 (will be clamped to 63)")
            
            print(f"  This will cause loss of information during reconstruction.")
        else:
            print(f"  Formula produces valid scaled values (1-63 range)")
        
        return formula
        
    except Exception as e:
        print(f"\nERROR: Failed to validate formula '{formula}'")
        print(f"Error details: {e}")
        print("\nSupported operations: +, -, *, /, ** (power), (), ln(), log(), log10(), exp(), sqrt(), abs()")
        print("Variable 'x' represents quality scores")
        raise