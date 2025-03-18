import numpy as np
import pandas as pd

def align_structures(df1: pd.DataFrame, df2: pd.DataFrame, agg_chrom=True, ret3d: bool = False) -> np.ndarray:
    """Compute the homogeneous transformation matrix or transformed coordinates to align structure 2 to structure 1.

    Args:
        df1 (pd.DataFrame): DataFrame with MultiIndex (chr, pos) and columns ['x', 'y', 'z']
        df2 (pd.DataFrame): DataFrame with MultiIndex (chr, pos) and columns ['x', 'y', 'z']
        agg_chrom (bool, optional): Whether to aggregate points by chromosome (chr) first. Defaults to True.
        ret3d (bool, optional): If True, return transformed coordinates instead of the transformation matrix. Defaults to False.

    Returns:
        np.ndarray: 
            - If ret3d=False: 4x4 homogeneous transformation matrix M
            - If ret3d=True: Transformed coordinates of structure 2 (n_points × 3)
    
    Raises:
        ValueError: If input DataFrames have mismatched indices or columns

    Example:
        # Return transformation matrix
        M = align_structures(df1, df2)
        
        # Return transformed coordinates (aggregated by chromosome)
        Q_aligned = align_structures(df1, df2, ret3d=True, agg_chrom=True)
        
        # Verify equivalence
        assert np.allclose(Q_aligned, (M @ np.hstack([df2.values, np.ones((len(df2),1))]).T).T[:,:3])
    """
    # Validate and sort by MultiIndex
    if agg_chrom:
        orig_df1 = df1.copy()
        orig_df2 = df2.copy()
        df1 = df1.groupby(level=[0]).mean()
        df2 = df2.groupby(level=[0]).mean()
        check_index_levels = ["chr"]
    else:
        check_index_levels = ["chr", "pos"]
    
    # Ensure valid rows exist in both DataFrames
    valid_rows = df1.index.intersection(df2.index)
    if len(valid_rows) == 0:
        raise ValueError("No common indices between input DataFrames after aggregation")
    
    df1 = df1.loc[valid_rows]
    df2 = df2.loc[valid_rows]
    
    df1_sorted = df1.sort_index().reset_index()
    df2_sorted = df2.sort_index().reset_index()
    
    # Check index match
    if not df1_sorted.set_index(check_index_levels).index.equals(
        df2_sorted.set_index(check_index_levels).index
    ):
        raise ValueError("Input DataFrames have mismatched indices after aggregation")
    
    # Extract coordinates for transformation calculation
    P = df1_sorted[['x', 'y', 'z']].values
    Q = df2_sorted[['x', 'y', 'z']].values
    
    # Compute transformation matrix (steps 1-7)
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    
    H = Q_centered.T @ P_centered
    U, S, Vh = np.linalg.svd(H)
    V = Vh.T
    d = np.linalg.det(V @ U.T)
    R = V @ np.diag([1, 1, d]) @ U.T
    t = centroid_P - R @ centroid_Q
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = t
    
    # Apply transformation if requested
    if ret3d:
        if agg_chrom:
            # Use original coordinates (before aggregation)
            Q_original = orig_df2[['x', 'y', 'z']].values
            Q_homo = np.hstack([Q_original, np.ones((len(Q_original), 1))])
            Q_transformed = (M @ Q_homo.T).T[:, :3]
            # Preserve original index and columns
            transformed_df = pd.DataFrame(
                Q_transformed,
                index=orig_df2.index,
                columns=['x', 'y', 'z']
            )
            return transformed_df
        else:
            # Use non-aggregated coordinates (original indices)
            Q_homo = np.hstack([Q, np.ones((len(Q), 1))])
            Q_transformed = (M @ Q_homo.T).T[:, :3]
            transformed_df = pd.DataFrame(
                Q_transformed,
                index=df2_sorted.set_index(check_index_levels).index,
                columns=['x', 'y', 'z']
            )
            return transformed_df
    else:
        return M