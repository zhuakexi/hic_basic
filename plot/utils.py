def filling_l2r_plotly(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Plotly_flavor, starts with 1.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i+1, j+1, k, feature
def filling_l2r_mpl(rows, cols, features):
    """
    Helper to iterate within row-cols. 
    Mpl_flavor, starts with 0.
    """
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            try:
                feature = features[k]
            except IndexError:
                feature = None
            yield i, j, k, feature