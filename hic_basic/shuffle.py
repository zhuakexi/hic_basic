import random
def get_shape(nested_list):
    """
    Get shape of a nested list with any depth
    """
    if isinstance(nested_list, list) or isinstance(nested_list, tuple):
        return [len(nested_list)] + get_shape(nested_list[0])
    else:
        return []
def random_group(group, size_of_each_group, N):
    """
    Random sample each group N times.
    Input:
        group: list of elements
        size_of_each_group: size of each group selected
        N: number of groups
    Output:
        list of random groups
    """
    # make soure the size of each group is even
    size_of_each_group = int(size_of_each_group)
    results = []
    for i in range(N):
        random_group = random.sample(group, size_of_each_group)
        results.append(random_group)
    return results
def exchange_group(groupA, groupB, fraction=0.5):
    """
    Swith fraction of elements of groupA with groupB.
    Note:
        1. groupA and groupB must have same length
        2. result groups have no intersection
    Input:
        groupA: list of elements
        groupB: list of elements
        fraction: fraction of elements to switch
    Output:
        mix_groupA, mix_groupB
    """
    assert len(groupA) == len(groupB), "groupA and groupB must have same length"
    n_switch = int(len(groupA) * fraction)
    mix_groupA = groupA.copy()
    mix_groupB = groupB.copy()
    idxs = random.sample(range(len(groupA)), n_switch)
    for idx in idxs:
        mix_groupA[idx], mix_groupB[idx] = mix_groupB[idx], mix_groupA[idx]
    return mix_groupA, mix_groupB
def AB_mix(obsA:list, obsB:list, size_of_each_group=100, fraction=0.5, N=1000, seed=None):
    """
    Shuffle two groups. Will first random sample each group N times, then exchange fraction of elements for each pair.
    Both group-pair and mix-pair are returned, all pairs are non-overlapping.
    This is used to calculate distribution of group-to-group and mix-to-mix differences.
    Input:
        obsA: list of elements
        obsB: list of elements
            Note: obsA and obsB better have same length(TODO: check this in code)
        fraction: fraction of elements to switch
        size_of_each_group: size of each group selected
        N: number of groups
        seed: random seed
    Output:
        groupA, groupB, mix1, mix2
    """
    if seed is not None:
        random.seed(seed)
    groupAs = random_group(obsA, size_of_each_group, N)
    groupBs = random_group(obsB, size_of_each_group, N)
    #mixed_groups = mix_select_group(groupA, groupB, size_of_each_group, N, seed=True)
    mixed_groups = [
        exchange_group(groupA, groupB, fraction=0.5)
        for groupA, groupB in zip(groupAs, groupBs)
    ]
    mix1s, mix2s = zip(*mixed_groups)
    return groupAs, groupBs, mix1s, mix2s