def find_longest_self_match(y, window=None, guaranteed_prefix_len=0, guarantee_loc=1):
    """
    Compute the length and location of the longest reoccurrence of a prefix of y.

    Inputs:
        y (str): A string.
        window (int, optional):  A positive integer for the window size for y.
                                      Effectively defaults to the length of y (all of y).                     
        guaranteed_prefix_len (int, optional): An integer for the prefix size already guaranteed.
                                        Used to exploit previous computation with smaller window.
                                        Defaults to 0 (no guarantee).
        guarantee_loc (int, optional): A positive integer for the location from where matches are sought.
                                       Used to exploit previous computation with smaller window.
                                       Defaults to 1 (sought anywhere).

    Outputs:
        tuple: the length and location of the reoccurrence.
    """
    if window== None:
        m = len(y)
    else:
        m = window

    attempted_len = guaranteed_prefix_len
    cursor = guarantee_loc
    
    prefix_len = 0
    prefix_loc = -1
    
    while cursor > 0:
        attempted_len = attempted_len + 1
        cursor = y.find(y[0:attempted_len], cursor, m)

        if cursor > -1:
            prefix_len = attempted_len
            prefix_loc = cursor
        else: break

    return [prefix_len, prefix_loc]

def find_longest_cross_match(y, x, window=None, guaranteed_prefix_len=0, guarantee_loc=0):
    """
    Compute the length and location of the longest reoccurrence of a prefix of y in x.

    Inputs:
        y (str): A string.
        x (str): A string.
        window (int, optional):  A positive integer for the window size for y.
                                      Effectively defaults to the length of y (all of y).                   
        guaranteed_prefix_len (int, optional): An integer for the minimum prefix sized to bee tested.
                                        Used to exploit previous computation with smaller window.
                                        Defaults to 0 (no guarantee).
        guarantee_loc (int, optional): A positive integer for the location from where matches are sought.
                                       Used to exploit previous computation with smaller window.
                                       Defaults to 0 (sought anywhere).

    Outputs:
        tuple: the length and location of the reoccurrence.
    """
    if window == None:
        m = len(y)
    else:
        m = window
    
    prefix_len = 0
    prefix_loc = -1

    attempted_len = guaranteed_prefix_len
    cursor = guarantee_loc
    
    while cursor > -1:
        attempted_len = attempted_len + 1
        cursor = y.find(x[0:attempted_len], cursor, m)

        if cursor > -1:
            prefix_len = attempted_len
            prefix_loc = cursor
        else: break

    return [prefix_len, prefix_loc]

def compute_match_ent_estimate(y, window=None, logbase=0):
    """
    Compute a longest-match estimate for the specific entropy of the law of Y from a sample y.

    Inputs:
        y (str): A string.
        window (int, optional):  A positive integer for the window size for y.
                                 Effectively defaults to the length of y (all of y).
        logbase (float, optional): Base for log operations. 
                                   Default is natural log.

    Outputs:
        float: longest-match estimation of the entropy.
    """

    if window == None:
        m = len(y)
    else:
        m = window

    import math

    if not logbase:
            logbase = math.e
    def log(a):
        return math.log(a) / math.log(logbase)

    return log(m) / (find_longest_self_match(y, m)[0])

def compute_match_cross_ent_estimate(y, x, window=None, logbase=0):
    """
    Compute longest-match estimator for the cross-entropy rate of the law of Y
    relative to the law of X

    Inputs:
        y (str): A sample y from Y.
        x (str): A sample x from X.
        window (int, optional):  A positive integer for the window size for y.
                                 Effectively defaults to the length of y (all of y).
        logbase (float, optional): Base for log operations. 
                                   Default is natural log.

    Outputs:
        float: longest-match estimation of the cross entropy.
    """

    if window == None:
        m = len(y)
    else:
        m = window

    import math

    if not logbase:
            logbase = math.e
    def log(a):
        return math.log(a) / math.log(logbase)

    return log(m) / (find_longest_cross_match(y, x, m)[0])

def compute_match_KL_estimate(y, x, window=None, logbase=0):
    """
    Compute longest-match estimator for the KL-divergence rate of the law of Y
    relative to the law of X

    Inputs:
        y (str): A sample y from Y.
        x (str): A sample x from X.
        window (int, optional):  A positive integer for the window size for y.
                                 Effectively defaults to the length of y (all of y).
        logbase (float, optional): Base for log operations. 
                                   Default is natural log.

    Outputs:
        float: longest-match estimation of the KL divergence.
    """
    return  compute_match_cross_ent_estimate(y, x, window, logbase) - compute_match_ent_estimate(y, window, logbase)

def compute_LZ_len(y):
    """
    Compute the number of words C in LZ-type parsing of y.

    Inputs:
        y (str): A string.

    Outputs:
        float: The number of words C.
    """

    dict = {}
    i = 0
    N = len(y)

    while i < N:
        cur_str = y[i]

        while cur_str in dict and i < N-1:
            i += 1
            cur_str += y[i]
        if cur_str not in dict:
            dict[cur_str] = None
        i += 1

    return len(dict)

def compute_ZM_len(y, x, mZM=True):
    """
    Compute the number of words C in ZM-type parsings of y using x.

    Inputs:
        y (str): The first string.
        x (str): The second string, possibly of different length.
        mZM (bool, optional): If True, proceeds with the modified ZM parsing. 
                              If False, proceeds with the original ZM parsing. 
                              Defaults to True.

    Outputs:
        float: The number of words C.
    """

    C, j, i = 1, 0, 0
    N = len(y)
    if mZM:
        N -=1
    
    while j < N:
        if y[i:j+1] in x:
            j += 1
        else:
            C += 1
            if mZM:
                j += 1
            else:
                if i == j: # to handle the case where supp y is not in supp x
                    j += 1
            i = j
    return C

def compute_LZ_ent_estimate(y, logbase=0):
    """
    Compute LZ estimate for the specific entropy of the law of Y from a sample y.

    Inputs:
        y (str): A finite string realization of Y.
        logbase (float, optional): Base for log operations. Default is natural log.

    Outputs:
        float: The LZ estimate of the specific entropy.
    """

    import math
    if not logbase:
            logbase = math.e
    def log(a):
        return math.log(a) / math.log(logbase)

    N = len(y)
    C = compute_LZ_len(y)

    return (C * log(C)) / N

def compute_ZM_cross_ent_estimate(y, x, mZM=True, logbase=0):
    """
    Compute ZM-type estimator Q for the cross-entropy rate of the law of Y
    relative to the law of X

    Inputs:
        y (str): The first string.
        x (str): The second string, possibly of different length.
        mZM (bool, optional): If True, computes the modified ZM estimator.
                              If False, computes the original ZM estimator.
                              Defaults to True.
        logbase (float, optional): Base for log operations. Default is natural log.

    Outputs:
        float: The estimator Q.
    """

    import math 

    if not logbase:
            logbase = math.e
    def log(a):
        return math.log(a) / math.log(logbase)

    N = len(y)

    C = compute_ZM_len(y, x, mZM)

    if mZM:
        N = N - C

    if N == 0:
        Q = math.inf # to handle division by 0
    else:
        Q = (C * log(len(x))) / N

    return Q

def compute_ZM_KL_estimate(y, x, mZM=True, logbase=0):
    """
    Compute an estimate of the KL-divergence of the law of Y relative to 
    the law of X using the ZM-type estimator Q for the cross-entropy rate 
    of the law of Y relative to the law of X and the LZ estimate for the 
    specific entropy of the law of Y

    Inputs:
        y (str): The first string.
        x (str): The second string, possibly of different length.
        mZM (bool, optional): If True, uses the modified ZM estimator.
                              If False, uses the original ZM estimator.
                              Defaults to True.
        logbase (float, optional): Base for log operations. Default is natural log.

    Outputs:
        float: The estimate of the KL-divergence.
    """

    return  compute_ZM_cross_ent_estimate(y, x, mZM, logbase) - compute_LZ_ent_estimate(y, logbase)

def compute_ZM_cross_ent_estimate_mult(y, x, Ns, Ms=[], mZM=True, verbose=True, logbase=0):
    """
    Compute ZM-type estimator Q for the cross-entropy rate of the law of Y
    relative to the law of X for different limits on the lengths

    Inputs:
        y (str): The first string.
        x (str): The second string.
        Ns (list of int): List of positive integers representing window sizes for y.
        Ms (list of int, optional): List of positive integers for window sizes for x.
                                    By default, or if given [], becomes Ns.
        mZM (bool, optional): If True, computes the modified ZM estimator. 
                              If False, computes the original ZM estimator. 
                              Defaults to True.
        logbase (float, optional): Base for log operations. Default is natural log.

    Outputs:
        list: A list of Q values.
    """

    Qs = []

    if len(Ms) == 0:
        Ms = Ns

    if not logbase:
        logbasestring = str(logbase)
    else:
        logbasestring = 'e'

    if verbose:
        if mZM: print(f'Computing mZM estimates (base '+logbasestring+').')
        else: print(f'Computing orignal ZM estimates (base '+logbasestring+').')

    L = min(len(Ms),len(Ns))

    if L != max(len(Ns),len(Ms)): print('Unequal lists of window sizes. Dropping some window sizes.')

    for ell in range(L):
        N = Ns[ell]
        M = Ms[ell]

        if N > len(y) or M > len(x):
            print(f'Window exceeds string length. Stopping early.')
            break
        
        Q = compute_ZM_cross_ent_estimate(y[:N], x[:M], mZM, logbase)

        if verbose: print(f'N: {N}, M: {M}, Q: {Q}')

        Qs.append(Q)

    if verbose: print('Finished.')
    
    return Qs

def compute_ZM_KL_estimate_mult(y, x, Ns, Ms=[], mZM=True, verbose=True, logbase=0):
    """
    Compute an estimate of the KL-divergence of the law of Y relative to 
    the law of X for different limits on the lengths. Uses the ZM-type 
    estimator Q for the cross-entropy rate of the law of Y relative to the
    law of X and the LZ estimate for the specific entropy of the law of Y

    Inputs:
        y (str): The first string.
        x (str): The second string.
        Ns (list of int): List of positive integers for window sizes for y.
        Ms (list of int, optional): List of positive integers for window sizes for x.
                          By default, or if given [], becomes Ns.
        mZM (bool, optional): If True, computes the modified ZM estimator. 
                              If False, computes the original ZM estimator. 
                              Defaults to True.
        logbase (float, optional): Base for log operations. Default is natural log.

    Outputs:
        list: A list of Q values.
    """

    KLs = []

    if len(Ms) == 0:
        Ms = Ns

    if not logbase:
        logbasestring = str(logbase)
    else:
        logbasestring = 'e'

    if verbose: print(f'Computing KL estimates (base '+logbasestring+').')

    L = min(len(Ms),len(Ns))

    if L != max(len(Ns),len(Ms)): print('Unequal lists of window sizes. Dropping some window sizes.')

    for ell in range(L):
        N = Ns[ell]
        M = Ms[ell]

        if N > len(y) or M > len(x):
            print(f'Window exceeds string length. Stopping early.')
            break
        
        KL = compute_ZM_KL_estimate(y[:N], x[:M], mZM, logbase)

        if verbose: print(f'N: {N}, M: {M}, KL: {KL}')

        KLs.append(KL)

    if verbose: print('Finished.')
    
    return KLs