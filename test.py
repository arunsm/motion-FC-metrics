def calculate_FD_J(in_file):
    """
    @ Krsna
    May 2013
    compute
    1) Jenkinson FD from 3dvolreg's *.affmat12.1D file from -1Dmatrix_save option
    input: subject ID, rest_number, name of 6 parameter motion correction file (an output of 3dvolreg)
    output: FD_J.1D file
    Assumptions:    1) subject is available in BASE_DIR
    2) 3dvolreg is already performed and the 1D motion parameter and 1D_matrix file file is present in sub?/rest_? called as --->'lfo_mc_affmat.1D'
    """
    import math
    import os
    import numpy as np
    # TODO: update docstrings
    # in_file = CPAC's "coordinate transformation" resource
    out_file = os.path.join(os.getcwd(), 'FD_J.1D')
    pm_ = np.genfromtxt(in_file)
    pm = np.zeros((pm_.shape[0],pm_.shape[1]+4))
    pm[:, :12] = pm_
    pm[:, 12:] = [0.0, 0.0, 0.0, 1.0]
    # The default radius (as in FSL) of a sphere represents the brain
    rmax = 80.0
    # rigid body transformation matrix
    T_rb_prev = np.matrix(np.eye(4))
    out_lines = []
    for i in range(0, pm.shape[0]):
        # making use of the fact that the order of aff12 matrix is
        # "row-by-row"
        T_rb = np.matrix(pm[i].reshape(4,4))
        if not out_lines:
            out_lines.append('0')
        else:
            M = np.dot(T_rb, T_rb_prev.I) - np.eye(4)
            A = M[0:3, 0:3]
            b = M[0:3, 3]
            FD_J = math.sqrt((rmax*rmax/5)*np.trace(np.dot(A.T, A)) + np.dot(b.T, b))
            out_lines.append('\n{0:.8f}'.format(FD_J))
        T_rb_prev = T_rb
    with open(out_file, "w") as f:
        for line in out_lines:
            f.write(line)
    return out_file

out_file = calculate_FD_J('100408_Affine_Parameters.txt')