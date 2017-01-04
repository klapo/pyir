def findXYZ(m, U, V, val, flag):
    #
    #  XYZ = findXYZ(m, UV, val, flag)
    #
    #  Routine to find new [x y z] location based on the image geometry, m,
    #  and image coordinates, UV.  The Direct Linear Transformation (DLT)
    #  equations, which transform the XYZ world coordinates into UV image
    #  coordinates:
    #
    #       U = (Ax + By + Cz + D)/(Ex + Fy + Gz + 1);
    #       V = (Hx + Jy + Kz + L)/(Ex + Fy + Gz + 1);
    #
    #  can be rearranged to solve for x, y, or z in terms of the image coordinates.
    #  Because the system is underdetermined, one of the coordinates must be
    #  specified.  The variable, flag, indicates which of the three coordinates
    #  is known: flag=1 for x, flag=2 for y, and flag=3 for z.  The value of the
    #  known coordinate is entered into the variable val, which can be a scalar or
    #  an array.
    #
    #  Inputs:
    #   m       - geometry vector (Walton m-vector = DLT coefficients)
    #   UV      - Nx2 matrix of [U V] coordinates
    #   val     - the specific value for the known x, y, or z coordinate
    #             If val is a scalar, it is used for all UV pairs.
    #             If val is an array, it must be as long as UV.
    #   flag    - flag indicating which of the x, y, or z coordinates is known,
    #             e.g., flag = 3 implies that z is known.
    #
    #  Output:
    #   XYZ     - Nx3 matrix of [x y z] values
    #
    # Copyright by Oregon State University, 2002
    # Developed through collaborative effort of the Argus Users Group
    # For official use by the Argus Users Group or other licensed activities.
    # Holman & Paden 04/28/04

    # Build val array if not already there.
    Nu = U.size
    Nv = V.size
    if not Nv == Nu:
        raise ValueError('U and V are not the same length')
    else:
        N = Nu
    if val.size == 1:
        val = np.ones_like(U) * val
    if not val.size == N:
        raise ValueError('val variable not same length as UV list')

    # Change from Walton m-vector notation to DLT notation so don't have to
    # use subscripts
    A, B, C, D, E, F, G, H, J, K, L = m

    # Assign variable names to coefficients derived in solving for x,y, or z
    M = E * U - A
    N = F * U - B
    O = G * U - C
    P = D - U
    Q = E * V - H
    R = F * V - J
    S = G * V - K
    T = L - V

    # Solve for unknown coordinates for given known coordinate

    if flag == 0:
        X = val
        Y = ((O * Q - S * M) * X + (S * P - O * T)) / (S * N - O * R)
        Z = ((N * Q - R * M) * X + (R * P - N * T)) / (R * O - N * S)

    elif flag == 1:
        Y = val
        X = ((O * R - S * N) * Y + (S * P - O * T)) / (S * M - O * Q)
        Z = ((M * R - Q * N) * Y + (Q * P - M * T)) / (Q * O - M * S)
    elif flag == 2:
        Z = val
        X = ((N * S - R * O) * Z + (R * P - N * T)) / (R * M - N * Q)
        Y = ((M * S - Q * O) * Z + (Q * P - M * T)) / (Q * N - M * R)
    else:
        raise ValueError('Invalid flag value (must be int of 0, 1, 2)')

    return(X, Y, Z)
