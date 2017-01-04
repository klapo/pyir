def infoPRY2M(angles, camera, ip,
              oPitch=None, oRoll=None, oYaw=None):
    '''
    infoPRY2M -- convert ptich, roll yaw angles and camera/IP data to M vector

    m = info2M( angles, camera, ip ) converts the angle data into
    a Walton M vector based on camera location and IP data.

    Angles is a 4x1 or 1x4 array of angles in the order [pitch
    yaw hfov roll], in radians

    Camera is a 3 element vector of camera location values in the order
    [x y z].

    IP is a 6 element vector in order of [lx ly uWidth vWidth u0 v0].
    Or IP is the usual ip structure format.
    '''

    # set parameters
    pitch, yaw, fov, roll = angles
    cx, cy, cz = camera
    lu = ip[0]
    lv = ip[1]
    width = ip[2]
    u0 = ip[4]
    v0 = ip[5]

    # use camera space
    f = width / 2 / np.tan(fov / 2)

    # INS rotation matrix
    a = np.array([np.cos(yaw), np.sin(yaw), 0,
                  -np.sin(yaw), np.cos(yaw), 0,
                  0, 0, 1]).reshape(3, 3)
    b = np.array([np.cos(roll), 0, -np.sin(roll),
                  0, 1, 0,
                  np.sin(roll), 0, np.cos(roll)]).reshape(3, 3)
    c = np.array([1, 0, 0,
                  0, np.cos(pitch), np.sin(pitch),
                  0, -np.sin(pitch), np.cos(pitch)]).reshape(3, 3)

    # offset rotation to the cameras
    if oPitch is not None or oRoll is not None or oYaw is not None:
        a2 = np.array([np.cos(oYaw), np.sin(oYaw), 0,
                       -np.sin(oYaw), np.cos(oYaw), 0,
                       0, 0, 1]).reshape(3, 3)
        b2 = np.array([np.cos(oRoll), 0, -np.sin(oRoll),
                       0, 1, 0,
                       np.sin(oRoll), 0, np.cos(oRoll)]).reshape(3, 3)
        c2 = np.array([1, 0, 0,
                       0, np.cos(oPitch), np.sin(oPitch),
                       0, -np.sin(oPitch), np.cos(oPitch)]).reshape(3, 3)
        Rg = np.dot(np.dot(a2, b2), c2)
        # INS to local (ENU) rotation, with offset to camera rotaion
        R = np.dot(np.dot(np.dot(a, b), c), Rg).T
    else:
        # INS to local (ENU) rotation
        R = np.dot(np.dot(b, c), a)
    L = -(cx * R[2, 0] + cy * R[2, 1] + cz * R[2, 2])

    # m vector
    m = np.zeros(11)
    m[0] = lu * (u0 * R[2, 0] + f * R[0, 0]) / L
    m[1] = lu * (u0 * R[2, 1] + f * R[0, 1]) / L
    m[2] = lu * (u0 * R[2, 2] + f * R[0, 2]) / L

    m[3] = -(m[0] * cx + m[1] * cy + m[2] * cz)

    m[4] = R[2, 0] / L
    m[5] = R[2, 1] / L
    m[6] = R[2, 2] / L

    m[7] = lv * (v0 * R[2, 0] + f * R[1, 0]) / L
    m[8] = lv * (v0 * R[2, 1] + f * R[1, 1]) / L
    m[9] = lv * (v0 * R[2, 2] + f * R[1, 2]) / L

    m[10] = -(m[7] * cx + m[8] * cy + m[9] * cz)

    return m
